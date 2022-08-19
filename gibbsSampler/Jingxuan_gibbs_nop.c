// BINF 8500 assignment 5
// gibbs sampler
// Jingxuan Chen, Nov 18, 2019
// `gcc Jingxuan_gibbs.c pssm.c -lm -o gibbs` to compile

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>
//#include <unistd.h> // sleep(1)
#include "pssm.h" // the code from pssm assignment

#define READNAMELEN 100 // read name in fasta
#define SEQLEN 500 // max sequence length
#define SEQNUM 50 // max sequence count
// #define NSEED 3
// #define NUPDATE 300
#define ADJUST 5 // number of adjustment frequency

typedef struct
{
    char readName[READNAMELEN];
    char seq[SEQLEN];
    int seqInt[SEQLEN];
    int posMotif;
    int newposMotif;
} rawSeq; // raw sequences directly from file

typedef struct
{
    double pssm;
    double prob;
} posSeq; // store the probability for each motif

typedef struct
{
    double score;
    int seed;
    int lMotif;
    int posMotif[SEQNUM];
} seedOut; // the output for each seed

int NSEED, NUPDATE; // nseed and nupdate were user specified
rawSeq *inputSeq;
time_t now; // this variable for ramdom seeding

char *getMotif(rawSeq inputSeq, int lMotif);
int checkUpdate(rawSeq *inputSeq, int cSeq);
int updateMotif(rawSeq updateSeq, char **otherSeqMotif, int lMotif, int lSeq, int cSeq);
double samplerScore(rawSeq *inputSeq, int cSeq, int lSeq, int lMotif);
seedOut sampler(int seed, int cSeq, int lSeq, int lMotif);
int checkLeft(int cSeq, int lSeq, int lMotif, rawSeq *inputSeq);
int checkRight(int cSeq, int lSeq, int lMotif, rawSeq *inputSeq);

int main(int argc, char *argv[])
{
    FILE *fp;
    char *fileName, c, *motif;
    int lMotif, lSeq, cSeq, i, j, seed;
    // rawSeq *inputSeq;
    seedOut *samplerOut;

    if (argc != 5)
    {
        printf("\n\nUsage: ./gibbs sequenceFile estimatedMotifLength seed update\n\n");
        printf("e.g. time ./gibbs E.coliRpoN-sequences-16-100nt.fasta 20 200 200\n\n");
        return -1;
    }

    fileName = argv[1];
    fp = fopen(fileName, "r");
    if (fp == NULL)
    {
        printf("\n\nUsage: ./gibbs sequenceFile estimatedMotifLength \n\n");
        printf("e.g. time ./gibbs E.coliRpoN-sequences-16-100nt.fasta 20\n\n");
        perror("Invalid input file!\n\n");
        return -1;
    }
    lMotif = atoi(argv[2]); // get the length of motif from user
    NSEED = atoi(argv[3]); // get seed and update number from command line
    NUPDATE = atoi(argv[4]); // to test parameter

    /* calculate the number of sequences in the input file */
    cSeq = 0;
    while (!feof(fp))
    {
        if ((c = fgetc(fp)) == '>')
            cSeq++;
    }
    // printf("The count for squences: %d\n", cSeq);
    /**** read input sequences ****/
    inputSeq = (rawSeq *)malloc(cSeq * sizeof(rawSeq));
    fseek(fp, 0, SEEK_SET);
    fgetc(fp); // get rid of the first '>'
    for (i = 0; i < cSeq; i++)
    {
        fscanf(fp, "%s", inputSeq[i].readName);
        // fgetc(fp);
        j = 0;
        while (((c = fgetc(fp)) != '>') && !feof(fp))
        {
            if (isalpha(c) == 0)
                continue;
            inputSeq[i].seq[j] = toupper(c);
            j++;
        }
        inputSeq[i].seq[j] = '\0';
    }
    lSeq = j;
    // debug
    // printf("lsq: %d\n", lSeq);
    // for (i = 0; i < cSeq; i++)
    // {
    //     // printf("%s  %s\n", inputSeq[i].readName, inputSeq[i].seq);
    //     for (j = 0; j < lSeq; j++)
    //         printf("[%d]%c", (j + 1), inputSeq[i].seq[j]);
    //     printf("\n");
    // }
    fclose(fp);
    // generate seqInt for each input seq
    for (i = 0; i < cSeq; i++)
    {
        for (j = 0; j < lSeq; j++)
        {
            inputSeq[i].seqInt[j] = getIndex(inputSeq[i].seq[j]);
            // printf("%d", inputSeq[i].seqInt[j]);
        }
        // printf("\n");
    }

    /******** Test different initializations to find the global max ********/
    samplerOut = (seedOut *)malloc(NSEED * sizeof(seedOut));
    now = time(NULL); // record the time when sampler starts
    for (seed = 0; seed < NSEED; seed++)
    {
        samplerOut[seed] = sampler(seed, cSeq, lSeq, lMotif);
        // printf("Seed %d finished with score %lf\n", seed, samplerOut[seed].score);
        // sleep(1); // in order to generate different seed (time()) for each test,
                    // use sleep to make sure the time interval is greater than 1 sec 
    }

    /***************** Find the max seed ********************/
    int finalLmotif; // the motif positions and length of motif that generate maximum score
    double maxScore = samplerOut[0].score; // initialize 
    finalLmotif = samplerOut[0].lMotif;
    for (j = 0; j < cSeq; j++)
        inputSeq[j].newposMotif = samplerOut[0].posMotif[j]; // store corresponding motif position

    for (i = 1; i < NSEED; i++)
    {
        if (samplerOut[i].score > maxScore)
        {
            maxScore = samplerOut[i].score; // store the max score
            for (j = 0; j < cSeq; j++)
                inputSeq[j].newposMotif = samplerOut[i].posMotif[j]; // store corresponding motif position
            finalLmotif = samplerOut[i].lMotif; // store the resulting length of motif
        }
    }
    /**** results output ****/
    printf("Global maximum score is %lf\n", maxScore);
    motif = (char *)malloc((finalLmotif + 1) * sizeof(char));
    printf("The final motif length is %d\n", finalLmotif);
    printf("Seq\tStartPos\tEndPos\tSeqName\n");
    for (i = 0; i < cSeq; i++)
    {
        strncpy(motif, getMotif(inputSeq[i], finalLmotif), finalLmotif);
        motif[finalLmotif] = '\0';
        printf("%s\t%d\t%d\t%s\n", motif, (inputSeq[i].newposMotif + 1), (inputSeq[i].newposMotif + finalLmotif), inputSeq[i].readName);
    }

    free(inputSeq);
    free(samplerOut);
    return 0;
}
char *getMotif(rawSeq inputSeq, int lMotif) // function to extract motif from input data structure
{
    char *seqMotif, *pSeq;
    pSeq = inputSeq.seq;
    seqMotif = (char *)malloc((lMotif + 1) * sizeof(char));
    seqMotif = pSeq + inputSeq.newposMotif; // get the motif seq for this seq
    seqMotif[lMotif] = '\0';                // set the last char of seq
    // printf("the motif is: %s\n", seqMotif);
    return seqMotif; // note that here it just returns an address
}

/*  if stop changing return 1, else -1 */
// int checkUpdate(rawSeq *inputSeq, int cSeq)
// {
//     int i;
//     for (i = 0; i < cSeq; i++)
//     {
//         if(inputSeq[i].posMotif != inputSeq[i].newposMotif)
//             return -1;
//     }
//     return 1;
// }

/*** function to update motif position ****/
int updateMotif(rawSeq updateSeq, char **otherSeqMotif, int lMotif, int lSeq, int cSeq)
{
    double *seqBg;
    double **freqMatrix, **probMatrix, **scoreMatrix;
    double score, p, sumP = 0.0, len = 0.0, randNum;
    int i, maxpos= lSeq - lMotif, *thisSeqInt, *pseqInt;
    posSeq *posUpdateSeq;

    /***** PSSM *****/
    seqBg = (double *)malloc(4 * sizeof(double));
    seqBg = bg(updateSeq.seqInt, lSeq);

    freqMatrix = (double **)malloc(lMotif * sizeof(double *));
    for (i = 0; i < lMotif; i++)
        freqMatrix[i] = (double *)malloc(4 * sizeof(double));
    freqMatrix = fMatrix(otherSeqMotif, lMotif, (cSeq - 1)); // calculate the freq matrix other than the update seq

    probMatrix = (double **)malloc(lMotif * sizeof(double *));
    for (i = 0; i < lMotif; i++)
        probMatrix[i] = (double *)malloc(4 * sizeof(double));
    probMatrix = freq2prob(freqMatrix, lMotif);

    scoreMatrix = (double **)malloc(lMotif * sizeof(double *));
    for (i = 0; i < lMotif; i++)
        scoreMatrix[i] = (double *)malloc(4 * sizeof(double));
    scoreMatrix = prob2score(probMatrix, lMotif, seqBg);

    /// debug
    // printf("PSSM:\npos\tA\tC\tG\tT\t\n");
    // for (i = 0; i < lMotif; i++)
    // {
    //     printf("%d\t", i);
    //     for (j = 0; j < 4; j++)
    //     {
    //         printf("%.3lf\t", scoreMatrix[i][j]);
    //     }
    //     printf("\n");
    // }

    /**** scan through each position ****/
    posUpdateSeq = (posSeq *)malloc((maxpos + 1) * sizeof(posSeq)); // maxpos is the index
    thisSeqInt = (int *)malloc(lMotif * sizeof(int));
    for (i = 0; i <= maxpos; i++) // maxpos is the index
    {
        pseqInt = updateSeq.seqInt;
        thisSeqInt = pseqInt + i;
        score = seqScore(scoreMatrix, lMotif, thisSeqInt);
        posUpdateSeq[i].pssm = score;
        p = pow(2, score); // covert score to probability to avoid negative scores
        posUpdateSeq[i].prob = p;

        /* calculate the sum of probabitity for proportional sampling */
        sumP += p;
        // printf("pssm, p, sumP: %lf, %lf, %lf ", score, p, sumP);
    }
    //printf("%lf ", sumP);
    // this is the old score
    score = posUpdateSeq[updateSeq.newposMotif].pssm;
    randNum = drand48(); // seed was specified in sampler() func
    for (i = 0; i <= maxpos; i++) 
    {
        len += (posUpdateSeq[i].prob / sumP);
        if (len > randNum)
            return i;
    }
    return -1;
}

/****** function to calculate final score for a specifc set of ******/
double samplerScore(rawSeq *inputSeq, int cSeq, int lSeq, int lMotif)
{
    char **pssmMotif, *motif;
    int i, j, k;
    double *seqBg, score = 0.0, tmp;
    double **freqMatrix, **probMatrix, **scoreMatrix;
    
    pssmMotif = (char **)malloc((cSeq - 1) * sizeof(char *));
    for (k = 0; k < (cSeq - 1); k++)
        pssmMotif[k] = (char *)malloc((lMotif + 1) * sizeof(char));
    seqBg = (double *)malloc(4 * sizeof(double));
    freqMatrix = (double **)malloc(lMotif * sizeof(double *));
    for (k = 0; k < lMotif; k++)
        freqMatrix[k] = (double *)malloc(4 * sizeof(double));
    scoreMatrix = (double **)malloc(lMotif * sizeof(double *));
    for (k = 0; k < lMotif; k++)
        scoreMatrix[k] = (double *)malloc(4 * sizeof(double));
    probMatrix = (double **)malloc(lMotif * sizeof(double *));
    for (k = 0; k < lMotif; k++)
        probMatrix[k] = (double *)malloc(4 * sizeof(double));
    
    for (i = 0; i < cSeq; i++)
    {
        // get the motif matrix to calculate pssm


        for (j = 0; j < cSeq; j++)
        {
            if (j < i)
            {
                strncpy(pssmMotif[j], getMotif(inputSeq[j], lMotif), lMotif);
                pssmMotif[j][lMotif] = '\0';
                // printf("motif %d is %s\n", j, pssmMotif[j]);
            }
            else if (j == i)
                continue;
            else if (j > i)
            {
                strncpy(pssmMotif[j - 1], getMotif(inputSeq[j], lMotif), lMotif);
                pssmMotif[j - 1][lMotif] = '\0';
                // printf("motif %d is %s\n", j, pssmMotif[j - 1]);
            }
        }
        /***** PSSM *****/

        seqBg = bg(inputSeq[i].seqInt, lSeq);

        freqMatrix = fMatrix(pssmMotif, lMotif, (cSeq - 1)); // calculate the freq matrix other than the update seq

        probMatrix = freq2prob(freqMatrix, lMotif);

        scoreMatrix = prob2score(probMatrix, lMotif, seqBg);

        // the score for this motif
        motif = (char *)malloc((lMotif + 1) * sizeof(char));
        strncpy(motif, getMotif(inputSeq[i], lMotif), lMotif);
        motif[lMotif] = '\0';
        // printf("motif for %d is %s\n", i, motif);
        tmp = motifScore(scoreMatrix, lMotif, motif);
        score += tmp;
        // printf("sum score is %lf, this score is %lf\n", score, tmp);
    }
    return score;
}

/***** sampling function for each seed *****/
seedOut sampler(int seed, int cSeq, int lSeq, int lMotif_ori)
{
    int i, j, upd, tmplMotif, lMotif, maxMotif[cSeq];
    char **pssmMotif;
    seedOut samplerOut;
    rawSeq *tmpSeq;
    double score, tmpscore, tmpscore_l, tmpscore_r;

    lMotif = lMotif_ori; // the originial len of motif from args

    pssmMotif = (char **)malloc((cSeq - 1) * sizeof(char *));
    for (i = 0; i < (cSeq - 1); i++)
        pssmMotif[i] = (char *)malloc((lMotif + 1) * sizeof(char));

    srand48((long)now + seed); // use this weird expression to make sure
                            // using different seed every time run the program
                            // and do not need to wait 1 sec
    /**** initialize motif locs ****/
    for (i = 0; i < cSeq; i++)
    {
        // the last possible start point of motif is lSeq - lMotif
        inputSeq[i].newposMotif = (int)(drand48() * (lSeq - lMotif));
        // printf("the pos is %d\n", inputSeq[i].posMotif);
    }

    // calculate the initialized motif score
    score = tmpscore = samplerScore(inputSeq, cSeq, lSeq, lMotif);

    /**** modify the position of motif until stop changing ****/
    for (upd = 0; upd < NUPDATE; upd++) // update until some limit time
    {
        /*** assign old motif pos to `new` variable ***/
        // printf("update start:\n");
        for (i = 0; i < cSeq; i++)
            inputSeq[i].posMotif = inputSeq[i].newposMotif;

        // motif = (char *)malloc((lMotif + 1) * sizeof(char)); // debug
        // for (i = 0; i < cSeq; i++)
        // {
        //     motif = getMotif(inputSeq[i], lMotif);
        //     motif[lMotif] = '\0';
        //     // strncpy(motif, getMotif(inputSeq[i], lMotif), lMotif);
        //     printf("%s\t%d\t%s\n", inputSeq[i].readName, inputSeq[i].newposMotif, motif);
        // }

        /**** update locs for each sequence ****/
        for (i = 0; i < cSeq; i++)
        {
            for (j = 0; j < cSeq; j++)
            {
                // generate the motifs for updating
                if (j < i)
                {
                    strncpy(pssmMotif[j], getMotif(inputSeq[j], lMotif), lMotif);
                    pssmMotif[j][lMotif] = '\0';
                    // printf("motif %d is %s\n", j, pssmMotif[j]);
                }
                else if (j == i)
                    continue;
                // printf("This is the seq being updated.\n");
                else if (j > i)
                {
                    strncpy(pssmMotif[j - 1], getMotif(inputSeq[j], lMotif), lMotif);
                    pssmMotif[j - 1][lMotif] = '\0';
                    // printf("motif %d is %s\n", j, pssmMotif[j - 1]);
                }
            }
            // for (j = 0; j < (cSeq - 1); j++)
            // {
            //     printf("motif %d is %s\n", j, pssmMotif[j]);
            // }
            inputSeq[i].newposMotif = updateMotif(inputSeq[i], pssmMotif, lMotif, lSeq, cSeq);
            // printf("seq%d: new: %d old: %d\n", i, inputSeq[i].newposMotif, inputSeq[i].posMotif);
        }
//        tmpscore = samplerScore(inputSeq, cSeq, lSeq, lMotif);
//        printf("%d %lf\n", upd, tmpscore);

        /****** adjust the border of motif ********/
        // if (0)
        if (upd % ADJUST == 0)
        {
            // printf("lMotif now is %d, lSeq is %d\t", lMotif, lSeq);
            tmpSeq = inputSeq;
            /***** firstly move the left border *****/
            // check left border
            if (checkLeft(cSeq, lSeq, lMotif, tmpSeq) == -1)
                tmpscore_l = -1000; // could not move to left
            else
            {
                // left border move to left
                for (i = 0; i < cSeq; i++)
                    tmpSeq[i].newposMotif--;
                tmplMotif = lMotif + 1;
                tmpscore_l = samplerScore(tmpSeq, cSeq, lSeq, tmplMotif);
            }
            // left border move to right
            tmpSeq = inputSeq;
            for (i = 0; i < cSeq; i++)
                tmpSeq[i].newposMotif++;
            tmplMotif = lMotif - 1;
            tmpscore_r = samplerScore(tmpSeq, cSeq, lSeq, tmplMotif);
            if ((tmpscore_l > tmpscore_r) && (tmpscore_l > score))
            {
                score = tmpscore_l;
                lMotif++;
                for (i = 0; i < cSeq; i++)
                    inputSeq[i].newposMotif--;
            }
            else if ((tmpscore_r > tmpscore_l) && (tmpscore_r > score))
            {
                score = tmpscore_r;
                lMotif--;
                for (i = 0; i < cSeq; i++)
                    inputSeq[i].newposMotif++;
            }
            // printf("left border l=%lf, r=%lf, score=%lf\n", tmpscore_l, tmpscore_r, score);
            /***** then move the right border *****/
            tmpSeq = inputSeq;
            // printf("can I move to right?: %d\n", checkRight(cSeq, lSeq, lMotif, tmpSeq));
            if (checkRight(cSeq, lSeq, lMotif, tmpSeq) == -1)
                tmpscore_r = -1000;
            else
            {
                // right border move to right
                tmpscore_r = samplerScore(tmpSeq, cSeq, lSeq, (lMotif + 1));
            }
            // right border move to left
            tmpscore_l = samplerScore(tmpSeq, cSeq, lSeq, (lMotif - 1));
            if ((tmpscore_l > tmpscore_r) && (tmpscore_l > score))
            {
                score = tmpscore_l;
                lMotif--;
            }
            else if ((tmpscore_r > tmpscore_l) && (tmpscore_r > score))
            {
                score = tmpscore_r;
                lMotif++;
            }
            // printf("right border l=%lf, r=%lf, score=%lf\n", tmpscore_l, tmpscore_r, score);
            // printf("lmotif after adjust: %d\n", lMotif);
        }

        /*** after one update, calculate the new score and 
         * check if it is larger than previous
         * if yes, save the new motif info,
         * if not, continue next update ***/
        tmpscore = samplerScore(inputSeq, cSeq, lSeq, lMotif);
        if (tmpscore > score)
        {
            score = tmpscore;
            for (i = 0; i < cSeq; i++)
                maxMotif[i] = inputSeq[i].newposMotif;
        }
        // motif = (char *)malloc((lMotif + 1) * sizeof(char));
        // for (i = 0; i < cSeq; i++)
        // {
        //     motif = getMotif(inputSeq[i], lMotif);
        //     motif[lMotif] = '\0';
        //     // strncpy(motif, getMotif(inputSeq[i], lMotif), lMotif);
        //     printf("%s\t%d\t%s\n", inputSeq[i].readName, inputSeq[i].newposMotif, motif);
        // }
        // printf("\n");
    }
    // printf("seed update OK\n");

    /******* At the end of each seed, calculate the score. *******/
    for (i = 0; i < cSeq; i++)
        inputSeq[i].newposMotif = maxMotif[i];

    score = samplerScore(inputSeq, cSeq, lSeq, lMotif);
    // printf("score for seed %d is %lf\n", seed, score);
    // store the out statistics for each seed iteration.
    samplerOut.score = score;
    samplerOut.seed = seed;
    samplerOut.lMotif = lMotif;
    for (i = 0; i < cSeq; i++)
        samplerOut.posMotif[i] = maxMotif[i];
    // motif = (char *)malloc((lMotif + 1) * sizeof(char));
    // for (i = 0; i < cSeq; i++)
    // {
    //     motif = getMotif(inputSeq[i], lMotif);
    //     motif[lMotif] = '\0';
    //     // strncpy(motif, getMotif(inputSeq[i], lMotif), lMotif);
    //     printf("%s\t%d\t%s\n", inputSeq[i].readName, inputSeq[i].newposMotif, motif);
    // }
    // printf("\n");
    return samplerOut;
}

/**** function to check if left or right border could be moved ***/
int checkLeft(int cSeq, int lSeq, int lMotif, rawSeq *inputSeq)
{
    int i;
    for (i = 0; i < cSeq; i++)
    {
        if (inputSeq[i].newposMotif <= 0) // could not move to left
            return -1;
    }
    return 1;
}

int checkRight(int cSeq, int lSeq, int lMotif, rawSeq *inputSeq)
{
    int i;
    for (i = 0; i < cSeq; i++)
    {
        if (inputSeq[i].newposMotif >= (lSeq - lMotif)) // could not move to right
            return -1;
    }
    return 1;
}
