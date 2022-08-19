#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <pthread.h>
#include <time.h>
#include "pssm.h" // the code from pssm assignment

#define READNAMELEN 100
#define SEQLEN 500
// #define NSEED 3
// #define NUPDATE 200
#define ADJUST 5
// #define STOPUPDATE 1
int NSEED, NUPDATE;
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
} seedOut; // the output for each seed

rawSeq *inputSeq;

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
    int lMotif, lSeq, cSeq, i, j, seed, upd;
    char **pssmMotif;
    double *bgProb, score;
    // rawSeq *inputSeq;
    seedOut *samplerOut, finalOut;

    if (argc != 5)
    {
        printf("\n\nUsage: ./gibbs sequenceFile estimatedMotifLength \n\n");
        printf("e.g. time ./gibbs E.coliRpoN-sequences-16-100nt.fasta 20\n\n");
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
    NSEED = atoi(argv[3]);
    NUPDATE = atoi(argv[4]);
    /* calculate the number of sequences in the input file */
    cSeq = 0;
    while (!feof(fp))
    {
        if ((c = fgetc(fp)) == '>')
            cSeq++;
    }
    // printf("The count for squences: %d\n", cSeq);
    inputSeq = (rawSeq *)malloc(cSeq * sizeof(rawSeq));
    fseek(fp, 0, SEEK_SET);
    fgetc(fp); // geet rid of the first '>'
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
    samplerOut = (seedOut *)malloc(NSEED * sizeof(seedOut));

    // /* run nseeds times gibs algorithem to find the local max*/
    // pssmMotif = (char**)malloc((cSeq - 1) * sizeof(char*));
    // for (i = 0; i < (cSeq - 1); i++)
    //     pssmMotif[i] = (char*)malloc((lMotif + 1) * sizeof(char));

    // // test samplerScore with jan's output
    // inputSeq[0].newposMotif = 31;
    // inputSeq[1].newposMotif = 32;
    // inputSeq[2].newposMotif = 42;
    // inputSeq[3].newposMotif = 36;
    // inputSeq[4].newposMotif = 53;
    // inputSeq[5].newposMotif = 3;
    // inputSeq[6].newposMotif = 52;
    // inputSeq[7].newposMotif = 48;
    // inputSeq[8].newposMotif = 32;
    // inputSeq[9].newposMotif = 48;
    // inputSeq[10].newposMotif = 43;
    // inputSeq[11].newposMotif = 45;
    // inputSeq[12].newposMotif = 72;
    // inputSeq[13].newposMotif = 11;
    // inputSeq[14].newposMotif = 39;
    // inputSeq[15].newposMotif = 29;
    // score = samplerScore(inputSeq, cSeq, lSeq, lMotif);
    // printf("sample score is: %lf\n", score); // should be 287

    /******** Test different initializations to find the global max ********/
    for (seed = 0; seed < NSEED; seed++)
    {
        samplerOut[seed] = sampler(seed, cSeq, lSeq, lMotif);
        // printf("score for seed %d is %lf\n", seed, samplerOut[seed].score);
    }

    /***************** Find the max seed ********************/
    int maxScore = samplerOut[0].score, maxSeed = samplerOut[0].seed;
    for (i = 1; i < NSEED; i++)
    {
        if (samplerOut[i].score > maxScore)
        {
            maxScore = samplerOut[i].score;
            maxSeed = samplerOut[i].seed;
        }
    }

    /***************** Rerun sampler for max seed and print out results ********************/
    finalOut = sampler(maxSeed, cSeq, lSeq, lMotif);
    printf("score for seed %d is %lf\n", maxSeed, finalOut.score);
    printf("the final len of motif is: %d\n", finalOut.lMotif);
    motif = (char *)malloc((finalOut.lMotif + 1) * sizeof(char));
    for (i = 0; i < cSeq; i++)
    {
        // motif = getMotif(inputSeq[i], finalOut.lMotif);
        strncpy(motif, getMotif(inputSeq[i], finalOut.lMotif), finalOut.lMotif);
        motif[finalOut.lMotif] = '\0';
        printf("%s\t%d\t%s\n", inputSeq[i].readName, inputSeq[i].newposMotif, motif);
    }

    return 0;
}
char *getMotif(rawSeq inputSeq, int lMotif)
{
    char *seqMotif, *pSeq;
    pSeq = inputSeq.seq;
    seqMotif = (char *)malloc((lMotif + 1) * sizeof(char));
    seqMotif = pSeq + inputSeq.newposMotif; // get the motif seq for this seq
    seqMotif[lMotif] = '\0';                // set the last char of seq
    // printf("the motif is: %s\n", seqMotif);
    return seqMotif;
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

int updateMotif(rawSeq updateSeq, char **otherSeqMotif, int lMotif, int lSeq, int cSeq)
{
    double *seqBg;
    double **freqMatrix, **probMatrix, **scoreMatrix;
    double score, tmps, p, sumP, len, randNum;
    int i, j, maxpos = lSeq - lMotif, *thisSeqInt, *pseqInt, count;
    posSeq *posUpdateSeq;
    int updatedMotif;

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
    posUpdateSeq = (posSeq *)malloc((maxpos) * sizeof(posSeq));
    thisSeqInt = (int *)malloc(lMotif * sizeof(int));
    sumP = 0.0;
    for (i = 0; i <= maxpos; i++) // maxpos is the index
    {
        pseqInt = updateSeq.seqInt;
        thisSeqInt = pseqInt + i;
        score = seqScore(scoreMatrix, lMotif, thisSeqInt);
        posUpdateSeq[i].pssm = score;
        p = pow(2, score); // covert score to probability
        posUpdateSeq[i].prob = p;

        /* calculate the sum of probabitity for proportional sampling */
        sumP += p;
        // printf("pssm, p, sumP: %lf, %lf, %lf ", score, p, sumP);
    }
    //printf("%lf ", sumP);
    // this is the old score
    score = posUpdateSeq[updateSeq.newposMotif].pssm;
    randNum = drand48();
    len = 0;
    for (i = 0; i <= maxpos; i++)
    {
        len += (posUpdateSeq[i].prob / sumP);
        if (len > randNum)
        {
            updatedMotif = i;
            return updatedMotif;
        }
    }
    // count = 0;
    // while (count <= STOPUPDATE)
    // {
    //     randNum = drand48();
    //     for (i = 0; i <= maxpos; i++)
    //     {
    //         len += (posUpdateSeq[i].prob / sumP);
    //         // printf("len = %lf; i = %d \n", len, i);
    //         if (len > randNum)
    //         {              
    //             // printf("old score: %lf, new: %lf\n", score, posUpdateSeq[i].pssm);
    //             if(count == 0)
    //             {
    //                 score = posUpdateSeq[i].pssm;
    //                 updatedMotif = i;
    //             } else
    //             {
    //                 if(score < posUpdateSeq[i].pssm)
    //                 {
    //                     score = posUpdateSeq[i].pssm;
    //                     updatedMotif = i;
    //                 }
    //             }
    //             break;
    //         }
    //     }
    //     count++;
    // }
    return updatedMotif;
}

double samplerScore(rawSeq *inputSeq, int cSeq, int lSeq, int lMotif)
{
    char **pssmMotif, *motif;
    int i, j, k, x, y;
    double *seqBg, score = 0.0, tmp;
    double **freqMatrix, **probMatrix, **scoreMatrix;

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
    pssmMotif = (char **)malloc((cSeq - 1) * sizeof(char *));
    for (i = 0; i < (cSeq - 1); i++)
        pssmMotif[i] = (char *)malloc((lMotif + 1) * sizeof(char));

    for (i = 0; i < cSeq; i++)
    {
        // get the motif matrix to calculate pssm
        for (j = 0; j < cSeq; j++)
        {
            if (j < i)
            {
                strncpy(pssmMotif[j], getMotif(inputSeq[j], lMotif), lMotif);
                // pssmMotif[j] = getMotif(inputSeq[j], lMotif);
                pssmMotif[j][lMotif] = '\0';
                // printf("motif %d is %s\n", j, pssmMotif[j]);
            }
            else if (j == i)
                continue;
            else if (j > i)
            {
                strncpy(pssmMotif[j - 1], getMotif(inputSeq[j], lMotif), lMotif);
                // pssmMotif[j - 1] = getMotif(inputSeq[j], lMotif);
                pssmMotif[j - 1][lMotif] = '\0';
                // printf("motif %d is %s\n", j, pssmMotif[j - 1]);
            }
        }
        // for (j = 0; j < ( cSeq - 1); j++)
        // {
        //     printf("%s\n",pssmMotif[j]);
        // }
        /***** PSSM *****/
        seqBg = bg(inputSeq[i].seqInt, lSeq);

        freqMatrix = fMatrix(pssmMotif, lMotif, (cSeq - 1)); // calculate the freq matrix other than the update seq

        probMatrix = freq2prob(freqMatrix, lMotif);

        scoreMatrix = prob2score(probMatrix, lMotif, seqBg);
        // debug
        // printf("PSSM:\npos\tA\tC\tG\tT\t\n");
        // for (x = 0; x < lMotif; x++)
        // {
        //     printf("%d\t", x);
        //     for (y = 0; y < 4; y++)
        //     {
        //         printf("%.3lf\t", freqMatrix[x][y]);
        //     }
        //     printf("\n");
        // }

        // the score for this motif
        motif = (char*)malloc((lMotif + 1) * sizeof(char));
        strncpy(motif, getMotif(inputSeq[i], lMotif), lMotif);
            // motif = getMotif(inputSeq[i], lMotif);
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
    int i, j, upd, finallMotif, tmplMotif, lMotif, maxMotif[cSeq];
    char **pssmMotif, *motif;
    seedOut samplerOut;
    rawSeq *tmpSeq;
    double score, tmpscore, tmpscore_l, tmpscore_r;

    lMotif = lMotif_ori; // the originial len of motif from args

    pssmMotif = (char **)malloc((cSeq - 1) * sizeof(char *));
    for (i = 0; i < (cSeq - 1); i++)
        pssmMotif[i] = (char *)malloc((lMotif + 1) * sizeof(char));

    srand48(seed);
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
                    // pssmMotif[j] = getMotif(inputSeq[j], lMotif);
                    strncpy(pssmMotif[j], getMotif(inputSeq[j], lMotif), lMotif);
                    pssmMotif[j][lMotif] = '\0';
                    // printf("pssm motif %d is %s\n", j, pssmMotif[j]);
                }
                else if (j == i)
                    continue;
                // printf("This is the seq being updated.\n");
                else if (j > i)
                {
                    // pssmMotif[j - 1] = getMotif(inputSeq[j], lMotif);
                    strncpy(pssmMotif[j - 1], getMotif(inputSeq[j], lMotif), lMotif);
                    pssmMotif[j - 1][lMotif] = '\0';
                    // printf("pssm motif %d is %s\n", j, pssmMotif[j - 1]);
                }
            }
            // for (j = 0; j < (cSeq - 1); j++)
            // {
            //     printf("motif %d is %s\n", j, pssmMotif[j]);
            // }
            inputSeq[i].newposMotif = updateMotif(inputSeq[i], pssmMotif, lMotif, lSeq, cSeq);
            // printf("seq%d: new: %d old: %d\n", i, inputSeq[i].newposMotif, inputSeq[i].posMotif);
        }


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
        // motif = (char *)malloc((lMotif + 1) * sizeof(char));
        // for (i = 0; i < cSeq; i++)
        // {
        //     motif = getMotif(inputSeq[i], lMotif);
        //     motif[lMotif] = '\0';
        //     // strncpy(motif, getMotif(inputSeq[i], lMotif), lMotif);
        //     printf("%s\t%d\t%s\n", inputSeq[i].readName, inputSeq[i].newposMotif, motif);
        // }
        // printf("\n");
        /*** after one update, calculate the new score and 
         * check if it is larger than previous
         * if yes, save the new motif info,
         * if not, just continue the process ***/
        tmpscore = samplerScore(inputSeq, cSeq, lSeq, lMotif);
        // printf("tmpscore is: %lf score is: %lf \n", tmpscore, score);
        // printf("%lf\n", tmpscore);

        if (tmpscore > score)
        {
            score = tmpscore;
            for (i = 0; i < cSeq; i++)
            {
                maxMotif[i] = inputSeq[i].newposMotif;
            }
        }
        // for (i = 0; i < cSeq; i++)
        // {
        //     printf("%s  %s\n", maxSeq[i].readName, maxSeq[i].seq);
        // }
    }
    // printf("seed update OK\n");

    /******* At the end of each seed, calculate the score. *******/
    for (i = 0; i < cSeq; i++)
    {
        inputSeq[i].newposMotif = maxMotif[i];
    }
    score = samplerScore(inputSeq, cSeq, lSeq, lMotif);
    // store the out statistics for each seed iteration.
    samplerOut.score = score;
    samplerOut.seed = seed;
    samplerOut.lMotif = lMotif;
    // printf("score for seed %d is %lf\n", seed, score);
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