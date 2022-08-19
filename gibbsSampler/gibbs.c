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
#define NSEED 100
#define NUPDATE 500
#define ADJUST 5
#define STOPUPDATE 1

typedef struct
{
    char readName[READNAMELEN];
    char seq[SEQLEN];
    int seqInt[SEQLEN];
    int posMotif;
    int newposMotif;
}rawSeq; // raw sequences directly from file

typedef struct 
{
    double pssm;
    double prob;
}posSeq; // store the probability for each motif

typedef struct 
{
    double score;
    int seed;
    int lMotif;
}seedOut; // the output for each seed

rawSeq *inputSeq;

char *getMotif(rawSeq inputSeq, int lMotif);
int checkUpdate(rawSeq *inputSeq, int cSeq);
int updateMotif(rawSeq updateSeq, char **otherSeqMotif, int lMotif, int lSeq, int cSeq);
double samplerScore(rawSeq *inputSeq, int cSeq, int lSeq, int lMotif);
seedOut sampler(int seed, int cSeq, int lSeq, int lMotif, rawSeq *inputSeq);
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

    if (argc != 3)
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
    /* calculate the number of sequences in the input file */
    cSeq = 0;
    while (!feof(fp))
    {
        if ((c = fgetc(fp)) == '>')
            cSeq++;
    }
    // printf("The count for squences: %d\n", cSeq);
    inputSeq = (rawSeq*) malloc( cSeq * sizeof(rawSeq));
    fseek(fp, 0, SEEK_SET);
    fgetc(fp); // geet rid of the first '>'
    for (i = 0; i < cSeq; i++)
    {
        fscanf(fp, "%s", inputSeq[i].readName);
        // fgetc(fp);
        j = 0;
        while (((c = fgetc(fp)) != '>') && !feof(fp))
        {
            if (isalpha(c) == 0) continue;
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
    //     printf("%s  %s\n", inputSeq[i].readName, inputSeq[i].seq);
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
    samplerOut = (seedOut*)malloc(NSEED * sizeof(seedOut));

    // /* run nseeds times gibs algorithem to find the local max*/
    // pssmMotif = (char**)malloc((cSeq - 1) * sizeof(char*));
    // for (i = 0; i < (cSeq - 1); i++)
    //     pssmMotif[i] = (char*)malloc((lMotif + 1) * sizeof(char));

    /******** Test different initializations to find the global max ********/
    for (seed = 0; seed < NSEED; seed++)
    {
        samplerOut[seed] = sampler(seed, cSeq, lSeq, lMotif, inputSeq);
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
    finalOut = sampler(maxSeed, cSeq, lSeq, lMotif, inputSeq);
    printf("score for seed %d is %lf\n", maxSeed, finalOut.score);
    printf("the final len of motif is: %d\n", finalOut.lMotif);
    motif = (char *)malloc((finalOut.lMotif + 1) * sizeof(char));
    for (i = 0; i < cSeq; i++)
    {
        motif = getMotif(inputSeq[i], finalOut.lMotif);
        motif[finalOut.lMotif] = '\0';
        // strncpy(motif, getMotif(inputSeq[i], lMotif), lMotif);
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
    seqMotif[lMotif] = '\0'; // set the last char of seq
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
    double score, p, sumP = 0.0, len = 0.0, randNum;
    int i, j, maxpos = lSeq - lMotif, *thisSeqInt, *pseqInt, count;
    posSeq *posUpdateSeq;
    int updatedMotif;

    /***** PSSM *****/
    seqBg = (double*)malloc(4 * sizeof(double));
    seqBg = bg(updateSeq.seqInt, lSeq);

    freqMatrix = (double**)malloc(lMotif * sizeof(double*));
    for(i = 0; i < lMotif; i++)
        freqMatrix[i] = (double*)malloc(4 * sizeof(double));
    freqMatrix = fMatrix(otherSeqMotif, lMotif, (cSeq - 1)); // calculate the freq matrix other than the update seq

    probMatrix = (double**)malloc(lMotif * sizeof(double*));
    for(i = 0; i < lMotif; i++)
        probMatrix[i] = (double*)malloc(4 * sizeof(double));
    probMatrix = freq2prob(freqMatrix, lMotif);

    scoreMatrix = (double**)malloc(lMotif * sizeof(double*));
    for(i = 0; i < lMotif; i++)
        scoreMatrix[i] = (double*)malloc(4 * sizeof(double));
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
    posUpdateSeq = (posSeq*)malloc((maxpos + 1) * sizeof(posSeq));
    thisSeqInt = (int *)malloc(lMotif * sizeof(int));
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
    count = 0;
    while (count <= STOPUPDATE)
    {
        randNum = drand48();
        for (i = 0; i <= maxpos; i++)
        {
            len += (posUpdateSeq[i].prob / sumP);
            // printf("len = %lf; i = %d \n", len, i);
            if (len > randNum)
            {
                updatedMotif = i;
                // printf("old score: %lf, new: %lf\n", score, posUpdateSeq[i].pssm);
                if (posUpdateSeq[i].pssm > score) 
                    return updatedMotif;
                break;
            }
        }
        if (count == STOPUPDATE)
        {
            return updateSeq.newposMotif;
        }
        count++;
    }
    return -1;
}

double samplerScore(rawSeq *inputSeq, int cSeq, int lSeq, int lMotif)
{
    double *bgProb, **freqMatrix, **probMatrix, **scoreMatrix, score = 0.0, tmp;
    int *allSeqInt, i, j, pos;
    char **pssmMotif, *allSeq;
    // printf("enter func OK\n");

    // all seq background
    // allSeq = (char*)malloc((cSeq * lMotif) * sizeof(char));
    // allSeqInt = (int*)malloc((cSeq * lMotif) * sizeof(int));
    allSeqInt = (int *)malloc((cSeq * lSeq) * sizeof(int));
    for (i = 0; i < cSeq; i++)
    {
        // for (j = 0; j < lMotif; j++)
        // {
        //     pos = inputSeq[i].newposMotif + j;
        //     allSeq[(i * lMotif) + j] = inputSeq[i].seq[pos];
        // }
        for (j = 0; j < lSeq; j++)
        {
            allSeqInt[(i * lSeq) + j] = inputSeq[i].seqInt[j];
        }
    }
    // for (i = 0; i < (cSeq * lMotif); i++)
    // {
    //     allSeqInt[i] = getIndex(allSeq[i]);
    // }
    bgProb = bg(allSeqInt, (cSeq * lSeq));
    // printf("bg OK\n");
    // Motif alignment
    pssmMotif = (char **)malloc(cSeq * sizeof(char *));
    for (i = 0; i < cSeq; i++)
        pssmMotif[i] = (char *)malloc((lMotif + 1) * sizeof(char));
    for (i = 0; i < cSeq; i++)
    {
        strncpy(pssmMotif[i], getMotif(inputSeq[i], lMotif), lMotif);
        pssmMotif[i][lMotif] = '\0';
        // printf("motif %d = %s\n", i, pssmMotif[i]);
    }
    // printf("alignment OK\n");

    // generate pssm score matrix
    freqMatrix = (double **)malloc(lMotif * sizeof(double *));
    for (i = 0; i < lMotif; i++)
        freqMatrix[i] = (double *)malloc(4 * sizeof(double));
    freqMatrix = fMatrix(pssmMotif, lMotif, cSeq); // calculate the freq matrix other than the update seq

    probMatrix = (double **)malloc(lMotif * sizeof(double *));
    for (i = 0; i < lMotif; i++)
        probMatrix[i] = (double *)malloc(4 * sizeof(double));
    probMatrix = freq2prob(freqMatrix, lMotif);

    scoreMatrix = (double **)malloc(lMotif * sizeof(double *));
    for (i = 0; i < lMotif; i++)
        scoreMatrix[i] = (double *)malloc(4 * sizeof(double));
    scoreMatrix = prob2score(probMatrix, lMotif, bgProb);
    // printf("matrix OK\n");

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

    // calculate the motif score for each motif and then sum
    for (i = 0; i < cSeq; i++)
    {
        tmp = (double)motifScore(scoreMatrix, lMotif, pssmMotif[i]);
        score += tmp;
        // printf("score for %d is %lf\n", i, tmp);
    }
    // printf("score OK\n");

    return score;
}

    /***** sampling function for each seed *****/
seedOut sampler(int seed, int cSeq, int lSeq, int lMotif_ori, rawSeq *inputSeq)
{
    int i, j, upd, finallMotif, tmplMotif, lMotif;
    char **pssmMotif, *motif;
    seedOut samplerOut;
    rawSeq* tmpSeq;
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
                if (j < i)
                {
                    pssmMotif[j] = getMotif(inputSeq[j], lMotif);
                    // strncpy(pssmMotif[j], getMotif(inputSeq[j], lMotif), lMotif);
                    pssmMotif[j][lMotif] = '\0';
                    // printf("motif %d is %s\n", j, pssmMotif[j]);
                }
                else if (j == i)
                    continue;
                // printf("This is the seq being updated.\n");
                else if (j > i)
                {
                    pssmMotif[j - 1] = getMotif(inputSeq[j], lMotif);
                    // strncpy(pssmMotif[j - 1], getMotif(inputSeq[j], lMotif), lMotif);
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
        /*** after one update, calculate the new score and 
         * check if it is larger than previous
         * if yes, save the new motif info,
         * if not, back to old motif info ***/
        tmpscore = samplerScore(inputSeq, cSeq, lSeq, lMotif);
        // printf("tmpscore is: %lf\n", tmpscore);
        /**** after each update, check if stop changing ****/
        if (tmpscore < score)
        {
            for (i = 0; i < cSeq; i++)
                inputSeq[i].newposMotif = inputSeq[i].posMotif;
        } else
        {
            score = tmpscore;
        }
        /****** adjust the border of motif ********/
        if (0)
        // if (upd % ADJUST == 0)
        {
            // printf("lMotif now is %d, score is %lf\n", lMotif, score);
            tmpSeq = inputSeq;
            /***** firstly move the left border *****/
            // check left border
            if (checkLeft(cSeq, lSeq, lMotif, tmpSeq) == -1)
                tmpscore_l = -1; // could not move to left
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
            } else if ((tmpscore_r > tmpscore_l) && (tmpscore_r > score))
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
                tmpscore_r = -1;
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
        // if (checkUpdate(inputSeq, cSeq) == 1)
        //     break;
    }
    // printf("seed update OK\n");

    /******* At the end of each seed, calculate the score. *******/
    score = samplerScore(inputSeq, cSeq, lSeq, lMotif);
    // printf("score for seed %d is %lf\n", seed, score);
    // store the out statistics for each seed iteration.
    samplerOut.score = score;
    samplerOut.seed = seed;
    samplerOut.lMotif = lMotif;
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
int checkLeft(int cSeq, int lSeq, int lMotif, rawSeq* inputSeq)
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