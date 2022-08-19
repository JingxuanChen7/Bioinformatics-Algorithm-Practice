#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <pthread.h>
#include <time.h>
#include "pssm.h"

double *bg(int *seqInt, int len)
{
    double *bgProb, lc, lg, la, lt; // back ground prob for ACGT
    int lseq, i; // 1:A 2:C 3:G 4:T
    lseq = lc = lg = la = lt = 0;
    for (i = 0; i < len; i++)
    {
        if (seqInt[i] >= 0)
        {
            lseq++;
            if (seqInt[i] == 0)
                la++;
            else if (seqInt[i] == 1)
                lc++;
            else if (seqInt[i] == 2)
                lg++;
            else if (seqInt[i] == 3)
                lt++;
            // if ((seqInt[i] == 0) || (seqInt[i] == 3)) // A or T
            //     lat++;
            // else if ((seqInt[i] == 1) || (seqInt[i] == 2)) // C or G
            //     lcg++;
        }
    }
    // pseudocouts for freq
    la += 0.25;
    lc += 0.25;
    lg += 0.25;
    lt += 0.25;
    // while (seqInt != '\0')
    // {
    //     if (isalpha(*pseq) != 0)
    //     {
    //         lseq++;
    //         if (toupper(*pseq) == 'A')
    //             l1++;
    //         else if (toupper(*pseq) == 'C')
    //             l2++;
    //         else if (toupper(*pseq) == 'G')
    //             l3++;
    //         else if (toupper(*pseq) == 'T')
    //             l4++;
    //     }
    //     pseq++;
    // }
    bgProb = (double *)malloc(4 * sizeof(double));
    bgProb[0] = la / (double)lseq; // A
    bgProb[1] = lc / (double)lseq; // C
    bgProb[2] = lg / (double)lseq; // G
    bgProb[3] = lt / (double)lseq; // T
    // bgProb[0] = 0.25; // A
    // bgProb[1] = 0.25; // C
    // bgProb[2] = 0.25; // G
    // bgProb[3] = 0.25; // T
    // printf("%lf %lf %lf %lf\n", bgProb[0], bgProb[1], bgProb[2], bgProb[3]);
    return bgProb;
}

double **fMatrix(char **align, int lMotif, int cMotif)
{
    double **freqMatrix;
    int countA, countC, countG, countT, i, j;

    freqMatrix = (double **)malloc(lMotif * sizeof(double *));
    for (i = 0; i < lMotif; i++)
        freqMatrix[i] = (double *)malloc(4 * sizeof(double));

    for (i = 0; i < lMotif; i++)
    {
        countA = countC = countG = countT = 0;
        for (j = 0; j < cMotif; j++)
        {
            if (align[j][i] == 'A')
                countA++;
            else if (align[j][i] == 'C')
                countC++;
            else if (align[j][i] == 'G')
                countG++;
            else if (align[j][i] == 'T')
                countT++;
        }
        freqMatrix[i][0] = (double)(countA + 0.25); // add qseudocount to avoid zeros
        freqMatrix[i][1] = (double)(countC + 0.25);
        freqMatrix[i][2] = (double)(countG + 0.25);
        freqMatrix[i][3] = (double)(countT + 0.25);
    }
    return freqMatrix;
}

double **freq2prob(double **freqMatrix, int lMotif)
{
    double **probMatrix, sumPos;
    int i, j;
    probMatrix = (double **)malloc(lMotif * sizeof(double *));
    for (i = 0; i < lMotif; i++)
        probMatrix[i] = (double *)malloc(4 * sizeof(double));

    for (i = 0; i < lMotif; i++)
    {
        sumPos = 0.0;
        for (j = 0; j < 4; j++)
            sumPos += freqMatrix[i][j];
        for (j = 0; j < 4; j++)
            probMatrix[i][j] = freqMatrix[i][j] / sumPos;
    }

    return probMatrix;
}

double **prob2score(double **probMatrix, int lMotif, double *bgProb)
{
    double **scoreMatrix;
    int i, j;

    scoreMatrix = (double **)malloc(lMotif * sizeof(double *));
    for (i = 0; i < lMotif; i++)
        scoreMatrix[i] = (double *)malloc(4 * sizeof(double));

    for (i = 0; i < lMotif; i++)
    {
        for (j = 0; j < 4; j++)
            scoreMatrix[i][j] = log2(probMatrix[i][j] / bgProb[j]);
    }

    return scoreMatrix;
}

double motifScore(double **scoreMatrix, int lMotif, char *motifSeq)
{
    double score = 0.0;
    int index, i;
    for (i = 0; i < lMotif; i++)
    {
        index = getIndex(motifSeq[i]);
        score += scoreMatrix[i][index];
    }
    return score;
}

double seqScore(double **scoreMatrix, int lMotif, int *motifSeqInt)
{
    double score = 0.0;
    int index, i;
    for (i = 0; i < lMotif; i++)
    {
        index = motifSeqInt[i];
        score += scoreMatrix[i][index];
    }
    return score;
}

int getIndex(char n)
{
    if (n == 'A' || n == 'a')
        return 0;
    else if (n == 'C' || n == 'c')
        return 1;
    else if (n == 'G' || n == 'g')
        return 2;
    else if (n == 'T' || n == 't')
        return 3;
    else
        return -1;
}

// char *int2seq(int *seqInt, int l)
// {
//     int i;
//     char *seq = (char*)malloc((l + 1) * sizeof(char));
//     seq[l] = '\0';
//     for (i = 0; i < l; i++)
//     {
//         if (seqInt[i] == 0)
//             seq[i] = 'A';
//         else if (seqInt[i] == 1)
//             seq[i] = 'C';
//         else if (seqInt[i] == 2)
//             seq[i] = 'G';
//         else if (seqInt[i] == 3)
//             seq[i] = 'T';
//     }
//     return seq;
// }
