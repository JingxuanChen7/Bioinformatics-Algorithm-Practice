/*
    Assignment 4-PSSM
    Jingxuan Chen, Oct 18, 2019
    math included in this code,
    include -lm while compiling
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>


// #define GENOME 100000000
double *bg(int *seqInt, long len); // background prob for GC and AT
double **fMatrix(char **align, long lMotif, long cMotif); // frequency matrix
double **freq2prob(double **freqMatrix, long lMotif); // probability matrix
double **prob2score(double **probMatrix, long lMotif, double *bgProb); // score matrix
long getIndex(char n); // change ACGT to 0123
// char *int2seq(int *seqInt, long l); // change int to seq for printing
double seqScore(double **scoreMatrix, long lMotif, int *motifSeqInt); // calculate int seq score
double motifScore(double **scoreMatrix, long lMotif, char *motifSeq); // calculate alpha seq score
double **complement(double **m, long nrow, long ncolumn); // get the complementary score matrix
void findMotif(int *seqInt, char *seqChar, double **scoreMatrix, double **cscoreMatrix, long lMotif, long len, double minScore); // screening through genome

int main(int argc, char *argv[])
{
    clock_t start_t, file_t, bg_t, freq_t, prob_t, pssm_t, minscore_t, c_t, print_t;
    FILE *fp1, *fp2;
    char *alignFile, *seqFile, *seq, **align, c;
    long end, end1, lseq, len, lMotif, cMotif, i, j;
    int *seqInt;
    double *bgProb, **freqMatrix, **probMatrix, **scoreMatrix, **cscoreMatrix, score, minScore;

    start_t = clock();

    if ((argc != 3) && (argc != 4))
    {
        printf("\n\nUsage:Usage: ./pssm AlignmentFile SequenceFile [ScoreCutoff] \n\n");
        printf("If ScoreCutoff is omitted the higher of zero and the lowest\nscore among the sequences in the training set is used\n\n\n");
        printf("e.g., ./pssm sigma54.txt ecoK12-MG1655.fasta\n\n\n");

        return -1;
    }

    alignFile = argv[1]; seqFile = argv[2];
    fp1 = fopen(alignFile, "r");
    fp2 = fopen(seqFile, "r");
    if (fp1 == NULL || fp2 == NULL)
    {
        printf("\n\nUsage:Usage: ./pssm AlignmentFile SequenceFile [ScoreCutoff] \n\n");
        printf("If ScoreCutoff is omitted the higher of zero and the lowest\nscore among the sequences in the training set is used\n\n\n");
        perror("File error.\n");
        return -1;
    }

    /* find the size of fp1 and read fp1 */
    fscanf(fp1, "%*[^\n]\n", NULL);
    lMotif = ftell(fp1); // actual length + '\n'
    fseek(fp1, 0, SEEK_END);
    end = ftell(fp1);
    cMotif = (long) (end / lMotif);
    lMotif--; // get rid of '\n', actual align length
    align = (char**) malloc(cMotif * sizeof(char*));
    for (i = 0; i < cMotif; i++)
        align[i] = (char*) malloc((lMotif + 1) * sizeof(char));
    // printf("Length for each motif: %ld, Motif count: %ld\n", lMotif, cMotif);
    fseek(fp1, 0, SEEK_SET);
    for(i = 0; i < cMotif; i++)
    {
        for(j = 0; j < lMotif; j++)
        {
            align[i][j] = fgetc(fp1);
        }
        fgetc(fp1); // skip '\n'
        align[i][lMotif] = '\0';
    }
    // for (i = 0; i < cMotif; i++) // debug
    // {
    //     for (j = 0; j < lMotif; j++)
    //     {
    //         printf("%c", align[i][j]);
    //     }
    //     printf("\n");
    // }
    fclose(fp1);

    /* file the size of fp2 */
    fscanf(fp2, ">%*[^\n]\n", NULL);  // skip the first line
    end1 = ftell(fp2);
    fseek(fp2, 0, SEEK_END);
    end = ftell(fp2);
    lseq = end - end1;  // estimate length of sequence
    // printf("The length for genome seq is: %ld\n", lseq);
    fseek(fp2, 0, SEEK_SET);
    fscanf(fp2, ">%*[^\n]\n", NULL);
    seq = (char*) malloc(lseq * sizeof(char));

    len = 0;
    // fread(seq, sizeof(char), lseq, fp2);
    while (!feof(fp2))
    {
        if (isalpha(seq[len] = fgetc(fp2)) == 0) // if not letter, do not fget to seq
            continue;
        len++;
    }
    // printf("%s\n", seq);
    fclose(fp2);
    seq[len] = '\0';
    // printf("%s %ld", seq, len);

    /* change the char seq into int seq */
    seqInt = (int*)malloc(len * sizeof(int));
    for (i = 0; i < len; i++)
        seqInt[i] = getIndex(seq[i]); 
    // free(seq);

    file_t = clock();
    printf("\nTime until read file: %lf\n", (double)(file_t - start_t) / CLOCKS_PER_SEC);

    /* calculate background probs */
    bgProb = (double *)malloc(4 * sizeof(double));
    bgProb = bg(seqInt, len);
    // printf("background prob: %lf, %lf, %lf, %lf.\n", bgProb[0], bgProb[1], bgProb[2], bgProb[3]);
    bg_t = clock();
    printf("\nTime until calculating background: %lf\n", (double)(bg_t - start_t) / CLOCKS_PER_SEC);

    /* Freq matrix with pseudocounts */
    freqMatrix = (double**)malloc(lMotif * sizeof(double*));
    for (i = 0; i < lMotif; i++)
        freqMatrix[i] = (double*)malloc(4 * sizeof(double));
    freqMatrix = fMatrix(align, lMotif, cMotif);
    // for (i = 0; i < lMotif; i++)
    // {
    //     for (j = 0; j < 4; j++)
    //     {
    //         printf("%lf\t", freqMatrix[i][j]);
    //     }
    //     printf("\n");
    // }
    freq_t = clock();
    printf("\nTime until calculate frequency: %lf\n", (double)(freq_t - start_t) / CLOCKS_PER_SEC);

    /* Freq matrix to prob matrix */
    probMatrix = (double**)malloc(lMotif * sizeof(double*));
    for (i = 0; i < lMotif; i++)
        probMatrix[i] = (double*)malloc(4 * sizeof(double));
    probMatrix = freq2prob(freqMatrix, lMotif);
    // for (i = 0; i < lMotif; i++)
    // {
    //     for (j = 0; j < 4; j++)
    //     {
    //         printf("%lf\t", probMatrix[i][j]);
    //     }
    //     printf("\n");
    // }
    prob_t = clock();
    printf("\nTime until calculate probability: %lf\n", (double)(prob_t - start_t) / CLOCKS_PER_SEC);

    /* prob matrix to score matrix */
    scoreMatrix = (double**)malloc(lMotif * sizeof(double*));
    for (i = 0; i < lMotif; i++)
        scoreMatrix[i] = (double*)malloc(4 * sizeof(double));
    scoreMatrix = prob2score(probMatrix, lMotif, bgProb);
    pssm_t = clock();
    printf("\nTime until calculate pssm: %lf\n", (double)(pssm_t - start_t) / CLOCKS_PER_SEC);

    printf("PSSM:\npos\tA\tC\tG\tT\t\n");
    for (i = 0; i < lMotif; i++)
    {
        printf("%ld\t", i);
        for (j = 0; j < 4; j++)
        {
            printf("%.3lf\t", scoreMatrix[i][j]);
        }
        printf("\n");
    }

    /* calculate score for each training motif */
    printf("\nThe scores for each training motif:\n");
    for (i = 0; i < cMotif; i++)
    {
        printf("%s\t", align[i]);
        score = motifScore(scoreMatrix, lMotif, align[i]);
        printf("%.3lf\n", score);
        if (i == 0)  // find the min score
            minScore = score;
        else if (minScore > score)
            minScore = score;
    }

    /* check whether user input the threshold */
    if (argc == 3)
        printf("\nNo threshold was entered. Use the min score in training motif set: %.3lf\n", minScore);
    else if (argc == 4)
    {
        minScore = atoi(argv[3]);
        printf("\nUse the entered threshold: %.2lf", minScore);
    }
    minscore_t = clock();
    printf("\nTime until got threshold: %lf\n", (double)(minscore_t - start_t) / CLOCKS_PER_SEC);

    /* find the complementary score matrix */
    cscoreMatrix = (double**)malloc(lMotif * sizeof(double*));
    for (i = 0; i < lMotif; i++)
        cscoreMatrix[i] = (double*)malloc(4 * sizeof(double));
    cscoreMatrix = complement(scoreMatrix, lMotif, 4);
    c_t = clock();
    printf("\nTime until filp the matrix: %lf\n", (double)(c_t - start_t) / CLOCKS_PER_SEC);

    /* find possible motif in genome */

    findMotif(seqInt, seq, scoreMatrix, cscoreMatrix, lMotif, len, minScore);

    print_t = clock();
    printf("\nTime until print results: %lf\n", (double)(print_t - start_t) / CLOCKS_PER_SEC);

    free(freqMatrix);
    free(scoreMatrix);
    free(probMatrix);
    return 0;
}

double* bg(int* seqInt, long len)
{
    double *bgProb; // back ground prob for ACGT
    long lseq, lcg, lat, i; // 1:A 2:C 3:G 4:T
    lseq = lcg = lat = 0;
    for (i = 0; i < len; i++)
    {
        if(seqInt[i] >= 0)
        {
            lseq++;
            if ((seqInt[i] == 0) || (seqInt[i] == 3))  // A or T
                lat++;
            else if ((seqInt[i] == 1) || (seqInt[i] == 2))  // C or G
                lcg++;
        }
    }
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
    bgProb[0] = ((double)lat / (double)lseq) / 2; // A
    bgProb[1] = ((double)lcg / (double)lseq) / 2; // C
    bgProb[2] = ((double)lcg / (double)lseq) / 2; // G
    bgProb[3] = ((double)lat / (double)lseq) / 2; // T
    return bgProb;
}

double **fMatrix(char **align, long lMotif, long cMotif)
{
    double **freqMatrix;
    long countA, countC, countG, countT, i, j;

    freqMatrix = (double**)malloc(lMotif * sizeof(double*));
    for (i = 0; i < lMotif; i++)
        freqMatrix[i] = (double*)malloc(4 * sizeof(double));
    
    for (i = 0; i < lMotif; i++)
    {
        countA = countC = countG = countT = 0;
        for (j = 0; j < cMotif; j++)
        {
            if(align[j][i] == 'A')
                countA++;
            else if(align[j][i] == 'C')
                countC++;
            else if(align[j][i] == 'G')
                countG++;
            else if(align[j][i] == 'T')
                countT++;
        }
        freqMatrix[i][0] = (double)(countA + 0.25);  // add qseudocount to avoid zeros
        freqMatrix[i][1] = (double)(countC + 0.25);
        freqMatrix[i][2] = (double)(countG + 0.25);
        freqMatrix[i][3] = (double)(countT + 0.25);
    }
    return freqMatrix;
}

double **freq2prob(double **freqMatrix, long lMotif)
{
    double **probMatrix, sumPos;
    long i, j;
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

double **prob2score(double **probMatrix, long lMotif, double *bgProb)
{
    double **scoreMatrix;
    long i, j;

    scoreMatrix = (double**)malloc(lMotif * sizeof(double*));
    for (i = 0; i < lMotif; i++)
        scoreMatrix[i] = (double*)malloc(4 * sizeof(double));

    for (i = 0; i < lMotif; i++)
    {
        for(j = 0; j < 4; j++)
            scoreMatrix[i][j] = log2(probMatrix[i][j] / bgProb[j]);
    }

    return scoreMatrix;
}

double motifScore(double **scoreMatrix, long lMotif, char *motifSeq)
{
    double score = 0.0;
    long index, i;
    for (i = 0; i < lMotif; i++)
    {
        index = getIndex(motifSeq[i]);
        score += scoreMatrix[i][index];
    }
    return score;
}

double seqScore(double **scoreMatrix, long lMotif, int *motifSeqInt)
{
    double score = 0.0;
    long index, i;
    for (i = 0; i < lMotif; i++)
    {
        index = motifSeqInt[i]; 
        score += scoreMatrix[i][index];
    }
    return score;
}

long getIndex(char n)
{
    if (n == 'A')
        return 0;
    else if (n == 'C')
        return 1;
    else if (n == 'G')
        return 2;
    else if (n == 'T')
        return 3;
    else 
        return -1;
}

// char *int2seq(int *seqInt, long l)
// {
//     long i;
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

void findMotif(int *seqInt, char *seqChar, double **scoreMatrix, double **cscoreMatrix, long lMotif, long len, double minScore)
// void *findMotif(void* arg)
{
    // args *aarg = arg;
    int *motifSeq, *pseq;
    char *seq, *pseqChar;
    long i, maxpos;
    double score, cscore;

    maxpos = len - lMotif; // max first position
    pseq = seqInt;                           // point to first element
    pseqChar = seqChar;
    motifSeq = (int *)malloc((lMotif) * sizeof(int));
    seq = (char*)malloc((lMotif + 1) * sizeof(char));
    // motifSeq[lMotif] = '\0';
    for ( i = 0; i < maxpos; i++)
    {
        // strncpy(motifSeq, pseq + i, (lMotif * sizeof(char)));
        motifSeq = pseq + i;
        score = seqScore(scoreMatrix, lMotif, motifSeq);
        if (score > minScore)
        {
            // seq = int2seq(motifSeq, lMotif);
            seq = pseqChar + i;
            seq[lMotif] = '\0';
            printf("%ld\t%ld\t%c\t%s\t%.3lf\n", i + 1, i + lMotif + 1, '+', seq, score);
        }
        cscore = seqScore(cscoreMatrix, lMotif, motifSeq);
        if (cscore > minScore)
        {
            // seq = int2seq(motifSeq, lMotif);
            seq = pseqChar + i;
            seq[lMotif] = '\0';
            printf("%ld\t%ld\t%c\t%s\t%.3lf\n", i + 1, i + lMotif + 1, '-', seq, cscore);
        }
    }
    // return NULL;
}

double **complement(double **m, long nrow, long ncolumn)
{
    double **cm;
    long i, j;
    cm = (double**)malloc(nrow * sizeof(double*));
    for (i = 0; i < nrow; i++)
        cm[i] = (double*)malloc(ncolumn * sizeof(double));
    for (i = 0; i < nrow; i++)
    {
        for (j = 0; j < ncolumn; j++)
            cm[i][j] = m[nrow-i-1][ncolumn-j-1];
    }
    return cm;
}