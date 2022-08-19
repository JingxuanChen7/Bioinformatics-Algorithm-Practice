#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
// #include <pthread.h>

// #define GENOME 100000000
double *bg(char *seq);
double **fMatrix(char **align, long lMotif, long cMotif);
double **freq2prob(double **freqMatrix, long lMotif);
double **prob2score(double **probMatrix, long lMotif, double *bgProb);
long getIndex(char n);
double motifScore(double **scoreMatrix, long lMotif, char *motifSeq);
double **complement(double **m, long nrow, long ncolumn);
void findMotif(char *seq, double **scoreMatrix, long lMotif, long len, double minScore, char strand);

int main(int argc, char *argv[])
{
    FILE *fp1, *fp2;
    char *alignFile, *seqFile, *seq, **align, c;
    long end, end1, lseq, len, lMotif, cMotif, i, j;
    double *bgProb, **freqMatrix, **probMatrix, **scoreMatrix, **cscoreMatrix, score, minScore;

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
    fscanf(fp2, ">%*[^\n]\n", NULL);
    end1 = ftell(fp2);
    fseek(fp2, 0, SEEK_END);
    end = ftell(fp2);
    lseq = end - end1;
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

    /* calculate background probs */
    bgProb = (double*)malloc(4 * sizeof(double));
    bgProb = bg(seq);
    // printf("background prob: %lf, %lf, %lf, %lf.\n", bgProb[0], bgProb[1], bgProb[2], bgProb[3]);

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

    /* prob matrix to score matrix */
    scoreMatrix = (double**)malloc(lMotif * sizeof(double*));
    for (i = 0; i < lMotif; i++)
        scoreMatrix[i] = (double*)malloc(4 * sizeof(double));
    scoreMatrix = prob2score(probMatrix, lMotif, bgProb);
    for (i = 0; i < lMotif; i++)
    {
        for (j = 0; j < 4; j++)
        {
            printf("%lf\t", scoreMatrix[i][j]);
        }
        printf("\n");
    }

    /* calculate score for each training motif */
    for (i = 0; i < cMotif; i++)
    {
        printf("%s\t", align[i]);
        score = motifScore(scoreMatrix, lMotif, align[i]);
        printf("%.3lf\n", score);
        if (i == 0)
            minScore = score;
        else if (minScore > score)
            minScore = score;
    }
    /* check whether user input the threshold */
    if (argc == 3)
        printf("No threshold was entered. Use the min score in training motif set: %.3lf\n", minScore);
    else if (argc == 4)
    {
        minScore = atoi(argv[3]);
        printf("Use the entered threshold: %.2lf", minScore);
    }

    /* find the complementary score matrix */
    cscoreMatrix = (double**)malloc(lMotif * sizeof(double*));
    for (i = 0; i < lMotif; i++)
        cscoreMatrix[i] = (double*)malloc(4 * sizeof(double));
    cscoreMatrix = complement(scoreMatrix, lMotif, 4);

    /* find possible motif in genome */
    // pthread_t thread_id;
    findMotif(seq, scoreMatrix, lMotif, len, minScore, '+');
    findMotif(seq, cscoreMatrix, lMotif, len, minScore, '-');

    return 0;
}

double* bg(char* seq)
{
    char* pseq;
    double *bgProb; // back ground prob for ACGT
    long lseq, l1, l2, l3, l4; // 1:A 2:C 3:G 4:T
    pseq = seq;
    lseq = l1 = l2 = l3 = l4 = 0;
    while (*pseq != '\0')
    {
        if (isalpha(*pseq) != 0)
        {
            lseq++;
            if (toupper(*pseq) == 'A')
                l1++;
            else if (toupper(*pseq) == 'C')
                l2++;
            else if (toupper(*pseq) == 'G')
                l3++;
            else if (toupper(*pseq) == 'T')
                l4++;
        }
        pseq++;
    }
    bgProb = (double *)malloc(4 * sizeof(double));
    bgProb[0] = (double)l1 / (double)lseq; // A
    bgProb[1] = (double)l2 / (double)lseq; // C
    bgProb[2] = (double)l3 / (double)lseq; // G
    bgProb[3] = (double)l4 / (double)lseq; // T
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
        freqMatrix[i][0] = (double)(countA + 0.25);
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

void findMotif(char *seq, double **scoreMatrix, long lMotif, long len, double minScore, char strand)
{
    char *motifSeq, *pseq;
    long i, maxpos;
    double score;

    maxpos = len - lMotif; // max first position
    pseq = seq; // point to first element
    motifSeq = (char*)malloc((lMotif + 1) * sizeof(char));
    motifSeq[lMotif] = '\0';
    for ( i = 0; i < maxpos; i++)
    {
        strncpy(motifSeq, pseq + i, (lMotif * sizeof(char)));
        score = motifScore(scoreMatrix, lMotif, motifSeq);
        if (score > minScore)
        {
            printf("%ld\t%ld\t%c\t%s\t%.3lf\n", i + 1, i + lMotif + 1, strand, motifSeq, score);
        }
    }
    
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