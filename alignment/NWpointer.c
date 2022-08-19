#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <time.h>

#define MAXSEQ 1000000 // the max length of sequences

int isLetter(char c);
long findMax(long a, long b, long c);
int **matrixMatch(char *seq1, char *seq2, long lseq1, long lseq2, int sMatch, int sMismatch);
long **matrixScore(int **mMatch, long lseq1, long lseq2, int sGap);
void printAlign(long **mScore, int **mMatch, char *seq1, char *seq2, long lseq1, long lseq2, int sMatch, int sMismatch, int sGap);

int main(int argc, char *argv[])
{
    clock_t start_t, file_t, match_t, score_t, print_t;
    char *seqFile1, *seqFile2, *seq1, *seq2;
    int sMatch, sMismatch, sGap, **mMatch;
    long lseq1 = 0, lseq2 = 0, nline = 0, i, j, **mScore;

    start_t = clock();
    if (argc != 6)
    {
        printf("\nUsage: ./NWalignment <seqFile1> <seqFile2> <MatchScore> <MisMatch> <Gap>\n\n");
        printf("e.g. ./NWalignment Test01.fasta Test02.fasta 1 0 -1\n\n");
        return -1;
    }

    seqFile1 = argv[1]; seqFile2 = argv[2];
    FILE *fp1, *fp2;
    fp1 = fopen(seqFile1, "r");
    fp2 = fopen(seqFile2, "r");
    if(fp1 == NULL || fp2 == NULL)
    {
        printf("\nUsage: ./NWalignment <seqFile1> <seqFile2> <MatchScore> <MisMatch> <Gap>\n\n");
        printf("e.g. ./NWalignment seq1.fasta seq2.fasta 1 0 -1\n\n");
        perror("File error.\n");
        return -1;
    }
    sMatch = atoi(argv[3]);
    sMismatch = atoi(argv[4]);
    sGap = atoi(argv[5]);

    seq1 = malloc(MAXSEQ * sizeof(char));
    seq2 = malloc(MAXSEQ * sizeof(char));

    /* skip the sequence name line*/
    fscanf(fp1, ">%*[^\n]\n", NULL);
    while(!feof(fp1))
    {
        if (isLetter(*(seq1 + lseq1) = fgetc(fp1)) == -1) // if not letter, do not fget to seq
            continue;
        lseq1++;
        
    }
    * (seq1 + lseq1) = '\0';
    seq1 = (char*)realloc(seq1, (lseq1 + 1) * sizeof(char));
    for(i = 0; i < lseq1; i++)
    {
        *(seq1 + i) = toupper(*(seq1 + i)); // standardize to uppercase
        //printf("%c", seq1[i]);
    }
    //printf("\nLength is %ld\n", lseq1);
    fclose(fp1);
    /* skip the sequence name line*/
    fscanf(fp2, ">%*[^\n]\n", NULL);
    while (!feof(fp2))
    {
        if (isLetter(*(seq2 + lseq2) = fgetc(fp2)) == -1)
            continue;
        lseq2++;
    }
    *(seq2 + lseq2) = '\0';
    seq2 = realloc(seq2, (lseq2 + 1) * sizeof(char));
    for (i = 0; i < lseq2; i++)
    {
        *(seq2 + i) = toupper(*(seq2 + i)); // standardize to uppercase
        //printf("%c", seq2[i]);
    }
    //printf("\nLength is %ld\n", lseq2);
    fclose(fp2);
    /* print out the input sequences */
    printf("\nNeedleman-Wunsch global alignment\n\n");
    printf("\nSequence 1 (%s):\n\n", seqFile1);
    for (i = 0; i < lseq1; i++)
        printf("%c", *(seq1 + i));
    printf("\n\nLength is %ld\n\n", lseq1);
    printf("\nSequence 2 (%s):\n\n", seqFile2);
    for (i = 0; i < lseq2; i++)
        printf("%c", *(seq2 + i));
    printf("\n\nLength is %ld\n\n", lseq2);

    file_t = clock();

    /* initialize matrix with match & mismatch score */
    mMatch = (int **)malloc((lseq2) * sizeof(int *));
    // for (i = 0; i <= lseq2; i++)
    //     mMatch[i] = (int *) malloc((lseq1) * sizeof(int));
    mMatch[0] = (int *)malloc((lseq1 * lseq2) * sizeof(int));
    mMatch = matrixMatch(seq1, seq2, lseq1, lseq2, sMatch, sMismatch);
    printf("match finished\n");
    // printf("\t");
    // for (j = 0; j < lseq1; j++)
    // {
    //     printf("%c\t", seq1[j]);
    // }
    // printf("\n");

    // for (i = 0; i < lseq2; i++) // debug initialization
    // {
    //     printf("%c\t", seq2[i]);
    //     for (j = 0; j < lseq1; j++)
    //     {
    //         printf("%d\t", mMatch[i][j]);
    //     }
    //     printf("\n");
    // }
    // printf("\n");
    match_t = clock();

    /* fill the matrix score */
    mScore = (long **)malloc((lseq2 + 1) * sizeof(long *));
    // for (i = 0; i <= lseq2; i++)
    // {
    //     mScore[i] = (long *)malloc((lseq1 + 1) * sizeof(long));
    // }
    mScore[0] = (long *)malloc((lseq1 + 1) * (lseq2 + 1) * sizeof(long));
    mScore = matrixScore(mMatch, lseq1, lseq2, sGap);
    printf("score finish\n");
    // for (i = 0; i < lseq2 + 1; i++)  // debug calculating score
    // {
    //     for (j = 0; j < lseq1 + 1; j++)
    //     {
    //         printf("%ld\t", mScore[i][j]);
    //     }
    //     printf("\n");
    // }

    score_t = clock();

    /* print out alignment */
    printAlign(mScore, mMatch, seq1, seq2, lseq1, lseq2, sMatch, sMismatch, sGap);
    print_t = clock();


    printf("\nTime until read file: %lf\n", (double)(file_t - start_t) / CLOCKS_PER_SEC);
    printf("Time until make match matrix: %lf\n", (double)(match_t - start_t) / CLOCKS_PER_SEC);
    printf("Time until make score matrix: %lf\n", (double)(score_t - start_t) / CLOCKS_PER_SEC);
    printf("Time until print output: %lf\n", (double)(print_t - start_t) / CLOCKS_PER_SEC);
    // free(mMatch); // report bug on teaching cluster
    // free(mScore);
    return 0;
}

int isLetter(char c)
{
    if ((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z'))
        return 1;
    else
        return -1;
}

/* initialize a matrix with match/mismatch score */
int **matrixMatch(char *seq1, char *seq2, long lseq1, long lseq2, int sMatch, int sMismatch)
{
    int **mMatch;
    long i, j;
    mMatch = (int**) malloc((lseq2) * sizeof(int*));
    // for (i = 0; i <= lseq2; i++)
    //     *mMatch[i] = (int*) malloc((lseq1) * sizeof(int));
    mMatch[0] = (int *)malloc((lseq1 * lseq2) * sizeof(int));
    printf("malloc\n");
    for (i = 0; i < lseq2; i++)
    {
        for (j = 0; j < lseq1; j++)
        {
            if(*(seq1 + j) == *(seq2 + i))          
                **(mMatch + i * lseq1 + j) = sMatch;
            else
                **(mMatch + i * lseq1 + j) = sMismatch;
            printf("match %ld %ld\n", i, j);
        }
    }

    return mMatch;
}

/* fill in a martix for score */
long **matrixScore (int **mMatch, long lseq1, long lseq2, int sGap)
{
    long **mScore, right, down, diag;
    long i, j;
    mScore = (long**) malloc((lseq2 + 1) * sizeof(long*));
    // for (i = 0; i < lseq2 + 1; i++)
    // {
    //     mScore[i] = (long*) malloc((lseq1 + 1) * sizeof(long));
    // }
    mScore[0] = (long *)malloc((lseq1 + 1) * (lseq2 + 1) * sizeof(long));

    /* start to fill in the score matrix */
    mScore[lseq2][lseq1] = 0;
    for (i = 1; i < lseq2 + 1; i++) // fill the last row and column
    {
        **(mScore + (lseq2-i) * lseq1 + lseq1) = i * sGap;
    }
    for (j = 1; j <= lseq1 + 1; j++)
    {
        **(mScore + lseq2 * lseq1 + lseq1 - j) = j * sGap;
    }
    for (i = 1; i < lseq2 + 1; i++)
    {
        for (j = 1; j < lseq1 + 1; j++)
        {
            right = **(mScore + (lseq2 - i) * lseq1 + lseq1 - j + 1)+ sGap;
            down = **(mScore + (lseq2 - i + 1) * lseq1 + lseq1 - j) + sGap;
            diag = **(mScore + (lseq2 - i + 1) * lseq1 + lseq1 - j + 1) + **(mMatch + (lseq2 - i) * lseq1 + lseq1 - j);
            **(mScore + (lseq2 - i) * lseq1 + lseq1 - j) = findMax(right, down, diag);
            // printf("[%ld,%ld]: %ld %ld %ld max:%ld\n", lseq2 - i, lseq1 - j, right, down, diag, mScore[lseq2 - i][lseq1 - j]);
            // printf("right:%d down:%d diag:%d sGap:%d match:%d\n", mScore[lseq2 - i][lseq1 - j + 1], mScore[lseq2 - i + 1][lseq1 - j], mScore[lseq2 - i + 1][lseq1 - i + 1], sGap, mMatch[lseq2 - i][lseq1 - j]);
        }
    }

    return mScore;
}

long findMax(long a, long b, long c)
{
    long max = a;
    if (max < b)
        max = b;
    if (max < c) 
        max = c;
    return max;
}

void printAlign(long **mScore, int **mMatch, char *seq1, char *seq2, long lseq1, long lseq2, int sMatch, int sMismatch, int sGap)
{
    long i = 0, j = 0, score, lalign = 0;
    long right, down, diag, pos = 0, nline, line, pos1, pos2; 
                                                // coordinates for seq1 and seq2
    char *alignSeq1, *alignSeq2;

    alignSeq1 = malloc((lseq1 + lseq2 + 1) * sizeof(char));
    alignSeq2 = malloc((lseq1 + lseq2 + 1) * sizeof(char));
    score = mScore[0][0]; // the final score is the top left elements of mScore

    /* track back and store the aligned sequences */
    while ((i < lseq2) || (j < lseq1))
    {
        // printf("lalign=%ld\n", lalign);

        // if comes from <direction>, then the value should be:
        if ((i < lseq2) && (j < lseq1)) // avoid out of boundary issue
        {
            diag = **(mScore + (i + 1) * lseq1 + j + 1) + **(mMatch + i * lseq1 + j);
            if (**(mScore + i * lseq1 + j) == diag)
            {
                
                *alignSeq1 = *(seq1 + j);
                *alignSeq2 = *(seq2 + i);
                alignSeq1++;
                alignSeq2++;
                i++;
                j++;
                lalign++;
                continue;
            }
        }
        if (j < lseq1)
        {
            right = **(mScore + i * lseq1 + j + 1) + sGap;
            if (**(mScore + i * lseq1 + j) == right)
            {

                *alignSeq1 = *(seq1 + j);
                *alignSeq2 = '-';
                alignSeq1++;
                alignSeq2++;
                j++;
                lalign++;
                continue;
            }
        }

        if (i < lseq2)
        {
            down = **(mScore + (i + 1) * lseq1 + j) + sGap;
            if (**(mScore + i * lseq1 + j) == down)
            {
                
                *alignSeq1 = '-';
                *alignSeq2 = *(seq2 + i);
                alignSeq1++;
                alignSeq2++;
                i++;
                lalign++;
                continue;
            }
        }
        i++;
        j++;
    }
    /* point to the beginning of aligned seq */
    alignSeq1 -= lalign;
    alignSeq2 -= lalign;
    alignSeq1 = realloc(alignSeq1, (lalign + 1) * sizeof(char));
    alignSeq2 = realloc(alignSeq2, (lalign + 1) * sizeof(char));
    *(alignSeq1 + lalign) = '\0';
    *(alignSeq2 + lalign) = '\0';

    /* print output */

    printf("\nScores:\n  Match: %d,  Mismatch: %d,  Gap: %d\n", sMatch, sMismatch, sGap);
    printf("\n\nAlignment score: %ld\nAlignment:\n", score);
    nline = (long)(lalign / 60) + 1;
    line = pos1 = pos2 = 0;
    while (line < nline)
    {
        /* seq 1 */
        printf("\n%ld:\t", pos1 + 1);
        for (pos = line * 60; (pos < line * 60 + 60) && (pos < lalign) ; pos++) // start position in aligned seq for each line
        {
            printf("%c", *(alignSeq1 + pos));
            if (*(alignSeq1 + pos) != '-') 
                pos1++;
        }
        printf("\n\t");
        /* print matching line */
        for (pos = line * 60; (pos < line * 60 + 60) && (pos < lalign) ; pos++)
        {
            if (*(alignSeq1 + pos) == *(alignSeq2 + pos))
                printf("*");
            else
                printf(" ");
        }
        /* seq 2 */
        printf("\n%ld:\t", pos2 + 1);
        for (pos = line * 60; (pos < line * 60 + 60) && (pos < lalign) ; pos++) // start position in aligned seq for each line
        {
            printf("%c", *(alignSeq2 + pos));
            if (*(alignSeq2 + pos) != '-')
                pos2++;
        }
        printf("\n");
        line++;
    }
    // printf("aligned seq 1: ");
    // for (j = 0; j < lalign; j++)
    // {
    //     printf("%c", *alignSeq1);
    //     alignSeq1++;
    // }
    // printf("\naligned seq 2: ");
    // for (i = 0; i < lalign; i++)
    // {
    //     printf("%c", *alignSeq2);
    //     alignSeq2++;
    // }
}