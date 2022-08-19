/*  Jingxuan Chen 
    09/03/2019
    Implement QuickSort to 
    sort fastq file by SEQ 
    in alphabetical order  
    To run the code, compile and type `./quicksort sample1M.fastq >> output` in terminal*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define LINELEN 100    /* This is useful because if you want to change the value you can 
                              do it here and you don't have to read through the ode and find
                              every place you use it. This is processed by the preprocessor,
                              which replaces LINELEN with 100 in the code before it is 
                              passed to the compiler. See, e.g., https://en.wikipedia.org/wiki/C_preprocessor */

/* structure to store the fastq*/
typedef struct 
{
    char readname[LINELEN+2];   // note that LINELEN is not a variable but a constant -- it is replaced with 100 before the compiler is called.
    char seq[LINELEN+2];        // I added +2 because you need one character for the terminal \0, adn one for \n so if the line contains 100 characters, you need 102 to store the string
    char plus[LINELEN+2];       //  you are still using far more memory than you need because most lines are shorter
    char qual[LINELEN+2];
//    char readname[100];
//    char seq[100]; 
//    char plus[100];
//    char qual[100];
}fastq;

void swap(fastq **p1, fastq **p2);
void QuickSort(fastq *preads[], long min, long max);
long Partition(fastq *preads[], long min, long max);

int main(int argc, char *argv[])
{
    char *pinputFile = argv[1];
    char line[LINELEN+2];
//    char line[100];

// I change all 'int' to 'long'; you cannot rely on int being interpreted as long, even if it is the case on 64-bit computers
    long nline = 0, nread, i, n = 1000001; // n is read count     // You could also use #define for that, too
// you could also find first how many reads are in the file and reserve only the amount of memory that you need

    /* read files into an array of structure */
    fastq *reads = malloc(n * sizeof(fastq));
//    fastq *preads[n];     //   I changed it to dynamic allocation to avoid making the big array in stack
    fastq **preads = malloc(n * sizeof(fastq*)); 
    
    FILE *fp = NULL;   //   I don't think you need to assign fp to NULL because fopen returns NULL upon error
    fp = fopen(pinputFile,"r");

    if (fp == NULL){

        perror("Error opening the file");
        return -1;

    } else if (argc == 2){


        nread=0;    // I thought rewriting this part could make it faster but it made almost no difference
        while (!feof(fp))
        {
            if (fgets(reads[nread].readname, LINELEN-2, fp) != NULL)   // I added the -2 to make sure there is room for the \n and \0
// this is still problematic because if there are more characters than LINELEN it will read the rest of the same line with the next fgets and everything will get messed up 
            {
              fgets(reads[nread].seq, LINELEN-2, fp); 
              fgets(reads[nread].plus, LINELEN-2, fp);
              fgets(reads[nread].qual, LINELEN-2, fp);
              ++nread;
            };
            
        }
/*        while (!feof(fp))
        {

            nline++; // count the line number
            nread = (nline - 1) / 4; // count the read number

            switch (nline % 4)
            {
            case 1:
                fgets(reads[nread].readname, sizeof(reads[nread].readname), fp);
                reads[nread].readname[strlen(reads[nread].readname) - 1] = '\0';
                break;
            case 2:
                fgets(reads[nread].seq, sizeof(reads[nread].seq), fp);
                reads[nread].seq[strlen(reads[nread].seq) - 1] = '\0';
                break;
            case 3:
                fgets(reads[nread].plus, sizeof(reads[nread].plus), fp);
                reads[nread].plus[strlen(reads[nread].plus) - 1] = '\0';
                break;
            case 0:
                fgets(reads[nread].qual, sizeof(reads[nread].qual), fp);
                reads[nread].qual[strlen(reads[nread].qual) - 1] = '\0';
                break;
            default:
                break;
            }
            
        }*/
    } else {
        perror("Error argv");    // better thing to do is print a brief instructuion for the user how to run the program, e.g., "quicksort <FastqFile>"
        return -1;
    }

    fclose(fp);
//    fp = NULL;    //   I don't think you need to do this because you never use fp again.

    /* assign pointers to the array */
    for (i = 0; i < nread; i++)
    {
        
        preads[i] = &reads[i];

    }

    /* apply sorting algorithm and print out result */
    QuickSort(preads, 0, nread - 1);

    for(i = 0; i < nread; i++){
//        printf("%s\n%s\n%s\n%s\n", preads[i]->readname, preads[i]->seq, preads[i]->plus, preads[i]->qual);
        printf("%s%s%s%s", preads[i]->readname, preads[i]->seq, preads[i]->plus, preads[i]->qual);  // I had to change this because I left the \n in the string
    }

    free(reads);
    return 0;
}

void swap(fastq **p1, fastq **p2)
{
    fastq *temp = *p1;
    *p1 = *p2;
    *p2 = temp;
}

void QuickSort(fastq *preads[], long min, long max)
{
    if (min < max)
    {
        long pos = Partition(preads, min, max);
        QuickSort(preads, min, pos - 1);
        QuickSort(preads, pos + 1, max);
    }
}

long Partition(fastq *preads[], long min, long max)
{
    fastq *pivot = preads[max];
    long x = min, y = max - 1;

    while(1)
    {
        while (strcmp(pivot->seq, preads[x]->seq) > 0) x++;   // it should be 0, not 1
        while (y>min && strcmp(pivot->seq, preads[y]->seq) < 0) y--;   // If you are using the modified Hoare partition the pivot in the last
                                                         // spot you need also the y>min in the condition. To understand why, imagine what happens
                                                         // if the smallest element becomes the pivot. I think this is why the program was crashing.
//        while (strcmp(pivot->seq, preads[x]->seq) > 1) x++;
//        while (strcmp(pivot->seq, preads[y]->seq) < 1) y--;
        if (x < y)
        {
            swap(&preads[x],&preads[y]);
            ++x;    // The Hoare partition is usually implemented with do-while loops; if you use while loop you need to add this
            --y;
        } else
        {
            swap(&preads[x], &preads[max]);
            return x;
        }
    }

}
