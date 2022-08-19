/*  Jingxuan Chen 
    09/03/2019
    Implement QuickSort to 
    sort fastq file by SEQ 
    in alphabetical order  
    To run the code, compile and type `./quicksort sample1M.fastq >> output` in terminal */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/* structure to store the fastq*/
typedef struct 
{
    char readname[100];
    char seq[100]; 
    char plus[100];
    char qual[100];
}fastq;

void swap(fastq **p1, fastq **p2);
void QuickSort(fastq *preads[], int min, int max);
int Partition(fastq *preads[], int min, int max);

int main(int argc, char *argv[])
{
    char *pinputFile = argv[1];
    char line[100];
    int nline = 0, nread, i, n = 1000001; // n is read count

    /* read files into an array of structure */
    fastq *reads = malloc(n * sizeof(fastq));
    fastq *preads[n]; 
    
    FILE *fp = NULL;
    fp = fopen(pinputFile,"r");

    if (fp == NULL){

        perror("Error opening the file");
        return -1;

    } else if (argc == 2){

        while (!feof(fp))
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
            
        }
    } else {
        perror("Error argv");
        return -1;
    }

    fclose(fp);
    fp = NULL;

    /* assign pointers to the array */
    for (i = 0; i < nread; i++)
    {
        
        preads[i] = &reads[i];

    }

    /* apply sorting algorithm and print out result */
    QuickSort(preads, 0, nread - 1);

    for(i = 0; i < nread; i++){
        printf("%s\n%s\n%s\n%s\n", preads[i]->readname, preads[i]->seq, preads[i]->plus, preads[i]->qual);
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

void QuickSort(fastq *preads[], int min, int max)
{
    if (min < max)
    {
        int pos = Partition(preads, min, max);
        QuickSort(preads, min, pos - 1);
        QuickSort(preads, pos + 1, max);
    }
}

int Partition(fastq *preads[], int min, int max)
{
    fastq *pivot = preads[max];
    int x = min, y = max - 1;

    while(1)
    {
        while (strcmp(pivot->seq, preads[x]->seq) > 1) x++;
        while (strcmp(pivot->seq, preads[y]->seq) < 1 && y > min) y--;
        if (x < y)
        {
            swap(&preads[x],&preads[y]);
        } else
        {
            swap(&preads[x], &preads[max]);
            return x;
        }
    }

}
