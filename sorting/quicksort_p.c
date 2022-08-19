#include <stdio.h>
#include <string.h>

int main(int argc, char *argv[])
{
    char *pinputFile = argv[1];
    char c, line[100];
    int nline = 0, nread, i;

    struct pfastq
    {
        FILE *preadname;
        FILE *pseq; // read length
        FILE *pplus;
        FILE *pqual;
    };
    struct pfastq reads[1001], *preads = reads;

    FILE *fp = NULL;
    fp = fopen(pinputFile, "r");
    printf("%p",fp);
    if(fp == NULL)
    {
        perror("Error opening the file");
        return -1;
    } else if (argc == 2)
    {
        reads[0].preadname = fp;
        while(!feof(fp))
        {
            fgets(line, sizeof(line), fp);
            line[strlen(line) - 1] = '\0';
            nline++;
            nread = (nline - 1) / 4;

            switch (nline % 4)
            {
            case 1:
                reads[nread].preadname = fp;
                break;
            case 2:
                reads[nread].pseq = fp;
                break;
            case 3:
                reads[nread].pplus = fp;
                break;
            case 0:
                reads[nread].pqual = fp;
                break;
            default:
                break;
            }

            // switch (nline % 4)
            // {
            // case 1:
            //     reads[nread].preadname = line;
            //     printf("preadname get\n%s\n", reads[nread].preadname);
            //     break;
            // case 2:
            //     reads[nread].pseq = line;
            //     break;
            // case 3:
            //     reads[nread].pplus = line;
            //     break;
            // case 0:
            //     reads[nread].pqual = line;
            //     printf("pqual get\n%s\n", reads[nread].pqual);
            //     break;
            // default:
            //     break;
            // }

            // while (*preads->preadname != '\0')
            // {
            //     printf("%c", *preads->preadname);
            //     preads->preadname++;
            // }

        }
    } else
    {
        perror("Error argv");
        return -1;
    }

    // do
    // {
    //     c = fgetc(reads[0].preadname);
    //     printf("%c", c);
    // } while (c != '\n');
    printf("%p", reads[0].preadname);
    printf("%p", reads[3].pqual);
    fclose(fp);
    fp = NULL;
    return 0;
}