/* Assignment 2 for BINF 9500
   Jingxuan Chen
   9/17/2019
   may need to include `-lm` 
   to link <math.h> when compiling
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


#define MAXD 100 // max 100 dimensions
#define MAXP 1000 // max 1000 points
#define MAXLEN 100 // max length of point information
#define FINDMIN(k) (int)(300 * k * sqrt(k)) // iteration times to find the minimum WCSS

typedef struct
{
    char info[MAXLEN + 2]; // the species information of point
    double data[MAXD]; // cordinates
    double cluster[MAXD]; // centroid cordinates
    double newCluster[MAXD];
    int clusterIndex; 
} point;

typedef struct {
    int seed;
    double WCSS;
} kmeansOut; // structure used to find the minimum WCSS

double dist(double ptData[], double ct[], int ndimension); // calculate distance between 2 points
int cmpCt(int npoint, int ndimension, point pt[]); // check if cluster changed
// int findMaxDistp(int npoint, int ndimension, point pt[]);
// implement kmeans
kmeansOut kmeans(int k, int ndimension, int npoint, point pt[], int seed, char *pinputFile, int checkMin);
double **kmeans_plus(int k, int ndimension, int npoint, point pt[]); // kmeans++ initialization
double findMinWCSS(int k, int ndimension, int npoint, point pt[], char *pinputFile);

int main(int argc, char *argv[])
{
    char *pinputFile = argv[1];
    int kmax , i, j, npoint = 0, ndimension = 0, len;
    double mean = 0, sd = 0, sum = 0, sqsum = 0; // for data normalization
    double WCSS, AIC, BIC;
    point *pt = malloc(MAXP * sizeof(point));    
    
    FILE *fp = NULL;
    fp = fopen(pinputFile, "r");
    
    if (argc == 1)
    {
        printf("\nUsage: ./program <file> <k>\n\n<k> represents the max k value\n\n");
        return 0;
    } else if (fp == NULL)
    {
        perror("File error!\n");
        return -1;
    } else if (argc != 3)
    {
        printf("\nUsage: ./program <file> <k>\n\n<k> represents the max k value\nPlease enter variables correctly.\n\n");
        return 0;
    } else
    {
        // alternatively, kmax = atoi(argv[2]);
        sscanf(argv[2], "%d", &kmax);
        /* skip the first line */
        fscanf(fp, "%*[^\n]\n", NULL);
        npoint = 0;
        while (!feof(fp))
        {

            fscanf(fp, "%[^\t]s", pt[npoint].info);
            i = 0;
            while(fgetc(fp) == '\t')
            {
                fscanf(fp, "%lf", &pt[npoint].data[i]);                
                i++;
                ndimension = i;
            }
            npoint++;
        }

        /* check if there is one empty point in the end
        (there should be one empty line in the end)
        */
        if(pt[npoint-1].info[0] == '\0') 
        {
            len = npoint - 1;
        } else
        {
            len = npoint;
        }
        // remove the empty point if exists
        pt = realloc(pt, len * sizeof(point));

        fclose(fp);


        /* find the mean and sd for each dimension and normalize */
        for (i = 0; i < ndimension; i++)
        {
            for (j = 0; j < len; j++) sum += pt[j].data[i];
            mean = sum / len;
            for (j = 0; j < len; j++) sqsum += (pt[j].data[i] - mean) * (pt[j].data[i] - mean);
            sd = sqrt(sqsum / (len - 1));
            for (j = 0; j < len; j++) pt[j].data[i] = (pt[j].data[i] - mean) / sd;
            sum = 0; sqsum = 0; // reset sums
            //printf("min[%d] is %lf, max[%d] is %lf\n", i, pmin[i], i, pmax[i]);
        }
        //printf("Start kmeans!\n");

        /* output of summary statistics */
        // char outFilename[50];
        // sprintf(outFilename, "%s_summary.txt", pinputFile);
        // fp = fopen(outFilename, "w+"); // output for plotting
        // fprintf(fp, "k\tMeanDistance\tWCSS\tAIC\tBIC\n");
        printf("k\tMeanDistance\tWCSS\tAIC\tBIC\n");
        for (i = 1; i <= kmax; i++)
        {
            WCSS = findMinWCSS(i, ndimension, len, pt, pinputFile); // find the min final WCSS
            AIC = (double)(2.0 * i * ndimension + WCSS);
            BIC = (double)(log((double) len) * i * ndimension + WCSS);
            // fprintf(fp, "%d\t%.5lf\t%.4lf\t%.4lf\t%.4lf\n", i, sqrt(WCSS / npoint), WCSS, AIC, BIC);
            printf("%d\t%.5lf\t%.4lf\t%.4lf\t%.4lf\n", i, sqrt(WCSS / npoint), WCSS, AIC, BIC);
            
        }
        // fclose(fp);

    }
    
    free(pt);
    return 0;
}

double findMinWCSS(int k, int ndimension, int npoint, point pt[], char *pinputFile)
{
    // int FINDMIN = (int)(300 * k * sqrt(k));
    kmeansOut minKmOut, kmOut[FINDMIN(k)];
    int i, j, minIndex;
    double minWCSS, pt2ct;

    for (i = 0; i < FINDMIN(k); i++)
    {
        // last parameter used to check if should print output
        kmOut[i] = kmeans(k, ndimension, npoint, pt, i, pinputFile, -1);
        // printf("%d finished\n", i);
        // printf("WCSS %d is %.4lf\t seed is %d\n", i, kmOut[i].WCSS, kmOut[i].seed);
    }
    minIndex = 0;
    minWCSS = kmOut[0].WCSS;
    for (i = 1; i < FINDMIN(k); i++)
    {
        if (minWCSS > kmOut[i].WCSS)
        {
            minWCSS = kmOut[i].WCSS;
            minIndex = i;
        }
    }
    // find the min WCSS run and print output
    kmeans(k, ndimension, npoint, pt, kmOut[minIndex].seed, pinputFile, 1);

    return minWCSS;
}

/* implement kmeans process and then return WCSS and seed */
kmeansOut kmeans(int k, int ndimension, int npoint, point pt[], int seed, char* pinputFile, int checkMin)
{
    kmeansOut kmOut;
    double** ct;
    // double ct[k][ndimension];
    int i, j, l, x, minIndex, count, maxIndex;
    double pt2ct, minDist, thisWCSS = 0.0, len, sumDist, randNum, sum, maxDist; // the distance of point to each centroid

    // srand48(time(0));

    kmOut.seed = seed;
    srand48(kmOut.seed);

    /* initialization for kmeans */ 

    // for (i = 0; i < k; i++)
    // {
    //     /* random x-th point as centroid*/
    //     x = (int) (drand48() * npoint);
    //     // printf("the x is %d\n",x);
    //     for (j = 0; j < ndimension; j++)
    //     {
    //         ct[i][j] = pt[x].data[j];
    //         // printf("the cent is %lf\n", ct[i][j]);
    //     }
    // }

    ct = malloc(k * sizeof(double*));
    // alternatively, ct[0] = (double *)malloc(k * ndimension * sizeof(double));
    for (i = 0; i <= k; i++)
    {
        ct[i] = malloc(ndimension * sizeof(double));
    }

    ct = kmeans_plus(k, ndimension, npoint, pt);
    // printf("kmeans++ finished\n");

    /* for debug use */
    // for (i = 0; i < k; i++)
    // {
    //     printf("Centroid %d\t", i);
    //     for (j = 0; j < ndimension; j++)
    //     {
    //         printf("%.2lf\t", ct[i][j]);
    //     }
    //     printf("\n");
    // }
    // printf("Initialization finished!\n");

    /* first assignment step */

    for (i = 0; i < npoint; i++)
    {
        minDist = dist(pt[i].data, ct[0], ndimension);
        minIndex = 0;
        for (j = 1; j < k; j++)
        {
            pt2ct = dist(pt[i].data, ct[j], ndimension);
            if (minDist > pt2ct)
            {
                minDist = pt2ct;
                minIndex = j;
            }
        }
        pt[i].clusterIndex = minIndex;
        for(j = 0; j < ndimension; j++)
        {
            pt[i].newCluster[j] = ct[minIndex][j];
        }
    }
    /* iteration */
    //int cc = 0;
    do
    {
        /* Assign new cluster to cluster */
        for(i = 0; i < npoint; i++)
        {
            for (j = 0; j < ndimension; j++)
            {
                pt[i].cluster[j] = pt[i].newCluster[j];
            }
        }

        /* update centroid - find the mean of all points in one cluster */

        for (i = 0; i < k; i++)
        {
            for (x = 0; x < ndimension; x++)
            {
                sum = 0.0;
                count = 0;
                for (j = 0; j < npoint; j++)
                {
                    if (pt[j].clusterIndex == i)
                    {
                        sum += pt[j].data[x];
                        count++;
                    }
                }
                ct[i][x] = (double)(sum / count);
            }
        }

        /* update clustering */
        for (i = 0; i < npoint; i++)
        {
            minDist = dist(pt[i].data, ct[0], ndimension);
            minIndex = 0;
            for (j = 1; j < k; j++)
            {
                pt2ct = dist(pt[i].data, ct[j], ndimension);
                if (minDist > pt2ct)
                {
                    minDist = pt2ct;
                    minIndex = j;
                }
            }
            pt[i].clusterIndex = minIndex;
            for (j = 0; j < ndimension; j++)
            {
                pt[i].newCluster[j] = ct[minIndex][j];
            }
        }

        /* Check for empty and 1-point clusters 
        reassign points if it exists */  // not necessary if run FUNMIN times
        // for (i = 0; i < k; i++)
        // {
        //     count = 0;
        //     for (j = 0; j < npoint; j++)
        //     {
        //         if (pt[j].clusterIndex == i)
        //             count++;
        //     }
        //     if (count <= 1)
        //     {
        //         x = findMaxDistp(npoint, ndimension, pt);
        //         pt[x].clusterIndex = i;
        //         for (l = 0; l < ndimension; l++)
        //         {
        //             pt[x].newCluster[l] = ct[i][l];
        //         }
        //     }
        // }

    } while (cmpCt(npoint, ndimension, pt) == -1); // check if cluster = new cluster

    //printf("kmeans finished!\n");

    /* Check for empty and 1-point clusters 
    rerun kmeans process if it exists */ // not necessary if run FUNMIN times
    // for(i = 0; i < k; i++)
    // {
    //     count = 0;
    //     for (j = 0; j < npoint; j++)
    //     {
    //         if(pt[j].clusterIndex == i) count++;
    //     }
    //     if (count <= 1) goto repeat;
    // }
    // // for (i = 0; i < k; i++)
    // {
    //     printf("Cluster %d:", i + 1);
    //     for (j = 0; j < ndimension; j++)
    //         printf("\t%.2lf", ct[i][j]);
    //     printf("\n");
    // }


    for (i = 0; i < npoint; i++)
    {
        pt2ct = dist(pt[i].data, pt[i].cluster, ndimension);
        thisWCSS += pt2ct;
    }

    kmOut.WCSS = thisWCSS;

    /* will not output to a file if this is not the min WCSS run */
    if (checkMin == -1)
        goto noPrint;

    /* output to a file */
    FILE *fp;
    char outFilename[50];
    sprintf(outFilename, "%s_k_%d.txt", pinputFile, k);
    fp = fopen(outFilename, "w+");

    fprintf(fp, "Detailed output for k=%d\n", k);
    fprintf(fp, "\nList of data points in each cluster preceded by distance to centroid:\n");
    for (i = 0; i < k; i++)
    {
        fprintf(fp, "\nCluster %d:\n", i + 1);
        for (j = 0; j < npoint; j++)
        {
            if (pt[j].clusterIndex == i)
            {
                /* distance of each point to assigned centroid */
                pt2ct = dist(pt[j].data, pt[j].cluster, ndimension);
                fprintf(fp, "\t%.2lf\t%s\n", sqrt(pt2ct), pt[j].info);
            }
        }
    }
    fprintf(fp, "\n\nMutual pairwise distances among centroids:\n");
    for (i = 0; i < k; i++)
    {
        fprintf(fp, "Cluster %d\t", i + 1);
    }
    fprintf(fp, "\n");
    for (i = 0; i < k; i++)
    {
        for (j = 0; j < k; j++)
        {
            if (i <= j)
                fprintf(fp, "\t%.4lf\t", dist(ct[i], ct[j], k));
            else
                fprintf(fp, "\t\t");
        }
        fprintf(fp, "Cluster %d\n", i + 1);
    }

    fprintf(fp, "\n\nCentroid coordinates:\n");
    for (i = 0; i < k; i++)
    {
        fprintf(fp, "Cluster %d:", i + 1);
        for (j = 0; j < ndimension; j++)
            fprintf(fp, "\t%.2lf", ct[i][j]);
        fprintf(fp, "\n");
    }
    fclose(fp);

    noPrint:
    return kmOut;
}

/* calculate the squared distance of 2 vecters */
double dist(double ptData[], double ct[], int ndimension)
{
    int i;
    double sum = 0;
    for (i = 0; i < ndimension; i++)
    {
        sum += ((ptData[i] - ct[i]) * (ptData[i] - ct[i]));
    }
    return sum;
}

/* check if the cluster changed */
int cmpCt(int npoint, int ndimension, point pt[])
{
    int i, j;
    for (i = 0; i < npoint; i++)
    {
        for (j = 0; j < ndimension; j++)
        {
            if (pt[i].cluster[j] != pt[i].newCluster[j])
                return -1;
        }
    }
    return 1;
}

// int findMaxDistp(int npoint, int ndimension, point pt[]) // not necessary if not avoid empty cluster
// {
//     int p, i, maxIndex;
//     double maxDist, d;
//     maxDist = dist(pt[0].data, pt[0].cluster, ndimension);
//     maxIndex = 0;
//     for (p = 1; p < npoint; p++)
//     {
//         d = dist(pt[p].data, pt[p].cluster, ndimension);
//         if (maxDist < d)
//         {
//             maxDist = d;
//             maxIndex = p;
//         }
//     }
//     return maxIndex;
// }

double** kmeans_plus(int k, int ndimension, int npoint, point pt[])
{
    int i, j, x, l;
    double len, sumDist, randNum, pt2ct;
    double **ct;
    ct = malloc(k * sizeof(double*));
    for (i = 0; i <= k; i++)
    {
        ct[i] = malloc(ndimension * sizeof(double));
    }
    /* kmeans++ initialization */
    for (i = 0; i < k; i++)
    {
        if (i == 0) // initialize 1st centroid
        {
            x = (int)(drand48() * npoint);
            // printf("the x is %d\n",x);
            for (j = 0; j < ndimension; j++)
            {
                ct[i][j] = pt[x].data[j];
                //printf("the cent is %lf\n", ct[i][j]);
            }
        }
        else // initialize next few centroids
        {
            len = 0.0;
            sumDist = 0.0;
            randNum = drand48();
            // printf("%d\trand=%lf\n", i,randNum);
            for (l = 0; l < npoint; l++)
            {
                pt2ct = dist(pt[l].data, ct[i - 1], ndimension);
                sumDist += pt2ct;
            }

            for (l = 0; l < npoint; l++)
            {
                pt2ct = dist(pt[l].data, ct[i - 1], ndimension);
                len += (pt2ct / sumDist);
                if (len > randNum)
                {
                    x = l;
                    for (j = 0; j < ndimension; j++)
                    {
                        ct[i][j] = pt[x].data[j];
                        //printf("the cent is %lf\n", ct[i][j]);
                    }
                    break;
                }
            }
        }
    }
    return ct;
}