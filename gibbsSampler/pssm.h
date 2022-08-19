#ifndef PSSM_H_ 
#define PSSM_H_

double *bg(int *seqInt, int len); // background prob for GC and AT
double **fMatrix(char **align, int lMotif, int cMotif); // frequency matrix
double **freq2prob(double **freqMatrix, int lMotif); // probability matrix
double **prob2score(double **probMatrix, int lMotif, double *bgProb); // score matrix
int getIndex(char n); // change ACGT to 0123
// char *int2seq(int *seqInt, int l); // change int to seq for printing
double seqScore(double **scoreMatrix, int lMotif, int *motifSeqInt); // calculate int seq score
double motifScore(double **scoreMatrix, int lMotif, char *motifSeq); // calculate alpha seq score

#endif