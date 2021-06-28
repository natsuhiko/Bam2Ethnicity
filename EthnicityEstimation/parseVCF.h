#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <string.h>
#include <math.h>
#include <zlib.h>

#define FORMAT_GT 1
#define FORMAT_GL 2
#define FORMAT_AP 3
#define FORMAT_PP 4
#define FORMAT_AS 5
#define FORMAT_RD 6
#define FORMAT_BF 7
#define FORMAT_DS 8
#define FORMAT_OTHER 0

#define FORMAT_AP_ORIG 0

#define VT_SNP 1
#define VT_INDEL 2
#define VT_SV 3
#define VT_OTHER 4

#define BUFS 30
#define CELLS 10000

typedef struct{
        int VT;
        double RSQ;
        double AF;
        double CER;
}VCF_info;

int LOG10;
void Log10();

double getRsq(int* dip1, int* dip2, long N);
double getAF(double* gen, double* w, long N);
double getHWE(int* dip, double* w, long N);
double getIA(double* gl, double* w, long N);
double getCR(double* gl, double* w, long N);
int init();
int parseHeader(int argc, char** argv);
int parseFormat(char* str, int* formatID);
int parseInfo(char* str, VCF_info* vinfo);
int parseCell(char* cell, int* dip, double* gl, double* ap, long* ase, double* rd, int* formatID);

int Fread_parseVCF(char* buf, FILE* fp);
int parseLine0(char* chr, long* pos, char* rs, char* al, VCF_info* vinfo, int* dip, double* gen, double* gl, double* ap, long* ase, double* rd, long N, int genType, FILE* fp);
int parseLine(char* chr, long* pos, char* rs, char* al, VCF_info* vinfo, int* dip, double* gen, double* gl, double* ap, long* ase, double* rd, long N, int genType);
int parseHap(char* hap);
int parseGL(char* gl);
void em(double* g, double* h, int* dip);


/*
char* format=NULL;
char* buf=NULL;
char* cell=NULL;
int ret=0;
int i0=0;
int* formatID;
*/
