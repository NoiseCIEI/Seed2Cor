#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include "mysac64.h"
#include "64_sac_db.h"
using namespace std;

//global control
extern   SAC_DB *sdb;
extern   int NTHRDS;
extern   int imonth;
extern   int currevn;
extern   int fdprompt;
extern   struct NOTE *reports;
extern   pthread_attr_t attr_j;
extern   pthread_mutex_t cevlock, fftlock, fiolock;
struct NOTE {
   char *head;
   char *tail;
};
//parameters
extern   char rdsexe[200];
extern   char evrexe[200];
extern   char stalst[200];
extern   char seedlst[200];
extern   char ch[4];
extern   int sps;
extern   float gapfrac;
extern   float t1;
extern   float tlen;
extern   float perl, perh;
extern   int tnorm_flag;
extern   float Eperl, Eperh;
extern   float timehlen;
extern   float frechlen;
extern   char fwname[200];
extern   int ftlen;
extern   int fprcs;
extern   float memomax;
extern   int lagtime;
extern   int mintlen;
extern   int fdel1;
extern   int fdel2;
extern   int fskip1;
extern   int fskip2;
extern   int fskip3;
extern   int fskip4;
extern   int CorOutflag;

/* Finction prorotypes */ 
void fRemove (const char *fname);

int dRemove(char *dirname);

void Move (const char *oldname, const char *newname);

void MKDir(const char *dirname);

char * List(const char *dir, const char *pattern, int type, int *nfile);

SAC_HD *read_sac (char *fname, float **sig, SAC_HD *SHD);

void write_sac (const char *fname, float *sig, SAC_HD *SHD);

SAC_HD *read_shd (char *fname, SAC_HD* shd);

void UpdateTime(SAC_HD *shd);

int read_rec(int rec_flag, char *fname, int len, int *rec_b, int *rec_e, int *nrec);

void Filter (double f1, double f2, double f3, double f4, double dt, int n, float *seis_in, float *seis_out);

/*
void PrintParams() {
   cout<<sdb<<endl;
   cout<<imonth<<endl;
   cout<<rdsexe<<endl;
   cout<<sacexe<<endl;
   cout<<stalst<<endl;
   cout<<seedlst<<endl;
   cout<<ch<<endl;
   cout<<sps<<endl;
   cout<<gapfrac<<endl;
   cout<<t1<<endl;
   cout<<tlen<<endl;
   cout<<perl<<" "<<perh<<endl;
   cout<<tnorm_flag<<endl;
   cout<<Eperl<<" "<<Eperh<<endl;
   cout<<timehlen<<endl;
   cout<<frechlen<<endl;
   cout<<fwname<<endl;
   cout<<ftlen<<endl;
   cout<<fprcs<<endl;
   cout<<memomax<<endl;
   cout<<lagtime<<endl;
   cout<<mintlen<<endl;
   cout<<fdel1<<endl;
   cout<<fdel2<<endl;
   cout<<fskip1<<endl;
   cout<<fskip2<<endl;
   cout<<fskip3<<endl;
   cout<<fskip4<<endl;
}
*/
