#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include "mysac64.h"
#include <pthread.h>
using namespace std;

#define max(a, b) ( ((a) > (b)) ? (a) : (b) )
#define min(a, b) ( ((a) < (b)) ? (a) : (b) )

pthread_mutex_t fiolock;

int
jday(int y, int m, int d)
{
    int i, jd = 0;

    for (i = 1; i < m; i++) {
        if ( (i == 1) || (i == 3) || (i == 5) || (i == 7) || (i == 8) || (i == 10) ) { jd += 31; } else if (i == 2) {
            if ( (y % 400 == 0) || (y % 100 != 0 && y % 4 == 0) ) jd += 29;
            else jd += 28;
        } else { jd += 30; }
    }
    return jd + d;
}

double
abs_time(int yy, int jday, int hh, int mm, int ss, int ms)
{
    // computes time in s relative to 1900
    int nyday = 0, i;

    for (i = 1901; i < yy; i++) {
        if ( (i % 400 == 0) || (i % 100 != 0 && i % 4 == 0) ) nyday += 366;
        else nyday += 365;
    }
    return 24. * 3600. * (nyday + jday) + 3600. * hh + 60. * mm + ss + 0.001 * ms;
}

void
UpdateTime(SAC_HD * shd)
{
    if (shd->nzmsec < 1000 && shd->nzmsec >= 0) return;

    int i = (int) floor(shd->nzmsec / 1000);
    shd->nzmsec -= i * 1000;
    shd->nzsec  += i;
    if (shd->nzsec < 60 && shd->nzsec >= 0) return;

    i = (int) floor(shd->nzsec / 60);
    shd->nzsec -= 60 * i;
    shd->nzmin += i;
    if (shd->nzmin < 60 && shd->nzmin >= 0) return;

    i = (int) floor(shd->nzmin / 60);
    shd->nzmin  -= i * 60;
    shd->nzhour += i;
    if (shd->nzhour < 24) return;

    shd->nzhour -= 24;
    shd->nzjday++;
    if (shd->nzjday < 366) return;

    if ( ((shd->nzyear % 400 == 0) || (shd->nzyear % 100 != 0 && shd->nzyear % 4 == 0)) && shd->nzjday < 367) return;

    shd->nzjday = 1;
    shd->nzyear++;
}

bool
read_shd(char * fname, SAC_HD * SHD)
{
    if (!SHD) return false;  // SHD = &SAC_HEADER;

    FILE * fsac;
    // SAC_HD *SHD;
    if ((fsac = fopen(fname, "r")) == NULL) return false;

    pthread_mutex_lock(&fiolock);
    fread(SHD, sizeof(SAC_HD), 1, fsac);
    fclose(fsac);
    pthread_mutex_unlock(&fiolock);

    return true;
}

SAC_HD *
read_sac(char * fname, float ** sig, SAC_HD * SHD)
{
    FILE * fsac;

    if ((fsac = fopen(fname, "rb")) == NULL) return NULL;

    pthread_mutex_lock(&fiolock);
    if (!SHD) SHD = &SAC_HEADER;
    fread(SHD, sizeof(SAC_HD), 1, fsac);
    if (*sig == NULL) *sig = (float *) malloc(SHD->npts * sizeof(float));
    // else *sig = (float *) realloc (*sig, SHD->npts * sizeof(float));
    fread(*sig, sizeof(float), SHD->npts, fsac);
    fclose(fsac);
    pthread_mutex_unlock(&fiolock);

    /*-------------  calcule de t0  ----------------*/
    {
        int eh, em, i;
        float fes;
        char koo[9];

        for (i = 0; i < 8; i++) koo[i] = SHD->ko[i];
        koo[8] = 0;

        SHD->o = SHD->b + SHD->nzhour * 3600. + SHD->nzmin * 60
          + SHD->nzsec + SHD->nzmsec * .001;

        sscanf(koo, "%d%*[^0123456789]%d%*[^.0123456789]%g", &eh, &em, &fes);

        SHD->o -= (eh * 3600. + em * 60. + fes);
        /*-------------------------------------------*/ }
    return SHD;
}

void
write_sac(const char * fname, float * sig, SAC_HD * SHD)
{
    FILE * fsac;

    if ( (fsac = fopen(fname, "wb")) == NULL) {
        cerr << "ERROR(write_sac): Cannot open file " << fname << endl;
        return;
    }

    if (!SHD) SHD = &SAC_HEADER;
    SHD->iftype    = (int) ITIME;
    SHD->leven     = (int) TRUE;
    SHD->lovrok    = (int) TRUE;
    SHD->internal4 = 6L;

    /*+++++++++++++++++++++++++++++++++++++++++*/
    SHD->depmin = sig[0];
    SHD->depmax = sig[0];
    int i;
    for (i = 0; i < SHD->npts; i++) {
        if (SHD->depmin > sig[i]) SHD->depmin = sig[i];
        else if (SHD->depmax < sig[i]) SHD->depmax = sig[i];
    }

    // memset(SHD, 0, sizeof(SAC_HD));
    pthread_mutex_lock(&fiolock);
    fwrite(SHD, sizeof(SAC_HD), 1, fsac);
    fwrite(sig, sizeof(float), SHD->npts, fsac);
    pthread_mutex_unlock(&fiolock);

    fclose(fsac);
}

int
read_rec(int rec_flag, char * fname, int len, int * rec_b, int * rec_e, int * nrec)
{
    FILE * frec;
    int irec;

    if (rec_flag) {
        if ((frec = fopen(fname, "r")) == NULL) return 0;

        pthread_mutex_lock(&fiolock);
        for (irec = 0;; irec++)
            if (fscanf(frec, "%d %d", &rec_b[irec], &rec_e[irec]) != 2) break;
        *nrec = irec;
        fclose(frec);
        pthread_mutex_unlock(&fiolock);
        if (irec == 0) return 0;
    } else   {
        rec_b[0] = 0;
        rec_e[0] = len - 1;
        *nrec    = 1;
    }
    return 1;
}
