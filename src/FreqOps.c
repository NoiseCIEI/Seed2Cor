#include <stdlib.h>
#include <fftw3.h>
#include <math.h>
#include <string.h>
#include <iostream>
#include <pthread.h>
using namespace std;

#define PI 3.14159265358979323846

extern pthread_mutex_t fftlock;

void
FFTW_F(fftw_plan plan, fftw_complex * out, int ns, float * seis, int n)
{
    fftw_execute(plan);
    pthread_mutex_lock(&fftlock);
    fftw_destroy_plan(plan);
    pthread_mutex_unlock(&fftlock);
    int k;
    for (k = 0; k < n; k++) seis[k] = out[k][0];
    // for(k=0; k<n; k+=1000) if(seis[k] != 0.) printf("%d %f\n", k, seis[k]);
}

void
FFTW_B(int type, float * seis, int n, fftw_complex ** in, fftw_complex ** out, int * nso, fftw_plan * planF, int Fflag)
{
    int ns = (int) (log((double) n) / log(2.)) + 1;

    if (ns < 13) ns = 13;
    ns   = (int) pow(2, ns);
    *nso = ns;
    *in  = (fftw_complex *) fftw_malloc(ns * sizeof(fftw_complex) );// fftw_alloc_complex(ns);
    *out = (fftw_complex *) fftw_malloc(ns * sizeof(fftw_complex) ); // fftw_alloc_complex(ns);

    // measure plan using the allocated in/out blocks
    pthread_mutex_lock(&fftlock);
    fftw_plan plan = fftw_plan_dft_1d(ns, *in, *out, FFTW_BACKWARD, type); // FFTW_ESTIMATE / FFTW_MEASURE
    if (Fflag == 1) *planF = fftw_plan_dft_1d(ns, *out, *in, FFTW_FORWARD, type);
    if (plan == NULL || (Fflag == 1 && *planF == NULL) ) {
        fprintf(stderr, "Error(FFTW_B): fftw_plan creation failed!!\n");
        pthread_mutex_unlock(&fftlock);
        exit(0);
    }
    pthread_mutex_unlock(&fftlock);
    // initialize input array and excute
    memset(*in, 0, ns * sizeof(fftw_complex));
    int k;
    for (k = 1; k < n; k++) (*in)[k][0] = seis[k];
    fftw_execute(plan);
    // cleanup
    pthread_mutex_lock(&fftlock);
    fftw_destroy_plan(plan);
    pthread_mutex_unlock(&fftlock);
    // if( Fflag==0 ) fftw_free(*in);

    // kill half spectrum and correct ends
    int nk = ns / 2 + 1;
    for (k = nk; k < ns; k++) {
        (*out)[k][0] = 0.;
        (*out)[k][1] = 0.;
    }
    (*out)[0][0]     *= 0.5;
    (*out)[0][1]     *= 0.5;
    (*out)[nk - 1][1] = 0.;
} /* FFTW_B */

void
TaperL(double f3, double f4, double dom, int nk, fftw_complex * sf)
{
    double f, ss;

    int i = (int) ceil(f3 / dom);

    for (f = i * dom; f < f4; i++, f += dom) {
        ss        = ( 1. + cos(PI * (f3 - f) / (f4 - f3)) ) * 0.5;
        sf[i][0] *= ss;
        sf[i][1] *= ss;
    }
    for (; i < nk; i++) {
        sf[i][0] = 0.;
        sf[i][1] = 0.;
    }
}

void
TaperB(double f1, double f2, double f3, double f4, double dom, int nk, fftw_complex * sf)
{
    int i;
    double f, ss;

    for (i = 0, f = 0.; f < f1; i++, f += dom) {
        sf[i][0] = 0.;
        sf[i][1] = 0.;
    }
    for (; f < f2; i++, f += dom) {
        ss        = ( 1. - cos(PI * (f1 - f) / (f2 - f1)) ) * 0.5;
        sf[i][0] *= ss;
        sf[i][1] *= ss;
    }
    i = (int) ceil(f3 / dom);
    for (f = i * dom; f < f4; i++, f += dom) {
        ss        = ( 1. + cos(PI * (f3 - f) / (f4 - f3)) ) * 0.5;
        sf[i][0] *= ss;
        sf[i][1] *= ss;
    }
    for (; i < nk; i++) {
        sf[i][0] = 0.;
        sf[i][1] = 0.;
    }
}

void
Filter(double f1, double f2, double f3, double f4, double dt, int n, float * seis_in, float * seis_out)
{
    if (f4 > 0.5 / dt) {
        fprintf(stdout, "\n*** Warning(Filter): filter band out of range! ***\n");
        f4 = 0.49999 / dt;
    }
    fftw_plan planF = NULL;
    // backward FFT: s ==> sf
    int ns;
    fftw_complex * s, * sf;
    FFTW_B(FFTW_ESTIMATE, seis_in, n, &s, &sf, &ns, &planF, 1);
    // make tapering
    int nk     = ns / 2 + 1;
    double dom = 1. / dt / ns;
    if (f2 == -1.) TaperL(f3, f4, dom, nk, sf);
    else TaperB(f1, f2, f3, f4, dom, nk, sf);

    // forward FFT: sf ==> s
    FFTW_F(planF, s, ns, seis_out, n);
    fftw_free(s);
    fftw_free(sf);

    // forming final result
    int k;
    float ftmp = 2. / ns;
    for (k = 0; k < n; k++) {
        if (seis_in[k] == 0) seis_out[k] = 0.;
        else seis_out[k] *= ftmp;
    }
}

void
IFFT(int nlen, float * amp, float * pha, float * seis_out, int * nsig)
{
    int i, ns = (nlen - 1) * 2;
    fftw_complex * in, * out;

    // measure the plan
    in  = (fftw_complex *) fftw_malloc(ns * sizeof(fftw_complex) );
    out = (fftw_complex *) fftw_malloc(ns * sizeof(fftw_complex) );
    pthread_mutex_lock(&fftlock);
    fftw_plan planF = fftw_plan_dft_1d(ns, in, out, FFTW_FORWARD, FFTW_MEASURE); // FFTW_ESTIMATE / FFTW_MEASURE / FFTW_PATIENT
    if (planF == NULL) {
        fprintf(stderr, "Error(IFFT): fftw_plan creation failed!!\n");
        pthread_mutex_unlock(&fftlock);
        exit(0);
    }
    pthread_mutex_unlock(&fftlock);
    // initialize and excute
    memset(in, 0, ns * sizeof(fftw_complex));
    for (i = 0; i < nlen; i++) {
        in[i][0] = amp[i] * cos(pha[i]);
        in[i][1] = -amp[i] * sin(pha[i]);
    }
    FFTW_F(planF, out, ns, seis_out, ns);

    fftw_free(in);
    fftw_free(out);
    *nsig = ns;
}

void
fDiv(double dom, int nsig, fftw_complex * sf, double * freq, double * amp, double * pha, int ntra)
{
    int isig, itra = 1;
    double f, ampcur, phacur, realsig, imagsig;
    double sintmp, costmp;

    for (isig = (int) ceil(freq[0] / dom); isig < nsig; isig++) {
        f = isig * dom;
        while (f > freq[itra]) {
            itra++;
            if (itra >= ntra) break;
        }
        // interpolate to get current amp and pha
        sintmp = (f - freq[itra - 1]) / (freq[itra] - freq[itra - 1]);
        ampcur = amp[itra - 1] + (amp[itra] - amp[itra - 1]) * sintmp;
        phacur = pha[itra - 1] + (pha[itra] - pha[itra - 1]) * sintmp;
        // divide sf by (ampcur, phacur)
        realsig     = sf[isig][0];
        imagsig     = sf[isig][1];
        sintmp      = sin(phacur);
        costmp      = cos(phacur);
        sf[isig][0] = (realsig * costmp - imagsig * sintmp) / ampcur;
        sf[isig][1] = (imagsig * costmp + realsig * sintmp) / ampcur;
    }
}

void
FDivide(double f1, double f2, double f3, double f4, double dt, int n, float * seis_in, float * seis_out, double * freq,
  double * amp, double * pha, int nf)
{
    if (f4 > 0.5 / dt) {
        fprintf(stdout, "\n*** Warning(FDivide): filter band out of range! ***\n");
        f4 = 0.4999 / dt;
        if (f3 >= f4) f3 = f4 / 1.1;
    }
    fftw_plan planF = NULL;
    // backward FFT: s ==> sf
    int ns;
    fftw_complex * s, * sf;
    FFTW_B(FFTW_ESTIMATE, seis_in, n, &s, &sf, &ns, &planF, 1);
    // make tapering and divide
    int nk     = ns / 2 + 1;
    double dom = 1. / dt / ns;
    fDiv(dom, nk, sf, freq, amp, pha, nf);
    TaperB(f1, f2, f3, f4, dom, nk, sf);

    // forward FFT: sf ==> s
    FFTW_F(planF, s, ns, seis_out, n);
    fftw_free(s);
    fftw_free(sf);

    // forming final result
    int k;
    float ftmp = 2. / ns;
    for (k = 0; k < n; k++) {
        // if( seis_in[k]==0 ) seis_out[k] = 0.;
        // else
        seis_out[k] *= ftmp;
    }
}
