#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <fftw3.h>
#include <math.h>
#include <string.h>
#include <pthread.h>

#define PI 3.14159265358979323846

// extern pthread_mutex_t fftlock;

void
FFTW_B(int type, float * seis, int n, fftw_complex ** in, fftw_complex ** out, int * nso, fftw_plan * planF, int Fflag);

void
TaperB(double f1, double f2, double f3, double f4, double dom, int nk, fftw_complex * sf);

int
Smooth2(double f1, double f4, double dom, int nk, fftw_complex * sf, float * seissm)
{
    return 1;
}

int
Smooth1(double f1, double f4, double dom, int nk, fftw_complex * sf, int num)
{
    int i, j, nb, ne, wb, we;
    double * amp = new double[nk];
    double wsum, wavg;

    for (i = 0; i < nk; i++) amp[i] = sqrt(pow(sf[i][0], 2.) + pow(sf[i][1], 2.) );

    if (num * 2 > nk - 1) num = (nk - 1) / 2;
    nb = (int) ceil(f1 / dom);
    ne = (int) floor(f4 / dom) + 1;
    for (i = nb, wsum = 0.; i <= nb + num; i++) wsum += amp[i];
    wb = nb;
    we = i;
    for (i = nb + 1; i <= nb + num; i++, we++) {
        wavg          = wsum / (we - wb);
        sf[i - 1][0] /= wavg;
        sf[i - 1][1] /= wavg;
        wsum         += amp[we];
    }
    for (j = 2 * num + 1; i < ne - num; i++, we++, wb++) {
        wavg          = wsum / j;
        sf[i - 1][0] /= wavg;
        sf[i - 1][1] /= wavg;
        wsum         += amp[we] - amp[wb];
    }
    for (; i <= ne; i++, wb++) {
        wavg          = wsum / (we - wb);
        sf[i - 1][0] /= wavg;
        sf[i - 1][1] /= wavg;
        wsum         -= amp[wb];
    }


    delete [] amp;
    amp = NULL;
    return 1;
} /* Smooth1 */

int
Whiten(double f1, double f2, double f3, double f4, double dt, int n, float hlen, float * seis_in, float * seissm,
  float ** outam, float ** outph, int * nkout, double * domout)
{
    if (f4 > 0.5 / dt) {
        fprintf(stdout, "*** Warning: whiten band upper end out of range! Corrected to spectrum end ***");
        f4 = 0.5 / dt;
        if (f3 >= f4) {
            f3 = f1 / f2 * f4;
            if (f3 < f2) f3 = f2;
        }
    }

    int ns;
    fftw_plan planF = NULL;
    fftw_complex * s, * sf;

    // backward FFT: s ==> sf
    FFTW_B(FFTW_ESTIMATE, seis_in, n, &s, &sf, &ns, &planF, 0); // 0 tells FFTW_B to not create planF
    int nk = ns / 2 + 1;
    *nkout = nk;
    double dom = 1. / dt / ns;
    *domout = dom;
    // free(s);

    // Smoothing
    if (hlen == -1.) {
        if (!Smooth2(f1, f4, dom, nk, sf, seissm) ) {
            fftw_free(sf);
            return 0;
        }
    } else if (hlen != 0.) {
        int num = (int) floor(hlen / dom + 0.5);
        Smooth1(f1, f4, dom, nk, sf, num);
    }

    // make tapering
    TaperB(f1, f2, f3, f4, dom, nk, sf);

    // final am & ph results
    *outam = new float[nk];
    *outph = new float[nk];
    int k;
    for (k = 0; k < nk; k++) {
        (*outam)[k] = sqrt(pow(sf[k][0], 2) + pow(sf[k][1], 2) );
        (*outph)[k] = atan2(sf[k][1], sf[k][0]);
    }
    fftw_free(s);
    fftw_free(sf);

    return 1;
} /* Whiten */
