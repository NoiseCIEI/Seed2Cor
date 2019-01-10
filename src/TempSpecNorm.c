#include "Param.h"
#include <vector>
#include <algorithm>

int
Whiten(double f1, double f2, double f3, double f4, double dt, int n, float hlen, float * seis_in, float * seissm,
  float ** outam, float ** outph, int * nk, double * dom);

void
UpdateRec(char * name, int * rec_b, int * rec_e, int nrec, int ithread)
{
    int rec_b1[1000], rec_e1[1000], nrec1;
    FILE * frec;
    int irec, irec1;

    if (!read_rec(1, name, 0, rec_b1, rec_e1, &nrec1) ) {
        reports[ithread].tail += sprintf(reports[ithread].tail, "*** Warning: cannot open record file %s ***", name);
        frec = fopen(name, "w");
        for (irec = 0; irec < nrec; irec++)
            fprintf(frec, "%d %d\n", rec_b[irec], rec_e[irec]);
        fclose(frec);
    }
    int recB, recE;
    char name2[100];
    sprintf(name2, "%s2", name);
    frec = fopen(name2, "w");
    for (irec = 0; irec < nrec; irec++)
        for (irec1 = 0; irec1 < nrec1; irec1++) {
            if (rec_b[irec] >= rec_e1[irec1]) continue;
            if (rec_e[irec] <= rec_b1[irec1]) break;
            recB = max(rec_b[irec], rec_b1[irec1]);
            recE = min(rec_e[irec], rec_e1[irec1]);
            fprintf(frec, "%d %d\n", recB, recE);
        }
    fclose(frec);
}

void
OneBit(float * sig, SAC_HD * shd)
{
    int i;

    for (i = 0; i < shd->npts; i++) {
        if (sig[i] > 0.) sig[i] = 1.;
        else if (sig[i] < 0.) sig[i] = -1.;
    }
}

void
RunAvgNorm(float * sig, SAC_HD * shd, float * sigw)
{
    int i, j, wb, we, n = shd->npts;
    float wsum, dt = shd->delta;
    int half_l = (int) floor(timehlen / dt + 0.5);

    if (half_l * 2 > n - 1) half_l = (n - 1) / 2;
    for (i = 0, wsum = 0.; i <= half_l; i++) wsum += fabs(sigw[i]);
    wb = 0;
    we = i;
    for (i = 1; i <= half_l; i++, we++) {
        if (wsum > 1.e-15) sig[i - 1] *= ((double) we / wsum);
        wsum += fabs(sigw[we]);
    }
    for (j = we; i < n - half_l; i++, wb++, we++) {
        // if( i>80000/dt && i<82000/dt ) std::cerr<<(i-1)*dt<<" "<<sig[i-1]<<" "<<wsum<<" / "<<j<<std::endl;
        if (wsum > 1.e-15) sig[i - 1] *= ((double) j / wsum);
        wsum += ( fabs(sigw[we]) - fabs(sigw[wb]) );
    }
    for (; i < n; i++, wb++) {
        if (wsum > 1.e-15) sig[i - 1] *= ((double) (we - wb) / wsum);
        wsum -= fabs(sigw[wb]);
    }
    if (wsum > 1.e-15) sig[n - 1] *= ((double) (we - wb) / wsum);
}

void
RunAvg(float * sig, SAC_HD * shd)
{
    int n = shd->npts;
    float * sigw, dt = shd->delta;

    sigw = new float[n];
    double f2 = 1. / Eperh, f1 = f2 * 0.6, f3 = 1. / Eperl, f4 = f3 * 1.4;
    //   double f1 = 1./60., f2 = 1./50., f3 = 1./10., f4 = 1./7.5;
    if (Eperl == -1) memcpy(sigw, sig, n * sizeof(float));
    else Filter(f1, f2, f3, f4, (double) dt, n, sig, sigw);

    // normalize sig by a running average of sigw
    RunAvgNorm(sig, shd, sigw);

    delete [] sigw;
    sigw = NULL;
    // write_sac("temp1.sac", sig, shd );
    // exit(0);
}

inline void
TaperCos(float * sig, int winb, int wine)
{
    int step = 1;

    if (wine < winb) step = -1;
    float pi = 4.0 * atan(1.0), dtmp = pi / (wine - winb);
    for (int i = winb; i != wine; i += step) sig[i] *= ( 0.5 * ( cos( (i - winb) * dtmp) + 1 ) );
}

int
EqkCut(float * sig, SAC_HD * shd, char * recname, int ithread)
{
    int i, n = shd->npts, ninc = n / 1000, npole = 0;

    for (i = 0; i < n; i += ninc) if (fabs(sig[i]) < 1.e-20) npole++;
    if (npole > 600) {
        reports[ithread].tail += sprintf(reports[ithread].tail, "*** Warning: Signal time length not long enough. ***");
        return 0;
    }

    float * sigw = new float[n];
    double f2 = 1. / Eperh, f1 = f2 * 0.8, f3 = 1. / Eperl, f4 = f3 * 1.2;
    double dt = (double) (shd->delta);
    if (Eperl == -1) memcpy(sigw, sig, shd->npts * sizeof(float));
    else Filter(f1, f2, f3, f4, dt, n, sig, sigw);

    /* compute noise level */
    int ii, is;
    int s1k   = (int) floor(1000. / dt + 0.5); // npts of a 1000 sec window
    int nos1k = (int) (n / s1k);

    // for (each of) the ith 1000 sec window, search for maximum amplitude and store into win_max[i]
    double win_max[nos1k];
    memset(win_max, 0, nos1k * sizeof(double));
    for (i = 0; i < n; i++) sigw[i] = fabs(sigw[i]);
    for (ii = 0, is = 0; is < nos1k; is++) {
        for (; ii < (is + 1) * s1k; ii++)
            if (win_max[is] < sigw[ii]) win_max[is] = sigw[ii];
    }
    for (; ii < n; ii++)
        if (win_max[is] < sigw[ii]) win_max[is] = sigw[ii];

    // sort win_max
    std::vector < double > win_max_sorted(win_max, win_max + nos1k);
    std::sort(win_max_sorted.begin(), win_max_sorted.end() );

    // and define max noise level as 2 times the smallest 10 average max
    double noisemax = 0.;
    std::vector < double > ::iterator iter, itermin;
    for (itermin = win_max_sorted.begin(); itermin < win_max_sorted.end(); itermin++) if (*itermin > 1.e-20) break;
    if (itermin < win_max_sorted.end() ) {
        for (iter = itermin; iter < win_max_sorted.end() && iter < itermin + 30; iter++) noisemax += *iter;
        noisemax *= 2. / (iter - itermin);
    }

    /* compute noise average and noise std between windows */
    double window_avg = 0.;
    for (iter = itermin; iter < win_max_sorted.end() && *iter < noisemax; iter++) window_avg += *iter;
    ii          = iter - itermin;
    window_avg /= ii;
    double window_std = 0., dtmp;
    for (iter = itermin; iter < itermin + ii; iter++) {
        dtmp        = window_avg - *iter;
        window_std += dtmp * dtmp;
    }
    window_std = sqrt(window_std / (ii - 1));

    /* mark windows with a max amp > window_avg+2.0*window_std to be 'zero' */
    dtmp = window_avg + 2.0 * window_std;
    short keep[nos1k];
    for (i = 0; i < nos1k; i++) keep[i] = win_max[i] > dtmp ? 0 : 1;

    /* and zero out invalidated windows */
    for (i = 0; i < nos1k; i++)
        if (keep[i] == 0) for (ii = i * s1k; ii < (i + 1) * s1k; ii++) sig[ii] = 0.;

    /* locate contigious valid windows and apply a cosine taper */
    int rec_b[1000], rec_e[1000], rec_i = 0;
    int winlen_min = (int) ceil(2500. / dt);

    /* locate all rec_begin and rec_end pairs
     * zero out the windows that are shorter than winlen_min */
    rec_b[0] = 0;
    for (i = 1; i < nos1k; i++) {
        if (keep[i] - keep[i - 1] == 1) {
            rec_b[rec_i] = i * s1k; // a new window begins
        } else if (keep[i] - keep[i - 1] ==
          -1) // the current window ends
        {
            rec_e[rec_i] = i * s1k;
            /* invalidate the current window if it is shorter than winlen_min */
            if ((rec_e[rec_i] - rec_b[rec_i]) < winlen_min)
                for (ii = rec_b[rec_i]; ii < rec_e[rec_i]; ii++) sig[ii] = 0.;
            else rec_i++;
        }
    }
    /* mark the last rec_end and check its window length */
    if (keep[nos1k - 1] == 1) {
        if ( (n - rec_b[rec_i]) < winlen_min * 0.6) {
            for (ii = rec_b[rec_i]; ii < n; ii++) sig[ii] = 0;
        } else { rec_e[rec_i] = n; rec_i++; }
    }

    /* check if there's enough data left */
    for (ii = 0, i = 0; i < rec_i; i++) ii += rec_e[i] - rec_b[i];
    if (ii < 0.2 * n) {
        reports[ithread].tail += sprintf(reports[ithread].tail,
            "*** Warning: Time length not enough after removing earthquakes. Skipped. ***");
        return 0;
    }

    /* taper 300 sec of data on each side of each window just to be safe */
    int taperhl = (int) ceil(150. / dt);
    for (i = 0; i < rec_i; i++) {
        if (rec_b[i] != 0) {
            TaperCos(sig, rec_b[i] + 2 * taperhl, rec_b[i] - 1);
            rec_b[i] += taperhl;
        }
        if (rec_e[i] != n) {
            TaperCos(sig, rec_e[i] - 2 * taperhl, rec_e[i] + 1);
            rec_e[i] -= taperhl;
        }
    }

    /* produce a new rec file named ft_name_rec2 */
    UpdateRec(recname, rec_b, rec_e, rec_i, ithread);
    /* norm by running average if required */
    if (tnorm_flag == 4) RunAvgNorm(sig, shd, sigw);

    delete [] sigw;
    sigw = NULL;
    return 1;
} /* EqkCut */

int
TemperalNorm(char * fname, float ** sig, SAC_HD * shd, int ithread)
{
    *sig = NULL;
    if (read_sac(fname, sig, shd) == NULL) {
        reports[ithread].tail += sprintf(reports[ithread].tail, "*** Warning: Cannot open file %s ***", fname);
        return 0;
    }
    char recname[300];
    sprintf(recname, "%s_rec", fname);
    if (tnorm_flag == 1) { OneBit(*sig, shd); } else if (tnorm_flag == 2) {
        RunAvg(*sig, shd);
    } else if (tnorm_flag == 3 ||
      tnorm_flag ==
      4) { if (!EqkCut(*sig, shd, recname, ithread)) return 0; } else if (tnorm_flag != 0) {
        cerr << "ERROR(TemperalNorm): Undefined normalization method!" << endl;
        exit(0);
    }

    /*
     * //check norm results
     * std::string outname( (shd->kstnm) );
     * outname.replace(outname.begin()+4, outname.end(), "_check_eqkcut.SAC");
     * write_sac(outname.c_str(), *sig, shd);
     * std::cerr<<outname<<std::endl;
     * exit(0);
     */

    return 1;
}

int
SpectralNorm(char * fname, float * sig, SAC_HD shd, int ithread)
{
    int nk, flag_whiten;
    double dom;
    double dt = (double) shd.delta;
    int n     = shd.npts;
    // double f1 = 1./100., f2 = 1./80., f3 = 1./4., f4 = 1./3.;
    double f2 = 1. / perh, f1 = f2 * 0.8, f3 = 1. / perl, f4 = f3 * 1.2;
    float * sigw = NULL;

    // if(strcmp(fwname,"0") != 0) {
    if (frechlen == -1.) {
        SAC_HD shdw;
        if (read_sac(fwname, &sigw, &shdw) == NULL) {
            cerr << "ERROR(SpectralNorm): Cannot open file " << fwname << endl;
            free(sig);
            exit(0);
        }
        if (shdw.npts != shd.npts || fabs(shdw.delta - shd.delta) > 1.e-3) {
            cerr << "ERROR(SpectralNorm): Smoothing spectrum is incompatibale with input signals!" << endl;
            exit(0);
        }
    }
    // float seis_out[n];
    // int nf = (int)(log((double)n)/log(2.))+1;
    // if(nf<13) nf = 13;
    // nf = (int)pow(2,nf);
    float * outam = NULL, * outph = NULL; // [nf];
    // memset (seis_outamp,0,nf*sizeof(float));
    // memset (seis_outph,0,nf*sizeof(float));
    // pthread_mutex_lock(&fftlock);
    // whiten_(&f1,&f2,&f3,&f4,&npow,&dt,&n,&frechlen,sig,sigw,seis_out,outam,outph,&ns,&dom,&flag_whiten);
    flag_whiten = Whiten(f1, f2, f3, f4, dt, n, frechlen, sig, sigw, &outam, &outph, &nk, &dom);
    // pthread_mutex_unlock(&fftlock);
    if (flag_whiten == 0) {
        reports[ithread].tail += sprintf(reports[ithread].tail, "*** Warning: Skipped due to probamatic spectrum. ***");
    } else {
        char nameamp[200], nameph[200];
        sprintf(nameamp, "%s.am", fname);
        sprintf(nameph, "%s.ph", fname);
        shd.npts   = nk;
        shd.delta  = dom;
        shd.b      = 0;
        shd.iftype = IXY;
        write_sac(nameamp, outam, &shd);
        write_sac(nameph, outph, &shd);
    }

    delete [] outam;
    outam = NULL;
    delete [] outph;
    outph = NULL;
    free(sig);
    free(sigw);

    return flag_whiten;
} /* SpectralNorm */

void *
TSNormEntrance(void * tid)
{
    int ithread = *((int *) tid);
    int iev, ist, nst;
    char amname[200], phname[200];
    float * sig;
    SAC_HD shd;

    setvbuf(stdout, NULL, _IOLBF, 0);
    for (;;) {
        // get/update current event number
        pthread_mutex_lock(&cevlock);
        iev = currevn;
        currevn++;
        pthread_mutex_unlock(&cevlock);
        if (iev >= NEVENTS) break;
        // Normalizations
        if (strcmp(sdb->mo[imonth].seedf[iev], "0") == 0) continue;
        nst = 0;
        reports[ithread].tail += sprintf(reports[ithread].tail,
            "### Am & ph records produced for event %s from thread %d: ",
            sdb->ev[iev].name, ithread);
        for (ist = 0; ist < sdb->nst; ist++) {
            sprintf(amname, "%s.am", sdb->rec[iev][ist].ft_fname);
            sprintf(phname, "%s.ph", sdb->rec[iev][ist].ft_fname);
            if (fskip3 == 2 || (fskip3 == 1 && access(amname, R_OK) != -1 && access(phname, R_OK) != -1) ) continue;
            if (sdb->rec[iev][ist].n <= 0) continue;
            /* report */
            if (nst % 20 == 0) reports[ithread].tail += sprintf(reports[ithread].tail, "\n   ");
            reports[ithread].tail += sprintf(reports[ithread].tail, "%s ", sdb->st[ist].name);
            /* normalizing */
            if (!TemperalNorm(sdb->rec[iev][ist].ft_fname, &sig, &shd, ithread) )
            { sdb->rec[iev][ist].n = 0; continue; }
            if (!SpectralNorm(sdb->rec[iev][ist].ft_fname, sig, shd, ithread) )
            { sdb->rec[iev][ist].n = 0; continue; }
            nst++;
        }
        reports[ithread].tail += sprintf(reports[ithread].tail, "\n   %d stations processed. ###\n", nst);
        cout << reports[ithread].head;
        reports[ithread].tail = reports[ithread].head;
    }

    pthread_exit(NULL);
} /* TSNormEntrance */

void
TempSpecNorm()
{
    int ithread;

    // initialize report arrays
    reports = (struct NOTE *) malloc(NTHRDS * sizeof(struct NOTE));
    for (ithread = 0; ithread < NTHRDS; ithread++) {
        reports[ithread].head = (char *) malloc( (sdb->nst + 1) * 100 * sizeof(char) );
        reports[ithread].tail = reports[ithread].head;
    }

    currevn = 0;
    // Thread ids and lock
    pthread_t tid[NTHRDS];

    // Create threads to produce normalized am/ph files
    int rc, targs[NTHRDS];
    for (ithread = 0; ithread < NTHRDS; ithread++) {
        targs[ithread] = ithread;
        rc = pthread_create(&tid[ithread], &attr_j, TSNormEntrance, (void *) (&(targs[ithread])) );
        if (rc) {
            cerr << " Thread creation failed!  ERR: " << strerror(rc) << endl;
            exit(0);
        }
    }
    // Wait for all threads to finish
    for (ithread = 0; ithread < NTHRDS; ithread++) pthread_join(tid[ithread], NULL);

    // pthread_mutex_destroy(&rdslock);

    // free report arrays
    for (ithread = 0; ithread < NTHRDS; ithread++) free(reports[ithread].head);
    free(reports);
} /* TempSpecNorm */
