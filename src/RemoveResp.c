#include "Param.h"

pthread_mutex_t evrlock;

// Transfer(): Remove instrument response
// CutRec(): cut the signal based on event origin time

void
FDivide(double f1, double f2, double f3, double f4, double dt, int n, float * seis_in, float * seis_out, double * freq,
  double * amp, double * pha, int nf);

int
Resampling(char * sacname, float ** sig2, SAC_HD * sd, int ithread);


int
CheckExistenceft(int ne, int ns)
{
    SAC_HD shd;

    if (!read_shd(sdb->rec[ne][ns].ft_fname, &shd) ) return 0;

    sdb->rec[ne][ns].n  = shd.npts;
    sdb->rec[ne][ns].dt = shd.delta;
    return 1;
}

void
UpdateRecCut(char * name, int nstart, int rec_b, int rec_e)
{
    int rec_b1[1000], rec_e1[1000], nrec1;
    FILE * frec;
    char recname1[200], recname[200];

    sprintf(recname1, "%s_rec1", name);
    sprintf(recname, "%s_rec", name);
    if (!read_rec(1, recname1, 0, rec_b1, rec_e1, &nrec1) ) {
        frec = fopen(recname, "w");
        fprintf(frec, "%d	%d\n", rec_b, rec_e);
        fclose(frec);
        return;
    }
    int i;
    frec = fopen(recname, "w");
    for (i = 0; i < nrec1; i++) {
        rec_b1[i] -= nstart;
        rec_e1[i] -= nstart;
        if (rec_e1[i] <= rec_b || rec_b1[i] >= rec_e) continue;
        if (rec_b1[i] < rec_b) rec_b1[i] = rec_b;
        if (rec_e1[i] > rec_e) rec_e1[i] = rec_e;
        fprintf(frec, "%d	%d\n", rec_b1[i], rec_e1[i]);
    }
    fclose(frec);
}

int
CutRec(int ne, int ns, float * sig1, SAC_HD shd1, int ithread)
{
    // cerr<<"CutRec: dt "<<sdb->rec[ne][ns].dt<<" rect0 "<<sdb->rec[ne][ns].t0<<" b "<<shd1.b<<" evt0 "<<sdb->ev[ne].t0<<endl;
    int n, nstt, nend;
    float t1b; // t2, t1e;

    float dt = sdb->rec[ne][ns].dt;

    n = (int) floor(tlen / dt + 0.5) + 1;
    // t2 = t1 + (n-1)*dt;

    //   if( read_sac(sdb->rec[ne][ns].ft_fname, sig1, shd1) == NULL ) {
    //      //cerr<<"Cannot open file "<<sdb->rec[ne][ns].ft_fname<<endl;
    //      sdb->rec[ne][ns].n = 0;
    //      return 0;
    //   }

    t1b  = sdb->rec[ne][ns].t0 - sdb->ev[ne].t0;
    t1b += shd1.b;
    // t1e = t1b + (sdb->rec[ne][ns].n-1)*dt;

    nstt = (int) floor((t1 - t1b) / dt + 0.5);
    nend = nstt + n;
    float tshift = fabs(nstt * dt - (t1 - t1b));
    // cerr<<"t1: "<<t1<<" t1b: "<<t1b<<" dt: "<<dt<<" nstt: "<<nstt<<endl;
    if (tshift > 1.e-3) {
        reports[ithread].tail += sprintf(reports[ithread].tail,
            "*** Warning: signal shifted by %fsec when cutting! ***", tshift);
        // cerr<<"Shifted "<<tshift<<endl;
    }

    int flag = 0;
    float * sig2 = NULL;
    int rec_b, rec_e;
    rec_b = 0;
    rec_e = n;
    if ( (nstt < 0) || (nend > shd1.npts) ) {
        int i = 0;
        if (nstt < 0) i += -nstt;
        if (nend > shd1.npts) i += nend - shd1.npts;
        if ((float) i / n > gapfrac) {
            sdb->rec[ne][ns].n = 0;
            fRemove(sdb->rec[ne][ns].ft_fname);
            return 0;
        }
        flag = 1;
        reports[ithread].tail += sprintf(reports[ithread].tail,
            "*** Warning: cut range isn't fully covered. zeros padded ***");
        // fprintf(stderr, "*** Warning: cut range isn't fully covered. zeros padded ***");
        sig2 = (float *) malloc(n * sizeof(float) );
        for (i = nstt; i < nend; i++) {
            if (i < 0 || i > shd1.npts - 1) sig2[i - nstt] = 0.;
            else sig2[i - nstt] = sig1[i];
        }
        if (nstt < 0) rec_b = -nstt;
        if (nend > shd1.npts) rec_e = shd1.npts - 1 - nstt;
    }
    UpdateRecCut(sdb->rec[ne][ns].ft_fname, nstt, rec_b, rec_e);

    shd1.npts = n;
    sdb->rec[ne][ns].n  = n;
    sdb->rec[ne][ns].dt = shd1.delta;
    // shd1.nzmsec += (int)floor((t1-t1b)*1000+0.5);
    shd1.b      = 0.;
    shd1.nzhour = 0;
    shd1.nzmin  = 0;
    shd1.nzsec  = 0;
    shd1.nzmsec = (int) floor(t1 * 1000 + 0.5);
    UpdateTime(&shd1);
    if (flag) {
        write_sac((const char *) (sdb->rec[ne][ns].ft_fname), sig2, &shd1);
        free(sig2);
    } else { write_sac((const char *) (sdb->rec[ne][ns].ft_fname), &(sig1[nstt]), &shd1); }
    free(sig1);

    return 1;
} /* CutRec */

/*
 * int TransferSac(int ne, int ns) {
 *
 * if(ne>=NEVENTS || ns>=sdb->nst) {
 *    fprintf(stderr, "*** Warning: event/station # out of range! ***");
 *    //sdb->rec[ne][ns].n = 0;
 *    return 0;
 * }
 *
 * FILE *ff;
 * float fl2 = 1./perh*0.7, fl1 = fl2*0.8, fl3 = 1./perl*1.3, fl4 = fl3*1.2;
 * //float fl1 = 1./100., fl2 = 1./80., fl3 = 1./4., fl4 = 1./3.;
 *
 * ff = fopen("sac_bp_respcor","w");
 * fprintf(ff, "%s << END\n", sacexe);
 * fprintf(ff, "r %s\n", sdb->rec[ne][ns].fname);
 * fprintf(ff,"rmean\n");
 * fprintf(ff,"rtrend\n");
 * fprintf(ff,"transfer from  EVALRESP FNAME  %s to vel freqlimits %f %f %f %f\n", sdb->rec[ne][ns].resp_fname, fl1, fl2, fl3, fl4 );
 * fprintf(ff,"w %s\n", sdb->rec[ne][ns].ft_fname);
 * fprintf(ff,"quit\n");
 * fprintf(ff,"END\n");
 * fclose(ff);
 *
 * System("sh sac_bp_respcor >& /dev/null");
 * if( fdel1 ) fRemove(sdb->rec[ne][ns].fname);
 * return 1;
 * }
 */

void
RTrend(float * sig, SAC_HD * shd)
{
    // fit a*x+b
    int i, npts = shd->npts;
    float X = 0., Y = 0., X2 = 0., Y2 = 0., XY = 0.;

    for (i = 0; i < npts; i++) {
        X  += i;
        Y  += sig[i];
        X2 += i * i;
        Y2 += sig[i] * sig[i];
        XY += i * sig[i];
    }
    float a = (npts * XY - X * Y) / (npts * X2 - X * X);
    float b = (-X * XY + X2 * Y) / (npts * X2 - X * X);
    // correct sig and DEPMEN
    float mean = 0., max = sig[0], min = sig[0];
    float shift = b;
    for (i = 0; i < npts; i++, shift += a) {
        sig[i] -= shift;
        mean   += sig[i];
        if (min > sig[i]) min = sig[i];
        else if (max < sig[i]) max = sig[i];
    }
    shd->depmin = min;
    shd->depmax = max;
    shd->depmen = mean / npts;
}

bool
CheckSPS(int ne, int ns, int ithread)
{
    char * fname = sdb->rec[ne][ns].fname;
    SAC_HD shd;

    if (!read_shd(fname, &shd) ) return false;

    if (fabs(shd.delta * sps - 1.) > 1.0e-3) {
        float * sig;
        if (!Resampling(fname, &sig, &shd, ithread) ) return false;

        write_sac(fname, sig, &shd);
        sdb->rec[ne][ns].dt = shd.delta;
        sdb->rec[ne][ns].n  = shd.npts;
    }
    return true;
}

int
TransferEvr(int ne, int ns, float ** sig, SAC_HD * sd, int ithread)
{
    // read in sac file
    *sig = NULL;
    if ( (read_sac(sdb->rec[ne][ns].fname, sig, sd)) == NULL) {
        reports[ithread].tail += sprintf(reports[ithread].tail, "*** Warning: cannot read sac file %s ***",
            sdb->rec[ne][ns].fname);
        // fprintf(stderr, "*** Warning: cannot read sac file %s ***", sdb->rec[ne][ns].fname);
        return 0;
    }
    // running evalresp
    int nf = 100;
    // std::cerr<<"channel from param: "<<ch<<std::endl;
    char buff[300], sta[8], ch[8], net[8];
    float f2 = 1. / perh * 0.7, f1 = f2 * 0.8, f3 = 1. / perl * 1.3, f4 = f3 * 1.2;
    float fmax = 0.499 / (sd->delta);
    if (f4 > fmax) {
        f4 = fmax;
        if (f3 >= fmax) f3 = fmax / 1.1;
        if (f2 >= f3) f2 = f3;
        // reports[ithread].tail += sprintf(reports[ithread].tail, "*** Warning: Filter band out of range. Filter reshaped! ***");
    }
    double pi = 4 * atan(1.0), pio180 = pi / 180.;
    sscanf(sd->kstnm, "%s", sta);
    sscanf(sd->kcmpnm, "%s", ch);
    // std::cerr<<"channel from sacheader: "<<ch<<std::endl;
    sscanf(sd->knetwk, "%s", net);
    sprintf(buff, "%s %s %s %4d %3d %f %f %d -f %s -v >& /dev/null", evrexe, sta, ch, sd->nzyear, sd->nzjday, f1, f4,
      nf, sdb->rec[ne][ns].resp_fname);
    pthread_mutex_lock(&evrlock); // lock
    system(buff);
    char nameam[50], nameph[50];
    // sprintf(nameam, "AMP.%s.%s.*.%s", net, sta, ch);
    // sprintf(nameph,"PHASE.%s.%s.*.%s", net, sta, ch);
    sprintf(nameam, "AMP.*.%s*.%s", sta, ch);
    sprintf(nameph, "PHASE.*.%s*.%s", sta, ch);

    // find am file
    FILE * fam = NULL, * fph = NULL;
    int nlist;
    char * list = List(".", nameam, 0, &nlist);
    if (nlist != 1) {
        cerr << "ERROR(TransferEvr): " << nlist << " AMP file(s) found with pattern " << nameam << endl;
        pthread_mutex_unlock(&evrlock);
        free(*sig);
        return 0; // exit(0);
    }
    sscanf(list, "%s", nameam);
    free(list);
    if ( (fam = fopen(nameam, "r")) == NULL) {
        cerr << "ERROR(TransferEvr): Cannot open file " << nameam << endl;
        pthread_mutex_unlock(&evrlock);
        free(*sig);
        return 0; // exit(0);
    }
    // find ph file
    list = List(".", nameph, 0, &nlist);
    if (nlist != 1) {
        cerr << "ERROR(TransferEvr): " << nlist << " PHASE file(s) found with pattern !" << nameph << endl;
        pthread_mutex_unlock(&evrlock);
        free(*sig);
        return 0; // exit(0);
    }
    sscanf(list, "%s", nameph);
    free(list);
    if ( (fph = fopen(nameph, "r")) == NULL) {
        cerr << "ERROR(TransferEvr): Cannot open file " << nameph << endl;
        pthread_mutex_unlock(&evrlock);
        free(*sig);
        return 0; // exit(0);
    }
    // read in am and ph data
    double freq[nf], dtmp, amp[nf], pha[nf];
    int i = 0;
    while (i < nf) {
        if (fgets(buff, 300, fam) == NULL) break;
        sscanf(buff, "%lf %lf", &freq[i], &amp[i]);
        if (fgets(buff, 300, fph) == NULL) break;
        sscanf(buff, "%lf %lf", &dtmp, &pha[i]);
        if (dtmp != freq[i]) {
            reports[ithread].tail += sprintf(reports[ithread].tail, "*** Warning: incompatible AMP - PHASE pair! ***");
            continue;
        }
        amp[i] *= 0.000000001;
        pha[i] *= pio180;
        i++;
    }
    fclose(fam);
    fclose(fph);
    fRemove(nameam);
    fRemove(nameph);
    pthread_mutex_unlock(&evrlock); // unlock
    // remove trend ( and mean )
    RTrend(*sig, sd);
    // run rmresponse
    FDivide(f1, f2, f3, f4, (double) (sd->delta), sd->npts, *sig, *sig, freq, amp, pha, nf);
    if (fdel1) fRemove(sdb->rec[ne][ns].fname);
    // cerr<<"Transfer: dt "<<sdb->rec[ne][ns].dt<<" rect0 "<<sdb->rec[ne][ns].t0<<" b "<<sd->b<<" evt0"<<sdb->ev[ne].t0<<endl;
    return 1;
} /* TransferEvr */

void *
RmRESPEntrance(void * tid)
{
    int ithread = *((int *) tid);
    int ne, ns, nst, flag;
    float * sig = NULL; // pointer to signal produced by TransferEvr
    SAC_HD shd;

    setvbuf(stdout, NULL, _IOLBF, 0);
    for (;;) {
        // get/update current event number
        pthread_mutex_lock(&cevlock);
        ne = currevn;
        currevn++;
        pthread_mutex_unlock(&cevlock);
        if (ne >= NEVENTS) break;
        if (strcmp(sdb->mo[imonth].seedf[ne], "0") == 0) continue;
        nst = 0;
        reports[ithread].tail += sprintf(reports[ithread].tail, "### Response removed for event %s from thread %d: ",
            sdb->ev[ne].name, ithread);
        // remove RESP for each station
        for (ns = 0; ns < sdb->nst; ns++) {
            if (fskip2 == 2 || fskip2 == 1) {
                flag = CheckExistenceft(ne, ns);
                if (fskip2 == 2) {
                    if (!flag) sdb->rec[ne][ns].n = 0;
                    continue;
                } else if (flag) { continue; }
            }
            if (sdb->rec[ne][ns].n <= 0) continue;
            if (!CheckSPS(ne, ns, ithread) ) continue;                 // check sampling rate and resample if necessary
            if (!TransferEvr(ne, ns, &sig, &shd, ithread) ) continue;  // sig allocated here
            if (!CutRec(ne, ns, sig, shd, ithread) ) continue;         // sig freed here
            if (nst % 20 == 0) reports[ithread].tail += sprintf(reports[ithread].tail, "\n   ");
            reports[ithread].tail += sprintf(reports[ithread].tail, "%s ", sdb->st[ns].name);
            nst++;
        }
        for (ns = 0; ns < sdb->nst; ns++) {
            if ( (sdb->rec[ne][ns].resp_fname) == NULL) continue;
            if (fdel1) fRemove(sdb->rec[ne][ns].resp_fname);
            delete [] sdb->rec[ne][ns].resp_fname;
        }
        reports[ithread].tail += sprintf(reports[ithread].tail, "\n   %d stations processed. ###\n", nst);
        cout << reports[ithread].head;
        reports[ithread].tail = reports[ithread].head;
    }
    fRemove("sac_bp_respcor");
    fRemove("RESP_tmp");

    pthread_exit(NULL);
} /* RmRESPEntrance */

void
RmRESP()
{
    int ithread;

    // initialize report arrays
    reports = (struct NOTE *) malloc(NTHRDS * sizeof(struct NOTE));
    for (ithread = 0; ithread < NTHRDS; ithread++) {
        // reports[ithread].head = (char *) malloc ( (sdb->nst+1) * 100 * sizeof(char) );
        reports[ithread].head = new char[(sdb->nst + 2) * 200];
        reports[ithread].tail = reports[ithread].head;
    }

    currevn = 0;
    // Thread ids and lock
    pthread_t tid[NTHRDS];
    pthread_mutex_init(&evrlock, NULL);

    // Create threads to produce event-station sac files
    int rc, targs[NTHRDS];
    for (ithread = 0; ithread < NTHRDS; ithread++) {
        targs[ithread] = ithread;
        rc = pthread_create(&tid[ithread], &attr_j, RmRESPEntrance, (void *) (&(targs[ithread])) );
        if (rc) {
            cerr << "ERROR(RmRESP): Thread creation failed!  ERR: " << strerror(rc) << endl;
            exit(0);
        }
    }
    // Wait for all threads to finish
    for (ithread = 0; ithread < NTHRDS; ithread++) pthread_join(tid[ithread], NULL);

    pthread_mutex_destroy(&evrlock);

    // free report arrays
    // for(ithread=0;ithread<NTHRDS;ithread++) free(reports[ithread].head);
    for (ithread = 0; ithread < NTHRDS; ithread++) delete [] reports[ithread].head;
    free(reports);
} /* RmRESP */
