// #include <errno.h>
#include "Param.h"


int currevn;
struct NOTE * reports;
pthread_attr_t attr_j;
pthread_mutex_t cevlock, fftlock, rdslock;// , syslock;

double
abs_time(int yy, int jday, int hh, int mm, int ss, int ms);

char *
wMove(const char * odir, const char * pattern, const char * tdir, int retlst, int * nfile);

void
SetName(int ne, int ns)
{
    sprintf(sdb->rec[ne][ns].fname, "%s/%s/%s.%s.%s.SAC", sdb->mo[imonth].name, sdb->ev[ne].name, sdb->ev[ne].name,
      sdb->st[ns].name, ch);
    sprintf(sdb->rec[ne][ns].ft_fname, "%s/%s/ft_%s.%s.%s.SAC", sdb->mo[imonth].name, sdb->ev[ne].name,
      sdb->ev[ne].name, sdb->st[ns].name, ch);
    sprintf(sdb->rec[ne][ns].chan, "%s", ch);
}

int
CheckExistence(int ne, int ns, int ithread)
{
    // check for sac file
    SAC_HD shd;

    if (!read_shd(sdb->rec[ne][ns].fname, &shd) ) return 0;

    sdb->rec[ne][ns].n  = shd.npts;
    sdb->rec[ne][ns].t0 = abs_time(shd.nzyear, shd.nzjday, shd.nzhour, shd.nzmin, shd.nzsec, shd.nzmsec);
    sdb->rec[ne][ns].dt = shd.delta;
    int nlist;
    // check for RESP file
    char dir[150], respname[150];
    sprintf(dir, "%s/%s", sdb->mo[imonth].name, sdb->ev[ne].name);
    sprintf(respname, "RESP.*.%s.*.%s", sdb->st[ns].name, ch);
    char * list = List(dir, respname, 0, &nlist);
    // exit(0);
    if (list == NULL) return 0;

    if (nlist > 1) reports[ithread].tail += sprintf(reports[ithread].tail,
            "*** Warning: more than one RESP files found ***");
    sdb->rec[ne][ns].resp_fname = new char[150];
    sscanf(list, "%s", sdb->rec[ne][ns].resp_fname);
    free(list);
    return 1;
}

int
Resampling(char * sacname, float ** sig2, SAC_HD * sd, int ithread)
{
    SAC_HD shd   = sac_null;
    float * sig1 = NULL;

    if (read_sac(sacname, &sig1, &shd) == NULL) return 0;

    float dt = 1. / sps;
    int iinc = (int) floor(dt / shd.delta + 0.5);
    // shd.delta = 1./(int)floor(1./shd.delta+0.5);

    /*
     * if(fabs(iinc*shd.delta-dt)>1e-7) {
     *    cerr<<"Error: "<<sps<<" is not a factor of "<<int(1/shd.delta)<<endl;
     *    exit(0);
     * }
     */
    if (iinc != 1) {
        double f1 = -1., f2 = -1., f3 = sps / 2.2, f4 = sps / 2.01;
        Filter(f1, f2, f3, f4, (double) shd.delta, shd.npts, sig1, sig1);
    }

    int i, j;
    int nptst = (int) floor((shd.npts - 1) * shd.delta * sps + 0.5) + 10;
    float nb;
    if ( (*sig2 = (float *) malloc(nptst * sizeof(float))) == NULL) perror("malloc sig2");
    long double fra1, fra2;
    nb = ceil((shd.nzmsec * 0.001 + shd.b) * sps);
    i  = (int) floor((nb * dt - shd.nzmsec * 0.001 - shd.b) / shd.delta);
    if (fabs(iinc * shd.delta - dt) < 1.e-7) { // sps is a factor of 1/delta
        fra2 = (nb * dt - i * shd.delta - shd.nzmsec * 0.001 - shd.b) / shd.delta;
        fra1 = 1. - fra2;
        if (fra2 == 0) {
            for (j = 0; i < shd.npts; j++) {
                (*sig2)[j] = sig1[i];
                i += iinc;
            }
        } else {
            for (j = 0; i < shd.npts - 1; j++) {
                (*sig2)[j] = sig1[i] * fra1 + sig1[i + 1] * fra2;
                i += iinc;
            }
        }
    } else { // sps isn't a factor, slower way
        reports[ithread].tail += sprintf(reports[ithread].tail,
            "*** Warning: sps isn't a factor of %d, watch out for rounding error! ***",
            (int) floor(1 / shd.delta + 0.5));
        long double ti, tj;
        iinc = (int) floor(dt / shd.delta);
        ti   = i * shd.delta + shd.nzmsec * 0.001 + shd.b;
        tj   = nb * dt;
        for (j = 0; i < shd.npts - 1; j++) {
            fra2       = tj - ti;
            (*sig2)[j] = sig1[i] + (sig1[i + 1] - sig1[i]) * fra2;
            tj        += dt;
            i  += iinc;
            ti += iinc * shd.delta;
            if (ti + shd.delta <= tj) { ti += shd.delta; i++; }// if(j%1000==0)cerr<<i<<" "<<ti<<j<<" "<<tj<<" "<<endl;}
        }
    }
    free(sig1);
    shd.nzmsec = (int) (nb * dt * 1000 + 0.5);
    shd.b      = 0.;
    if (shd.nzmsec >= 1000) UpdateTime(&shd);
    shd.delta = dt;
    shd.npts  = j;
    *sd       = shd;
    //   sprintf(sacname,"temp.sac");
    //   write_sac(sacname,*sig2,&shd);
    // exit(0);
    return 1;
} /* Resampling */

char *
Seed2Sac(int ne, int ns, char * nseed, int * nfile, int ithread)
{
    char tdir[100];

    sprintf(tdir, "./Working_Thread_%d", ithread);

    char fname[100], str[300];
    sprintf(fname, "%s/from_seed", tdir);
    FILE * ff = fopen(fname, "w");
    fprintf(ff, "%s <<END\n", rdsexe);
    fprintf(ff, "%s\n", nseed);
    fprintf(ff, "\n");                     /* out file */
    fprintf(ff, "\n");                     /* volume */
    fprintf(ff, "d\n");                    /* option */
    fprintf(ff, "\n");                     /* summary file */
    fprintf(ff, "%s\n", sdb->st[ns].name); /* station list */
    fprintf(ff, "%s\n", ch);               /* channel list */
    fprintf(ff, "\n");                     /* network list */
    fprintf(ff, "\n");                     /* Loc Ids */
    fprintf(ff, "1\n");                    /* out format */
    fprintf(ff, "N\n");                    // new version!!!!!!!!!!
    fprintf(ff, "N\n");                    /* Output poles & zeroes */
    fprintf(ff, "0\n");                    /* Check Reversal */
    fprintf(ff, "\n");                     /* Select Data Type */
    fprintf(ff, "\n");                     /* Start Time */
    fprintf(ff, "\n");                     /* End Time */
    fprintf(ff, "\n");                     /* Sample Buffer Length  */
    fprintf(ff, "Y\n");                    /* Extract Responses */
    fprintf(ff, "quit\n");
    fprintf(ff, "END\n");
    fclose(ff);

    // extract SAC&RESP and mv into thread working directory
    pthread_mutex_lock(&rdslock); // lock for rdseed and shell operations
    sprintf(str, "sh %s >& /dev/null", fname);
    system(str);

    /*---------- mv response file -----------*/
    sprintf(str, "RESP.*.%s.*.%s", sdb->st[ns].name, ch);
    int nlist   = 0;
    char * list = List(".", str, 0, &nlist); // list RESP files in the current depth
    if (list == NULL) {
        pthread_mutex_unlock(&rdslock);
        return NULL;
    }
    char list_name[150];
    int offset, curp = 0;
    while ( (sscanf(&list[curp], "%s%n", list_name, &offset)) == 1) {
        if (curp == 0) {
            sdb->rec[ne][ns].resp_fname = new char[150];
            sprintf(sdb->rec[ne][ns].resp_fname, "%s/%s/%s", sdb->mo[imonth].name, sdb->ev[ne].name, list_name);
            Move(list_name, sdb->rec[ne][ns].resp_fname);
        } else { fRemove(list_name); }
        curp += offset;
    }
    free(list);
    /*---------mv sac files and produce saclst---------*/
    // list SAC files
    sprintf(str, "*%s*%s*SAC", sdb->st[ns].name, ch);
    char * filelst = wMove(".", str, tdir, 1, nfile);
    // check if IRIS went nuts
    if (*nfile > 100) { // ignore this station if they did
        char list_name[150];
        int offset, curp = 0;
        while ( (sscanf(&filelst[curp], "%s%n", list_name, &offset)) == 1) {
            fRemove(list_name);
            curp += offset;
        }
        free(filelst);
        *nfile  = 0;
        filelst = NULL;
    }
    // cerr<<filelst<<endl;

    /*
     * list = List(".", str, 0, &nlist); //list SAC files in the current depth
     * if(list==NULL) {
     *    pthread_mutex_unlock(&rdslock);
     *    return NULL;
     * }
     * //move and rename;
     * int sleng = 0, lstlen = strlen(list) + nlist*(strlen(tdir)+3);
     * char *filelst = new char[lstlen];
     * curp = 0;
     * while( (sscanf(&list[curp], "%s%n", list_name, &offset)) == 1 ) {
     *    sprintf(str, "%s/%s", tdir, list_name);
     *    Move(list_name, str);
     *    sleng += sprintf(&filelst[sleng], "%s\n", str);
     *    curp += offset;
     * }
     * nfile = nlist;
     * free(list);
     */
    pthread_mutex_unlock(&rdslock); // unlock
    return filelst;
} /* Seed2Sac */

float
av_sig(float * sig, int i, int N, int nwin)
{
    int n1, n2, j, nav = 0;
    float av = 0.;

    if (nwin > N) nwin = N;
    n1 = i - nwin / 2;
    if (n1 < 0) n1 = 0;
    n2 = n1 + nwin - 1;
    if (n2 > N - 1) n2 = N - 1;
    n1 = n2 - nwin + 1;

    for (j = n1; j <= n2; j++)
        if (sig[j] < 1.e29) {
            av += sig[j];
            nav++;
        }

    if (nav < 1) av = 1.e30;
    else av = av / (float) nav;

    return av;
}

int
merge_sac(float * sig[], SAC_HD * sd, int nfile, int ne, int ns, int ithread)
{
    if (nfile == 1) {
        sdb->rec[ne][ns].n  = sd[0].npts;
        sdb->rec[ne][ns].t0 =
          abs_time(sd[0].nzyear, sd[0].nzjday, sd[0].nzhour, sd[0].nzmin, sd[0].nzsec, sd[0].nzmsec);
        sdb->rec[ne][ns].dt = sd[0].delta;
        write_sac((sdb->rec[ne][ns].fname), sig[0], &sd[0]);
        return 1;
    }

    int i, nb, j, jj, N, nfirst, Nholes;
    float * sig0;
    double t1[1000], t2[1000], T1 = 1.e25, T2 = -100.;
    SAC_HD s0;

    for (i = 0; i < nfile; i++) {
        t1[i] = abs_time(sd[i].nzyear, sd[i].nzjday, sd[i].nzhour, sd[i].nzmin, sd[i].nzsec, sd[i].nzmsec);
        t2[i] = t1[i] + (sd[i].npts - 1) * sd[i].delta;
        if (t1[i] < T1) {
            T1     = t1[i];
            nfirst = i;
        }
        if (t2[i] > T2) T2 = t2[i];
    }

    memcpy(&s0, &(sd[nfirst]), sizeof(SAC_HD) );
    double dt = (int) floor(s0.delta * 1e8 + 0.5) / 1e8, tshift;
    N       = (int) floor((T2 - T1) / s0.delta + 0.5) + 1;
    s0.npts = N;
    sdb->rec[ne][ns].n  = N;
    sdb->rec[ne][ns].t0 = T1;
    sdb->rec[ne][ns].dt = dt;

    sig0 = (float *) malloc(N * sizeof(float));
    for (j = 0; j < N; j++) sig0[j] = 1.e30;

    for (i = 0; i < nfile; i++) {
        if (fabs(sd[i].delta - s0.delta) > .0001) {
            reports[ithread].tail += sprintf(reports[ithread].tail, "*** Warning: sps mismatch! ***");
            continue;
        }
        nb     = (int) floor((t1[i] - T1) / dt + 0.5);
        tshift = fabs((sd[i].b - s0.b) + (nb * dt - (t1[i] - T1)));
        if (tshift > 1.e-3) {
            reports[ithread].tail += sprintf(reports[ithread].tail,
                "*** Warning: signal shifted by %fsec when merging! ***", tshift);
        }
        for (j = 0, jj = nb; j < sd[i].npts; j++, jj++) if (sig0[jj] > 1.e29) sig0[jj] = sig[i][j];
    }

    for (Nholes = 0, j = 0; j < N; j++) {
        if (sig0[j] > 1.e29) Nholes++;
    }
    if ( (float) Nholes / (float) N > gapfrac) {
        sdb->rec[ne][ns].n = -1;
        free(sig0);
        return 0;
    }

    int rec_b[1000], rec_e[1000];
    char recname[200];
    rec_b[0] = 0;
    for (i = 1, j = 0; i < N; i++) {
        if (sig0[i - 1] > 1.e29) { if (sig0[i] < 1.e29) rec_b[j] = i; } else if (sig0[i] > 1.e29) { rec_e[j++] = i; }
    }
    if (sig0[N - 1] < 1.e29) rec_e[j++] = N;
    sprintf(recname, "%s_rec1", sdb->rec[ne][ns].ft_fname);
    FILE * frec = fopen(recname, "w");
    for (i = 0; i < j; i++) fprintf(frec, "%d %d\n", rec_b[i], rec_e[i]);
    fclose(frec);

    float av;
    int npart;
    for (j = 0; j < N; j++) if (sig0[j] > 1.e29) {
            for (npart = 16; npart != 1; npart /= 2) {
                av = av_sig(sig0, j, N, N / npart);
                if (av < 1.e29) break;
            }
            if (npart == 1) av = 0.;
            sig0[j] = av;
        }

    write_sac((sdb->rec[ne][ns].fname), sig0, &s0);

    free(sig0);

    return 1;
} /* merge_sac */

int
MakeRecord(int ne, int ns, int ithread)
{
    // if ( sdb->rec[ne][ns].n > 0 ) return 0;
    int i, nfile;

    // create working directory for the current thread
    // char tdir[100];
    // sprintf(tdir, "Working_Thread_%d", ithread);

    char * filelst = Seed2Sac(ne, ns, sdb->mo[imonth].seedf[ne], &nfile, ithread);

    if (filelst == NULL) return 0;

    // read sacfile name from filelst
    char sacname[150];
    int offset, curp = 0;
    float * sigrspd[nfile];
    SAC_HD sd[nfile];
    i = 0;
    while ( (sscanf(&filelst[curp], "%s%n", sacname, &offset)) == 1) {
        curp += offset;
        Resampling(sacname, &sigrspd[i], &sd[i], ithread);
        i++;
        fRemove(sacname);
    }
    delete [] filelst;
    int rc = merge_sac(sigrspd, sd, nfile, ne, ns, ithread);
    // cleanup
    for (i = 0; i < nfile; i++) free(sigrspd[i]);
    if (!rc) return 0;

    return 1;
}

void *
ExtractSacEntrance(void * tid)
{
    int ithread = *((int *) tid);
    int iev, ist, nst, flag;

    char tdir[100];

    sprintf(tdir, "./Working_Thread_%d", ithread);
    MKDir(tdir);
    setvbuf(stdout, NULL, _IOLBF, 0);

    for (;;) {
        // get/update current event number
        pthread_mutex_lock(&cevlock);
        iev = currevn;
        currevn++;
        pthread_mutex_unlock(&cevlock);
        if (iev >= NEVENTS) break;
        // produce sac files
        for (ist = 0; ist < sdb->nst; ist++) {
            sdb->rec[iev][ist].n = 0;
            sdb->rec[iev][ist].resp_fname = NULL;
        }
        if (strcmp(sdb->mo[imonth].seedf[iev], "0") == 0) continue;
        nst = 0;
        reports[ithread].tail += sprintf(reports[ithread].tail, "### Sac files extracted for event %s from thread %d: ",
            sdb->ev[iev].name, ithread);
        for (ist = 0; ist < sdb->nst; ist++) {
            SetName(iev, ist);
            // sdb->rec[iev][ist].n = 0;
            if (fskip1 == 2 || fskip1 == 1) {
                flag = CheckExistence(iev, ist, ithread);
                if (fskip1 == 2) continue;
                else if (flag) continue;
            }
            // pthread_mutex_lock(&cevlock);
            if (MakeRecord(iev, ist, ithread) ) {
                if (nst % 20 == 0) reports[ithread].tail += sprintf(reports[ithread].tail, "\n   ");
                reports[ithread].tail += sprintf(reports[ithread].tail, "%s ", sdb->st[ist].name);
                nst++;
            }
            // pthread_mutex_unlock(&cevlock);
        }
        reports[ithread].tail += sprintf(reports[ithread].tail, "\n   %d stations processed. ###\n", nst);
        cout << reports[ithread].head;
        reports[ithread].tail = reports[ithread].head;
    }

    dRemove(tdir);
    pthread_exit(NULL);
} /* ExtractSacEntrance */

void
ExtractSac()
{
    int itmp;

    MKDir("old_sac_files");
    wMove(".", "*.SAC", "old_sac_files", 0, &itmp);
    wMove(".", "RESP.*", "old_sac_files", 0, &itmp);
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
    pthread_mutex_init(&rdslock, NULL);

    // Create threads to produce event-station sac files
    int rc, targs[NTHRDS];
    for (ithread = 0; ithread < NTHRDS; ithread++) {
        targs[ithread] = ithread;
        rc = pthread_create(&tid[ithread], &attr_j, ExtractSacEntrance, (void *) (&(targs[ithread])) );
        if (rc) {
            cerr << "ERROR(ExtractSac): Thread creation failed!  ERR: " << strerror(rc) << endl;
            exit(0);
        }
    }
    // Wait for all threads to finish
    for (ithread = 0; ithread < NTHRDS; ithread++) pthread_join(tid[ithread], NULL);

    pthread_mutex_destroy(&rdslock);

    // free report arrays
    for (ithread = 0; ithread < NTHRDS; ithread++) free(reports[ithread].head);
    free(reports);

    wMove("old_sac_files", "*.SAC", ".", 0, &itmp);
    wMove("old_sac_files", "RESP.*", ".", 0, &itmp);
    fRemove("old_sac_files");
    // pthread_exit(NULL);
} /* ExtractSac */
