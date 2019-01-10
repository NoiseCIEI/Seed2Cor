// #include <sys/sysinfo.h>
#include "Param.h"
#include "DisAzi.h"

#define max(a, b) ( ((a) > (b)) ? (a) : (b) )
#define min(a, b) ( ((a) < (b)) ? (a) : (b) )

/* Global parameters */
SAC_DB * sdbnew;
int cortype;
int nrec;
int i1, i2, is1, is2, ig1, iread;
float * sigcor;
SAC_HD shdcor;
short ** dnum, *** dflag;
struct mstarec * stag, * starec2;
pthread_mutex_t addlock;

// station record structure
struct starec {
    char    name[100];
    SAC_HD  shd;
    float * amp, * pha, dtrec;
    int     nr, rec[2][1000];
};

struct mstarec {
    struct starec dayrec[NEVENTS];
};

/* Finction prorotypes */

void
Copy(char * oldname, char * newname);

void
IFFT(int nlen, float * amp, float * pha, float * seis_out, int * nsig);

// int calc_dist(double lati1, double long1, double lati2, double long2, double *dist);

int
StaValid(int ista)
{
    char fname[100];

    for (int iev = 0; iev < NEVENTS; iev++) {
        sprintf(fname, "%s.am", sdbnew->rec[iev][ista].ft_fname);
        if (access(fname, R_OK) == -1) continue;
        sprintf(fname, "%s.ph", sdbnew->rec[iev][ista].ft_fname);
        if (access(fname, R_OK) != -1) return 1;
    }
    return 0;
}

void
UpdateStaLst()
{
    int i, ival, iev;

    for (i = 0, ival = 0; i < sdbnew->nst; i++) {
        if (StaValid(i) ) {
            if (ival < i) {
                sdbnew->st[ival] = sdbnew->st[i];
                for (iev = 0; iev < NEVENTS; iev++) sdbnew->rec[iev][ival] = sdbnew->rec[iev][i];
            }
            ival++;
        }
    }
    sdbnew->nst = ival;
}

void
ReadGroup(int size, int ig, struct mstarec * stag)
{
    // float *am=NULL, *ph=NULL;
    char am_sac[200], ph_sac[200], recname[200];
    int ftype = 1;

    if (size < 0) { ftype = 0; size = -size; }
    int i, ist = ig * size, nst, iev;
    SAC_HD shdam, shdph;
    if (ftype) iread = 0;  // fprintf(stdout, "### Reading in records for station group %d: ", ig+1);
    else fprintf(stdout, "\n  *Reading in records for the %d(/%d)th station: ", ig + 1, sdbnew->nst);
    for (i = 0, nst = 0; i < size; i++, ist++) {
        // cerr<<"ReadGroup:  reading ist="<<ist<<endl;
        if (ftype) iread = i;
        else fprintf(stdout, "\n   %s: ", sdbnew->st[ist].name);
        for (iev = 0; iev < NEVENTS; iev++) {
            stag[i].dayrec[iev].nr = 0;
            if (ist >= sdbnew->nst) continue;
            if (!ftype) fprintf(stdout, " %d", iev + 1);
            if (sdbnew->rec[iev][ist].n <= 0) {
                if (!ftype) fprintf(stdout, "*** No record ***");
                continue;
            }
            sprintf(stag[i].dayrec[iev].name, "%s", sdbnew->rec[iev][ist].ft_fname);
            sprintf(am_sac, "%s.am", sdbnew->rec[iev][ist].ft_fname);
            if (read_sac(am_sac, &(stag[i].dayrec[iev].amp), &shdam) == NULL) {
                if (!ftype) fprintf(stdout, "*** Warning: Cannot open file %s ***", am_sac);
                continue;
            }
            sprintf(ph_sac, "%s.ph", sdbnew->rec[iev][ist].ft_fname);
            if (read_sac(ph_sac, &(stag[i].dayrec[iev].pha), &shdph) == NULL) {
                if (!ftype) fprintf(stdout, "*** Warning: Cannot open file %s ***", ph_sac);
                continue;
            }
            if (shdam.npts != shdph.npts || shdam.delta != shdph.delta) {
                if (!ftype) fprintf(stdout, "*** Warning: ph-am mismatch ***");
                continue;
            }
            stag[i].dayrec[iev].shd   = shdam;
            stag[i].dayrec[iev].dtrec = sdbnew->rec[iev][ist].dt;
            if (tnorm_flag == 3 || tnorm_flag == 4) sprintf(recname, "%s_rec2", sdbnew->rec[iev][ist].ft_fname);
            else sprintf(recname, "%s_rec", sdbnew->rec[iev][ist].ft_fname);
            if (!read_rec(ftlen, recname, sdbnew->rec[iev][ist].n, stag[i].dayrec[iev].rec[0],
              stag[i].dayrec[iev].rec[1], &(stag[i].dayrec[iev].nr)) )
            {
                if (!ftype) fprintf(stdout, "*** Warning: Cannot open file %s ***", recname);
                continue;
            }

            /*
             *  for(j=0; j<shdam.npts; j++) {
             *     stag[i].dayrec[iev].amp[j] = am[j];//am[j]*cos(ph[j]);
             *     stag[i].dayrec[iev].pha[j] = ph[j];//am[j]*sin(ph[j]);
             *  }
             *  free(am);
             *  free(ph);
             */
        }
        nst++;
    }
    if (ftype) iread = i;  // fprintf(stdout, "\n");
    else fprintf(stdout, "\n   %d station(s) read in. ###", nst);
} /* ReadGroup */

void *
ReadGroupEntrance(void * tid)
{
    int size = *((int *) tid);

    iread = 0;
    ReadGroup(size, ig1, stag);
    pthread_exit(NULL);
}

void
DeleteGroup(int size)
{
    int i, iev;
    char fname[150];

    cout << "### Deleting am&ph files for the current station group... ###" << endl;
    for (i = 0; i < size; i++) {
        for (iev = 0; iev < NEVENTS; iev++) {
            if (stag[i].dayrec[iev].nr == 0) continue;
            sprintf(fname, "%s.am", stag[i].dayrec[iev].name);
            fRemove(fname);
            sprintf(fname, "%s.ph", stag[i].dayrec[iev].name);
            fRemove(fname);
        }
    }
}

int
CalcRecCor(struct starec stag1, struct starec stag2, float * cor_rec, int lag, float dt)
{
    int t, irec1, irec2;
    int recB, recE;

    for (t = 0; t <= lag; t++) {
        cor_rec[lag + t] = 0;
        cor_rec[lag - t] = 0;
        for (irec1 = 0; irec1 < stag1.nr; irec1++) {
            for (irec2 = 0; irec2 < stag2.nr; irec2++) {
                if (stag1.rec[0][irec1] >= stag2.rec[1][irec2] - t) continue;
                if (stag1.rec[1][irec1] <= stag2.rec[0][irec2] - t) break;
                recB = max(stag1.rec[0][irec1], stag2.rec[0][irec2] - t);
                recE = min(stag1.rec[1][irec1], stag2.rec[1][irec2] - t);
                cor_rec[lag + t] += recE - recB;
            }
            for (irec2 = 0; irec2 < stag2.nr; irec2++) {
                if (stag1.rec[0][irec1] >= stag2.rec[1][irec2] + t) continue;
                if (stag1.rec[1][irec1] <= stag2.rec[0][irec2] + t) break;
                recB = max(stag1.rec[0][irec1], stag2.rec[0][irec2] + t);
                recE = min(stag1.rec[1][irec1], stag2.rec[1][irec2] + t);
                cor_rec[lag - t] += recE - recB;
            }
        }
    }
    cor_rec[lag] /= 2;

    if (cor_rec[0] < mintlen / dt || cor_rec[lag * 2] < mintlen / dt) {
        fprintf(stdout, "*** cor time less than %d sec. Skipped! ***", mintlen);
        return 0;
    }
    return 1;
}

int
CheckPrecNoise(int lag, float dt, float * cor, float dist)
{
    float noise = 0., cor_pre = 0.;
    int i, ndis, nb;

    ndis = (int) floor((dist / 0.8 + 50.) / dt + 0.5);
    if (ndis > lag) return 1;

    nb = lag * 4 / 5;
    if (ndis > nb) nb = ndis;
    for (i = lag + nb; i < 2 * lag; i++) noise += pow(cor[i], 2);
    noise = sqrt(noise / (lag - nb));

    ndis = (int) floor((dist / 4.5 - 50.) / dt + 0.5);
    if (ndis < 10. / dt) return 1;

    nb = int(100. / dt);
    if (ndis < nb) nb = ndis;
    for (i = lag - nb; i <= lag + nb; i++) cor_pre += pow(cor[i], 2);
    cor_pre = sqrt(cor_pre / (2. * nb));
    if (cor_pre > noise * 5) {
        fprintf(stdout, "*** Warning: Large precursor signal. Skipped! ***");
        return 0;
    }
    return 1;
}

int
ComprRec(struct starec * stag1, struct starec * stag2)
{
    if (stag1->dtrec != stag2->dtrec) return 0;

    if (stag1->shd.delta != stag2->shd.delta) return 0;

    if (stag1->shd.nzyear != stag2->shd.nzyear) return 0;

    if (stag1->shd.nzjday != stag2->shd.nzjday) return 0;

    if (stag1->shd.nzhour != stag2->shd.nzhour) return 0;

    if (stag1->shd.nzmin != stag2->shd.nzmin) return 0;

    if (stag1->shd.nzsec != stag2->shd.nzsec) return 0;

    if (abs(stag1->shd.nzmsec - stag2->shd.nzmsec) > 1) return 0;

    return 1;
}

int
InitCor(int is1, int is2, struct starec * dayrec2)
{
    int iev;

    for (iev = 0; iev < NEVENTS && dayrec2[iev].nr == 0; iev++);
    if (iev == NEVENTS) return 0;

    float dt = dayrec2[iev].dtrec;
    int i, lag = (int) floor(lagtime / dt + 0.5);
    sigcor = new float[2 * lag + 1];
    for (i = 0; i < 2 * lag + 1; i++) sigcor[i] = 0.;

    shdcor       = dayrec2[iev].shd;
    shdcor.delta = dt;
    strcpy(shdcor.kevnm, sdbnew->st[is1].name);
    shdcor.evla = sdbnew->st[is1].lat;
    shdcor.evlo = sdbnew->st[is1].lon;
    shdcor.stla = sdbnew->st[is2].lat;
    shdcor.stlo = sdbnew->st[is2].lon;
    // calc_dist((double)(shdcor.evla), (double)(shdcor.evlo), (double)(shdcor.stla), (double)(shdcor.stlo), &distmp);
    shdcor.dist   = Path < float > (shdcor.evlo, shdcor.evla, shdcor.stlo, shdcor.stla).Dist();
    shdcor.npts   = 2 * lag + 1;
    shdcor.e      = lag * dt;
    shdcor.b      = -shdcor.e;
    shdcor.user0  = 0;
    shdcor.nzjday = 1;
    return 1;
}

int
DoCor(struct starec * stag1, struct starec * stag2, short * dnum, short * dflag, char * outname)
{
    // fprintf(stderr, "%s-%s ", sdbnew->st[is1].name, sdbnew->st[is2].name);
    float dt1 = stag1->dtrec;

    if (!ComprRec(stag1, stag2) ) {
        fprintf(stdout, "*** Warning: Incompatible record! Skipped. ***");
        return 0;
    }
    int i, ns;
    int len = stag1->shd.npts, lag = (int) floor(lagtime / dt1 + 0.5);

    float * atmp, * ptmp, * seis_out;
    atmp     = new float[len];
    ptmp     = new float[len];
    seis_out = new float[2 * len];
    if (cortype == 0) {
        for (i = 0; i < len; i++) atmp[i] = stag1->amp[i] * stag2->amp[i];
    } else if (cortype == 1) {
        float ampmin = 0.;
        for (i = 0; i < len; i++) ampmin += stag2->amp[i];
        ampmin /= len * 20.;
        for (i = 0; i < len; i++) {
            if (stag2->amp[i] > ampmin) atmp[i] = stag1->amp[i] / stag2->amp[i];
            else atmp[i] = stag1->amp[i] / ampmin;
        }
    } else if (cortype == 2) {
        float ampmin = 0.;
        for (i = 0; i < len; i++) ampmin += stag1->amp[i];
        ampmin /= len * 20.;
        for (i = 0; i < len; i++) {
            if (stag1->amp[i] > ampmin) atmp[i] = stag2->amp[i] / stag1->amp[i];
            else atmp[i] = stag2->amp[i] / ampmin;
        }
    } else { cerr << "Error(DoCor): Unknown Cor type!" << endl; exit(-1); }
    for (i = 0; i < len; i++) ptmp[i] = stag2->pha[i] - stag1->pha[i];
    IFFT(len, atmp, ptmp, seis_out, &ns);
    delete [] atmp;
    atmp = NULL;
    delete [] ptmp;
    ptmp = NULL;

    if (lag > ns / 2) {
        fprintf(stdout, "*** Warning: lagtime overflow, corrected back to max len ***");
        lag = ns / 2;
    }
    // Compute and correct for cc time-length
    float cor_rec[2 * lag + 1], cor[2 * lag + 1];
    if (!CalcRecCor(*stag1, *stag2, cor_rec, lag, dt1) ) return 0;

    for (i = 1; i < (lag + 1); i++) {
        cor[lag - i] = seis_out[i] / cor_rec[lag - i];
        cor[lag + i] = seis_out[ns - i] / cor_rec[lag + i];
    }
    cor[lag] = seis_out[0] / cor_rec[lag];
    delete [] seis_out;
    seis_out = NULL;

    // output daily cross-correlations
    if (CorOutflag > 0) {
        SAC_HD shdcorD = shdcor;
        shdcorD.user0 = 1;
        write_sac(outname, cor, &shdcorD);
    }

    // Update sigcor & shdcor
    if (fprcs) if (!CheckPrecNoise(lag, dt1, cor, shdcor.dist) ) return 0;

    if (dt1 != shdcor.delta) {
        fprintf(stdout, "*** Warning: Incompatible sampling rate ( %f %f )! Skipped. ***", dt1, shdcor.delta);
        return 0;
    }
    pthread_mutex_lock(&addlock); // lock
    *dnum        = *dnum + 1;
    shdcor.user0 = *dnum;
    *dflag       = 1;
    for (i = 0; i < (2 * lag + 1); i++) sigcor[i] += cor[i];
    pthread_mutex_unlock(&addlock); // unlock

    return 1;
} /* DoCor */

extern long MemAvail;
int
GroupSize(int * spnpts)
{
    int iev, ist, nmax = 0, npts;

    for (iev = 0; iev < NEVENTS; iev++)
        for (ist = 0; ist < sdbnew->nst; ist++)
            if (sdbnew->rec[iev][ist].n > nmax) nmax = sdbnew->rec[iev][ist].n;
    npts    = (int) pow(2., floor(log(nmax) / log(2.))) + 10;
    *spnpts = npts;
    int size =
      (int) floor((memomax * MemAvail * 0.8 - 2 * sizeof(SAC_DB) - (10000000. + 18. * npts) * sizeof(float)
        - (0.5 * (sdbnew->nst) * (sdbnew->nst - 1) * (NEVENTS + 1) * sizeof(short))) / NEVENTS
        / (sizeof(starec) + 2. * npts * sizeof(float))) - 1;
    if (size < 1) {
        cerr
            <<
          "ERROR(GroupSize): No enough memory for Crosscorrelating!! Increase max memmory% or switch to CrossCor_old.c "
            << endl;
        exit(0);
    }
    if (size > sdbnew->nst) size = sdbnew->nst;
    // if(size > 0.5*sdbnew->nst) size = sdbnew->nst;
    return size;
}

int
ExtractFlag(char * dirname, char * str, int * is1, int * is2, char * pflagout)
{
    char * name = strtok(str, " "), fname[200];

    sprintf(fname, "%s/%s", dirname, name);
    if (access(fname, R_OK) != 0) return -1;

    char * pflag = strtok(NULL, " ");
    strtok(name, "/");
    strtok(NULL, "_");
    char * sta1 = strtok(NULL, "_");
    char * sta2 = strtok(NULL, ".");
    int i, ii;
    for (i = 0; i < sdbnew->nst; i++) {
        ii = *is1 - i;
        if (ii >= 0) if (strcmp(sdbnew->st[ii].name, sta1) == 0) { *is1 = ii; break; }
        ii = *is1 + i;
        if (ii < sdbnew->nst) if (strcmp(sdbnew->st[ii].name, sta1) == 0) { *is1 = ii; break; }
    }
    if (i >= sdbnew->nst) return 0;

    for (i = 0; i < sdbnew->nst; i++) {
        ii = *is2 - i;
        if (ii >= 0) if (strcmp(sdbnew->st[ii].name, sta2) == 0) { *is2 = ii; break; }
        ii = *is2 + i;
        if (ii < sdbnew->nst) if (strcmp(sdbnew->st[ii].name, sta2) == 0) { *is2 = ii; break; }
    }
    if (i >= sdbnew->nst) return 0;

    sprintf(pflagout, "%s", pflag);
    return 1;
}

void
PrepareCor(char * dirname, char * dirnameD, short *** dnum, short **** dflag)
{
    char str[300];
    FILE * fd;
    int is1, is2, iev;
    int dex = 1;

    // short dnum[sdbnew->nst][sdbnew->nst], dflag[sdbnew->nst][sdbnew->nst][NEVENTS];
    *dnum  = (short **) malloc(sdbnew->nst * sizeof(short *));
    *dflag = (short ***) malloc(sdbnew->nst * sizeof(short **));
    for (is1 = 0; is1 < sdbnew->nst; is1++) {
        (*dnum)[is1]  = (short *) malloc(sdbnew->nst * sizeof(short));
        (*dflag)[is1] = (short **) malloc(sdbnew->nst * sizeof(short *));
        for (is2 = is1; is2 < sdbnew->nst; is2++) {
            (*dnum)[is1][is2] = 0; // mark dnum 0 for the station pairs to be processed.
            if (sdbnew->st[is1].flag == sdbnew->st[is2].flag) { // skip if both in group 0
                if (sdbnew->st[is1].flag == 0) (*dnum)[is1][is2] = -1;
            } else if (sdbnew->st[is1].flag * sdbnew->st[is2].flag != 0) {
                (*dnum)[is1][is2] = -1; // skip if in different groups that are both not group 0
            }
            (*dflag)[is1][is2] = (short *) malloc(NEVENTS * sizeof(short));
            for (iev = 0; iev < NEVENTS; iev++) (*dflag)[is1][is2][iev] = 0;
        }
    }

    // make dir for daily correlation
    if (CorOutflag > 0) {
        if (access(dirnameD, 0) != 0) MKDir(dirnameD);
        for (is1 = 0; is1 < sdbnew->nst; is1++) {
            sprintf(str, "%s/%s", dirnameD, sdbnew->st[is1].name);
            MKDir(str);
        }
    }

    if (access(dirname, 0) != 0) {
        dex = 0;
        MKDir(dirname);
        if (fskip4)
            cout << "   Cannot access " << dirname << ". No Correlation will be skipped!" << endl;
    }
    for (is1 = 0; is1 < sdbnew->nst; is1++) {
        sprintf(str, "%s/%s", dirname, sdbnew->st[is1].name);
        MKDir(str);
    }
    if (!dex) return;

    sprintf(str, "%s/Cor_dayflag.lst", dirname);
    if ((fd = fopen(str, "r")) == NULL) {
        if (fskip4)
            cout << "   Cannot access " << str << ". No Correlation will be skipped!" << endl;
        return;
    }
    char str2[300];
    sprintf(str, "%s/Cor_dayflag.lst", dirname);
    sprintf(str2, "%s/Cor_dayflag.lst_old", dirname);
    Copy(str, str2);
    char pflag[NEVENTS + 1];
    if (fskip4) {
        for (is1 = 0, is2 = 0; fgets(str, 300, fd) != NULL;) {
            if (ExtractFlag(dirname, str, &is1, &is2, pflag) <= 0) continue;
            // for(iev=0;iev<NEVENTS;iev++) (*dflag)[is1][is2][iev] = pflag[iev]-'0';
            (*dnum)[is1][is2] = -1;
        }
        sprintf(str, "%s/Cor_dayflag.lst_tmp", dirname);
        fRemove(str);
    } else {
        int sflag;
        FILE * fd2;
        fclose(fd);
        sprintf(str, "%s/Cor_dayflag.lst_old", dirname);
        fd = fopen(str, "r");
        sprintf(str, "%s/Cor_dayflag.lst_tmp", dirname);
        fd2 = fopen(str, "w");
        for (is1 = 0, is2 = 0; fgets(str, 300, fd) != NULL;) {
            sprintf(str2, "%s", str);
            sflag = ExtractFlag(dirname, str2, &is1, &is2, pflag);
            if ( (sflag == 1 && (*dnum)[is1][is2] == -1) || sflag == 0)
                fprintf(fd2, "%s", str);
        }
        fclose(fd2);
    }
    fclose(fd);
} /* PrepareCor */

void *
DoCorEntrance(void * outnameD)
{
    // int ithread = *((int *)tid);
    int iev;
    char outname[300];

    for (;;) {
        // get/update current event number
        pthread_mutex_lock(&cevlock);
        iev = currevn;
        currevn++;
        pthread_mutex_unlock(&cevlock);
        if (iev >= NEVENTS) break;
        // Do Correlation
        if (stag[i1].dayrec[iev].nr == 0 || starec2[0].dayrec[iev].nr == 0) continue; fprintf(stdout, " %d", iev + 1);
        sprintf(outname, "%s_%d.SAC", (char *) outnameD, iev + 1);
        if (DoCor(&(stag[i1].dayrec[iev]), &(starec2[0].dayrec[iev]),
          &(dnum[is1][is2]), &(dflag[is1][is2][iev]), outname) )
        {
            nrec++;
            // if( npair%10 == 0) fprintf(stderr, "\n   ");
        }
    }
    pthread_exit(NULL);
}

void
CrossCorr()
{
    if (fskip4 == 2) return;

    cortype = 0;
    int iev, npair;
    char dirnameM[200], dirnameD[200];

    // Update SAC_DB for target month (remove empty stations and records)
    sdbnew  = new SAC_DB;
    *sdbnew = *sdb;
    UpdateStaLst();
    // Check memory size and allocate memories for station groups
    int spnpts, size = GroupSize(&spnpts);
    int ng = (int) ceil(sdbnew->nst / (float) size);
    stag = new mstarec[size];
    for (i1 = 0; i1 < size; i1++) {
        // stag[i1] = (mstarec *) malloc ( sizeof(mstarec) );
        for (iev = 0; iev < NEVENTS; iev++) {
            stag[i1].dayrec[iev].amp = (float *) malloc(spnpts * sizeof(float));
            stag[i1].dayrec[iev].pha = (float *) malloc(spnpts * sizeof(float));
        }
    }
    starec2 = new mstarec[1];
    for (iev = 0; iev < NEVENTS; iev++) {
        starec2[0].dayrec[iev].amp = (float *) malloc(spnpts * sizeof(float));
        starec2[0].dayrec[iev].pha = (float *) malloc(spnpts * sizeof(float));
    }

    // Initialize day-num-flag list
    sprintf(dirnameM, "%s/COR", sdbnew->mo[imonth].name);
    sprintf(dirnameD, "%s/COR_D", sdbnew->mo[imonth].name);
    PrepareCor(dirnameM, dirnameD, &dnum, &dflag);

    // Thread id and locks
    pthread_t tid[NTHRDS + 1];
    pthread_mutex_init(&(addlock), NULL);
    int ithread, rc;
    int targs[NTHRDS + 1];
    char filename[200], outnameD[200];
    // do cross-correlatings for each group pairs
    setvbuf(stdout, NULL, _IOLBF, 0);
    for (ig1 = 0; ig1 < ng; ig1++) {
        // ReadGroup(size, ig1, stag);
        targs[NTHRDS] = NTHRDS;
        rc = pthread_create(&tid[NTHRDS], &attr_j, ReadGroupEntrance, (void *) (&(size)) );
        if (rc) {
            cerr << "ERROR(CrossCorr): Creation failed! ERR: " << strerror(rc) << endl;
            exit(0);
        }
        if (cortype == 0) fprintf(stdout, "### Cross-Correlating ");  // cout<<"### Cross-Correlating ";
        else if (cortype == 1) fprintf(stdout, "### Deconvolution(1o2) ");
        else if (cortype == 2) fprintf(stdout, "### Deconvolution(2o1) ");
        fprintf(stdout, "for month %s between group %d(/%d) stations (%d - %d) and other stations: ",
          sdbnew->mo[imonth].name, ig1 + 1, ng, ig1 * size, ig1 * size + size - 1);
        npair = 0;
        nrec  = 0;
        for (is2 = ig1 * size; is2 < sdbnew->nst; is2++) {
            i2 = is2 - ig1 * size;
            // if(i2<size) starec2[0] = stag[i2]; need modification here!
            ReadGroup(-1, is2, starec2);
            while (iread < size && i2 > iread) {
                fprintf(stdout, "\n   Waiting for file readin\n");
                sleep(1);
            }
            for (i1 = 0; i1 < size; i1++) {
                // initialize cross-correlation info
                if (i1 > i2) break;
                is1 = i1 + ig1 * size;
                if (dnum[is1][is2] == -1) continue;
                if (!InitCor(is1, is2, starec2[0].dayrec) ) continue;
                sprintf(outnameD, "%s/%s/COR_%s_%s", dirnameD, sdbnew->st[is1].name, sdbnew->st[is1].name,
                  sdbnew->st[is2].name);
                sprintf(filename, "%s/%s/COR_%s_%s.SAC", dirnameM, sdbnew->st[is1].name, sdbnew->st[is1].name,
                  sdbnew->st[is2].name);
                npair++;
                fprintf(stdout, "\n   %s-%s: ", sdbnew->st[is1].name, sdbnew->st[is2].name);
                // create threads to work on each event
                currevn = 0;
                for (ithread = 0; ithread < NTHRDS; ithread++) {
                    targs[ithread] = ithread;
                    rc = pthread_create(&tid[ithread], &attr_j, DoCorEntrance, (void *) (outnameD)); // (void *) (&(targs[ithread])) );
                    if (rc) {
                        cerr << "ERROR(CrossCorr): Creation failed! ERR: " << strerror(rc) << endl;
                        exit(0);
                    }
                }
                for (ithread = 0; ithread < NTHRDS; ithread++) pthread_join(tid[ithread], NULL);
                if (CorOutflag != 1 && shdcor.user0 > 0) write_sac(filename, sigcor, &shdcor);
                delete [] sigcor;
                sigcor = NULL;
            }
        }
        pthread_join(tid[NTHRDS], NULL);
        // if( npair==0 || npair%10 != 0 ) fprintf(stderr, "\n   ");
        fprintf(stdout, "\n");
        if (fdel2) DeleteGroup(size);
        fprintf(stdout, "%d station pairs (%d records) processed. ###\n", npair, nrec);
    }

    pthread_mutex_destroy(&addlock);

    // clean up group memories
    for (is1 = 0; is1 < size; is1++) {
        for (iev = 0; iev < NEVENTS; iev++) {
            free(stag[is1].dayrec[iev].amp);
            free(stag[is1].dayrec[iev].pha);
        }
        // free(stag[is1]);
    }
    delete [] stag;
    stag = NULL;
    for (iev = 0; iev < NEVENTS; iev++) {
        free(starec2[0].dayrec[iev].amp);
        free(starec2[0].dayrec[iev].pha);
    }
    delete [] starec2;
    starec2 = NULL;

    // output/cleanup day-num-flag list
    char str[300];
    FILE * frec;
    sprintf(filename, "%s/Cor_dayflag.lst", dirnameM);
    sprintf(str, "%s_tmp", filename);
    Move(str, filename);
    frec = fopen(filename, "a");
    for (is1 = 0; is1 < sdbnew->nst; is1++) {
        for (is2 = is1; is2 < sdbnew->nst; is2++) {
            if (dnum[is1][is2] > 0) {
                sprintf(filename, "%s/COR_%s_%s.SAC", sdbnew->st[is1].name, sdbnew->st[is1].name, sdbnew->st[is2].name);
                fprintf(frec, "%s   ", filename);
                for (iev = 0; iev < NEVENTS; iev++) fprintf(frec, "%d", dflag[is1][is2][iev]);
                fprintf(frec, "\n");
            }
            free(dflag[is1][is2]);
        }
        free(dnum[is1]);
        free(dflag[is1]);
    }
    fclose(frec);

    free(dnum);
    free(dflag);
    delete sdbnew;
    sdbnew = NULL;
} /* CrossCorr */
