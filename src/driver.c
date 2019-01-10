/***********************************************************************************************************
*  This is an organized & parallelized version of the ambient noise code, combining all the codes involved in
*  the process of ambient noise data processing, including sac_from_seed(), cut_trans_RESP(), whiten_rej_phamp(),
*  filter4(), set_sacdb() and justCOR().
*
*  To make the code easier to use, all the old-version codes have been modified and re-organized. All important
*  parameters have been pulled out and can be set easily by modifying the input parameter file.
*
*  The code calls system() at two places. The first to rdseed in ProduceSac.c and the second to evalresp in
*  RemoveResp.c. The old sac calls for removing instrument responses have been replaced by evalresp plus a
*  c subroutine for both safety and performance purpose. All other system() calls have been replaced by c
*  codes and can be found in SysTools.c.
*
*  The process was multi-threaded in a way that it handles daily records parallelly. Using on multiple nodes
*  will need further modifications
*
*  The cross-correlation code has been modified to 1) read in multiple station records at a time (The number of
*  records read-in depends on the amount of memory available), which boosted the speed by about 20 percent. and
*  2) read/write only once for each monthly cross-correlation (instead of once per day) to allow multiple threads
*  running at the same time. (otherwise the writing speed is gonna be the bottleneck) Switch back to
*  CrossCorr_old.c when there's no enough memory.
*
*  Some new features have been added to the code, including 1) the ability to re-sample and deal with BH data,
*  2) the ability to choose between different normalize methods or to skip normalization processes, 3) an option
*  to skip exist target files, 4) an option to control the cc of which station-pairs to be processed by adding
*  flags to the input station list, and 3) output time-length records of each event and an option to correct for
*  it.
*
*  Attention:
*  1) The temperal normalization method of 'earthquake cutting' has NOT been well tested. Use with caution!
*  2) The spectrum reshape function is currently removed from the code. (Smooth2() in Whiten.c does nothing)
*
*  Copyright C 2013 Tian's Corporation, All rights reserved.
***********************************************************************************************************/

#define MAIN
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <iostream>
#include <string.h>
#include <math.h>
#include "mysac64.h"
#include "64_sac_db.h"
using namespace std;

SAC_DB * sdb;
int imonth;
// -------------------------------------------------------------parameters----------------------------------------------------------//
char rdsexe[200];   // rdseed excutable
char evrexe[200];   // evalresp excutable
char stalst[200];   // station list (J23A -129.6825 +44.844 optional_flag) (flag controls which sta-pairs to be ccd)
                    // (skip cc between sta-pairs that are 1> both flagged 0 or 2> both not 0 but flagged with different numbers)
char seedlst[200];  // SEED file list (down_US_ALL/OBS_BHZ_2012.FEB.29.203762.seed 2012 2 29)
char ch[4];         // channel name
int sps;            // target sampling rate
float gapfrac;      // maximum allowed gap fraction in input sac record
float t1;           // cutting begining time in sec
float tlen;         // time-record length in sec
float perl, perh;   // signal period band
int tnorm_flag;     // temperal normalization method. (1 for onebit, 2 for running average, 3 for earthquake cutting)
float Eperl, Eperh; // estimated earthquake period band (no effect when tnorm_flag==1; set Eperl=-1 to turn off the filter)
float timehlen;     // half len of time window for running average (temperal normalization) (no effect when tnorm_flag==1)
float frechlen;     // half len of frec window for whitenning. recommended: 0.0002; 0: no norm; -1: use input-smoothing file.
char fwname[200];   // input spectrum reshaping signal file (takes effect only when frechlen==-1)
int ftlen;          // turn on/off (1/0) cross-correlation-time-length correction for amplitude
int fprcs;          // turn on/off (1/0) precursor signal checking
float memomax;      // maximum memory fraction to be used. (set according to number of threads)
int lagtime;        // cross-correlation signal half length in sec
int mintlen;        // allowed minimum time length for cross-correlation (takes effect only when ftlen = 1)
int fdel1;          // delete original sac files after removing response when set to 1
int fdel2;          // delete am&ph files after cross-correlation when set to 1
int fskip1;         // skip ExtractSac() when set to 2, skip upon existence of target file when set to 1
int fskip2;         // skip RmRESP() when set to 2, skip upon existence of target file when set to 1
int fskip3;         // skip TempSpecNorm() when set to 2, skip upon existence of target file when set to 1
int fskip4;         // skip CrossCorr() if target file exists when set to 1
int CorOutflag;     // controls the output of cross-corellation results: 0=monthly, 1=daily, 2=both
// ---------------------------------------------------------------------------------------------------------------------------------//

void
FillStations();

void
FillMonths();

void
FillEvents();

void
InitialPthread(int);

void
CleanupPthread();

void
ExtractSac();

void
RmRESP();

void
TempSpecNorm();

void
CrossCorr();

int
isTermi(int deno)
{
    while (deno % 2 == 0) deno /= 2;
    while (deno % 5 == 0) deno /= 5;
    return deno == 1;
}

int
GetParameters(char * fname)
{
    FILE * fparam, * filetmp;
    char buff[300], ctmp[200];
    int ERR = 0, itmp;
    float ftmp;

    if ((fparam = fopen(fname, "r")) == NULL) {
        cerr << "Error(GetParam): Cannot open parameter file " << fname << endl;
        perror("Error");
        exit(0);
    }
    cout << "---------------------------Checking input parameters---------------------------" << endl;
    // rdsexe
    if ( (fgets(buff, 300, fparam)) == NULL) return 0;

    sscanf(buff, "%s", rdsexe);
    cout << "rdseed excutable:\t" << rdsexe << endl;
    if (!strstr(rdsexe, "rdseed") ) cout << "   Warning: Are you sure this is an rdseed excutable?" << endl;
    if (access(rdsexe, R_OK) != 0) {
        cerr << "   Error: cannot access rdseed through " << rdsexe << endl;
        ERR = 1;
    }
    // evrexe
    if ( (fgets(buff, 300, fparam)) == NULL) return 0;

    sscanf(buff, "%s", evrexe);
    cout << "evalresp excutable:\t" << evrexe << endl;
    if (!strstr(evrexe, "evalresp") ) cout << "   Warning: Are you sure this is an evalresp excutable?" << endl;
    if (access(evrexe, R_OK) != 0) {
        cerr << "   Error: cannot access evalresp through " << evrexe << endl;
        ERR = 1;
    }
    // stalst
    if ( (fgets(buff, 300, fparam)) == NULL) return 0;

    sscanf(buff, "%s", stalst);
    cout << "station list:\t\t" << stalst << endl;
    if ((filetmp = fopen(stalst, "r")) == NULL) {
        cerr << "   Error: cannot access file " << stalst << endl;
        ERR = 1;
    } else {
        fgets(buff, 300, filetmp);
        if (sscanf(buff, "%s %f %f", ctmp, &ftmp, &ftmp) != 3) {
            cerr << "   Error: incorrect format in file " << stalst << "! Shoud be (sta lon lat)" << endl;
            ERR = 1;
        }
        fclose(filetmp);
    }
    // seedlst
    if ( (fgets(buff, 300, fparam)) == NULL) return 0;

    sscanf(buff, "%s", seedlst);
    cout << "seed list:\t\t" << seedlst << endl;
    if ((filetmp = fopen(seedlst, "r")) == NULL) {
        cerr << "   Error: cannot access file " << seedlst << endl;
        ERR = 1;
    } else {
        fgets(buff, 300, filetmp);
        if (sscanf(buff, "%s %d %d %d", ctmp, &itmp, &itmp, &itmp) != 4) {
            cerr << "   Error: incorrect format in file " << seedlst << "! (path/seedfile yyyy mm dd) is expected"
                 << endl;
            ERR = 1;
        }
        fclose(filetmp);
    }
    // ch
    if ( (fgets(buff, 300, fparam)) == NULL) return 0;

    sscanf(buff, "%s", ch);
    cout << "channel:\t\t" << ch << endl;
    // sps
    if ( (fgets(buff, 300, fparam)) == NULL) return 0;

    sscanf(buff, "%f", &ftmp);
    sps = (int) ftmp;
    cout << "target sampling rate:\t" << ftmp << endl;
    if (sps != ftmp || sps <= 0) {
        cerr << "   Error: a positive integer is expected!" << endl;
        ERR = 1;
    } else if (!isTermi(sps) ) {
        cerr << "   Error: 1/" << sps << " isn't a terminating decimal, will cause rounding error!" << endl;
        ERR = 1;
    }
    // gapfrac
    if ( (fgets(buff, 300, fparam)) == NULL) return 0;

    sscanf(buff, "%f", &gapfrac);
    cout << "max gap fraction:\t" << gapfrac * 100 << "%" << endl;
    if (gapfrac < 0 || gapfrac > 1) {
        cerr << "   Error: a number between 0. - 1. is expected!" << endl;
        ERR = 1;
    }
    // t1
    if ( (fgets(buff, 300, fparam)) == NULL) return 0;

    sscanf(buff, "%f", &t1);
    cout << "cutting begining:\t" << t1 << "sec" << endl;
    if (fabs(t1) > 86400) cout << "   Warning: " << t1 << "sec exceeds one day." << endl;
    // tlen
    if ( (fgets(buff, 300, fparam)) == NULL) return 0;

    sscanf(buff, "%f", &tlen);
    cout << "time-rec length:\t" << tlen << "sec" << endl;
    ftmp = t1 + tlen;
    if (ftmp > 86400 || ftmp < 0) cout << "   Warning: ending time '" << ftmp << "sec' out of range." << endl;
    // perl perh
    if ( (fgets(buff, 300, fparam)) == NULL) return 0;

    sscanf(buff, "%f", &perl);
    if ( (fgets(buff, 300, fparam)) == NULL) return 0;

    sscanf(buff, "%f", &perh);
    cout << "signal per band:\t" << perl << " - " << perh << "sec" << endl;
    if (perl < 0 || perh < 0) {
        cerr << "   Error: period band can not go below 0." << endl;
        ERR = 1;
    }
    if (perl < 2. / sps) {
        cerr << "   Error: " << perl << "sec is out of the lower limit at sps " << sps << endl;
        ERR = 1;
    }
    if (perl <
      5. / sps)
        cout << "   Warning: signal at " << perl << "sec might be lost at a sampling rate of " << sps << endl;
    // tnorm_flag
    if ( (fgets(buff, 300, fparam)) == NULL) return 0;

    sscanf(buff, "%d", &tnorm_flag);
    cout << "t-norm method:\t\t";
    if (tnorm_flag == 0) { cout << "none" << endl; } else if (tnorm_flag == 1) {
        cout << "One-bit" << endl;
    } else if (tnorm_flag == 2) { cout << "Running average" << endl; } else if (tnorm_flag == 3) {
        cout << "Earthquake cutting" << endl;
    } else if (tnorm_flag == 4) { cout << "Earthquake cutting + Running average" << endl; } else {
        cerr << endl << "   Error: Unknow method. integer between 0 - 4 is expected" << endl;
        ERR = 1;
    }
    // Eperl Eperh
    if ( (fgets(buff, 300, fparam)) == NULL) return 0;

    sscanf(buff, "%f", &Eperl);
    if ( (fgets(buff, 300, fparam)) == NULL) return 0;

    sscanf(buff, "%f", &Eperh);
    if (tnorm_flag != 1) {
        if (Eperl == -1.) { cout << "Eqk filter:\t\toff" << endl; } else {
            cout << "Eqk filter:\t\t" << Eperl << " - " << Eperh << "sec" << endl;
            if (Eperl < 0 || Eperh < 0) {
                cerr << "   Error: period band can not go below 0." << endl;
                ERR = 1;
            }
        }
    }
    // timehlen
    if ( (fgets(buff, 300, fparam)) == NULL) return 0;

    sscanf(buff, "%f", &timehlen);
    if (tnorm_flag != 1) {
        cout << "t-len for run-avg:\t" << timehlen << endl;
        if (timehlen < 0) {
            cerr << "   Error: positive number is expected" << endl;
            ERR = 1;
        }
    }
    // frechlen
    if ( (fgets(buff, 300, fparam)) == NULL) return 0;

    sscanf(buff, "%f", &frechlen);
    cout << "t-len for whitening:\t" << frechlen << "  ";
    if (frechlen == -1) { cout << "input smoothing file will be used"; } else if (frechlen < 0) {
        cerr << endl << "   Error: non-negative number is expected";
        ERR = 1;
    }
    cout << endl;
    // fwname
    if ( (fgets(buff, 300, fparam)) == NULL) return 0;

    sscanf(buff, "%s", fwname);
    if (frechlen == -1) {
        cout << "spec reshaping file:\t" << fwname << endl;
        if (access(fwname, R_OK) != 0) {
            cerr << "   Error: cannot access file " << fwname << endl;
            ERR = 1;
        }
    }
    // ftlen
    if ( (fgets(buff, 300, fparam)) == NULL) return 0;

    sscanf(buff, "%d", &ftlen);
    if (ftlen <= 0) cout << "cor-time-len correction\toff" << endl;
    else cout << "cor-time-len correction\ton" << endl; /*if(tnorm_flag==3) ftlen = 2;*/ // fprcs
    if ( (fgets(buff, 300, fparam)) == NULL) return 0;

    sscanf(buff, "%d", &fprcs);
    if (fprcs <= 0) cout << "prcsr-signal checking\toff" << endl;
    else cout << "prcsr-signal checking\ton" << endl;
    // memomax
    if ( (fgets(buff, 300, fparam)) == NULL) return 0;

    sscanf(buff, "%f", &memomax);
    cout << "memory fraction:\t" << memomax * 100 << "%" << endl;
    if (memomax < 0 || memomax > 1) {
        cerr << "   Error: a number between 0. - 1. is expected!" << endl;
        ERR = 1;
    }
    // lagtime
    if ( (fgets(buff, 300, fparam)) == NULL) return 0;

    sscanf(buff, "%d", &lagtime);
    cout << "cor lag time:\t\t" << lagtime << "sec" << endl;
    if (lagtime < 0) {
        cerr << "   Error: negative lagtime!" << endl;
        ERR = 1;
    } else if (lagtime > 524288) { cout << "   Warning: lag time exceeds the maximum" << endl; }
    // mintlen
    if ( (fgets(buff, 300, fparam)) == NULL) return 0;

    sscanf(buff, "%d", &mintlen);
    cout << "min time len:\t\t" << mintlen << "sec" << endl;
    if (mintlen > 86400) cout << "   Warning: allowed minimum time length larger than a day." << endl;
    // fdel1
    if ( (fgets(buff, 300, fparam)) == NULL) return 0;

    sscanf(buff, "%d", &fdel1);
    cout << "Delete orig sacs?\t";
    if (fdel1 == 1) cout << "Yes" << endl;
    else cout << "No" << endl;
    // fdel2
    if ( (fgets(buff, 300, fparam)) == NULL) return 0;

    sscanf(buff, "%d", &fdel2);
    cout << "Delete am&ph files?\t";
    if (fdel2 == 1) cout << "Yes" << endl;
    else cout << "No" << endl;
    // fskip1
    if ( (fgets(buff, 300, fparam)) == NULL) return 0;

    sscanf(buff, "%d", &fskip1);
    cout << "ExtractSac()\t\t";
    if (fskip1 == 2) cout << "skip" << endl;
    else if (fskip1 == 1) cout << "skip if target file exists" << endl;
    else cout << "overwrite if target file exists" << endl;
    // fskip2
    if ( (fgets(buff, 300, fparam)) == NULL) return 0;

    sscanf(buff, "%d", &fskip2);
    cout << "RmRESP()\t\t";
    if (fskip2 == 2) cout << "skip" << endl;
    else if (fskip2 == 1) cout << "skip if target file exists" << endl;
    else cout << "overwrite if target file exists" << endl;
    // fskip3
    if ( (fgets(buff, 300, fparam)) == NULL) return 0;

    sscanf(buff, "%d", &fskip3);
    cout << "TempSpecNorm()\t\t";
    if (fskip3 == 2) cout << "skip" << endl;
    else if (fskip3 == 1) cout << "skip if target file exists" << endl;
    else cout << "overwrite if target file exists" << endl;
    // fskip4
    if ( (fgets(buff, 300, fparam)) == NULL) return 0;

    sscanf(buff, "%d", &fskip4);
    cout << "CrossCorr()\t\t";
    if (fskip4 == 2) cout << "skip" << endl;
    else if (fskip4 == 1) cout << "skip if target file exists" << endl;
    else cout << "overwrite if target file exists" << endl;
    // CorOutflag
    if ( (fgets(buff, 300, fparam)) == NULL) return 0;

    sscanf(buff, "%d", &CorOutflag);
    cout << "CC Output: \t\t";
    switch (CorOutflag) {
        case 0:
            cout << "monthly" << endl;
            break;
        case 1:
            cout << "daily" << endl;
            break;
        case 2:
            cout << "monthly + daily" << endl;
            break;
        default:
            cerr << endl << "   Error: unknow outflag(" << CorOutflag << ")!" << endl;
            ERR = 1;
    }
    // check for EOF
    if ( (fgets(buff, 300, fparam)) != NULL) cout << "   Warning: End of file not reached!" << endl;

    cout << "-------------------------------Checking completed-------------------------------" << endl;
    fclose(fparam);
    if (ERR == 1) exit(0);
    char cin;
    cout << "Continue?  ";
    fgets(buff, 300, stdin);
    sscanf(buff, "%c", &cin);
    if (cin != 'Y' && cin != 'y') exit(0);
    return 1;
} /* GetParameters */

int
main(int argc, char * argv[])
{
    if (argc != 2 && argc != 3) {
        cerr << "Usage: " << argv[0] << " [Parameters_file] [num_of_threads (optinal)]" << endl;
        return 0;
    }

    sdb = new SAC_DB;
    if (!GetParameters(argv[1]) ) {
        cerr << "Error(main): Parameter loading failed. Incorrect format from file " << argv[1] << endl;
        exit(-1);
    }
    FillStations();
    FillMonths();
    int num_threads = -1;
    if (argc == 3) num_threads = atoi(argv[2]);
    InitialPthread(num_threads);
    for (imonth = 0; imonth < sdb->nmo; imonth++) {
        if (access(sdb->mo[imonth].name, R_OK) == -1) {
            cerr << "Warning(main): Skipping month " << sdb->mo[imonth].name << ". (Month directory deleted)" << endl;
            continue;
        }
        cout << "Working on month " << sdb->mo[imonth].name << endl;
        FillEvents();
        ExtractSac();
        RmRESP();
        TempSpecNorm();
        CrossCorr();
    }
    CleanupPthread();
    delete sdb;

    return 1;
} /* main */
