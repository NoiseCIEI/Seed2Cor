#include "Param.h"

///////////////
char buff[300];
///////////////

int
jday(int y, int m, int d);

double
abs_time(int yy, int jday, int hh, int mm, int ss, int ms);

void
SortSeed(int * year, int * month, int * day, char ** name, int n)
{
    int i, j, mtmp, ytmp, dtmp;
    char ctmp[200];

    for (i = 0; i < n; i++) {
        dtmp = day[i];
        mtmp = month[i];
        ytmp = year[i];
        sprintf(ctmp, "%s", name[i]);
        for (j = i;
          j > 0 &&
          (ytmp < year[j - 1] || (ytmp == year[j - 1] && mtmp < month[j - 1]) ||
          (ytmp == year[j - 1] && mtmp == month[j - 1] && dtmp < day[j - 1]));
          j--)
        {
            day[j]   = day[j - 1];
            month[j] = month[j - 1];
            year[j]  = year[j - 1];
            sprintf(name[j], "%s", name[j - 1]);
        }
        day[j]   = dtmp;
        month[j] = mtmp;
        year[j]  = ytmp;
        sprintf(name[j], "%s", ctmp);
    }
}

void
FillMonths()
{
    char month_name[12][6];

    sprintf(month_name[0], "JAN");
    sprintf(month_name[1], "FEB");
    sprintf(month_name[2], "MAR");
    sprintf(month_name[3], "APR");
    sprintf(month_name[4], "MAY");
    sprintf(month_name[5], "JUN");
    sprintf(month_name[6], "JUL");
    sprintf(month_name[7], "AUG");
    sprintf(month_name[8], "SEP");
    sprintf(month_name[9], "OCT");
    sprintf(month_name[10], "NOV");
    sprintf(month_name[11], "DEC");
    FILE * fseed;
    char * seedn[NMONTH * 31];
    int i, j, nseed;
    int year[NMONTH * 31], month[NMONTH * 31], day[NMONTH * 31];
    if ((fseed = fopen(seedlst, "r")) == NULL) {
        cerr << "ERROR(FillMonths): Cannot open file " << seedlst << endl;
        exit(-1);
    }
    for (i = 0; i < NMONTH * 31; i++) seedn[i] = (char *) malloc(200 * sizeof(char));
    for (i = 0; i < NMONTH * 31; i++) {
        if ((fgets(buff, 300, fseed)) == NULL) break;
        sscanf(buff, "%s %d %d %d", seedn[i], &(year[i]), &(month[i]), &(day[i]));
        // cerr<<i<<" "<<seedn[i]<<" "<<year[i]<<" "<<month[i]<<" "<<day[i]<<endl;
    }
    fclose(fseed);
    if (i == NMONTH * 31) {
        cerr << "ERROR(FillMonths): Max #month reached. Increase NMONTH!" << endl;
        exit(-1);
    }
    nseed = i;
    SortSeed(&(year[0]), &(month[0]), &(day[0]), &(seedn[0]), nseed);
    sdb->nmo = -1;
    for (i = 0; i < nseed; i++) {
        if (i == 0 || month[i] != month[i - 1]) {
            sdb->nmo = sdb->nmo + 1;
            sprintf(sdb->mo[sdb->nmo].name, "%d.%3s", year[i], month_name[month[i] - 1]);
            MKDir(sdb->mo[sdb->nmo].name);
            sdb->mo[sdb->nmo].y = year[i];
            sdb->mo[sdb->nmo].m = month[i];
            for (j = 0; j < 31; j++) sprintf(sdb->mo[sdb->nmo].seedf[j], "0");
            if (sdb->nmo >= NMONTH) {
                cerr << "ERROR(FillMonths): Max #month reached. Increase NMONTH!" << endl;
                exit(-1);
            }
        }
        sprintf(sdb->mo[sdb->nmo].seedf[day[i] - 1], "%s", seedn[i]);
    }
    for (i = 0; i < NMONTH * 31; i++) free(seedn[i]);
    sdb->nmo = sdb->nmo + 1;
    cout << sdb->nmo << " months filled" << endl;
} /* FillMonths */

void
FillStations()
{
    FILE * fsta;
    int ist;

    if ((fsta = fopen(stalst, "r")) == NULL) {
        cerr << "ERROR(FillStations): Cannot open file " << stalst << endl;
        exit(0);
    }
    for (ist = 0;; ist++) {
        if (!fgets(buff, 300, fsta) ) break;
        sdb->st[ist].flag = 1;
        sscanf(buff, "%s %f %f %d", sdb->st[ist].name, &(sdb->st[ist].lon), &(sdb->st[ist].lat), &(sdb->st[ist].flag) );
        // fprintf(stderr,"Station %s filled\n", sdb->st[ist].name,  sdb->st[ist].lon, sdb->st[ist].lat );
    }
    fclose(fsta);
    sdb->nst = ist;
    cout << ist << " stations filled" << endl;
}

void
FillEvents()
{
    int i, iev;

    for (i = 0, iev = 0; i < NEVENTS; i++) {
        if (strcmp(sdb->mo[imonth].seedf[i], "0") == 0) continue;
        sdb->ev[i].yy   = sdb->mo[imonth].y;
        sdb->ev[i].mm   = sdb->mo[imonth].m;
        sdb->ev[i].dd   = i + 1;
        sdb->ev[i].h    = 0;
        sdb->ev[i].m    = 0;
        sdb->ev[i].s    = 0;
        sdb->ev[i].ms   = 0;
        sdb->ev[i].jday = jday(sdb->ev[i].yy, sdb->ev[i].mm, sdb->ev[i].dd);
        sdb->ev[i].t0   = abs_time(sdb->ev[i].yy, sdb->ev[i].jday, sdb->ev[i].h, sdb->ev[i].m, sdb->ev[i].s,
            sdb->ev[i].ms);
        sprintf(sdb->ev[i].name, "%s.%d", sdb->mo[imonth].name, sdb->ev[i].dd);
        sprintf(buff, "%s/%s", sdb->mo[imonth].name, sdb->ev[i].name);
        MKDir(buff);
        iev++;
    }
    sdb->nev = iev;
    cout << "   " << iev << " events filled" << endl;
}
