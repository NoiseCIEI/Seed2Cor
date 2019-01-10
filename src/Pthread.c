#include "Param.h"

int NTHRDS;
long MemAvail;

void
EstimateMemAvail(long &MemAvail);

void
InitialPthread(int num_threads)
{
    size_t stacksize = 16777216;

    if (num_threads == -1) {
        NTHRDS  = sysconf(_SC_NPROCESSORS_ONLN);
        NTHRDS += NTHRDS / 10; // NTHRDS = 1;
        NTHRDS -= 1;
    } else {
        NTHRDS = num_threads;
    }
    cout << "*** a total of " << NTHRDS << " threads will be created. ***" << endl;
    pthread_attr_init(&attr_j);
    pthread_attr_setdetachstate(&attr_j, PTHREAD_CREATE_JOINABLE);
    pthread_attr_setstacksize(&attr_j, stacksize);
    pthread_mutex_init(&cevlock, NULL);
    pthread_mutex_init(&fftlock, NULL);
    pthread_mutex_init(&fiolock, NULL);
    EstimateMemAvail(MemAvail);
    cout << "*** Estimated total available memory = " << MemAvail / (1024. * 1024.) << "Mb. ***" << endl;
}

void
CleanupPthread()
{
    pthread_attr_destroy(&attr_j);
    pthread_mutex_destroy(&cevlock);
    pthread_mutex_destroy(&fftlock);
    pthread_mutex_destroy(&fiolock);
}
