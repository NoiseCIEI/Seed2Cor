#ifndef MYOMP_H
#define MYOMP_H

#ifdef _OPENMP
#include <omp.h>
#else
typedef int omp_int_t;
inline omp_int_t omp_get_num_threads() { return 1; }
inline omp_int_t omp_get_max_threads() { return 1; }
inline omp_int_t omp_get_thread_num() { return 0; }
#endif

#endif
