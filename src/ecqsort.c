#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <signal.h>
#include "intercept_alloc.h"
#include "raise.h"

/* ecqsort_() : Fortran-callable quick-sort */

/* 
   by Sami Saarinen, ECMWF, 7/7/2005 : Interface derived from rsort32.c & rsort64.c
                            4/9/2006 : Dr.Hook call for kwiksort_u64_index
           - " -            7/2/2007 : Intercepting alloc (IBM & NEC SX) + NEC SX vectorization
           - " -            3/7/2007 : Rewritten to use qsort() standard library routine
*/

/* 
   Methods:

   0 : Unsigned 32-bit ints
   1 :   Signed 32-bit ints
   2 :          64-bit doubles (IEEE) : signbit + 11-bit exp + 52-bits mantissa
   3 :          32-bit floats  (IEEE) : signbit +  8-bit exp + 23-bits mantissa
   4 :   Signed 64-bit ints
   5 : Unsigned 64-bit ints

*/

typedef unsigned long long int  Uint64;
typedef          long long int  Sint64;
typedef unsigned int            Uint32;
typedef          int            Sint32;

typedef    short int            Sint16;
typedef   signed char           Sint8;

#ifdef __uxppx__
#ifndef VPP
#define VPP
#endif
#endif

#ifdef VPP
#pragma global noalias
#pragma global novrec
#elif defined(NECSX)
#pragma cdir options -pvctl,nodep
static const int Sint16_limit = 1; /* To prevent generating non-vectorizable constructs */
static const int Sint8_limit  = 1; /* To prevent generating non-vectorizable constructs */
#elif defined(NECSX)
static const int Sint16_limit = 1; /* To prevent generating non-vectorizable constructs */
static const int Sint8_limit  = 1; /* To prevent generating non-vectorizable constructs */
#pragma cdir options -pvctl,nodep
#else
static const int Sint16_limit = 32768;
static const int Sint8_limit  = 128;
#endif

/* Offset adjustment for 64-bit double handling on big/little endian machines */

#define  ALLOC(x,size)    \
 { int bytes = sizeof(*x) * (size); \
   bytes = (bytes < 1) ? 1 : bytes; \
   x = THEmalloc(bytes); \
   if (!x) { fprintf(stderr, \
		     "malloc() of %s (%d bytes) failed in file=%s, line=%d\n", \
		     #x, bytes, __FILE__, __LINE__); RAISE(SIGABRT); } }

#define FREE(x)           if (x) { THEfree(x); x = NULL; }

#if defined(NO_TRUNC) || defined(VPP) || defined(NECSX)
/* For systems without trunc() -function [an extension of ANSI-C, but usually available] */
#define trunc(x) ((x) - fmod((x),1))
#else
extern double trunc(double d);
#endif

#define MakeKwikSort(T, IdxT) \
typedef struct { \
  T value; \
  IdxT j; \
  IdxT idx; \
} T##IdxT##Str_t; \
\
static int \
T##IdxT##cmp(const void *A, const void *B) { \
  const T##IdxT##Str_t *a = A; \
  const T##IdxT##Str_t *b = B; \
  if      ( a->value > b->value ) return  1; \
  else if ( a->value < b->value ) return -1; \
  else { /* a->value == b->value */ \
    /* the next line is the key for the stable qsort() */ \
    return (a->j > b->j) ? 1 : -1; \
  } \
} \
\
static void \
kwiksort_##T##IdxT(const T v[], int n, int index[], int inc, \
                   int index_adj, int mode) \
{ \
  int j; \
  T##IdxT##Str_t *x = NULL; \
  ALLOC(x, n); \
  if (mode < 10) { \
    /* index[] needs to be initialized */ \
    if (inc == 1) { \
      for (j=0; j<n; j++) { x[j].value = v[j]; x[j].j = x[j].idx = j; } \
    } \
    else { \
      for (j=0; j<n; j++) { x[j].value = v[j * inc]; x[j].j = x[j].idx = j; } \
    } \
  } \
  else { \
    if (inc == 1) { \
      for (j=0; j<n; j++) { \
        int tmpidx = index[j] - index_adj; /* Fortran -> C */ \
        x[j].value = v[tmpidx]; \
        x[j].j = j; \
        x[j].idx = tmpidx; \
      } \
    } \
    else { \
      for (j=0; j<n; j++) { \
        int tmpidx = index[j] - index_adj; /* Fortran -> C */ \
        x[j].value = v[tmpidx * inc]; \
        x[j].j = j; \
        x[j].idx = tmpidx; \
      } \
    } \
  } \
  qsort(x, n, sizeof(*x), T##IdxT##cmp); \
  for (j=0; j<n; j++) index[j] = x[j].idx + index_adj; /* C -> Fortran */ \
  FREE(x); \
}

#define kwiksort(T) \
MakeKwikSort(T, Sint8)  \
MakeKwikSort(T, Sint16) \
MakeKwikSort(T, Sint32)

#define SORT_UINT 0
kwiksort(Uint32)

#define SORT_INT  1
kwiksort(Sint32)

#define SORT_R64  2
kwiksort(double)

#define SORT_R32  3
kwiksort(float)

#define SORT_I64  4
kwiksort(Sint64)

#define SORT_U64  5
kwiksort(Uint64)

#define DoSort(T) { \
  T *data = Data; \
  if (n < Sint8_limit) { \
    kwiksort_##T##Sint8(&data[addr], n, index, inc, index_adj, mode); \
  } \
  else if (n < Sint16_limit) { \
    kwiksort_##T##Sint16(&data[addr], n, index, inc, index_adj, mode); \
  } \
  else { \
    kwiksort_##T##Sint32(&data[addr], n, index, inc, index_adj, mode); \
  } \
}

void 
ecqsort_(const    int *Mode,
	 const    int *N,
	 const    int *Inc,
	 const    int *Start_addr,
	         void *Data,
	          int  index[],
	 const    int *Index_adj,
	          int *retc)
{
  int mode = *Mode;
  int method = mode%10;
  int n = *N;
  int rc = n;
  int inc = *Inc;
  int index_adj = *Index_adj;
  int addr = (*Start_addr) - 1; /* Fortran to C */
  Sint32 *ivec = NULL;

  if (method != SORT_UINT   &&
      method != SORT_INT    &&
      method != SORT_R64    &&
      method != SORT_R32    &&
      method != SORT_I64    &&
      method != SORT_U64 ) {
    rc = -1;
    goto finish;
  }

  if (n <= 0) {
    if (n < 0) rc = -2;
    goto finish;
  }

  if (inc < 1) {
    rc = -3;
    goto finish;
  }

  if (method == SORT_R64) {
    /* Check if all values are in fact integers --> method = SORT_INT */
    int j, jj = addr;
    double *data = Data;
    int cnt = 0;
    for (j=0; j<n; j++) {
      double d = data[jj];
      double trc = trunc(d);
      if (d >= -INT_MAX && d <= INT_MAX && d == trc) ++cnt;
#if !defined(NECSX) && !defined(VPP)
      else break;
#endif
      jj += inc;
    }

    if (cnt == n) {
      ALLOC(ivec, n);
      jj = addr;
      for (j=0; j<n; j++) {
	double d = data[jj];
	ivec[j] = (Sint32)d;
	jj += inc;
      }
      method = SORT_INT;
      addr = 0;
      inc = 1;
      Data = ivec; /* voila!! Who says that pointers are no good ? */
    }
  }

  switch (method) {
  case SORT_UINT:
    DoSort(Uint32);
    break;
  case SORT_INT:
    DoSort(Sint32);
    break;
  case SORT_R64:
    DoSort(double);
    break;
  case SORT_R32:
    DoSort(float);
    break;
  case SORT_I64:
    DoSort(Sint64);
    break;
  case SORT_U64:
    DoSort(Uint64);
    break;
  }

 finish:

  FREE(ivec); /* In case ivec[] was allocated */

  *retc = rc;
}
