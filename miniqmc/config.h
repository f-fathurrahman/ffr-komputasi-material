/* src/ohmms-config.h.in.  Generated from configure.in by autoheader.  */
// -*- c++  -*-
//
//Ohmms Configuration Header. Automatically Generated
//
//See the LICENSE file in the top-level directory for copyright notices
//
#ifndef QMCPLUSPLUS_CONFIGURATION_H
#define QMCPLUSPLUS_CONFIGURATION_H

// clang-format off

/* define the major version */
#define QMCPACK_VERSION_MAJOR  3

/* define the minor version */
#define QMCPACK_VERSION_MINOR  1

/* define the patch version */
#define QMCPACK_VERSION_PATCH  0

/* define the release version */
/* #undef QMCPACK_RELEASE */

/* define the git last commit date */
/* #undef QMCPLUSPLUS_LAST_CHANGED_DATE */

/* Enable OpenMP parallelization. */
#define ENABLE_OPENMP 1

/* Define to 1 if you have the `hdf5' library (-lhdf5). */
/* #undef HAVE_LIBHDF5 */

/* Define to 1 if you want to use parallel hdf5 for frequent output */
/* #undef ENABLE_PHDF5 */

/* Define to 1 if you have MPI library */
/* #undef HAVE_MPI */

/* Define the physical dimension of appliation. */
#define OHMMS_DIM 3

/* Define the index type: int, long */
#define OHMMS_INDEXTYPE int

/* Define the base precision: float, double */
#define OHMMS_PRECISION double

/* Define the full precision: double, long double */
#define OHMMS_PRECISION_FULL double

/* Define to 1 if precision is mixed, only for the CPU code */
/* #undef MIXED_PRECISION */

/* Define to 1 if complex wavefunctions are used */
/* #undef QMC_COMPLEX */

/* Define if the code is specialized for orthorhombic supercell */
#define OHMMS_ORTHO 0

/* Define if sincos function exists */
#define HAVE_SINCOS 1

/* Define if posix_memalign function exists */
#define HAVE_POSIX_MEMALIGN 1

/* Find essl library */
/* #undef HAVE_ESSL */

/* Fund acml library */
/* #undef HAVE_ACML */

/* For AFQMC compilation  */
/* #undef BUILD_AFQMC */

/* For FCIQMC compilation  */
/* #undef BUILD_FCIQMC */

#if defined(__INTEL_COMPILER)
  #if defined(__AVX512F__)
    #define QMC_CLINE 64
    #define QMC_ALIGNAS alignas(64)
    #define ASSUME_ALIGNED(x) __assume_aligned(x,64)
  #else
    #define QMC_CLINE 32
    #define QMC_ALIGNAS alignas(32)
    #define ASSUME_ALIGNED(x) __assume_aligned(x,32)
  #endif
#elif defined(__GNUC__) && !defined(__ibmxl__)
  #if defined(__AVX512F__)
    #define QMC_CLINE 64
    #define QMC_ALIGNAS alignas(64)
    #define ASSUME_ALIGNED(x) (x) = (__typeof__(x)) __builtin_assume_aligned(x,64)
  #else
    #define QMC_CLINE 32
    #define QMC_ALIGNAS alignas(32)
    #define ASSUME_ALIGNED(x) (x) = (__typeof__(x)) __builtin_assume_aligned(x,32)
  #endif
#else
  #define QMC_CLINE 32
  #define QMC_ALIGNAS alignas(32)
  #define ASSUME_ALIGNED(x)
#endif

/* Internal timers */
#define ENABLE_TIMERS 1

/* Use VTune Task API with timers */
/* #undef USE_VTUNE_TASKS */

// clang-format on

#endif // QMCPLUSPLUS_CONFIGURATION_H

