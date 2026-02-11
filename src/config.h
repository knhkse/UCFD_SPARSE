/**
 * @file        config.h
 * @brief       Header file for solver configuration
 * 
 */
#ifndef CONFIG
#define CONFIG
#include <stdint.h>
#include <float.h>
#include <inttypes.h>

// TODO : Auto-generation by Makefile

/**
 * Integer type designation
 */
#if defined(UCFD_INT64)
    typedef int64_t UCFD_INT;
#else
    typedef int32_t UCFD_INT;
#endif

/**
 * Float type designation
 */
#if defined(UCFD_FLOAT32)
    typedef float UCFD_FLOAT;
#else
    typedef double UCFD_FLOAT;
#endif

/**
 * ? If below is written from bespoke generation,
 * ? no need to compile with -DNVARS=7... ? => YES!
 */
// #define NVARS 7
// #define NFVARS 5
// #define NTURBVARS 2
// #define NDIMS 3
// #define BLOCK 5

#endif
