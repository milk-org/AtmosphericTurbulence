/**
 * @file    make_wavefront_sequence.h
 * @brief   Create wavefront sequence
 *
 *
 */

#ifndef _ATMOSPHERICTURBULENCE_MAKE_WAVEFRONT_SEQUENCE_H
#define _ATMOSPHERICTURBULENCE_MAKE_WAVEFRONT_SEQUENCE_H

errno_t AtmosphericTurbulence_make_wavefront_sequence_addCLIcmd();

errno_t AtmosphericTurbulence_make_wavefront_sequence(float slambdaum,
        long  WFprecision,
        int   compmode);

#endif
