/**
 * @file    make_AtmosphericTurbulence_vonKarmanWind.h
 *
 */

#ifndef _MAKE_ATMOSPHERICTURBULENCE_VONKARMANWIND_H
#define _MAKE_ATMOSPHERICTURBULENCE_VONKARMANWIND_H

imageID
make_AtmosphericTurbulence_vonKarmanWind(uint32_t vKsize,
                                         float    pixscale,
                                         float    sigmawind,
                                         float    Lwind,
                                         long     size,
                                         const char *restrict IDout_name);

#endif
