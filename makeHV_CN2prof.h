/**
 * @file    makeHV_CN2prof.h
 * @brief   Create HV CN2 turbulence profile
 *
 *
 */

#ifndef _MAKE_HV_CN2PROF_H
#define _MAKE_HV_CN2PROF_H

errno_t AtmosphericTurbulence_makeHV_CN2prof(double      wspeed,
                                             double      r0,
                                             double      sitealt,
                                             long        NBlayer,
                                             const char *outfile);

#endif
