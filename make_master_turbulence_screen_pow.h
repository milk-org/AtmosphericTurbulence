/**
 * @file    make_master_turbulence_screen_pow.h
 *
 */

#ifndef _MAKE_MASTER_TURBULENCE_SCREEN_POW_H
#define _MAKE_MASTER_TURBULENCE_SCREEN_POW_H

errno_t AtmosphericTurbulence_make_master_turbulence_screen_pow_addCLIcmd();

errno_t AtmosphericTurbulence_make_master_turbulence_screen_pow(const char *ID_name1, const char *ID_name2,
                                                                uint32_t size, float power);

#endif
