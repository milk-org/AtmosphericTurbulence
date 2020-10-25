/**
 * @file    make_master_turbulence_screen.h
 *
 */

#ifndef _MAKE_MASTER_TURBULENCE_SCREEN_H
#define _MAKE_MASTER_TURBULENCE_SCREEN_H


errno_t AtmosphericTurbulence_make_master_turbulence_screen_addCLIcmd();

errno_t AtmosphericTurbulence_make_master_turbulence_screen(
    const char *ID_name1,
    const char *ID_name2,
    uint32_t size,
    float outerscale,
    float innerscale,
    long WFprecision
);

#endif
