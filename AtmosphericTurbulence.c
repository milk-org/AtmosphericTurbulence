/**
 * @file    AtmosphericTurbulence.c
 * @brief   Atmospheric turbulence simulator
 *
 * Generates wavefronts\n
 */


#define _GNU_SOURCE

/* ================================================================== */
/* ================================================================== */
/*  MODULE INFO                                                       */
/* ================================================================== */
/* ================================================================== */

// module default short name
// all CLI calls to this module functions will be <shortname>.<funcname>
// if set to "", then calls use <funcname>
#define MODULE_SHORTNAME_DEFAULT "atmturb"

// Module short description
#define MODULE_DESCRIPTION       "Atmospheric turbulence simulator"






/* ================================================================== */
/* ================================================================== */
/*  HEADER FILES                                                      */
/* ================================================================== */
/* ================================================================== */


#include "CommandLineInterface/CLIcore.h"



//
// Forward declarations are required to connect CLI calls to functions
// If functions are in separate .c files, include here the corresponding .h files
//
#include "make_master_turbulence_screen.h"
#include "make_master_turbulence_screen_pow.h"
#include "make_wavefront_sequence.h"





/* ================================================================== */
/* ================================================================== */
/*  INITIALIZE LIBRARY                                                */
/* ================================================================== */
/* ================================================================== */

// Module initialization macro in CLIcore.h
// macro argument defines module name for bindings
//
INIT_MODULE_LIB(AtmosphericTurbulence)



/**
 * @brief Initialize module CLI
 *
 * CLI entries are registered: CLI call names are connected to CLI functions.\n
 * Any other initialization is performed\n
 *
 */
static errno_t init_module_CLI()
{

    //CLI_CMD_CONNECT("func1", "create_image_with_value");
    
    AtmosphericTurbulence_make_master_turbulence_screen_addCLIcmd();
	AtmosphericTurbulence_make_master_turbulence_screen_pow_addCLIcmd();

	AtmosphericTurbulence_make_wavefront_sequence_addCLIcmd();
    // optional: add atexit functions here

    return RETURN_SUCCESS;
}








