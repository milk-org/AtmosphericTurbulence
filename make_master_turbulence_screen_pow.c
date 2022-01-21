/**
 * @file    make_master_turbulence_screen_pow.c
 * @brief   Create master turbulence screen following power law index
 *
 *
 */

#include <math.h>

#include "CommandLineInterface/CLIcore.h"

// required for create_2Dimage_ID
#include "COREMOD_arith/COREMOD_arith.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "fft/fft.h"
#include "image_gen/image_gen.h"

// ==========================================
// Forward declaration(s)
// ==========================================

errno_t AtmosphericTurbulence_make_master_turbulence_screen_pow(
    const char *ID_name1, const char *ID_name2, uint32_t size, float power);

// ==========================================
// Command line interface wrapper function(s)
// ==========================================

/** @brief Example CLI function
 *
 * Command Line Interface (CLI) wrapper to function\n
 * This is referred to as a "CLI function", written to connect a command
 * on the CLI prompt to a function.\n
 * A CLI function will check arguments entered on the prompt, and pass
 * them to the function.
 *
 * Naming conventions:
 * - CLI function <modulename>__<function>__cli()
 * - Execution function <modulename>__<function>()
 *
 *
 *
 * ### Checking if arguments are valid
 *
 * Each argument is checked by calling CLI_checkarg(i, argtype), which
 * checks if argument number <i> conforms to type <argtype>.\n
 *
 * Types are defined in CLIcore.h. Common types are:
 * - CLIARG_FLOAT             floating point number
 * - CLIARG_LONG              integer (int or long)
 * - CLIARG_STR_NOT_IMG       string, not existing image
 * - CLIARG_IMG               existing image
 * - CLIARG_STR               string
 */
static errno_t AtmosphericTurbulence_make_master_turbulence_screen_pow__cli()
{
    if (0 + CLI_checkarg(1, CLIARG_STR_NOT_IMG) +
            CLI_checkarg(2, CLIARG_STR_NOT_IMG) + CLI_checkarg(3, CLIARG_LONG) +
            CLI_checkarg(4, CLIARG_FLOAT) ==
        0)
    {
        // If arguments meet requirements, command is executed
        //
        AtmosphericTurbulence_make_master_turbulence_screen_pow(
            data.cmdargtoken[1].val.string,
            data.cmdargtoken[2].val.string,
            data.cmdargtoken[3].val.numl,
            data.cmdargtoken[4].val.numf);

        return CLICMD_SUCCESS;
    }
    else
    {
        // If arguments do not pass test, errror code returned
        return CLICMD_INVALID_ARG;
    }
}

// ==========================================
// Register CLI command(s)
// ==========================================

errno_t AtmosphericTurbulence_make_master_turbulence_screen_pow_addCLIcmd()
{

    RegisterCLIcommand(
        "mkturbscreenpow",
        __FILE__,
        AtmosphericTurbulence_make_master_turbulence_screen_pow__cli,
        "creates master turbulence screens with power law index",
        "<screen1> <screen2> <size> <powindex>",
        "mkturbscreenpow scr1 scr2 512 1.7",
        "AtmosphericTurbulence_make_master_turbulence_screen_pow(const char "
        "*ID_name1, const char "
        "*ID_name2, uint32_t size, float power)");

    return RETURN_SUCCESS;
}

// By convention, function name starts with <modulename>__
//
errno_t AtmosphericTurbulence_make_master_turbulence_screen_pow(
    const char *ID_name1, const char *ID_name2, uint32_t size, float power)
{
    imageID  ID;
    uint32_t ii;
    uint32_t jj;
    float    value, C1, C2;
    long     cnt;
    long     Dlim = 3;

    make_rnd("tmppha", size, size, "");
    arith_image_cstmult("tmppha", 2.0 * PI, "tmppha1");
    delete_image_ID("tmppha", DELETE_IMAGE_ERRMODE_WARNING);
    make_dist("tmpd", size, size, size / 2, size / 2);
    make_rnd("tmpg", size, size, "-gauss");

    arith_image_cstpow("tmpd", power, "tmpd1");
    delete_image_ID("tmpd", DELETE_IMAGE_ERRMODE_WARNING);
    arith_image_div("tmpg", "tmpd1", "tmpamp");
    delete_image_ID("tmpg", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("tmpd1", DELETE_IMAGE_ERRMODE_WARNING);
    arith_set_pixel("tmpamp", 0.0, size / 2, size / 2);
    mk_complex_from_amph("tmpamp", "tmppha1", "tmpc", 0);
    delete_image_ID("tmpamp", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("tmppha1", DELETE_IMAGE_ERRMODE_WARNING);
    permut("tmpc");
    do2dfft("tmpc", "tmpcf");
    delete_image_ID("tmpc", DELETE_IMAGE_ERRMODE_WARNING);
    mk_reim_from_complex("tmpcf", "tmpo1", "tmpo2", 0);
    delete_image_ID("tmpcf", DELETE_IMAGE_ERRMODE_WARNING);

    /* compute the scaling factor in the power law of the structure function */
    fft_structure_function("tmpo1", "strf");
    ID    = image_ID("strf");
    value = 0.0;
    cnt   = 0;
    for (ii = 1; ii < Dlim; ii++)
        for (jj = 1; jj < Dlim; jj++)
        {
            value += log10(data.image[ID].array.F[jj * size + ii]) -
                     power * log10(sqrt(ii * ii + jj * jj));
            /*	printf("%ld %ld %f\n",ii,jj,log10(data.image[ID].array.F[jj*size+ii])-5.0/3.0*log10(sqrt(ii*ii+jj*jj)));*/
            cnt++;
        }
    delete_image_ID("strf", DELETE_IMAGE_ERRMODE_WARNING);
    C1 = pow(10.0, value / cnt);

    fft_structure_function("tmpo2", "strf");
    ID    = image_ID("strf");
    value = 0.0;
    cnt   = 0;
    for (ii = 1; ii < Dlim; ii++)
        for (jj = 1; jj < Dlim; jj++)
        {
            value += log10(data.image[ID].array.F[jj * size + ii]) -
                     power * log10(sqrt(ii * ii + jj * jj));
            cnt++;
        }
    delete_image_ID("strf", DELETE_IMAGE_ERRMODE_WARNING);
    C2 = pow(10.0, value / cnt);
    /*  printf("%f %f\n",C1,C2);*/
    arith_image_cstmult("tmpo1", 1.0 / sqrt(C1), ID_name1);
    arith_image_cstmult("tmpo2", 1.0 / sqrt(C2), ID_name2);
    delete_image_ID("tmpo1", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("tmpo2", DELETE_IMAGE_ERRMODE_WARNING);

    return RETURN_SUCCESS;
}
