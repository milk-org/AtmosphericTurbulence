/**
 * @file    make_AtmosphericTurbulence_vonKarmanWind.c
 *
 */

#include <math.h>

#include "CommandLineInterface/CLIcore.h"

#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_arith/COREMOD_arith.h"
#include "image_gen/image_gen.h"

#include "fft/dofft.h"
#include "fft/permut.h"



//
// pixscale [m/pix]
// sigmawind [m/s]
// Lwind [m]
//
imageID make_AtmosphericTurbulence_vonKarmanWind(
    uint32_t vKsize,
    float pixscale,
    float sigmawind,
    float Lwind,
    long size,
    const char *restrict IDout_name
)
{
    imageID ID, IDc;
    double dx, r;
    double rms = 0.0;

    double sigmau;
    double sigmav;
    //double sigmaw;

    create_3Dimage_ID(IDout_name, vKsize, 1, 3, &IDc);
    sigmau = sigmawind;
    sigmav = sigmawind;
    //sigmaw = sigmawind;

    // longitudinal (u)
    make_rnd("tmppha0", vKsize, 1, "");
    arith_image_cstmult("tmppha0", 2.0 * PI, "tmppha");
    delete_image_ID("tmppha0", DELETE_IMAGE_ERRMODE_WARNING);

    printf("vK wind outer scale = %f m\n", Lwind);
    printf("pixscale            = %f m\n", pixscale);
    printf("Image size          = %f m\n", vKsize * pixscale);


    create_2Dimage_ID("tmpamp0", vKsize, 1, &ID);
    for(long ii = 0; ii < vKsize; ii++)
    {
        dx = 1.0 * ii - vKsize / 2;
        r = sqrt(dx * dx); // period = (size*pixscale)/r
        // spatial frequency = 2 PI / period
        data.image[ID].array.F[ii] = sqrt(1.0 / pow(1.0 + pow(1.339 * 2.0 * M_PI * r /
                                          (vKsize * pixscale) * Lwind, 2.0), 5.0 / 6.0));
    }


    make_rnd("tmpg", vKsize, 1, "-gauss");
    arith_image_mult("tmpg", "tmpamp0", "tmpamp");
    save_fits("tmpamp0", "vKwind_tmpamp0.fits");
    delete_image_ID("tmpamp0", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("tmpg", DELETE_IMAGE_ERRMODE_WARNING);
    arith_set_pixel("tmpamp", 0.0, vKsize / 2, 0);
    mk_complex_from_amph("tmpamp", "tmppha", "tmpc", 0);
    delete_image_ID("tmpamp", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("tmppha", DELETE_IMAGE_ERRMODE_WARNING);
    permut("tmpc");
    do2dfft("tmpc", "tmpcf");
    delete_image_ID("tmpc", DELETE_IMAGE_ERRMODE_WARNING);
    mk_reim_from_complex("tmpcf", "tmpo1", "tmpo2", 0);
    delete_image_ID("tmpcf", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("tmpo2", DELETE_IMAGE_ERRMODE_WARNING);
    ID = image_ID("tmpo1");
    rms = 0.0;
    for(uint32_t ii = 0; ii < vKsize; ii++)
    {
        rms += data.image[ID].array.F[ii] * data.image[ID].array.F[ii];
    }
    rms = sqrt(rms / vKsize);

    for(uint32_t ii = 0; ii < vKsize; ii++)
    {
        data.image[IDc].array.F[ii] = data.image[ID].array.F[ii] / rms * sigmau;
    }
    delete_image_ID("tmpo1", DELETE_IMAGE_ERRMODE_WARNING);


    // tangential (v)
    make_rnd("tmppha0", vKsize, 1, "");
    arith_image_cstmult("tmppha0", 2.0 * PI, "tmppha");
    delete_image_ID("tmppha0", DELETE_IMAGE_ERRMODE_WARNING);

    create_2Dimage_ID("tmpamp0", vKsize, 1, &ID);
    for(uint32_t ii = 0; ii < vKsize; ii++)
    {
        dx = 1.0 * ii - vKsize / 2;
        r = sqrt(dx * dx);
        data.image[ID].array.F[ii] = sqrt((1.0 + 8.0 / 3.0 * pow(
                                               2.678 * 2.0 * M_PI * r / (vKsize * pixscale) * Lwind,
                                               2.0)) / pow(1.0 + pow(2.678 * 2.0 * M_PI * r / (vKsize * pixscale) * Lwind,
                                                       2.0), 11.0 / 6.0));
    }
    make_rnd("tmpg", vKsize, 1, "-gauss");
    arith_image_mult("tmpg", "tmpamp0", "tmpamp");
    delete_image_ID("tmpamp0", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("tmpg", DELETE_IMAGE_ERRMODE_WARNING);
    arith_set_pixel("tmpamp", 0.0, vKsize / 2, 0);
    mk_complex_from_amph("tmpamp", "tmppha", "tmpc", 0);
    delete_image_ID("tmpamp", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("tmppha", DELETE_IMAGE_ERRMODE_WARNING);
    permut("tmpc");
    do2dfft("tmpc", "tmpcf");
    delete_image_ID("tmpc", DELETE_IMAGE_ERRMODE_WARNING);
    mk_reim_from_complex("tmpcf", "tmpo1", "tmpo2", 0);
    delete_image_ID("tmpcf", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("tmpo2", DELETE_IMAGE_ERRMODE_WARNING);
    ID = image_ID("tmpo1");
    rms = 0.0;
    for(uint32_t ii = 0; ii < vKsize; ii++)
    {
        rms += data.image[ID].array.F[ii] * data.image[ID].array.F[ii];
    }
    rms = sqrt(rms / (vKsize));

    for(uint32_t ii = 0; ii < vKsize; ii++)
    {
        data.image[IDc].array.F[vKsize + ii] = data.image[ID].array.F[ii] / rms *
                                               sigmav;
    }
    delete_image_ID("tmpo1", DELETE_IMAGE_ERRMODE_WARNING);




    // vertical (w)
    make_rnd("tmppha0", vKsize, 1, "");
    arith_image_cstmult("tmppha0", 2.0 * PI, "tmppha");
    delete_image_ID("tmppha0", DELETE_IMAGE_ERRMODE_WARNING);

    create_2Dimage_ID("tmpamp0", vKsize, 1, &ID);
    for(uint32_t ii = 0; ii < vKsize; ii++)
    {
        dx = 1.0 * ii - size / 2;
        r = sqrt(dx * dx);
        data.image[ID].array.F[ii] = sqrt((1.0 + 8.0 / 3.0 * pow(
                                               2.678 * 2.0 * M_PI * r / (vKsize * pixscale) * Lwind,
                                               2.0)) / pow(1.0 + pow(2.678 * 2.0 * M_PI * r / (vKsize * pixscale) * Lwind,
                                                       2.0), 11.0 / 6.0));
    }
    make_rnd("tmpg", vKsize, 1, "-gauss");
    arith_image_mult("tmpg", "tmpamp0", "tmpamp");
    delete_image_ID("tmpamp0", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("tmpg", DELETE_IMAGE_ERRMODE_WARNING);
    arith_set_pixel("tmpamp", 0.0, vKsize / 2, 0);
    mk_complex_from_amph("tmpamp", "tmppha", "tmpc", 0);
    delete_image_ID("tmpamp", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("tmppha", DELETE_IMAGE_ERRMODE_WARNING);
    permut("tmpc");
    do2dfft("tmpc", "tmpcf");
    delete_image_ID("tmpc", DELETE_IMAGE_ERRMODE_WARNING);
    mk_reim_from_complex("tmpcf", "tmpo1", "tmpo2", 0);
    delete_image_ID("tmpcf", DELETE_IMAGE_ERRMODE_WARNING);
    delete_image_ID("tmpo2", DELETE_IMAGE_ERRMODE_WARNING);
    ID = image_ID("tmpo1");
    rms = 0.0;
    for(uint32_t ii = 0; ii < vKsize; ii++)
    {
        rms += data.image[ID].array.F[ii] * data.image[ID].array.F[ii];
    }
    rms = sqrt(rms / vKsize);

    for(uint32_t ii = 0; ii < size; ii++)
    {
        data.image[IDc].array.F[vKsize * 2 + ii] = data.image[ID].array.F[ii] / rms *
                sigmav;
    }
    delete_image_ID("tmpo1", DELETE_IMAGE_ERRMODE_WARNING);




    return(IDc);
}

