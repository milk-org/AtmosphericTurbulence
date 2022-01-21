/**
 * @file    AtmosphericTurbulence_prediction.c
 *
 */

#include <math.h>

#include "CommandLineInterface/CLIcore.h"

#include "COREMOD_iofits/COREMOD_iofits.h"
#include "COREMOD_memory/COREMOD_memory.h"

#include "linopt_imtools/linopt_imtools.h"

#include "AtmosphericTurbulence_conf.h"
#include "ReadConf.h"

#include "COREMOD_tools/COREMOD_tools.h"
#include "statistic/statistic.h"

//
// build full predictor (all pixels of WF)
//
int AtmosphericTurbulence_Build_LinPredictor_Full(
    const char *restrict WFin_name,
    const char *restrict WFmask_name,
    int    PForder,
    float  PFlag,
    double SVDeps,
    double RegLambda)
{
    imageID ID_WFin;
    long    NBmvec; // number of entries in data matrix
    imageID ID_WFmask;

    imageID IDmatA; // data matrix
    long    NBpix;
    long    pix;
    long    PFpix;

    long  mvecsize;
    long *pixarray_x;
    long *pixarray_y;
    long *pixarray_xy;

    imageID IDmatC;

    int  Save = 1;
    char filtname[200];
    char filtfname[200];

    long    k0;
    imageID ID_Pfilt;
    float   alpha;
    imageID IDfiltC;

    int  REG = 0; // 1 if regularization
    long NBmvec1;

    if (RegLambda > 1e-20)
    {
        REG = 1;
    }

    float *valfarray;
    long   ind1;

    ID_WFin         = image_ID(WFin_name);
    uint32_t xsize  = data.image[ID_WFin].md[0].size[0];
    uint32_t ysize  = data.image[ID_WFin].md[0].size[1];
    uint32_t zsize  = data.image[ID_WFin].md[0].size[2];
    uint64_t xysize = xsize * ysize;

    ID_WFmask = image_ID(WFmask_name);
    NBpix     = 0;
    for (uint64_t ii = 0; ii < xsize * ysize; ii++)
        if (data.image[ID_WFmask].array.F[ii] > 0.5)
        {
            NBpix++;
        }
    pixarray_x  = (long *) malloc(sizeof(long) * NBpix);
    pixarray_y  = (long *) malloc(sizeof(long) * NBpix);
    pixarray_xy = (long *) malloc(sizeof(long) * NBpix);

    // PRE_PROCESS WAVEFRONTS : REMOVE PISTON TERM

    {
        double totm = 0.0;
        for (uint32_t ii = 0; ii < xsize; ii++)
            for (uint32_t jj = 0; jj < ysize; jj++)
                if (data.image[ID_WFmask].array.F[jj * xsize + ii] > 0.5)
                {
                    totm += 1.0;
                }

        for (uint32_t kk = 0; kk < zsize; kk++)
        {
            double tot = 0.0;
            for (uint32_t ii = 0; ii < xsize; ii++)
            {
                for (uint32_t jj = 0; jj < ysize; jj++)
                {
                    data.image[ID_WFin]
                        .array.F[kk * xysize + jj * xsize + ii] *=
                        data.image[ID_WFmask].array.F[jj * xsize + ii];
                    tot += data.image[ID_WFin]
                               .array.F[kk * xysize + jj * xsize + ii];
                }
                for (uint32_t ii = 0; ii < xsize; ii++)
                {
                    for (uint32_t jj = 0; jj < ysize; jj++)
                    {
                        if (data.image[ID_WFmask].array.F[jj * xsize + ii] >
                            0.5)
                        {
                            data.image[ID_WFin]
                                .array.F[kk * xysize + jj * xsize + ii] -=
                                tot / totm;
                        }
                    }
                }
            }
        }
    }
    if (Save == 1)
    {
        save_fits(WFin_name, "wfinm.fits");
    }

    // LOAD PIXELS COORDINATES INTO ARRAYS

    NBpix = 0;
    for (uint32_t ii = 0; ii < xsize; ii++)
        for (uint32_t jj = 0; jj < ysize; jj++)
            if (data.image[ID_WFmask].array.F[jj * xsize + ii] > 0.5)
            {
                pixarray_x[NBpix]  = ii;
                pixarray_y[NBpix]  = jj;
                pixarray_xy[NBpix] = jj * xsize + ii;
                NBpix++;
            }
    printf("NBpix = %ld\n", NBpix);

    EXECUTE_SYSTEM_COMMAND("mkdir -p pixfilters");

    // build data matrix
    NBmvec   = zsize - PForder - (int) (PFlag) -1;
    mvecsize = NBpix * PForder;

    if (REG == 0)
    {
        create_2Dimage_ID("PFmatD", NBmvec, mvecsize, &IDmatA);
    }
    else
    {
        create_2Dimage_ID("PFmatD", NBmvec + mvecsize, mvecsize, &IDmatA);
    }

    // each column is a measurement
    // m index is measurement
    // dt*NBpix+pix index is pixel

    // CREATE DATA MATRIX

    if (REG == 0)
    {
        printf("NBmvec   = %ld  -> %ld \n", NBmvec, NBmvec);
        NBmvec1 = NBmvec;
    }
    else
    {
        printf("NBmvec   = %ld  -> %ld \n", NBmvec, NBmvec + mvecsize);
        NBmvec1 = NBmvec + mvecsize;
    }

    printf("mvecsize = %ld  (%d x %ld)\n", mvecsize, PForder, NBpix);

    for (int m = 0; m < NBmvec; m++)
    {
        k0 = m + PForder - 1; // dt=0 index
        for (pix = 0; pix < NBpix; pix++)
            for (int dt = 0; dt < PForder; dt++)
            {
                data.image[IDmatA].array.F[(NBpix * dt + pix) * NBmvec1 + m] =
                    data.image[ID_WFin]
                        .array.F[(k0 - dt) * xysize + pixarray_xy[pix]];
            }
    }
    if (REG == 1)
    {
        for (int m = 0; m < mvecsize; m++)
        {
            //m1 = NBmvec + m;
            data.image[IDmatA].array.F[(m) *NBmvec1 + (NBmvec + m)] = RegLambda;
        }
    }

    if (Save == 1)
    {
        save_fits("PFmatD", "PFmatD.fits");
    }
    list_image_ID();

    printf("Compute reconstruction matrix\n");
    fflush(stdout);

#ifdef HAVE_MAGMA
    CUDACOMP_magma_compute_SVDpseudoInverse("PFmatD",
                                            "PFmatC",
                                            SVDeps,
                                            100000,
                                            "PF_VTmat",
                                            0);
#else
    linopt_compute_SVDpseudoInverse("PFmatD",
                                    "PFmatC",
                                    SVDeps,
                                    100000,
                                    "PF_VTmat",
                                    NULL);
#endif

    if (0)
    {
        save_fits("PFmatD", "test_PFmatD.fits");
        save_fits("PFmatC", "test_PFmatC.fits");
        save_fits("PF_VTmat", "test_PF_VTmat.fits");
#ifdef HAVE_MAGMA
        CUDACOMP_magma_compute_SVDpseudoInverse("PFmatD",
                                                "PFmatC_magma",
                                                SVDeps,
                                                100000,
                                                "PF_VTmat_magma",
                                                0);
#else
        linopt_compute_SVDpseudoInverse("PFmatD",
                                        "PFmatC_magma",
                                        SVDeps,
                                        100000,
                                        "PF_VTmat_magma",
                                        NULL);
#endif

        list_image_ID();
        save_fits("PFmatC_magma", "test_PFmatC_magma.fits");
        save_fits("PF_VTmat_magma", "test_PF_VTmat_magma.fits");
        exit(0);
    }

    if (Save == 1)
    {
        save_fits("PFmatC", "PFmatC.fits");
    }
    IDmatC = image_ID("PFmatC");

    printf("Compute filters\n");
    fflush(stdout);

    create_3Dimage_ID("filtC", NBpix, NBpix, PForder, &IDfiltC);

    valfarray = (float *) malloc(sizeof(float) * NBmvec);

    alpha = PFlag - ((long) PFlag);

    for (PFpix = 0; PFpix < NBpix; PFpix++)
    {
        // PFpix is the pixel for which the filter is created

        sprintf(filtname,
                "PFfilt_%06ld_%03ld_%03ld",
                pixarray_xy[PFpix],
                pixarray_x[PFpix],
                pixarray_y[PFpix]);
        sprintf(filtfname,
                "./pixfilters/PFfilt_%06ld_%03ld_%03ld.fits",
                pixarray_xy[PFpix],
                pixarray_x[PFpix],
                pixarray_y[PFpix]);
        create_3Dimage_ID(filtname, xsize, ysize, PForder, &ID_Pfilt);

        // fill in valfarray

        for (int m = 0; m < NBmvec; m++)
        {
            k0 = m + PForder - 1;
            k0 += (long) PFlag;

            valfarray[m] =
                (1.0 - alpha) *
                    data.image[ID_WFin]
                        .array.F[(k0) *xysize + pixarray_xy[PFpix]] +
                alpha * data.image[ID_WFin]
                            .array.F[(k0 + 1) * xysize + pixarray_xy[PFpix]];
        }

        for (pix = 0; pix < NBpix; pix++)
        {
            for (int dt = 0; dt < PForder; dt++)
            {
                double val = 0.0;
                ind1       = (NBpix * dt + pix) * NBmvec1;
                for (int m = 0; m < NBmvec; m++)
                {
                    val += data.image[IDmatC].array.F[ind1 + m] * valfarray[m];
                }

                data.image[ID_Pfilt].array.F[xysize * dt + pixarray_xy[pix]] =
                    val;
                data.image[IDfiltC]
                    .array.F[dt * NBpix * NBpix + PFpix * NBpix + pix] = val;
            }
        }
        save_fits(filtname, filtfname);
    }

    free(valfarray);
    free(pixarray_x);
    free(pixarray_y);

    return (0);
}

//
// extract translation-invariant kernel from predictive AR filter and expand into individual filters
//
imageID AtmosphericTurbulence_LinPredictor_filt_2DKernelExtract(
    const char *restrict IDfilt_name,
    const char *restrict IDmask_name,
    long krad,
    const char *restrict IDkern_name)
{
    imageID IDkern;
    imageID IDfilt;
    imageID IDmask;
    imageID IDkern_cnt;
    char    filtname[200];
    char    filtfname[200];

    long PForder;
    long PFpix, NBpix, pix;

    long dt;

    imageID     IDfiltC1, IDfiltC1n, IDfiltC2, IDfiltC2n;
    imageID     ID_Pfilt, IDfiltC1cnt;
    long double tmp1, tmp2;
    double      gain;

    IDfilt = image_ID(IDfilt_name);

    IDmask          = image_ID(IDmask_name);
    uint32_t xsize  = data.image[IDmask].md[0].size[0];
    uint32_t ysize  = data.image[IDmask].md[0].size[1];
    uint64_t xysize = xsize * ysize;
    NBpix           = data.image[IDfilt].md[0].size[0];
    PForder         = data.image[IDfilt].md[0].size[2];

    uint32_t *pixarray_x;
    uint32_t *pixarray_y;
    uint64_t *pixarray_xy;

    pixarray_x  = (uint32_t *) malloc(sizeof(uint32_t) * NBpix);
    pixarray_y  = (uint32_t *) malloc(sizeof(uint32_t) * NBpix);
    pixarray_xy = (uint64_t *) malloc(sizeof(uint64_t) * NBpix);

    long NBpix1 = 0;
    {
        for (uint32_t ii = 0; ii < xsize; ii++)
            for (uint32_t jj = 0; jj < ysize; jj++)
                if (data.image[IDmask].array.F[jj * xsize + ii] > 0.5)
                {
                    pixarray_x[NBpix1]  = ii;
                    pixarray_y[NBpix1]  = jj;
                    pixarray_xy[NBpix1] = jj * xsize + ii;
                    NBpix1++;
                }
        printf("NBpix1 = %ld / %ld\n", NBpix1, NBpix);
    }

    uint32_t xksize = 2 * krad + 1;
    uint32_t yksize = 2 * krad + 1;
    create_3Dimage_ID(IDkern_name, xksize, yksize, PForder, &IDkern);
    create_3Dimage_ID("kerncnt", xksize, yksize, PForder, &IDkern_cnt);

    // extract 2D kernel
    for (PFpix = 0; PFpix < NBpix; PFpix++)
    {
        // PFpix is the pixel for which the filter is created

        uint32_t iifilt = pixarray_x[PFpix];
        uint32_t jjfilt = pixarray_y[PFpix];

        for (dt = 0; dt < PForder; dt++)
        {
            for (pix = 0; pix < NBpix; pix++)
            {
                uint32_t ii = pixarray_x[pix];
                uint32_t jj = pixarray_y[pix];

                long dii = ii - iifilt;
                long djj = jj - jjfilt;

                if (dii * dii + djj * djj < krad * krad)
                {
                    data.image[IDkern]
                        .array.F[dt * xksize * yksize + (djj + krad) * xksize +
                                 dii + krad] +=
                        data.image[IDfilt]
                            .array.F[dt * NBpix * NBpix + PFpix * NBpix + pix];
                    data.image[IDkern_cnt]
                        .array.F[dt * xksize * yksize + (djj + krad) * xksize +
                                 dii + krad] += 1.0;
                }
            }
        }
    }

    for (uint32_t ii = 0; ii < xksize; ii++)
        for (uint32_t jj = 0; jj < yksize; jj++)
            for (dt = 0; dt < PForder; dt++)
            {
                data.image[IDkern]
                    .array.F[dt * xksize * yksize + jj * xksize + ii] /=
                    (data.image[IDkern_cnt]
                         .array.F[dt * xksize * yksize + jj * xksize + ii] +
                     1.0e-8);
            }

    // expand 2D kernel into new filter
    create_3Dimage_ID("filtCk", NBpix, NBpix, PForder, &IDfiltC1);
    create_3Dimage_ID("filtCkcnt", NBpix, NBpix, PForder, &IDfiltC1cnt);

    create_2Dimage_ID("filtCkn", xsize, ysize, &IDfiltC1n);
    for (PFpix = 0; PFpix < NBpix; PFpix++)
    {
        tmp1            = 0.0;
        uint32_t iifilt = pixarray_x[PFpix];
        uint32_t jjfilt = pixarray_y[PFpix];
        for (dt = 0; dt < PForder; dt++)
        {
            for (pix = 0; pix < NBpix; pix++)
            {
                uint32_t ii = pixarray_x[pix];
                uint32_t jj = pixarray_y[pix];

                long dii = ii - iifilt;
                long djj = jj - jjfilt;

                if ((dii * dii + djj * djj) < (krad * krad))
                {
                    data.image[IDfiltC1]
                        .array.F[dt * NBpix * NBpix + PFpix * NBpix + pix] =
                        data.image[IDkern]
                            .array.F[dt * xksize * yksize +
                                     (djj + krad) * xksize + dii + krad];
                    data.image[IDfiltC1cnt]
                        .array.F[dt * NBpix * NBpix + PFpix * NBpix + pix] =
                        1.0;
                }
            }
        }
    }

    // remove offset
    create_3Dimage_ID("filtCk2", NBpix, NBpix, PForder, &IDfiltC2);
    create_2Dimage_ID("filtC2n", xsize, ysize, &IDfiltC2n);
    for (dt = 0; dt < PForder; dt++)
    {
        for (PFpix = 0; PFpix < NBpix; PFpix++)
        {
            tmp1 = 0.0;
            tmp2 = 0.0;

            for (pix = 0; pix < NBpix; pix++)
            {
                tmp1 += data.image[IDfiltC1cnt]
                            .array.F[dt * NBpix * NBpix + PFpix * NBpix + pix];
                tmp2 += data.image[IDfiltC1]
                            .array.F[dt * NBpix * NBpix + PFpix * NBpix + pix];
            }
            tmp1 = 1.0 * NBpix - tmp1;

            for (pix = 0; pix < NBpix; pix++)
            {
                data.image[IDfiltC1]
                    .array.F[dt * NBpix * NBpix + PFpix * NBpix + pix] -=
                    (1.0 -
                     data.image[IDfiltC1cnt]
                         .array.F[dt * NBpix * NBpix + PFpix * NBpix + pix]) *
                    (tmp2 / tmp1);

                // non piston-compensated
                data.image[IDfiltC2]
                    .array.F[dt * NBpix * NBpix + PFpix * NBpix + pix] =
                    data.image[IDfiltC1]
                        .array.F[dt * NBpix * NBpix + PFpix * NBpix + pix] +
                    tmp2 / tmp1;
            }
        }
    }

    for (PFpix = 0; PFpix < NBpix; PFpix++)
    {
        tmp1 = 0.0;
        for (dt = 0; dt < PForder; dt++)
            for (pix = 0; pix < NBpix; pix++)
            {
                tmp1 += data.image[IDfiltC2]
                            .array.F[dt * NBpix * NBpix + PFpix * NBpix + pix];
            }
        data.image[IDfiltC2n].array.F[pixarray_xy[PFpix]] = tmp1;
    }

    // rescale edge gain effect
    for (dt = 0; dt < PForder; dt++)
    {
        for (PFpix = 0; PFpix < NBpix; PFpix++)
        {
            if (data.image[IDfiltC2n].array.F[pixarray_xy[PFpix]] > 0.01)
            {
                gain = 1.0 / data.image[IDfiltC2n].array.F[pixarray_xy[PFpix]];
            }
            tmp1 = 0.0;
            tmp2 = 0.0;
            for (pix = 0; pix < NBpix; pix++)
            {
                data.image[IDfiltC2]
                    .array.F[dt * NBpix * NBpix + PFpix * NBpix + pix] *= gain;
                tmp1 += 1.0;
                tmp2 += data.image[IDfiltC2]
                            .array.F[dt * NBpix * NBpix + PFpix * NBpix + pix];
            }

            for (pix = 0; pix < NBpix; pix++)
            {
                data.image[IDfiltC2]
                    .array.F[dt * NBpix * NBpix + PFpix * NBpix + pix] -=
                    (tmp2 / tmp1);
            }
        }
    }

    // expand into individual filters
    for (PFpix = 0; PFpix < NBpix; PFpix++)
    {
        sprintf(filtname,
                "PFfilt_%06lu_%03u_%03u",
                pixarray_xy[PFpix],
                pixarray_x[PFpix],
                pixarray_y[PFpix]);
        sprintf(filtfname,
                "./pixfilters/PFfilt_%06lu_%03u_%03u.fits",
                pixarray_xy[PFpix],
                pixarray_x[PFpix],
                pixarray_y[PFpix]);
        create_3Dimage_ID(filtname, xsize, ysize, PForder, &ID_Pfilt);
        tmp1 = 0.0;
        for (dt = 0; dt < PForder; dt++)
            for (pix = 0; pix < NBpix; pix++)
            {
                data.image[ID_Pfilt].array.F[xysize * dt + pixarray_xy[pix]] =
                    data.image[IDfiltC2]
                        .array.F[dt * NBpix * NBpix + PFpix * NBpix + pix];
                tmp1 += data.image[ID_Pfilt]
                            .array.F[xysize * dt + pixarray_xy[pix]];
            }
        save_fits(filtname, filtfname);
        data.image[IDfiltC1n].array.F[pixarray_xy[PFpix]] = tmp1;
    }

    free(pixarray_x);
    free(pixarray_y);
    free(pixarray_xy);

    return (IDkern);
}

//
// expand into individual filters
// also provides some statistical analysis on filter
//
long AtmosphericTurbulence_LinPredictor_filt_Expand(
    const char *restrict IDfilt_name, const char *restrict IDmask_name)
{
    char filtname[200];
    char filtfname[200];

    imageID IDfilt = image_ID(IDfilt_name);

    imageID  IDmask  = image_ID(IDmask_name);
    uint32_t xsize   = data.image[IDmask].md[0].size[0];
    uint32_t ysize   = data.image[IDmask].md[0].size[1];
    uint64_t xysize  = xsize * ysize;
    uint32_t NBpix   = data.image[IDfilt].md[0].size[0];
    uint32_t PForder = data.image[IDfilt].md[0].size[2];

    uint32_t *pixarray_x  = (uint32_t *) malloc(sizeof(uint32_t) * NBpix);
    uint32_t *pixarray_y  = (uint32_t *) malloc(sizeof(uint32_t) * NBpix);
    uint64_t *pixarray_xy = (uint64_t *) malloc(sizeof(uint64_t) * NBpix);

    uint32_t NBpix1 = 0;
    for (uint32_t ii = 0; ii < xsize; ii++)
        for (uint32_t jj = 0; jj < ysize; jj++)
            if (data.image[IDmask].array.F[jj * xsize + ii] > 0.5)
            {
                pixarray_x[NBpix1]  = ii;
                pixarray_y[NBpix1]  = jj;
                pixarray_xy[NBpix1] = jj * xsize + ii;
                NBpix1++;
            }
    printf("NBpix1 = %u / %u\n", NBpix1, NBpix);

    imageID IDnorm1;
    create_2Dimage_ID("filtmap_norm1", xsize, ysize, &IDnorm1);

    imageID IDnorm2;
    create_2Dimage_ID("filtmap_norm2", xsize, ysize, &IDnorm2);

    imageID IDtau;
    create_2Dimage_ID("filtmap_tau", xsize, ysize, &IDtau);

    // expand into individual filters
    for (uint32_t PFpix = 0; PFpix < NBpix; PFpix++)
    {
        sprintf(filtname,
                "PFfilt_%06lu_%03u_%03u",
                pixarray_xy[PFpix],
                pixarray_x[PFpix],
                pixarray_y[PFpix]);
        sprintf(filtfname,
                "./pixfilters/PFfilt_%06lu_%03u_%03u.fits",
                pixarray_xy[PFpix],
                pixarray_x[PFpix],
                pixarray_y[PFpix]);
        imageID ID_Pfilt;
        create_3Dimage_ID(filtname, xsize, ysize, PForder, &ID_Pfilt);

        double tau   = 0.0; // effective time averaging
        double norm1 = 0.0;
        double norm2 = 0.0;

        for (uint32_t dt = 0; dt < PForder; dt++)
            for (uint32_t pix = 0; pix < NBpix; pix++)
            {
                data.image[ID_Pfilt].array.F[xysize * dt + pixarray_xy[pix]] =
                    data.image[IDfilt]
                        .array.F[dt * NBpix * NBpix + PFpix * NBpix + pix];

                norm1 += fabs(data.image[ID_Pfilt]
                                  .array.F[xysize * dt + pixarray_xy[pix]]);
                norm2 += data.image[ID_Pfilt]
                             .array.F[xysize * dt + pixarray_xy[pix]] *
                         data.image[ID_Pfilt]
                             .array.F[xysize * dt + pixarray_xy[pix]];
                tau += data.image[ID_Pfilt]
                           .array.F[xysize * dt + pixarray_xy[pix]] *
                       data.image[ID_Pfilt]
                           .array.F[xysize * dt + pixarray_xy[pix]] *
                       dt;
            }
        tau /= norm2;
        norm2 = sqrt(norm2);

        data.image[IDnorm1].array.F[pixarray_xy[PFpix]] = norm1;
        data.image[IDnorm2].array.F[pixarray_xy[PFpix]] = norm2;
        data.image[IDtau].array.F[pixarray_xy[PFpix]]   = tau;

        save_fits(filtname, filtfname);
    }

    free(pixarray_x);
    free(pixarray_y);
    free(pixarray_xy);

    return (0);
}

//
// Apply full predictor (all pixels of WF)
//
// MODE = 0: apply pixfilter filters
// MODE = 1: simple WF averaging (temporal only)
//
// outp : prediction
// outf : time-shifted measurement
// outft : actual future value
// reft : true wavefront (no WFS noise)
//
// -> outp_res  (prediction residual)
// -> outf_res  (future measurement residual)
// -> outl_res  (last measurement residual)
//
int AtmosphericTurbulence_Apply_LinPredictor_Full(
    int MODE,
    const char *__restrict WFin_name,
    const char *__restrict WFmask_name,
    int   PForder,
    float PFlag,
    const char *__restrict WFoutp_name,
    const char *__restrict WFoutf_name)
{
    imageID ID_WFmask;
    imageID ID_Pfilt;
    imageID IDoutft; // truth (if exists)
    imageID IDoutp_res;
    imageID IDoutf_res;
    imageID IDoutl_res;

    imageID  ID_WFin = image_ID(WFin_name);
    uint32_t xsize   = data.image[ID_WFin].md[0].size[0];
    uint32_t ysize   = data.image[ID_WFin].md[0].size[1];
    uint32_t zsize   = data.image[ID_WFin].md[0].size[2];
    uint64_t xysize  = xsize * ysize;

    imageID IDoutp;
    create_3Dimage_ID(WFoutp_name, xsize, ysize, zsize, &IDoutp);

    imageID IDoutf;
    create_3Dimage_ID(WFoutf_name, xsize, ysize, zsize, &IDoutf);

    imageID IDreft = image_ID("reft");
    if (IDreft != -1)
    {
        create_3Dimage_ID("outft", xsize, ysize, zsize, &IDoutft);
        create_3Dimage_ID("outp_res", xsize, ysize, zsize, &IDoutp_res);
        create_3Dimage_ID("outf_res", xsize, ysize, zsize, &IDoutf_res);
        create_3Dimage_ID("outl_res", xsize, ysize, zsize, &IDoutl_res);
    }

    ID_WFmask = image_ID(WFmask_name);

    long NBpix = 0;
    {
        for (uint32_t ii = 0; ii < xsize * ysize; ii++)
        {
            if (data.image[ID_WFmask].array.F[ii] > 0.5)
            {
                NBpix++;
            }
        }
    }

    uint32_t *pixarray_x;
    uint32_t *pixarray_y;
    uint64_t *pixarray_xy;

    pixarray_x  = (uint32_t *) malloc(sizeof(uint32_t) * NBpix);
    pixarray_y  = (uint32_t *) malloc(sizeof(uint32_t) * NBpix);
    pixarray_xy = (uint64_t *) malloc(sizeof(uint64_t) * NBpix);

    double totm = 0.0;
    {
        for (uint32_t ii = 0; ii < xsize; ii++)
        {
            for (uint32_t jj = 0; jj < ysize; jj++)
            {
                if (data.image[ID_WFmask].array.F[jj * xsize + ii] > 0.5)
                {
                    totm += 1.0;
                }
            }
        }
    }

    NBpix = 0;
    for (uint32_t ii = 0; ii < xsize; ii++)
    {
        for (uint32_t jj = 0; jj < ysize; jj++)
        {
            if (data.image[ID_WFmask].array.F[jj * xsize + ii] > 0.5)
            {
                pixarray_x[NBpix]  = ii;
                pixarray_y[NBpix]  = jj;
                pixarray_xy[NBpix] = jj * xsize + ii;
                NBpix++;
            }
        }
    }
    printf("NBpix = %ld\n", NBpix);

    double alpha  = PFlag - ((long) PFlag);
    long   PFlagl = (long) PFlag;

    printf("Read and Apply filters\n");
    fflush(stdout);
    for (long PFpix = 0; PFpix < NBpix; PFpix++)
    {
        char filtname[STRINGMAXLEN_IMGNAME];
        WRITE_IMAGENAME(filtname,
                        "PFfilt_%06lu_%03u_%03u",
                        pixarray_xy[PFpix],
                        pixarray_x[PFpix],
                        pixarray_y[PFpix]);

        if (MODE == 0)
        {
            char filtfname[STRINGMAXLEN_FULLFILENAME];
            WRITE_FULLFILENAME(filtfname,
                               "./pixfilters/PFfilt_%06lu_%03u_%03u.fits",
                               pixarray_xy[PFpix],
                               pixarray_x[PFpix],
                               pixarray_y[PFpix]);
            load_fits(filtfname, filtname, LOADFITS_ERRMODE_WARNING, &ID_Pfilt);
        }
        else
        {
            create_3Dimage_ID(filtname, xsize, ysize, PForder, &ID_Pfilt);

            for (long step = 0; step < PForder; step++)
                data.image[ID_Pfilt]
                    .array.F[xysize * step + pixarray_y[PFpix] * xsize +
                             pixarray_x[PFpix]] = 1.0 / PForder;
        }

        for (long kk = PForder; kk < zsize; kk++)
        {
            double valp = 0.0;
            for (long step = 0; step < PForder; step++)
            {
                for (uint64_t ii = 0; ii < xsize * ysize; ii++)
                    valp +=
                        data.image[ID_Pfilt].array.F[xysize * step + ii] *
                        data.image[ID_WFin].array.F[(kk - step) * xysize + ii];
            }

            double valf = 0.0;
            if (kk + PFlag + 1 < zsize)
            {
                valf =
                    (1.0 - alpha) *
                        data.image[ID_WFin].array.F[(kk + PFlagl) * xysize +
                                                    pixarray_xy[PFpix]] +
                    alpha *
                        data.image[ID_WFin].array.F[(kk + PFlagl + 1) * xysize +
                                                    pixarray_xy[PFpix]];
            }

            double valft = 0.0;
            if (kk + PFlag + 1 < zsize)
            {
                valft =
                    (1.0 - alpha) *
                        data.image[IDreft].array.F[(kk + PFlagl) * xysize +
                                                   pixarray_xy[PFpix]] +
                    alpha *
                        data.image[IDreft].array.F[(kk + PFlagl + 1) * xysize +
                                                   pixarray_xy[PFpix]];
            }

            data.image[IDoutp].array.F[kk * xysize + pixarray_xy[PFpix]] = valp;
            data.image[IDoutf].array.F[kk * xysize + pixarray_xy[PFpix]] = valf;

            if (IDreft != -1)
            {
                valft = 0.0;
                if (kk + PFlag + 1 < zsize)
                    valft =
                        (1.0 - alpha) *
                            data.image[IDreft].array.F[(kk + PFlagl) * xysize +
                                                       pixarray_xy[PFpix]] +
                        alpha * data.image[IDreft]
                                    .array.F[(kk + PFlagl + 1) * xysize +
                                             pixarray_xy[PFpix]];
                data.image[IDoutft].array.F[kk * xysize + pixarray_xy[PFpix]] =
                    valft;

                data.image[IDoutp_res]
                    .array.F[kk * xysize + pixarray_xy[PFpix]] = valp - valft;
                data.image[IDoutf_res]
                    .array.F[kk * xysize + pixarray_xy[PFpix]] = valf - valft;
                data.image[IDoutl_res]
                    .array.F[kk * xysize + pixarray_xy[PFpix]] =
                    data.image[ID_WFin]
                        .array.F[kk * xysize + pixarray_xy[PFpix]] -
                    valft;
            }
        }
        delete_image_ID(filtname, DELETE_IMAGE_ERRMODE_WARNING);
    }

    free(pixarray_x);
    free(pixarray_y);
    free(pixarray_xy);

    return (0);
}

// use past and near pixels to predict current pixel value ( single pixel )

int AtmosphericTurbulence_Build_LinPredictor(long   NB_WFstep,
                                             double WFphaNoise,
                                             long   WFPlag,
                                             long   WFP_NBstep,
                                             long   WFP_xyrad,
                                             long   WFPiipix,
                                             long   WFPjjpix,
                                             float  slambdaum)
{
    int GHA = 0;

    imageID IDpha_measured;
    imageID IDpupamp;

    long frame;
    long NBFRAMES;
    //    long WFPxsize, WFPysize;
    double   pha;
    long     tspan;
    imageID  IDpha;
    long     k;
    uint32_t naxes[2];

    long  mvecsize; // measurement vector size
    long *mvecdx;
    long *mvecdy;
    long *mvecdz;
    long  NBmvec; // number of measurement vectors

    imageID IDmatA;
    long    l, m;
    long    k0;
    imageID IDmatC;
    double  val;
    imageID ID_WFPfilt;

    int Save = 1;

    // General Hessian Algorithm
    double GHA_eta = 1e-7; // learning rate

    imageID ID_GHA_x; // input vector [ n ]
    imageID ID_GHA_y; // output vector [ m ]
    imageID ID_GHA_z;
    imageID ID_GHA_UT;   // V [n x n]
    imageID ID_GHA_NT;   // V S+ [m x n]
    imageID ID_GHA_zzT;  // [m x m]
    imageID ID_GHA_V;    // [n x n]
    imageID ID_GHA_sval; // [m]
    imageID ID_GHA_Mest; // [m x n]
    imageID ID_GHA_WFPfilt;

    long        GHA_n; // size of input vector = n
    long        GHA_m; // size of output vector = m
    double      dval, dval0;
    long        ll;
    double      sval;
    long        offset;
    long        GHAiter;
    long        GHA_NBiter = 50;
    long double errval     = 0.0;
    imageID     ID_GHA_matA;
    long double cntval;
    long long   cnt;
    long        k1;

    long   IDpupmask;
    double maxPixSpeed = 0.3; // max wind speed in pix / frame

    double SVDeps = 1e-8;
    long   vID;

    double SLAMBDA;

    SLAMBDA = 1.0e-6 * slambdaum;

    if ((vID = variable_ID("SVDeps")) != -1)
    {
        SVDeps = data.variable[vID].value.f;
        printf("SVDeps = %f\n", SVDeps);
    }

    printf("WFP lag    = %ld\n", WFPlag);
    printf("WFP rad    = %ld\n", WFP_xyrad);
    printf("WFP NBstep = %ld\n", WFP_NBstep);

    printf("NOISE      = %f\n", WFphaNoise);
    fflush(stdout);

    ATMTURBCONF atmturbconf;
    AtmosphericTurbulence_ReadConf("WFsim.conf", &atmturbconf);

    // select center pixel
    long ii0 = WFPiipix;
    long jj0 = WFPjjpix;

    uint32_t WFPxsize = 1 + 2 * WFP_xyrad;
    uint32_t WFPysize = 1 + 2 * WFP_xyrad;
    printf("WFP_xyrad = %ld\n", WFP_xyrad);
    create_3Dimage_ID("WFP_pham",
                      WFPxsize,
                      WFPysize,
                      NB_WFstep,
                      &IDpha_measured);

    IDpupamp = image_ID("ST_pa");
    if (IDpupamp == -1)
    {
        printf("ERROR: pupil amplitude map not loaded");
        exit(0);
    }

    naxes[0] = data.image[IDpupamp].md[0].size[0];
    naxes[1] = data.image[IDpupamp].md[0].size[1];
    NBFRAMES = (long) (atmturbconf.TimeSpanCube / atmturbconf.TimeStep);

    printf("NBFRAMES = %ld\n", NBFRAMES);

    IDpupmask = image_ID("pupmask");
    if (IDpupmask == -1)
    {
        printf("ERROR: pupil mask not loaded");
        exit(0);
    }

    if (GHA == 1)
    {
        // prepare GHA
        // matrix convention
        // km x kn matrix
        // size[0] = km
        // size[1] = kn
        // m<n
        GHA_m = 1;
        GHA_n = WFPxsize * WFPysize * (WFP_NBstep - WFPlag);
        create_2Dimage_ID("GHA_x", GHA_n, 1, &ID_GHA_x); // input vector
        create_2Dimage_ID("GHA_y", GHA_m, 1, &ID_GHA_y); // output vector
        create_2Dimage_ID("GHA_z", GHA_m, 1, &ID_GHA_z); // z = UT y
        create_2Dimage_ID("GHA_UT", GHA_m, GHA_m, &ID_GHA_UT);
        create_2Dimage_ID("GHA_NT", GHA_m, GHA_n, &ID_GHA_NT);   // m x n matrix
        create_2Dimage_ID("GHA_zzT", GHA_m, GHA_m, &ID_GHA_zzT); // m x m matrix

        // initialization: set GHA_UT to Identity square matrix
        for (uint32_t ii = 0; ii < GHA_m; ii++)
            for (uint32_t jj = 0; jj < GHA_m; jj++)
            {
                data.image[ID_GHA_UT].array.F[jj * GHA_m + ii] = 0.0;
            }
        for (uint32_t ii = 0; ii < GHA_m; ii++)
        {
            data.image[ID_GHA_UT].array.F[ii * GHA_m + ii] = 1.0;
        }

        // set NT elements
        for (uint32_t ii = 0; ii < GHA_m; ii++)
            for (uint32_t jj = 0; jj < GHA_n; jj++)
            {
                data.image[ID_GHA_NT].array.F[jj * GHA_m + ii] = 0.0;
            }
        for (uint32_t ii = 0; ii < GHA_m; ii++)
        {
            data.image[ID_GHA_NT].array.F[ii * GHA_m + ii] = 1.0;
        }
        //data.image[ID_GHA_NT].array.F[10] = 1.0;
    }

    tspan = 0;
    k     = 0;
    printf("\n\n");
    cnt    = 0;
    cntval = 0.0;
    while (k < NB_WFstep)
    {
        char fnamepha[STRINGMAXLEN_FILENAME];
        char fnameamp[STRINGMAXLEN_FILENAME];

        printf("\r %4ld/%4ld   tspan = %4ld   ", k, NB_WFstep, tspan);
        fflush(stdout);
        //        printf("%ld/%ld\n", tspan, CONF_NB_TSPAN);

        WRITE_FILENAME(fnamepha,
                       "%s%08ld.%09ld.pha.fits",
                       atmturbconf.SWFfileprefix,
                       tspan,
                       (long) (1.0e12 * SLAMBDA + 0.5));

        WRITE_FILENAME(fnameamp,
                       "%s%08ld.%09ld.amp.fits",
                       atmturbconf.SWFfileprefix,
                       tspan,
                       (long) (1.0e12 * SLAMBDA + 0.5));

        load_fits(fnamepha, "wfpha", LOADFITS_ERRMODE_WARNING, &IDpha);

        for (frame = 0; frame < NBFRAMES; frame++)
        {
            printf("\r      %4ld/%4ld   tspan = %4ld   ", k, NB_WFstep, tspan);
            fflush(stdout);
            cnt    = 0;
            cntval = 0.0;
            if (k < NB_WFstep)
            {
                for (uint32_t ii = 0; ii < WFPxsize; ii++)
                    for (uint32_t jj = 0; jj < WFPysize; jj++)
                    {
                        long ii1 = ii0 + (ii - WFP_xyrad);
                        long jj1 = jj0 + (jj - WFP_xyrad);
                        if ((ii1 > 0) && (ii1 < atmturbconf.WFsize) &&
                            (jj1 > 0) && (jj1 < atmturbconf.WFsize))
                        {
                            pha = data.image[IDpha]
                                      .array.F[frame * naxes[0] * naxes[1] +
                                               jj1 * naxes[0] + ii1];
                        }
                        else
                        {
                            pha = 0;
                        }
                        data.image[IDpha_measured]
                            .array
                            .F[k * WFPxsize * WFPysize + jj * WFPxsize + ii] =
                            gauss() * WFphaNoise + pha;
                        cnt++;
                        cntval += data.image[IDpha_measured]
                                      .array.F[k * WFPxsize * WFPysize +
                                               jj * WFPxsize + ii];
                    }
                //            for(ii=0; ii<WFPxsize*WFPysize; ii++)
                //                  data.image[IDpha_measured].array.F[k*WFPxsize*WFPysize+ii] -= cntval/cnt;
            }
            k++;
        }
        delete_image_ID("wfpha", DELETE_IMAGE_ERRMODE_WARNING);
        tspan++;
    }

    //    for(ii=0; ii<WFPxsize*WFPysize*NB_WFstep; ii++)
    //      data.image[IDpha_measured].array.F[ii] -= cntval/cnt;

    if (Save == 1)
    {
        save_fits("WFP_pham", "WFP_pham.fits");
    }

    mvecsize = WFPxsize * WFPysize * (WFP_NBstep - WFPlag);
    mvecdx   = (long *) malloc(sizeof(long) * mvecsize);
    mvecdy   = (long *) malloc(sizeof(long) * mvecsize);
    mvecdz   = (long *) malloc(sizeof(long) * mvecsize);
    NBmvec   = NB_WFstep - WFP_NBstep;

    // indexes
    // m = measurement number
    // l = pixel index

    l = 0;
    for (k = WFPlag; k < WFP_NBstep; k++)
    {

        for (uint32_t ii = 0; ii < WFPxsize; ii++)
            for (uint32_t jj = 0; jj < WFPysize; jj++)
            {
                mvecdx[l] = ii - WFP_xyrad;
                mvecdy[l] = jj - WFP_xyrad;
                mvecdz[l] = k;

                if ((data.image[IDpupmask]
                         .array.F[(jj0 + mvecdy[l]) * naxes[0] +
                                  (ii0 + mvecdx[l])] > 0.1) &&
                    (sqrt(mvecdx[l] * mvecdx[l] + mvecdy[l] * mvecdy[l]) <
                     2.0 + maxPixSpeed * (k + 1)))
                {
                    l++;
                }
            }
    }

    printf("lmax = %ld / %ld\n", l, mvecsize);
    mvecsize = l;

    create_2Dimage_ID("WFPmatA", NBmvec, mvecsize, &IDmatA);
    // each column is a measurement
    // m index is measurement
    // l index is pixel
    for (m = 0; m < NBmvec; m++)
    {
        k0 = m + WFP_NBstep;
        for (l = 0; l < mvecsize; l++)
        {
            data.image[IDmatA].array.F[l * NBmvec + m] =
                data.image[IDpha_measured]
                    .array.F[(k0 - mvecdz[l]) * WFPxsize * WFPysize +
                             (mvecdy[l] + WFP_xyrad) * WFPxsize +
                             (mvecdx[l] + WFP_xyrad)];
        }
    }

    if (GHA == 1)
    {
        // for GHA: compute differences
        create_2Dimage_ID("GHAmatA", NBmvec, mvecsize, &ID_GHA_matA);
        for (k = 0; k < NB_WFstep; k++)
        {
            k1 = k + 50; // 50 frames offset
            if (k1 > NB_WFstep - 1)
            {
                k1 -= NB_WFstep;
            }
            for (l = 0; l < mvecsize; l++)
            {
                data.image[ID_GHA_matA].array.F[l * NBmvec + k] =
                    data.image[IDmatA].array.F[l * NBmvec + k] -
                    data.image[IDmatA].array.F[l * NBmvec + k1];
            }
        }
    }

    if (Save == 1)
    {
        save_fits("WFPmatA", "WFPmatA.fits");
        if (GHA == 1)
        {
            save_fits("GHAmatA", "GHAmatA.fits");
        }
    }

    if (GHA == 1)
    {
        // run GHA
        printf("\n");
        printf("RUNNING GHA  %ld x %ld  [%ld]... \n", GHA_m, GHA_n, NBmvec);
        fflush(stdout);
        for (GHAiter = 0; GHAiter < GHA_NBiter; GHAiter++)
        {
            errval = 0.0;
            for (k = 0; k < NBmvec; k++)
            {
                //printf("\n %ld\n", k);
                //fflush(stdout);

                // initialize input vector x
                for (uint32_t ii = 0; ii < GHA_n; ii++)
                {
                    data.image[ID_GHA_x].array.F[ii] =
                        data.image[ID_GHA_matA].array.F[ii * NBmvec + k];
                }

                // initialize output vector y
                for (uint32_t ii = 0; ii < GHA_m; ii++)
                {
                    data.image[ID_GHA_y].array.F[ii] =
                        data.image[ID_GHA_x].array.F[60];
                }
                //data.image[IDpha_measured].array.F[(k+WFP_NBstep)*WFPxsize*WFPysize+WFP_xyrad*WFPxsize+WFP_xyrad];
                //data.image[ID_GHA_x].array.F[24];

                // Compute vector z = UT y
                for (uint32_t ii = 0; ii < GHA_m; ii++)
                {
                    data.image[ID_GHA_z].array.F[ii] = 0.0;
                }

                for (uint32_t ii = 0; ii < GHA_m; ii++)
                    for (uint32_t jj = 0; jj < GHA_m; jj++)
                    {
                        data.image[ID_GHA_z].array.F[ii] +=
                            data.image[ID_GHA_UT].array.F[jj * GHA_m + ii] *
                            data.image[ID_GHA_y].array.F[jj];
                    }

                // compute LT[zzT]
                for (uint32_t ii = 0; ii < GHA_m; ii++)
                    for (uint32_t jj = 0; jj < GHA_m; jj++)
                        if (jj <= ii)
                        {
                            data.image[ID_GHA_zzT].array.F[jj * GHA_m + ii] =
                                data.image[ID_GHA_z].array.F[ii] *
                                data.image[ID_GHA_z].array.F[jj];
                        }

                // update UT
                for (uint32_t ii = 0; ii < GHA_m; ii++)
                    for (uint32_t jj = 0; jj < GHA_m; jj++)
                    {
                        dval = 0.0;
                        dval = data.image[ID_GHA_z].array.F[ii] *
                               data.image[ID_GHA_y].array.F[jj]; // z yT
                        dval0 = 0.0;
                        for (ll = 0; ll < GHA_m; ll++)
                        {
                            dval0 +=
                                data.image[ID_GHA_zzT]
                                    .array.F[ll * GHA_m + ii] *
                                data.image[ID_GHA_UT].array.F[jj * GHA_m + ll];
                        }
                        dval -= dval0;

                        data.image[ID_GHA_UT].array.F[jj * GHA_m + ii] +=
                            GHA_eta * dval;
                    }

                // update NT
                errval = 0.0;
                for (uint32_t ii = 0; ii < GHA_m; ii++)
                    for (uint32_t jj = 0; jj < GHA_n; jj++)
                    {
                        dval = 0.0;
                        dval = data.image[ID_GHA_z].array.F[ii] *
                               data.image[ID_GHA_x].array.F[jj]; // z xT
                        dval0 = 0.0;
                        for (ll = 0; ll < GHA_m; ll++)
                        {
                            dval0 +=
                                data.image[ID_GHA_zzT]
                                    .array.F[ll * GHA_m + ii] *
                                data.image[ID_GHA_NT].array.F[jj * GHA_m + ll];
                        }
                        dval -= dval0;

                        errval += dval * dval;

                        data.image[ID_GHA_NT].array.F[jj * GHA_m + ii] +=
                            GHA_eta * dval;
                    }
                //  printf("%05ld   z = %g    U = %g     NT0 = %g\n", k, data.image[ID_GHA_z].array.F[0], data.image[ID_GHA_UT].array.F[0], data.image[ID_GHA_NT].array.F[0]);
            }
            printf("%3ld  errval = %.18lf    %.10f\n",
                   GHAiter,
                   (double) errval,
                   data.image[ID_GHA_NT].array.F[0]);
            fflush(stdout);
        }
        printf("done\n");
        fflush(stdout);

        if (Save == 1)
        {
            save_fits("GHA_UT", "GHA_UT.fits");
            save_fits("GHA_NT", "GHA_NT.fits");
        }

        printf("Computing Matrix M\n");
        fflush(stdout);

        // separate vectors V and singular values
        create_2Dimage_ID("GHA_V", GHA_n, GHA_m, &ID_GHA_V);
        create_2Dimage_ID("GHA_sval", GHA_m, 1, &ID_GHA_sval);
        for (uint32_t jj = 0; jj < GHA_m; jj++)
        {
            val = 0.0;
            for (uint32_t ii = 0; ii < GHA_n; ii++)
            {
                dval = data.image[ID_GHA_NT].array.F[ii * GHA_m + jj];
                val += dval * dval;
                data.image[ID_GHA_V].array.F[jj * GHA_n + ii] = dval;
            }
            val = sqrt(val);
            for (uint32_t ii = 0; ii < GHA_n; ii++)
            {
                data.image[ID_GHA_V].array.F[jj * GHA_n + ii] /= val;
            }
            printf("Singular value %3u = %g\n", jj, 1.0 / val);
            data.image[ID_GHA_sval].array.F[jj] = 1.0 / val;
        }
        save_fits("GHA_V", "GHA_V.fits");

        // compute Mest
        create_2Dimage_ID("GHA_Mest", GHA_m, GHA_n, &ID_GHA_Mest);

        for (uint32_t ll = 0; ll < GHA_m; ll++) // singular value index
        {
            sval = data.image[ID_GHA_sval].array.F[ll];
            for (uint32_t jj = 0; jj < GHA_n; jj++)
                for (uint32_t ii = 0; ii < GHA_m; ii++)
                {
                    data.image[ID_GHA_Mest].array.F[jj * GHA_m + ii] +=
                        data.image[ID_GHA_UT].array.F[ii * GHA_m + ll] * sval *
                        data.image[ID_GHA_V].array.F[ll * GHA_n + jj];
                }
        }

        create_3Dimage_ID("GHA_WFPfilt",
                          WFPxsize,
                          WFPysize,
                          WFP_NBstep,
                          &ID_GHA_WFPfilt);
        offset = WFPxsize * WFPysize * WFPlag;
        for (k = 0; k < WFPxsize * WFPysize * (WFP_NBstep - WFPlag); k++)
        {
            data.image[ID_GHA_WFPfilt].array.F[offset + k] =
                data.image[ID_GHA_Mest].array.F[k];
        }
        save_fits("GHA_WFPfilt", "GHA_WFPfilt.fits");
    }
    //exit(0);

    linopt_compute_SVDpseudoInverse("WFPmatA",
                                    "WFPmatC",
                                    SVDeps,
                                    10000,
                                    "WFP_VTmat",
                                    NULL);
    if (Save == 1)
    {
        save_fits("WFPmatC", "WFPmatC.fits");
    }
    IDmatC = image_ID("WFPmatC");

    create_3Dimage_ID("WFPfilt", WFPxsize, WFPysize, WFP_NBstep, &ID_WFPfilt);

    for (l = 0; l < mvecsize; l++)
    {
        val = 0.0;
        for (m = 0; m < NBmvec; m++)
        {
            val += data.image[IDmatC].array.F[l * NBmvec + m] *
                   data.image[IDpha_measured]
                       .array.F[(m + WFP_NBstep) * WFPxsize * WFPysize +
                                WFP_xyrad * WFPxsize + WFP_xyrad];
        }
        // printf("%5ld  ->  %5ld / %5ld     %5ld / %5ld    %5ld / %5ld\n", l, mvecdz[l], WFP_NBstep, mvecdy[l]+WFP_xyrad, WFPysize, mvecdx[l]+WFP_xyrad, WFPxsize);
        data.image[ID_WFPfilt].array.F[WFPxsize * WFPysize * mvecdz[l] +
                                       WFPxsize * (mvecdy[l] + WFP_xyrad) +
                                       (mvecdx[l] + WFP_xyrad)] = val;
    }

    {
        char fname[STRINGMAXLEN_FILENAME];

        sprintf(fname,
                "WFPfilt_lag%ld_rad%ld_%03ld_%03ld.fits",
                WFPlag,
                WFP_xyrad,
                WFPiipix,
                WFPjjpix);
        save_fits("WFPfilt", fname);
        save_fits("WFPfilt", "WFPfilt.fits");
    }

    for (k = WFPlag; k < WFP_NBstep; k++)
    {
        val = 0.0;
        for (uint64_t ii = 0; ii < WFPxsize * WFPysize; ii++)
        {
            val +=
                data.image[ID_WFPfilt].array.F[WFPxsize * WFPysize * k + ii] *
                data.image[ID_WFPfilt].array.F[WFPxsize * WFPysize * k + ii];
        }
        printf("%5ld  %.10f\n", k, val);
    }
    list_image_ID();

    free(mvecdx);
    free(mvecdy);
    free(mvecdz);

    return 0;
}
