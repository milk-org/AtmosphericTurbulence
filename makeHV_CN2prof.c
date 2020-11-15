/**
 * @file    makeHV_CN2prof.c
 * @brief   Create HV CN2 turbulence profile
 *
 *
 */


#include <math.h>

#include "CommandLineInterface/CLIcore.h"



errno_t AtmosphericTurbulence_makeHV_CN2prof(
    double wspeed,
    double r0,
    double sitealt,
    long NBlayer,
    const char *outfile
)
{
    FILE *fp;
    double h;
    double CN2;
    double hstep, hmax;
    double A, A0;
    double CN2sum = 0.0;
    double r0val;
    double lambda = 0.55e-6;
    double Acoeff = 1.0;
    long iter;
    double Astep;
    long k;
    double l0, L0;
    double *layerarray_h;
    double *layerarray_CN2frac;
    double *layerarray_Wspeed;
    double *layerarray_sigmaWindSpeed;
    double *layerarray_Lwind;

    hmax = 30000.0;
    hstep = 1.0;
    A0 = 1.7e-14;

    Astep = 1.0;
    for(iter = 0; iter < 30; iter++)
    {
        A = A0 * Acoeff;
        CN2sum = 0.0;
        for(h = sitealt; h < hmax; h += hstep)
        {
            CN2 = 5.94e-53 * pow(wspeed / 27.0, 2.0) * pow(h,
                    10.0) * exp(-h / 1000.0) + 2.7e-16 * exp(-h / 1500.0) + A * exp(-
                            (h - sitealt) / 100.0);
            CN2sum += CN2 / hstep;
        }
        r0val = 1.0 / pow(0.423 * pow(2.0 * M_PI / lambda, 2.0) * CN2sum, 3.0 / 5.0);

        //  printf("Acoeff = %12f -> r0 = %10f m -> seeing = %10f arcsec\n", Acoeff, r0val, (lambda/r0val)/M_PI*180.0*3600.0);

        if(r0val > r0)
        {
            Acoeff *= 1.0 + Astep;
        }
        else
        {
            Acoeff /= 1.0 + Astep;
        }
        Astep *= 0.8;
    }


    fp = fopen("conf_turb.txt", "w");
    fprintf(fp, "%f\n", (lambda / r0val) / M_PI * 180.0 * 3600.0);
    fclose(fp);

    layerarray_h = (double *) malloc(sizeof(double) * NBlayer);
    layerarray_CN2frac = (double *) malloc(sizeof(double) * NBlayer);
    layerarray_Wspeed = (double *) malloc(sizeof(double) * NBlayer);
    layerarray_sigmaWindSpeed = (double *) malloc(sizeof(double) * NBlayer);
    layerarray_Lwind = (double *) malloc(sizeof(double) * NBlayer);

    for(k = 0; k < NBlayer; k++)
    {
        layerarray_h[k] = sitealt + pow(1.0 * k / (NBlayer - 1),
                                        2.0) * (hmax - sitealt);
        layerarray_CN2frac[k] = 0.0;
    }

    hstep = 1.0;
    for(h = sitealt; h < hmax; h += hstep)
    {
        CN2 = 5.94e-53 * pow(wspeed / 27.0, 2.0) * pow(h,
                10.0) * exp(-h / 1000.0) + 2.7e-16 * exp(-h / 1500.0) + A * exp(-
                        (h - sitealt) / 100.0);
        k  = (long)(sqrt((h - sitealt) / (hmax - sitealt)) * (1.0 * NBlayer - 1.0) +
                    0.5);
        layerarray_CN2frac[k] += CN2;
    }

    fp = fopen(outfile, "w");
    fprintf(fp,
            "# altitude(m)   relativeCN2     speed(m/s)   direction(rad) outerscale[m] innerscale[m] sigmaWsp[m/s] Lwind[m]\n");
    fprintf(fp, "\n");
    for(k = 0; k < NBlayer; k++)
    {
        layerarray_CN2frac[k] /= CN2sum;

        l0 = 0.008 + 0.072 * pow(layerarray_h[k] / 20000.0, 1.6);

        if(layerarray_h[k] < 14000.0)
        {
            L0 = pow(10.0, 2.0 - 0.9 * (layerarray_h[k] / 14000.0));
        }
        else
        {
            L0 = pow(10.0, 1.1 + 0.3 * (layerarray_h[k] - 14000.0) / 6000.0);
        }

        layerarray_Wspeed[k] = wspeed * (0.3 + 0.8 * sqrt(1.0 * k / (1.0 + NBlayer)));
        layerarray_sigmaWindSpeed[k] = 0.1 * layerarray_Wspeed[k];
        layerarray_Lwind[k] = 500.0;

        fprintf(fp, "%12f  %12f  %12f  %12f  %12f  %12f  %12f  %12f\n", layerarray_h[k],
                layerarray_CN2frac[k], layerarray_Wspeed[k], 2.0 * M_PI * k / (1.0 + NBlayer),
                L0, l0, layerarray_sigmaWindSpeed[k], layerarray_Lwind[k]);
    }
    fclose(fp);

    free(layerarray_h);
    free(layerarray_CN2frac);
    free(layerarray_Wspeed);
    free(layerarray_sigmaWindSpeed);
    free(layerarray_Lwind);

    return RETURN_SUCCESS;
}


