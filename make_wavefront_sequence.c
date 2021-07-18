/**
 * @file    make_wavefront_sequence.c
 * @brief   Create wavefront sequence
 *
 *
 */

#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <assert.h>


#include "CommandLineInterface/CLIcore.h"
#include "COREMOD_memory/COREMOD_memory.h"
#include "COREMOD_iofits/COREMOD_iofits.h"

#include "WFpropagate/WFpropagate.h"
#include "WFpropagate/Fresnel_propagate.h"

#include "AtmosphereModel/AtmosphereModel.h"
#include "AtmosphereModel/AtmosphereModel_Create_from_CONF.h"
#include "AtmosphereModel/AirMixture_ria.h"
#include "AtmosphereModel/AtmosphereModel_stdAtmModel_ria.h"

#include "AtmosphericTurbulence_conf.h"
#include "ReadConf.h"

#include "Z_Air.h"
#include "make_master_turbulence_screen.h"
#include "make_AtmosphericTurbulence_vonKarmanWind.h"



static char ATMCONFFILE[STRINGMAXLEN_FILENAME] = "WFsim.conf";

//float SiteLat;
//float SiteLong;
//float SiteAlt;



// constants
double C_me = 9.10938291e-31; // electon mass [kg]
double C_e0 = 8.854187817620e-12; // Vacuum permittivity [F.m-1]
double C_Na = 6.0221413e23; // Avogadro number
double C_e = 1.60217657e-19; // electron charge [C]
double C_ls = 2.686777447e25; // Loschmidt constant



double rhocoeff = 1.0;








// ==========================================
// Forward declaration(s)
// ==========================================

errno_t AtmosphericTurbulence_make_wavefront_sequence(
    float  slambdaum,
    long   WFprecision,
    int    compmode
);


// ==========================================
// Command line interface wrapper function(s)
// ==========================================


static errno_t AtmosphericTurbulence_make_wavefront_sequence__cli()
{
    if(0
            + CLI_checkarg(1, CLIARG_FLOAT)
            + CLI_checkarg(2, CLIARG_LONG)
            + CLI_checkarg(3, CLIARG_LONG)
            == 0)
    {
        // If arguments meet requirements, command is executed
        //
        AtmosphericTurbulence_make_wavefront_sequence(
            data.cmdargtoken[1].val.numf,
            data.cmdargtoken[2].val.numl,
            data.cmdargtoken[3].val.numl);

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

errno_t AtmosphericTurbulence_make_wavefront_sequence_addCLIcmd()
{

    RegisterCLIcommand(
        "mkwfseq",
        __FILE__,
        AtmosphericTurbulence_make_wavefront_sequence__cli,
        "make atmospheric turbulence wavefront sequence",
        "<wavelength [nm]> <precision 0=single, 1=double> <computation mode>",
        "mkwfseq 1650.0 1 1",
        "AtmosphericTurbulence_make_wavefront_sequence(float slambdaum, long WFprecision, int compmode)");

    return RETURN_SUCCESS;
}















// compmode = 0 : compute atmosphere model only, no turbulence
// compmode = 1 : full computation

errno_t AtmosphericTurbulence_make_wavefront_sequence(
    float slambdaum,
    long WFprecision,
    int compmode
)
{
    long naxes_MASTER[2];

    //long master;
    long NBMASTERS;
    imageID *ID_TM;
    imageID *ID_TML;
    imageID IDout_array_pha;
    imageID IDout_array_amp;
    complex_float *array;
    complex_double *array_double;

    // phase only
    imageID ID_array1;
    imageID ID_sarray1;

    // phase + amplitude
    imageID ID_array2;
    imageID ID_sarray2;


    double SLAMBDA; // wagvelength [um]

    imageID IDout_sarray_pha;
    imageID IDout_sarray_amp;
    complex_float *sarray;
    complex_double *sarray_double;

    double Nlambda, Nslambda, l;


    // cone effect wavefront
    /*
        int make_cwavefront = 0;
        char CWF_FILE_PREFIX[100];
        complex_float *carray;
    */

    long NBLAYERS; /* number of layers */
    double *LAYER_ALT;
    double *LAYER_CN2;
    double *LAYER_SPD;
    double *LAYER_DIR;
    double *LAYER_OUTERSCALE;
    double *LAYER_INNERSCALE;
    double *LAYER_SIGMAWSPEED;
    double *LAYER_LWIND;

    uint32_t *naxes;
    uint32_t *naxesout;
    double *xpos;
    double *ypos;
    double *xpos0;
    double *ypos0;
    long *xposfcnt; // counter to keep track of modulo
    long *yposfcnt;
    double *vxpix;
    double *vypix;

    long NBFRAMES;

    //double coeff = 0.0;



    double *SLAYER_ALT;
    double *SLAYER_CN2;
    long NBSLAYERS;

    int pfactor = 1;
    //   float OuterScale;


    long start_cubeindex = 0;

    imageID IDshmpha, IDshmamp;
    imageID IDshmspha, IDshmsamp;

    // timing
    struct timespec tnow;
    double tnowdouble;



    // optimal phase unwrapping
    imageID IDpeakpha_re, IDpeakpha_im;
    imageID IDpeakpha_re_bin, IDpeakpha_im_bin, IDpeakpha_bin, IDpeakpha_bin_ch;
    uint32_t xsizepeakpha;
    uint32_t ysizepeakpha;


    double pcoeff2 = 0.15;





    int BICUBIC = 1; // 0 if bilinear


    FILE *fpxypos;










    SLAMBDA = 1.0e-6 * slambdaum;


    naxes = (uint32_t *) malloc(sizeof(uint32_t) * 3);
    naxesout = (uint32_t *) malloc(sizeof(uint32_t) * 3);


    ID_sarray1 = -1;
    sarray = NULL;
    sarray_double = NULL;


    printf("Making the wavefront series...\n");
    fflush(stdout);


    ATMTURBCONF atmturbconf;

    AtmosphericTurbulence_ReadConf( ATMCONFFILE, &atmturbconf );


    naxesout[0] = atmturbconf.WFsize;
    naxesout[1] = atmturbconf.WFsize;
    naxes[0] = atmturbconf.WFrawSize;
    naxes[1] = atmturbconf.WFrawSize;




    DEBUG_TRACEPOINT("Creating atmosphere model");

    // create atm.txt file with concentrations as function of altitude
    // load RIA (Refractive index and Absorption) files if available
    ATMOSPHERE_MODEL atm = AtmosphereModel_Create_from_CONF(ATMCONFFILE,
                           slambdaum * 1e-6);

    DEBUG_TRACEPOINT("Atmosphere model completed");




    if(0)
    {
        double dens_N2, dens_O2, dens_Ar, dens_H2O, dens_CO2, dens_Ne, dens_He,
               dens_CH4, dens_Kr, dens_H2, dens_O3, dens_N, dens_O, dens_H; // [cm-3]
        double lambda;
        double xN2, xO2, xAr, xH2O, xCO2, xNe, xHe, xCH4, xKr, xH2, xO3, xN, xO, xH;
        double P, T, TC, Pw, CO2ppm, denstot, xtot;
        double LoschmidtConstant =  2.6867805e25;

        double RH;
        double Z, Z0;


        // SOME TESTING
        //  AtmosphereModel_RefractionPath(1.5, zenithangle, 0); //79.999/180.0*M_PI);//zenithangle);


        /*   fp = fopen("Rprof.txt", "w");
           for(h=0; h<100000.0; h+=10.0)
               fprintf(fp, "%8g %.16f %.16f\n", h, AtmosphereModel_stdAtmModel_N(atm, h, 0.6e-6, 0), AtmosphereModel_stdAtmModel_N(atm, h, 1.6e-6, 0));
           fclose(fp);
        */



        // TESTING VALUES IN CIDDOR 1996


        if(0) // dry, table 1, 10 C, 100kPa -> Ciddor value = 1.000277747
        {
            P = 100000.0;
            Pw = 0.0;
            CO2ppm = 450.0;
            TC = 10.0;
        }
        else
        {
            // wet case #1
            P = 102993.0;
            Pw = 641.0;
            CO2ppm = 450.0;
            TC = 19.173;

            if(1)
            {
                // wet case #2
                P = 103006.0;
                Pw = 642.0;
                CO2ppm = 440.0;
                TC = 19.173;
            }
        }


        T = TC + 273.15;

        // dry air composition (Harisson 1965)
        dens_N2 = 780840e-6 * LoschmidtConstant * 1e-6;
        dens_O2 = 209844e-6 * LoschmidtConstant * 1e-6; // assuming CO2 -> O2
        dens_Ar = 9340e-6 * LoschmidtConstant * 1e-6;
        dens_Ne = 18.18e-6 * LoschmidtConstant * 1e-6;
        dens_He = 5.24e-6 * LoschmidtConstant * 1e-6;
        dens_CH4 = 1.774e-6 * LoschmidtConstant * 1e-6;
        dens_Kr = 1.14e-6 * LoschmidtConstant * 1e-6;
        dens_H2 = 0.56e-6 * LoschmidtConstant * 1e-6;
        dens_O3 = 0.0;
        dens_N = 0.0;
        dens_O = 0.0;
        dens_H = 0.0;
        dens_O2 -= CO2ppm * 1e-6 * LoschmidtConstant * 1e-6;
        dens_CO2 = CO2ppm * 1e-6 * LoschmidtConstant * 1e-6;
        denstot = dens_N2 + dens_O2 + dens_Ar + dens_Ne + dens_He + dens_CH4 + dens_Kr +
                  dens_H2 + dens_O3 + dens_N + dens_O + dens_H + dens_CO2;
        dens_H2O = denstot * (Pw / P) / (1.0 - (Pw / P)); // - 1e-6*CO2ppm);
        //dens_CO2 = CO2ppm*1e-6 * LoschmidtConstant*1e-6;


        //dens_CO2 = denstot * CO2ppm*1e-6 / (1.0 - (Pw/P) - 1e-6*CO2ppm);

        denstot = dens_N2 + dens_O2 + dens_Ar + dens_Ne + dens_He + dens_CH4 + dens_Kr +
                  dens_H2 + dens_O3 + dens_N + dens_O + dens_H + dens_H2O + dens_CO2;
        xN2 = dens_N2 / denstot;
        xO2 = dens_O2 / denstot;
        xAr = dens_Ar / denstot;
        xNe = dens_Ne / denstot;
        xHe = dens_H2 / denstot;
        xCH4 = dens_CH4 / denstot;
        xKr = dens_Kr / denstot;
        xH2 = dens_H2 / denstot;
        xO3 = dens_O3 / denstot;
        xN = dens_N / denstot;
        xO = dens_O / denstot;
        xH = dens_H / denstot;

        xCO2 = dens_CO2 / denstot;
        xH2O = dens_H2O / denstot;

        xtot = xN2 + xO2 + xAr + xNe + xHe + xCH4 + xKr + xH2 + xO3 + xN + xO + xH +
               xCO2 + xH2O;
        xN2 /= xtot;
        xO2 /= xtot;
        xAr /= xtot;
        xNe /= xtot;
        xHe /= xtot;
        xCH4 /= xtot;
        xKr /= xtot;
        xH2 /= xtot;
        xO3 /= xtot;
        xN /= xtot;
        xO /= xtot;
        xH /= xtot;
        xCO2 /= xtot;
        xH2O /= xtot;

        printf("CO2 ppm = %.6f\n", xCO2 * 1e6);
        printf("H2O ppm = %.6f (goal: %.6f)\n", xH2O * 1e6, Pw / P);

        RH = 0.0;
        Z0 = Z_Air(101325.0, 0.0, 0.0);
        Z = Z_Air(P, TC, RH); //
        //LL *= rhocoeff;
        printf("rhocoeff = %f %f\n", 1.0 / rhocoeff, P / 101325.0 * (273.15 / T));
        printf("Z = %f\n", Z);
        printf("Z0 = %f\n", Z0);


        //denstot = 1.0/rhocoeff * LoschmidtConstant*1e-6;
        denstot = (1.0 / Z) * P / 101325.0 * (273.15 / T) * LoschmidtConstant *
                  1e-6; // approximation - does not take into account CHANGE of Z with pressure, temperature

        dens_N2 = xN2 * denstot;
        dens_O2 = xO2 * denstot;
        dens_Ar = xAr * denstot;
        dens_Ne = xNe * denstot;
        dens_He = xHe * denstot;
        dens_CH4 = xCH4 * denstot;
        dens_Kr = xKr * denstot;
        dens_H2 = xH2 * denstot;
        dens_O3 = xO3 * denstot;
        dens_N = xN * denstot;
        dens_O = xO * denstot;
        dens_H = xH * denstot;
        dens_CO2 = xCO2 * denstot;
        dens_H2O = xH2O * denstot;

        double *densarray;
        densarray = (double *) malloc(sizeof(double) * atm.speciesRIA.NBspecies);
        densarray[speciesN2]  = dens_N2;
        densarray[speciesO2]  = dens_O2;
        densarray[speciesAr]  = dens_Ar;
        densarray[speciesNe]  = dens_Ne;
        densarray[speciesHe]  = dens_He;
        densarray[speciesCH4] = dens_CH4;
        densarray[speciesKr]  = dens_Kr;
        densarray[speciesH2]  = dens_H2;
        densarray[speciesO3]  = dens_O3;
        densarray[speciesN]   = dens_N;
        densarray[speciesO]   = dens_O;
        densarray[speciesH]   = dens_H;
        densarray[speciesCO2] = dens_CO2;
        densarray[speciesH2O] = dens_H2O;


        printf("denstot = %g part/cm3      %f x LScst\n", denstot,
               denstot * 1.0e6 / LoschmidtConstant);
        lambda = 0.633e-6;
        printf("633nm test values\n");
        RIAvalue riav = AirMixture_ria(atm.speciesRIA, lambda, densarray);
        printf("HARISSON MODEL:    n = %.12g\n", riav.rindex);

        free(densarray);

        riav = AtmosphereModel_stdAtmModel_ria(atm, 0.0, lambda, 1);
        printf("STD MODEL, 1 atm : n = %.12g\n", riav.rindex);

        {
            FILE *fpRlambda;

            fpRlambda = fopen("Rlambda.txt", "w");
            for(l = 1e-6; l < 5e-6; l *= 1.0 + 1e-5)
            {
                RIAvalue riav0 =  AtmosphereModel_stdAtmModel_ria(atm, 10, l, 0);
                RIAvalue riav1 =  AtmosphereModel_stdAtmModel_ria(atm, 1000, l, 0);
                RIAvalue riav2 =  AtmosphereModel_stdAtmModel_ria(atm, 4200, l, 0);
                fprintf(fpRlambda, "%.16f %.16f %.16f %.16f\n", l, riav0.rindex, riav1.rindex,
                        riav2.rindex);
            }
            fclose(fpRlambda);
        }


        {
            FILE *fpAtmRefrac;
            fpAtmRefrac = fopen("AtmRefrac.txt", "w");
            for(l = 0.4e-6; l < 2.0e-6; l += 0.01e-6)
            {
                RIAvalue riav0 =  AtmosphereModel_stdAtmModel_ria(atm, atm.SiteAlt, l, 0);
                fprintf(fpAtmRefrac, "%.16f %.16f\n", l, asin(sin(atmturbconf.zenithangle) / riav0.rindex));
            }
            fclose(fpAtmRefrac);
        }
    }













    if(compmode == 0)
    {
        free(naxes);
        free(naxesout);
        return 0;
    }

    printf("zenithangle = %f  alt = %f\n", atmturbconf.zenithangle, atm.SiteAlt);



    // for(Temp=173.0; Temp<373.0; Temp+=5.0)
    //   printf("T= %lf K    Ps(H2O) [Pa] = %g\n", Temp, AtmosphereModel_H2O_Saturation(Temp));



    //    Scoeff = LAMBDA/SLAMBDA;
    //   Nlambda = 0.0000834213+0.0240603/(130.0-1.0/pow(atmturbconf.lambda*1000000.0,2.0))+0.00015997/(38.9-1.0/pow(atmturbconf.lambda*1000000.0,2.0));
    //   Nslambda = 0.0000834213+0.0240603/(130.0-1.0/pow(SLAMBDA*1000000.0,2.0))+0.00015997/(38.9-1.0/pow(SLAMBDA*1000000.0,2.0));

    //printf("method 1 : %f %f\n", Nlambda, Nslambda);


    double Scoeff;
    {
        // compute Scoeff
        RIAvalue riav;

        riav = AtmosphereModel_stdAtmModel_ria(atm, 0.0, atmturbconf.lambda, 0);
        Nlambda = 1.0 - riav.rindex;
        riav = AtmosphereModel_stdAtmModel_ria(atm, 0.0, SLAMBDA, 0);
        Nslambda = 1.0 - riav.rindex;

        // multiplicative coefficient to go from reference lambda phase to science lambda phase
        Scoeff =  atmturbconf.lambda / SLAMBDA * Nslambda / Nlambda;
    }




    //  printf("Scoeff is %f (%f)\n",Scoeff,Nslambda/Nlambda);
    //   fflush(stdout);


    // printf("Zenith angle = %f rad\n", zenithangle);



    //  pfactor = naxes[0]/CONF_WFsize;
    //pfactor = 1;

    /*
        contraction_factor = 4;
        if(naxes[0] / CONF_WFsize == 1)
        {
            contraction_factor = 0;
        }
        if(naxes[0] / CONF_WFsize == 2)
        {
            contraction_factor = 1;
        }
        if(naxes[0] / CONF_WFsize == 4)
        {
            contraction_factor = 2;
        }
        if(naxes[0] / CONF_WFsize == 8)
        {
            contraction_factor = 3;
        }

        if(contraction_factor == 4)
        {
            printf("ERROR: unknown contraction factor\n");
            fflush(stdout);
            exit(0);
        }
    */


    /*  ID=image_ID("ST_pa");*/



    NBFRAMES = (long)(1.0 * atmturbconf.TimeSpanCube / atmturbconf.TimeStep + 0.5);
    printf("%.16f  %.16f  ->  %ld\n", atmturbconf.TimeSpanCube, atmturbconf.TimeStep, NBFRAMES);



    naxes[2] = NBFRAMES;
    naxesout[2] = NBFRAMES;

    printf("Allocating memory...\n");
    fflush(stdout);

    // OUTPUT ARRAYS
    if(WFprecision == 0)   // single precision
    {
        create_3Dimage_ID("outarraypha", naxesout[0], naxesout[1],
                          naxesout[2], &IDout_array_pha);
        create_3Dimage_ID("outarrayamp", naxesout[0], naxesout[1],
                          naxesout[2], &IDout_array_amp);
        create_3Dimage_ID("outsarraypha", naxesout[0], naxesout[1],
                          naxesout[2], &IDout_sarray_pha);
        create_3Dimage_ID("outsarrayamp", naxesout[0], naxesout[1],
                          naxesout[2], &IDout_sarray_amp);
        for(uint64_t ii = 0; ii < naxesout[0]*naxesout[1]*naxesout[2]; ii++)
        {
            data.image[IDout_array_amp].array.F[ii] = 1.0;
            data.image[IDout_array_pha].array.F[ii] = 0.0;

            data.image[IDout_sarray_amp].array.F[ii] = 1.0;
            data.image[IDout_sarray_pha].array.F[ii] = 0.0;
        }
    }
    else // double precision
    {
        create_3Dimage_ID_double("outarraypha", naxesout[0],
                                 naxesout[1], naxesout[2], &IDout_array_pha);
        create_3Dimage_ID_double("outarrayamp", naxesout[0],
                                 naxesout[1], naxesout[2], &IDout_array_amp);
        create_3Dimage_ID_double("outsarraypha", naxesout[0],
                                 naxesout[1], naxesout[2], &IDout_sarray_pha);
        create_3Dimage_ID_double("outsarrayamp", naxesout[0],
                                 naxesout[1], naxesout[2], &IDout_sarray_amp);
        for(uint64_t ii = 0; ii < naxesout[0]*naxesout[1]*naxesout[2]; ii++)
        {
            data.image[IDout_array_amp].array.D[ii] = 1.0;
            data.image[IDout_array_pha].array.D[ii] = 0.0;

            data.image[IDout_sarray_amp].array.D[ii] = 1.0;
            data.image[IDout_sarray_pha].array.D[ii] = 0.0;
        }
    }



    /*
     * =================================================================
     * Look for existing master turbulence files
     * =================================================================
     */

    naxes_MASTER[0] = atmturbconf.MasterSize;
    naxes_MASTER[1] = atmturbconf.MasterSize;

    {
        long master = 0;
        int stop = 1;
        while(stop)
        {
            char fnamemasterwf[STRINGMAXLEN_FILENAME];
            if(WFprecision == 0)
            {
                WRITE_FILENAME(fnamemasterwf, "t%03ld_%ld_f.fits", master, (long) atmturbconf.MasterSize);
            }
            else
            {
                WRITE_FILENAME(fnamemasterwf, "t%03ld_%ld_d.fits", master, (long) atmturbconf.MasterSize);
            }
            if(!file_exists(fnamemasterwf))
            {
                stop = 0;
            }
            else
            {
                master++;
            }
        }

        NBMASTERS = master;
        printf("%ld turbulence master files found\n", NBMASTERS);
        fflush(stdout);
    }





    /*
     * =================================================================
     * Load turbulence profile for disk
     * =================================================================
     */
    {
        char line[2000];
        char word[1000];

        FILE *fpturbprof;

        if((fpturbprof = fopen(atmturbconf.turbulenceprof_fname, "r")) == NULL)
        {
            printf("Cannot open turbulence profile file \"%s\"\n",
                   atmturbconf.turbulenceprof_fname);
            exit(1);
        }
        NBLAYERS = 0;
        while(fgets(line, 2000, fpturbprof) != NULL)
        {
            sscanf(line, "%s", word);
            if(isdigit(word[0]))
            {
                NBLAYERS += 1;
            }
        }
        fclose(fpturbprof);

        LAYER_ALT = (double *) malloc(NBLAYERS * sizeof(double));
        LAYER_CN2 = (double *) malloc(NBLAYERS * sizeof(double));
        LAYER_SPD = (double *) malloc(NBLAYERS * sizeof(double));
        LAYER_DIR = (double *) malloc(NBLAYERS * sizeof(double));
        LAYER_OUTERSCALE = (double *) malloc(NBLAYERS * sizeof(double));
        LAYER_INNERSCALE = (double *) malloc(NBLAYERS * sizeof(double));
        LAYER_SIGMAWSPEED = (double *) malloc(NBLAYERS * sizeof(double));
        LAYER_LWIND = (double *) malloc(NBLAYERS * sizeof(double));

        if((fpturbprof = fopen(atmturbconf.turbulenceprof_fname, "r")) == NULL)
        {
            printf("Cannot open turbulence profile file \"%s\"\n",
                   atmturbconf.turbulenceprof_fname);
            exit(1);
        }
        long layer = 0;
        while(fgets(line, 2000, fpturbprof) != NULL)
        {
            sscanf(line, "%s", word);
            if(isdigit(word[0]))
            {
                double fl1, fl2, fl3, fl4, fl5, fl6, fl7, fl8;

                sscanf(line, "%lf %lf %lf %lf %lf %lf %lf %lf", &fl1, &fl2, &fl3, &fl4, &fl5,
                       &fl6, &fl7, &fl8);
                if(fl1 > atm.SiteAlt - 0.1)
                {
                    LAYER_ALT[layer] = fl1;
                    LAYER_CN2[layer] = fl2;
                    LAYER_SPD[layer] = fl3;
                    LAYER_DIR[layer] = fl4;
                    LAYER_OUTERSCALE[layer] = fl5;
                    LAYER_INNERSCALE[layer] = fl6;
                    LAYER_SIGMAWSPEED[layer] = fl7;
                    LAYER_LWIND[layer] = fl8;
                    layer += 1;
                }
            }
        }
        fclose(fpturbprof);
    }


    {
        /* CN2 normalisation for 1024x1024 -> 256x256*/
        /* S<0.7 : x=(S-0.1)/0.06*2
        S>0.7 : x=(S-0.38)/0.033*2 */
        /*  if(seeing<0.7)
        CN2total = (seeing-0.1)/0.03;
        else
        CN2total = (seeing-0.38)/0.0165;*/
        /* S  = sqrt(CN2/16.666)*0.443 = 0.1085 sqrt(CN2) */

        double CN2total;
        double tmp;

        CN2total = 1.0; //84.926*seeing*seeing;

        tmp = 0;
        for(long layer = 0; layer < NBLAYERS; layer++)
        {
            tmp += LAYER_CN2[layer];
        }
        for(long layer = 0; layer < NBLAYERS; layer++)
        {
            LAYER_CN2[layer] *= CN2total / tmp;
        }


        for(long layer = 0; layer < NBLAYERS; layer++)
        {
            printf("Turbulence layer %ld : alt = %f m   CN2 = %f   V = %f m/s   Angle = %f rad   outerscale = %f m    innerscale = %f m  sigmaWind = %f m/s  Lwind = %f m\n",
                   layer, LAYER_ALT[layer], LAYER_CN2[layer], LAYER_SPD[layer], LAYER_DIR[layer],
                   LAYER_OUTERSCALE[layer], LAYER_INNERSCALE[layer], LAYER_SIGMAWSPEED[layer],
                   LAYER_LWIND[layer]);
        }
    }




    SLAYER_ALT = (double *) malloc(NBLAYERS * sizeof(double));
    SLAYER_CN2 = (double *) malloc(NBLAYERS * sizeof(double));
    NBSLAYERS = NBLAYERS;
    for(long layer = 0; layer < NBSLAYERS; layer++)
    {
        SLAYER_ALT[layer] = LAYER_ALT[layer];
        SLAYER_CN2[layer] = LAYER_CN2[layer];
    }


    /// temporary arrays for phase unwrapping
    if(WFprecision == 0)   // single precision
    {
        create_2Dimage_ID("peakphare", naxesout[0], naxesout[1], &IDpeakpha_re);
        create_2Dimage_ID("peakphaim", naxesout[0], naxesout[1], &IDpeakpha_im);
        //IDpeakpha = create_2Dimage_ID("peakpha", naxesout[0], naxesout[1]);

        xsizepeakpha = (uint32_t)(naxesout[0] / 20);
        ysizepeakpha = (uint32_t)(naxesout[1] / 20);
        create_2Dimage_ID("peakphare_bin", xsizepeakpha,
                          ysizepeakpha, &IDpeakpha_re_bin);
        create_2Dimage_ID("peakphaim_bin", xsizepeakpha,
                          ysizepeakpha, &IDpeakpha_im_bin);
        create_2Dimage_ID("peakpha_bin", xsizepeakpha, ysizepeakpha, &IDpeakpha_bin);
        create_2Dimage_ID("peakpha_bin_ch", xsizepeakpha,
                          ysizepeakpha, &IDpeakpha_bin_ch);
    }
    else
    {
        create_2Dimage_ID_double("peakphare", naxesout[0], naxesout[1], &IDpeakpha_re);
        create_2Dimage_ID_double("peakphaim", naxesout[0], naxesout[1], &IDpeakpha_im);
        //IDpeakpha = create_2Dimage_ID_double("peakpha", naxesout[0], naxesout[1]);

        xsizepeakpha = (uint32_t)(naxesout[0] / 20);
        ysizepeakpha = (uint32_t)(naxesout[1] / 20);
        create_2Dimage_ID_double("peakphare_bin", xsizepeakpha,
                                 ysizepeakpha, &IDpeakpha_re_bin);
        create_2Dimage_ID_double("peakphaim_bin", xsizepeakpha,
                                 ysizepeakpha, &IDpeakpha_im_bin);
        create_2Dimage_ID_double("peakpha_bin", xsizepeakpha,
                                 ysizepeakpha, &IDpeakpha_bin);
        create_2Dimage_ID_double("peakpha_bin_ch", xsizepeakpha,
                                 ysizepeakpha, &IDpeakpha_bin_ch);
    }






    /*
     * =================================================================
     * Group layers that are close in altitude
     * =================================================================
     */
    long *super_layer_index;
    {

        int OK = 0;
        while(OK == 0)
        {
            double minaltd;
            long layerindex = 0;

            printf("--------------------\n");
            for(long layer = 0; layer < NBSLAYERS; layer++)
            {
                printf("Super layer %ld/%ld  alt: %f  CN2: %f\n", layer, NBSLAYERS, SLAYER_ALT[layer],
                       SLAYER_CN2[layer]);
            }

            /* look for minimum altitude difference */
            minaltd = LAYER_ALT[NBLAYERS - 1];

            for(long layer = 0; layer < NBSLAYERS - 1; layer++)
            {
                double altdvalue = SLAYER_ALT[layer + 1] - SLAYER_ALT[layer];
                if(altdvalue < minaltd)
                {
                    minaltd = altdvalue;
                    layerindex = layer;
                }
            }
            printf("minimumum distance between layers: %.2f m  (%ld %ld)   (FresnelPropBin = %.2f)\n",
                   minaltd, layerindex, layerindex + 1, atmturbconf.FresnelPropBin);
            if((minaltd > atmturbconf.FresnelPropBin) || (NBSLAYERS == 1))
            {
                OK = 1;
            }
            else
            {
                /* group SLAYERs i and i+1 */
                printf("Group slayers %ld and %ld\n", layerindex, layerindex + 1);
                SLAYER_ALT[layerindex] = (SLAYER_CN2[layerindex] * SLAYER_ALT[layerindex] + SLAYER_CN2[layerindex +
                                          1] * SLAYER_ALT[layerindex + 1]) / (SLAYER_CN2[layerindex] + SLAYER_CN2[layerindex + 1]);
                SLAYER_CN2[layerindex] = SLAYER_CN2[layerindex] + SLAYER_CN2[layerindex + 1];
                for(long layer = layerindex + 1; layer < NBSLAYERS - 1; layer++)
                {
                    SLAYER_ALT[layer] = SLAYER_ALT[layer + 1];
                    SLAYER_CN2[layer] = SLAYER_CN2[layer + 1];
                }
                NBSLAYERS -= 1;
            }
        }
        for(long layer = 0; layer < NBSLAYERS; layer++)
        {
            printf("Super layer %ld  alt: %f  CN2: %g\n", layer, SLAYER_ALT[layer], SLAYER_CN2[layer]);
        }



        long NB_alt_bin_sep = NBSLAYERS - 1;
        double *alt_bin_sep;
        alt_bin_sep = (double *) malloc(NB_alt_bin_sep * sizeof(double));

        for(long i = 0; i < NB_alt_bin_sep; i++)
        {
            alt_bin_sep[i] = 0.5 * (SLAYER_ALT[i] + SLAYER_ALT[i + 1]);
            printf("Super layer %ld Altitude threshhold : %.2f m\n", i, alt_bin_sep[i]);
        }

        free(SLAYER_CN2);


        super_layer_index = (long *) malloc(NBLAYERS * sizeof(long));
        for(long layer = 0; layer < NBLAYERS; layer++)
        {
            long index = 0;
            for(long i = 0; i < NB_alt_bin_sep; i++)
                if((alt_bin_sep[i] < LAYER_ALT[layer])
                        && (alt_bin_sep[i + 1] > LAYER_ALT[layer]))
                {
                    index = i + 1;
                }
            if(LAYER_ALT[layer] > alt_bin_sep[NB_alt_bin_sep - 1])
            {
                index = NB_alt_bin_sep;
            }
            super_layer_index[layer] = index;
            printf("Layer %ld belongs to superlayer %ld/%ld\n", layer,
                   super_layer_index[layer], NBSLAYERS);
        }

        free(alt_bin_sep);
    }









    xpos = (double *) malloc(sizeof(double) * NBLAYERS);
    ypos = (double *) malloc(sizeof(double) * NBLAYERS);
    xpos0 = (double *) malloc(sizeof(double) * NBLAYERS);
    ypos0 = (double *) malloc(sizeof(double) * NBLAYERS);
    xposfcnt = (long *) malloc(sizeof(long) * NBLAYERS);
    yposfcnt = (long *) malloc(sizeof(long) * NBLAYERS);

    vxpix = (double *) malloc(sizeof(double) * NBLAYERS);
    vypix = (double *) malloc(sizeof(double) * NBLAYERS);



    /*
     * =================================================================
     * Compute chromatic translations due to atmospheric dispersion
     * =================================================================
     */
    {
        double h;


        h = 0.0;
        printf("\n\n");
        printf("Refractivity = %g -> %g     %g -> %g\n", atmturbconf.lambda,
               AtmosphereModel_stdAtmModel_N(atm, h, atmturbconf.lambda, 0), SLAMBDA,
               AtmosphereModel_stdAtmModel_N(atm, h, SLAMBDA, 0));


        printf("Computing refraction and position offset for zenithangle = %f\n",
               atmturbconf.zenithangle);
        for(long layer = 0; layer < NBLAYERS; layer++)
        {
            double Roffset;
            double RoffsetS;
            double Rindex;

            xposfcnt[layer] = 0;
            yposfcnt[layer] = 0;

            xpos[layer] = 0.05 * atmturbconf.MasterSize;
            ypos[layer] = 0.05 * atmturbconf.MasterSize;

            // layer shift due to atmospheric refraction
            // computes here what the offset is at a the reference wavelength
            Roffset = 0.0;
            for(h = atm.SiteAlt; h < LAYER_ALT[layer]; h += 1.0)
            {
                double tmpf;

                Rindex = 1.0 + AtmosphereModel_stdAtmModel_N(atm, h, atmturbconf.lambda, 0);
                tmpf = sin(atmturbconf.zenithangle) / sqrt(Rindex * Rindex - sin(atmturbconf.zenithangle) * sin(
                            atmturbconf.zenithangle));
                tmpf -= sin(atmturbconf.zenithangle) / sqrt(1.0 - sin(atmturbconf.zenithangle) * sin(atmturbconf.zenithangle));
                Roffset += tmpf;
                //           printf("h = %12f   Rindex = %12f   Roffset = %12f\n", h, Rindex, Roffset);
            }

            // we compute here the offset at the science wavelength
            RoffsetS = 0.0;
            for(h = atm.SiteAlt; h < LAYER_ALT[layer]; h += 1.0)
            {
                double tmpf;

                Rindex = 1.0 + AtmosphereModel_stdAtmModel_N(atm, h, SLAMBDA, 0);
                tmpf = sin(atmturbconf.zenithangle) / sqrt(Rindex * Rindex - sin(atmturbconf.zenithangle) * sin(
                            atmturbconf.zenithangle));
                tmpf -= sin(atmturbconf.zenithangle) / sqrt(1.0 - sin(atmturbconf.zenithangle) * sin(atmturbconf.zenithangle));
                RoffsetS += tmpf;
            }

            // the refractive offset is the difference between the reference and science wavelength offsets
            ypos[layer] += (RoffsetS - Roffset) / (atmturbconf.PupilScale / pfactor);

            // add here offset due to source position
            // SOURCE_Xpos and SOURCE_Ypos are in radian
            xpos[layer] += atmturbconf.sourceXpos * LAYER_ALT[layer] / (atmturbconf.PupilScale /
                           pfactor);
            ypos[layer] += atmturbconf.sourceYpos * LAYER_ALT[layer] / (atmturbconf.PupilScale /
                           pfactor);

            xpos0[layer] = atmturbconf.sourceXpos * LAYER_ALT[layer] / (atmturbconf.PupilScale /
                           pfactor); // for realtime mode
            ypos0[layer] = atmturbconf.sourceYpos * LAYER_ALT[layer] / (atmturbconf.PupilScale /
                           pfactor);

            vxpix[layer] = LAYER_SPD[layer] * cos(LAYER_DIR[layer]) /
                           (atmturbconf.PupilScale /
                            pfactor); /* pixel coordinate speed, in pixel per sec, x axis */
            vypix[layer] = LAYER_SPD[layer] * sin(LAYER_DIR[layer]) /
                           (atmturbconf.PupilScale /
                            pfactor); /* pixel coordinate speed, in pixel per sec, y axis */

            printf("------ layer %5ld, SPEED = %12f x %12f pix/step, offset = %12f m  [ %12f m  %12f m ] ----------\n",
                   layer, vxpix[layer], vypix[layer], RoffsetS - Roffset, RoffsetS, Roffset);
        }
    }






    /*
     * =================================================================
     * Compute or load master turbulence screens
     * =================================================================
     */

    printf("NBMASTERS = %ld\n", NBMASTERS);
    if(NBMASTERS < NBLAYERS)
    {
        NBMASTERS = NBLAYERS;
    }
    ID_TM = (long *) malloc(sizeof(long) * NBMASTERS);
    for(long i = 0; i < NBMASTERS; i++)
    {
        char fnamemasterwf[STRINGMAXLEN_FILENAME];

        if(WFprecision == 0)
        {
            WRITE_FILENAME(fnamemasterwf, "t%03ld_%ld_f.fits", i, (long) atmturbconf.MasterSize);
        }
        else
        {
            WRITE_FILENAME(fnamemasterwf, "t%03ld_%ld_d.fits", i, (long) atmturbconf.MasterSize);
        }

        char imnamemasterwf[STRINGMAXLEN_IMGNAME];
        WRITE_IMAGENAME(imnamemasterwf, "TM%ld", i);

        imageID IDt;
        load_fits(fnamemasterwf, imnamemasterwf, LOADFITS_ERRMODE_WARNING, &IDt);
        if(IDt == -1)
        {
            printf("CREATING %s   (%f - %f)\n", fnamemasterwf,
                   LAYER_OUTERSCALE[i] / atmturbconf.PupilScale, LAYER_INNERSCALE[i] / atmturbconf.PupilScale);
            AtmosphericTurbulence_make_master_turbulence_screen(imnamemasterwf, "tursctmp",
                    atmturbconf.MasterSize,
                    LAYER_OUTERSCALE[i] / atmturbconf.PupilScale, LAYER_INNERSCALE[i] / atmturbconf.PupilScale,
                    WFprecision);

            if(WFprecision == 0)
            {
                save_fl_fits(imnamemasterwf, fnamemasterwf);
            }
            else
            {
                save_db_fits(imnamemasterwf, fnamemasterwf);
            }
            delete_image_ID("tursctmp", DELETE_IMAGE_ERRMODE_WARNING);
        }
        ID_TM[i] = image_ID(imnamemasterwf);
    }
    ID_TML = (long *) malloc(sizeof(long) * NBLAYERS);

    long layer1 = 0;
    for(long layer = 0; layer < NBLAYERS; layer++)
    {
        ID_TML[layer] = ID_TM[layer1];
        if(layer1 == NBMASTERS)
        {
            printf("ERROR: number of master turbulence phase screens (%ld) is too small\n",
                   NBMASTERS);
            exit(0);
        }
        layer1++;
    }







    /*
     * =================================================================
     * Measure r0 (pix) for each master turbulence screen
     * =================================================================
     */
    {
        long dpix = 50;
        double r0tot = 0.0;

        double r0;

        long r0cnt = 0;

        assert(atmturbconf.MasterSize > dpix);

        for(long k = 0; k < NBMASTERS; k++)
        {
            long cnt = 0;
            double tot = 0.0;
            if(WFprecision == 0)
            {
                for(uint32_t ii = 0; ii < atmturbconf.MasterSize - dpix; ii++)
                    for(uint32_t jj = 0; jj < atmturbconf.MasterSize; jj++)
                    {
                        float p1 = data.image[ID_TM[k]].array.F[jj * atmturbconf.MasterSize + ii];
                        float p2 = data.image[ID_TM[k]].array.F[jj * atmturbconf.MasterSize + ii + dpix];
                        tot += (p1 - p2) * (p1 - p2);
                        cnt++;
                    }
            }
            else
            {
                for(uint32_t ii = 0; ii < atmturbconf.MasterSize - dpix; ii++)
                    for(uint32_t jj = 0; jj < atmturbconf.MasterSize; jj++)
                    {
                        double p1 = data.image[ID_TM[k]].array.D[jj * atmturbconf.MasterSize + ii];
                        double p2 = data.image[ID_TM[k]].array.D[jj * atmturbconf.MasterSize + ii + dpix];
                        tot += (p1 - p2) * (p1 - p2);
                        cnt++;
                    }
            }
            r0 = 1.0 * dpix * pow((tot / cnt) / 6.88, -3.0 / 5.0);
            printf("TURBULENCE MASTER %ld    r0 = %g pix\n", k, r0);
            r0tot += r0;
            r0cnt++;
        }
        r0 = r0tot / r0cnt;
        printf("r0 = %g pix -> %g pix\n", r0,
               atmturbconf.lambda / (atmturbconf.seeing / 3600.0 / 180.0 * PI) / atmturbconf.PupilScale * pfactor);

        // renormalize turbulence screens such that a single screen has the right r0
        double normcoeff = pow(r0 / (atmturbconf.lambda / (atmturbconf.seeing / 3600.0 / 180.0 * PI) /
                                     atmturbconf.PupilScale * pfactor), 5.0 / 6.0);
        if(WFprecision == 0)
        {
            for(long k = 0; k < NBMASTERS; k++)
            {
                for(uint64_t ii = 0; ii < atmturbconf.MasterSize * atmturbconf.MasterSize; ii++)
                {
                    data.image[ID_TM[k]].array.F[ii] *= normcoeff;
                }
            }
        }
        else
        {
            for(long k = 0; k < NBMASTERS; k++)
            {
                for(uint64_t ii = 0; ii < atmturbconf.MasterSize * atmturbconf.MasterSize; ii++)
                {
                    data.image[ID_TM[k]].array.D[ii] *= normcoeff;
                }
            }
        }

        r0tot = 0.0;
        r0cnt = 0;
        for(long k = 0; k < NBMASTERS; k++)
        {
            long cnt = 0;
            double tot = 0.0;

            if(WFprecision == 0)
            {
                for(uint32_t ii = 0; ii < atmturbconf.MasterSize - dpix; ii++)
                    for(uint32_t jj = 0; jj < atmturbconf.MasterSize; jj++)
                    {
                        float p1 = data.image[ID_TM[k]].array.F[jj * atmturbconf.MasterSize + ii];
                        float p2 = data.image[ID_TM[k]].array.F[jj * atmturbconf.MasterSize + ii + dpix];
                        tot += (p1 - p2) * (p1 - p2);
                        cnt++;
                    }
            }
            else
            {
                for(uint32_t ii = 0; ii < atmturbconf.MasterSize - dpix; ii++)
                    for(uint32_t jj = 0; jj < atmturbconf.MasterSize; jj++)
                    {
                        double p1 = data.image[ID_TM[k]].array.D[jj * atmturbconf.MasterSize + ii];
                        double p2 = data.image[ID_TM[k]].array.D[jj * atmturbconf.MasterSize + ii + dpix];
                        tot += (p1 - p2) * (p1 - p2);
                        cnt++;
                    }
            }
            r0 = 1.0 * dpix * pow((tot / cnt) / 6.88, -3.0 / 5.0);
            printf("TURBULENCE MASTER %ld    r0 = %g pix\n", k, r0);
            r0tot += r0;
            r0cnt++;
        }
        r0 = r0tot / r0cnt;
        printf("r0 = %g pix\n", r0);
    }


    // target seeing = seeing [arcsec]
    // ref lambda = atmturbconf.lambda [m]
    // if single screen:
    // r0[m] = atmturbconf.lambda[m]/seeing[rad]
    // r0[pix] = atmturbconf.lambda*1.0e-6/(seeing/3600.0/180.0*PI)/PupilScale
    // multiply by (r0/r0goal)^6/5


    // each layer coeff mult by sqrt(fracCN2/cos(zenithangle)))


    for(long layer = 0; layer < NBLAYERS; layer++)
    {
        //      coeff = 3.645*183.8115*pow(10.0,-12)/LAMBDA/LAMBDA*PupilScale*PUPIL_SCALE*sqrt(LAYER_CN2[i]);
        //if(pfactor==2)
        //	coeff *= 1.0/(pfactor*pfactor*pfactor*pfactor*pfactor*pfactor);
        if(WFprecision == 0)
        {
            for(uint64_t ii = 0; ii < atmturbconf.MasterSize * atmturbconf.MasterSize; ii++)
            {
                data.image[ID_TML[layer]].array.F[ii] *= sqrt(LAYER_CN2[layer] / cos(atmturbconf.zenithangle));
            }
        }
        else
        {
            for(uint64_t ii = 0; ii < atmturbconf.MasterSize * atmturbconf.MasterSize; ii++)
            {
                data.image[ID_TML[layer]].array.D[ii] *= sqrt(LAYER_CN2[layer] / cos(atmturbconf.zenithangle));
            }
        }
        printf("Layer %ld, coeff = %g\n", layer, sqrt(LAYER_CN2[layer] / cos(atmturbconf.zenithangle)));
    }






    /*
     * =================================================================
     * Create image arrays for internal memory
     * =================================================================
     */

    if(WFprecision == 0)
    {
        create_2Dimage_ID("array1", naxes[0], naxes[1], &ID_array1);
        if(atmturbconf.flag_SWF_make == 1)
        {
            create_2Dimage_ID("sarray1", naxes[0], naxes[1], &ID_sarray1);
        }
    }
    else
    {
        create_2Dimage_ID_double("array1", naxes[0], naxes[1], &ID_array1);
        if(atmturbconf.flag_SWF_make == 1)
        {
            create_2Dimage_ID_double("sarray1", naxes[0], naxes[1], &ID_sarray1);
        }
    }



    if(atmturbconf.flag_WFampl == 1) // includes sub pixel translation
    {
        if(WFprecision == 0)
        {
            create_2DCimage_ID("array2", naxes[0], naxes[1], &ID_array2);
            if(atmturbconf.flag_SWF_make == 1)
            {
                create_2DCimage_ID("sarray2", naxes[0], naxes[1], &ID_sarray2);
            }
        }
        else
        {
            create_2DCimage_ID_double("array2", naxes[0], naxes[1], &ID_array2);
            if(atmturbconf.flag_SWF_make == 1)
            {
                create_2DCimage_ID_double("sarray2", naxes[0], naxes[1], &ID_sarray2);
            }
        }
    }



    if(WFprecision == 0)
    {
        if((array = (complex_float *) malloc(NBFRAMES * naxes[0] * naxes[1] * sizeof(
                complex_float))) == NULL)
        {
            printf("Memory allocation error (\"array\" in make_AtmosphericTurbulence_wavefront_series)\n");
            printf("Decrease the size of the wavefront cube\n");
            exit(0);
        }

        if(atmturbconf.flag_SWF_make == 1)
        {
            if((sarray = (complex_float *) malloc(NBFRAMES * naxes[0] * naxes[1] * sizeof(
                    complex_float))) == NULL)
            {
                printf("Memory allocation error (\"sarray\" in make_AtmosphericTurbulence_wavefront_series)\n");
                printf("Decrease the size of the wavefront cube\n");
                exit(0);
            }
        }
    }
    else
    {
        if((array_double = (complex_double *) malloc(NBFRAMES * naxes[0] * naxes[1] *
                           sizeof(complex_double))) == NULL)
        {
            printf("Memory allocation error (\"array_double\" in make_AtmosphericTurbulence_wavefront_series)\n");
            printf("Decrease the size of the wavefront cube\n");
            exit(0);
        }

        if(atmturbconf.flag_SWF_make == 1)
        {
            if((sarray_double = (complex_double *) malloc(NBFRAMES * naxes[0] * naxes[1] *
                                sizeof(complex_double))) == NULL)
            {
                printf("Memory allocation error (\"sarray_double\" in make_AtmosphericTurbulence_wavefront_series)\n");
                printf("Decrease the size of the wavefront cube\n");
                exit(0);
            }
        }
    }





    if(atmturbconf.flag_SWF_SHMoutput == 1)
    {
        if(WFprecision == 0)
        {
            create_image_ID("atmwfpha", 2, naxesout, _DATATYPE_FLOAT, 1, 0, 0, &IDshmpha);
            create_image_ID("atmwfamp", 2, naxesout, _DATATYPE_FLOAT, 1, 0, 0, &IDshmamp);
        }
        else
        {
            create_image_ID("atmwfpha", 2, naxesout, _DATATYPE_DOUBLE, 1, 0, 0, &IDshmpha);
            create_image_ID("atmwfamp", 2, naxesout, _DATATYPE_DOUBLE, 1, 0, 0, &IDshmamp);
        }
    }


    if(atmturbconf.flag_SHMoutput == 1)
    {
        char imnameoutpha[STRINGMAXLEN_IMGNAME];
        WRITE_IMAGENAME(imnameoutpha, "%spha", atmturbconf.SWFSHMprefix);
        if(WFprecision == 0)
        {
            create_image_ID(imnameoutpha, 2, naxesout, _DATATYPE_FLOAT, 1, 0, 0, &IDshmspha);
        }
        else
        {
            create_image_ID(imnameoutpha, 2, naxesout, _DATATYPE_DOUBLE, 1, 0, 0, &IDshmspha);
        }

        char imnameoutamp[STRINGMAXLEN_IMGNAME];
        WRITE_IMAGENAME(imnameoutamp, "%samp", atmturbconf.SWFSHMprefix);
        if(WFprecision == 0)
        {
            create_image_ID(imnameoutamp, 2, naxesout, _DATATYPE_FLOAT, 1, 0, 0, &IDshmsamp);
        }
        else
        {
            create_image_ID(imnameoutamp, 2, naxesout, _DATATYPE_DOUBLE, 1, 0, 0, &IDshmsamp);
        }

        int kw = 0;

        strcpy(data.image[IDshmspha].kw[kw].name, "TIME");
        data.image[IDshmspha].kw[kw].type = 'D';
        data.image[IDshmspha].kw[kw].value.numf = 0.0;
        strcpy(data.image[IDshmspha].kw[kw].comment, "Physical time [sec]");

        strcpy(data.image[IDshmsamp].kw[kw].name, "TIME");
        data.image[IDshmsamp].kw[kw].type = 'D';
        data.image[IDshmsamp].kw[kw].value.numf = 0.0;
        strcpy(data.image[IDshmsamp].kw[kw].comment, "Physical time [sec]");
    }







    /*
     * =================================================================
     * Look for existing output files to set start time
     * =================================================================
     */

    printf("SKIP_EXISTING = %d\n", atmturbconf.flag_SkipExisting);
    if(atmturbconf.flag_SkipExisting == 1)
    {
        start_cubeindex = 0;
        int OK = 1;
        while(OK == 1)
        {
            //char fnamewf[STRINGMAXLEN_FILENAME];
            char fnameswf[STRINGMAXLEN_FILENAME];

            // lambda in pm
            //WRITE_FILENAME(fnamewf, "%s%08ld.%09ld.pha.fits", WFfileprefix, start_cubeindex, (long)(SLAMBDA * 1e12 + 0.5));
            WRITE_FILENAME(fnameswf, "%s%08ld.%09ld.pha.fits", atmturbconf.SWFfileprefix,
                           start_cubeindex, (long)(SLAMBDA * 1e12 + 0.5));
            printf("TESTING FILE %s ... ", fnameswf);
            if(file_exists(fnameswf) == 1)
            {
                start_cubeindex ++;
                printf("exists\n");
                OK = 1;
            }
            else
            {
                printf("does not exist\n");
                OK = 0;
            }
        }
    }
    printf("Start cubeindex = %ld\n", start_cubeindex);




    /*
     * =================================================================
     * Set wind speed for each layer
     * =================================================================
     */
    {
        long WSPEEDsize = naxes[0];
        double WSPEEDpixscale = atmturbconf.PupilScale;
        uint32_t vKwindsize = (uint32_t)(20000.0 / WSPEEDpixscale);  // 20 km

        for(long layer = 0; layer < NBLAYERS; layer++)
        {
            char fnamewspeed[STRINGMAXLEN_FILENAME];
            WRITE_FILENAME(fnamewspeed, "wspeed_%03ld.fits", layer);

            char imnamewspeed[STRINGMAXLEN_IMGNAME];
            WRITE_IMAGENAME(imnamewspeed, "wspeed_%03ld", layer);

            imageID ID;
            load_fits(fnamewspeed, imnamewspeed, LOADFITS_ERRMODE_WARNING, &ID);
            if(ID == -1)
            {
                printf("COMPUTE WIND SPEED SCATTER - LAYER %ld   sigma = %f m/s\n", layer,
                       LAYER_SIGMAWSPEED[layer]);
                make_AtmosphericTurbulence_vonKarmanWind(vKwindsize, WSPEEDpixscale,
                        LAYER_SIGMAWSPEED[layer], LAYER_LWIND[layer], WSPEEDsize, imnamewspeed);
                save_fits(imnamewspeed, fnamewspeed);


                // print wind speed as a function of time and position
                char fnamewspeedtxt[STRINGMAXLEN_FILENAME];
                WRITE_FILENAME(fnamewspeedtxt, "wspeed_%03ld.txt", layer);
                ID = image_ID(imnamewspeed);

                FILE *fpwsp;
                if((fpwsp = fopen(fnamewspeedtxt, "w")) == NULL)
                {
                    printf("ERROR: cannot create file \"%s\"\n", fnamewspeedtxt);
                    exit(0);
                }
                for(uint64_t ii = 0; ii < vKwindsize; ii++)
                {
                    fprintf(fpwsp, "%6ld   %20.16f   %20.16f  %.16g  %.16g  %.16g\n", ii,
                            WSPEEDpixscale * ii, WSPEEDpixscale * ii * LAYER_SPD[layer],
                            data.image[ID].array.F[ii], data.image[ID].array.F[vKwindsize + ii],
                            data.image[ID].array.F[2 * vKwindsize + ii]);
                }
                fclose(fpwsp);
            }
        }

    }








    printf("WAVEFRONT_AMPLITUDE = %d\n", atmturbconf.flag_WFampl);

    fpxypos = fopen("xypos.log", "w");
    fclose(fpxypos);

    fpxypos = fopen("xypos3.log", "w");
    fclose(fpxypos);


    /*
     * =================================================================
     * Main wavefront computation loop
     * =================================================================
     */

    for(long cubeindex = start_cubeindex; cubeindex < atmturbconf.NBTimeCube; cubeindex++)
    {
        for(long frame = 0; frame < NBFRAMES; frame++)
        {
            if(atmturbconf.flag_SWF_make == 1)
            {
                if(WFprecision == 0)
                {
                    for(uint32_t ii = 0; ii < naxes[0]; ii++)
                        for(uint32_t jj = 0; jj < naxes[1]; jj++)
                        {
                            data.image[ID_array1].array.F[jj * naxes[0] + ii] = 0.0;
                            data.image[ID_sarray1].array.F[jj * naxes[0] + ii] = 0.0;
                        }
                    if(atmturbconf.flag_WFampl == 1)
                        for(uint32_t ii = 0; ii < naxes[0]; ii++)
                            for(uint32_t jj = 0; jj < naxes[1]; jj++)
                            {
                                data.image[ID_array2].array.CF[jj * naxes[0] + ii].re = 1.0;
                                data.image[ID_array2].array.CF[jj * naxes[0] + ii].im = 0.0;
                                data.image[ID_sarray2].array.CF[jj * naxes[0] + ii].re = 1.0;
                                data.image[ID_sarray2].array.CF[jj * naxes[0] + ii].im = 0.0;
                            }
                }
                else
                {
                    for(uint32_t ii = 0; ii < naxes[0]; ii++)
                        for(uint32_t jj = 0; jj < naxes[1]; jj++)
                        {
                            data.image[ID_array1].array.D[jj * naxes[0] + ii] = 0.0;
                            data.image[ID_sarray1].array.D[jj * naxes[0] + ii] = 0.0;
                        }
                    if(atmturbconf.flag_WFampl == 1)
                        for(uint32_t ii = 0; ii < naxes[0]; ii++)
                            for(uint32_t jj = 0; jj < naxes[1]; jj++)
                            {
                                data.image[ID_array2].array.CD[jj * naxes[0] + ii].re = 1.0;
                                data.image[ID_array2].array.CD[jj * naxes[0] + ii].im = 0.0;
                                data.image[ID_sarray2].array.CD[jj * naxes[0] + ii].re = 1.0;
                                data.image[ID_sarray2].array.CD[jj * naxes[0] + ii].im = 0.0;
                            }
                }
            }
            else
            {
                if(WFprecision == 0)
                {
                    for(uint32_t ii = 0; ii < naxes[0]; ii++)
                        for(uint32_t jj = 0; jj < naxes[1]; jj++)
                        {
                            data.image[ID_array1].array.F[jj * naxes[0] + ii] = 0.0;
                        }
                    if(atmturbconf.flag_WFampl == 1)
                        for(uint32_t ii = 0; ii < naxes[0]; ii++)
                            for(uint32_t jj = 0; jj < naxes[1]; jj++)
                            {
                                data.image[ID_array2].array.CF[jj * naxes[0] + ii].re = 1.0;
                                data.image[ID_array2].array.CF[jj * naxes[0] + ii].im = 0.0;
                            }
                }
                else
                {
                    for(uint32_t ii = 0; ii < naxes[0]; ii++)
                        for(uint32_t jj = 0; jj < naxes[1]; jj++)
                        {
                            data.image[ID_array1].array.D[jj * naxes[0] + ii] = 0.0;
                        }
                    if(atmturbconf.flag_WFampl == 1)
                        for(uint32_t ii = 0; ii < naxes[0]; ii++)
                            for(uint32_t jj = 0; jj < naxes[1]; jj++)
                            {
                                data.image[ID_array2].array.CD[jj * naxes[0] + ii].re = 1.0;
                                data.image[ID_array2].array.CD[jj * naxes[0] + ii].im = 0.0;
                            }
                }
            }

            usleep(atmturbconf.TimeDelayus);


            if(atmturbconf.flag_WaitSem == 1) // wait for semaphore to advance to next WF step
            {
                printf("WAITING for semaphore #0 \"%s\" ...\n", atmturbconf.WaitSemName);
                COREMOD_MEMORY_image_set_semwait(atmturbconf.WaitSemName, 0);
                printf("Done\n");
            }


            clock_gettime(CLOCK_REALTIME, &tnow);
            tnowdouble = 1.0 * tnow.tv_sec + 1.0e-9 * tnow.tv_nsec;
            tnowdouble *= atmturbconf.RealTimeFactor;
            for(long layer = NBLAYERS - 1; layer != -1; layer--)
            {
                if(atmturbconf.flag_RealTime == 0)
                {
                    tnowdouble = (cubeindex * NBFRAMES + frame) * atmturbconf.TimeStep;
                    printf("\rLayer %2ld/%2ld, Frame %4ld/%4ld, File %6ld/%6ld  [TIME = %10.4f s]  ",
                           layer, NBLAYERS, frame, NBFRAMES, cubeindex, atmturbconf.NBTimeCube,
                           (cubeindex * NBFRAMES + frame)*atmturbconf.TimeStep);
                }
                else
                {
                    printf("\rLayer %2ld/%2ld, Frame %4ld/%4ld, File %6ld/%6ld  [PHYSICAL TIME = %.9lf s]  ",
                           layer, NBLAYERS, frame, NBFRAMES, cubeindex, atmturbconf.NBTimeCube, tnowdouble);
                }
                fflush(stdout);

                // recompute Scoeff for this layer
                Nlambda = AtmosphereModel_stdAtmModel_N(atm, LAYER_ALT[layer], atmturbconf.lambda, 0);
                Nslambda = AtmosphereModel_stdAtmModel_N(atm, LAYER_ALT[layer], SLAMBDA, 0);
                Scoeff =  atmturbconf.lambda / SLAMBDA * Nslambda /
                          Nlambda; // multiplicative coefficient to go from reference lambda phase to science lambda phase


                if(layer != NBLAYERS - 1)
                {
                    if(super_layer_index[layer + 1] != super_layer_index[layer])
                    {
                        if(atmturbconf.flag_FresnelProp == 1)
                        {
                            Fresnel_propagate_wavefront("array2", "array2p", atmturbconf.PupilScale / pfactor,
                                                        (SLAYER_ALT[super_layer_index[layer + 1]] -
                                                         SLAYER_ALT[super_layer_index[layer]]) / cos(atmturbconf.zenithangle), atmturbconf.lambda);
                            delete_image_ID("array2", DELETE_IMAGE_ERRMODE_WARNING);
                            chname_image_ID("array2p", "array2");
                        }

                        ID_array2 = image_ID("array2");
                        if(atmturbconf.flag_SWF_make == 1)
                        {
                            if(atmturbconf.flag_FresnelProp == 1)
                            {
                                //				printf("FRESNEL PROPAGATION\n");  //TEST
                                //			fflush(stdout);

                                Fresnel_propagate_wavefront("sarray2", "sarray2p", atmturbconf.PupilScale / pfactor,
                                                            (SLAYER_ALT[super_layer_index[layer + 1]] -
                                                             SLAYER_ALT[super_layer_index[layer]]) / cos(atmturbconf.zenithangle), SLAMBDA);
                                delete_image_ID("sarray2", DELETE_IMAGE_ERRMODE_WARNING);
                                chname_image_ID("sarray2p", "sarray2");
                            }
                            ID_sarray2 = image_ID("sarray2");
                        }
                    }
                }


                // layer_scale = (SODIUM_ALT-LAYER_ALT[layer])/SODIUM_ALT;

                //vpix = 0.0; //0.1*sin(11.0*vindex*(layer+3))*sqrt(vxpix[layer]*vxpix[layer]+vypix[layer]*vypix[layer]);
                //PA = sin(10.0 * vindex * (layer + 2));

                if(atmturbconf.flag_RealTime == 1) // real time
                {
                    xpos[layer] = xpos0[layer] + vxpix[layer] * tnowdouble + 1.0 * xposfcnt[layer] *
                                  naxes_MASTER[0];
                    ypos[layer] = ypos0[layer] + vypix[layer] * tnowdouble + 1.0 * yposfcnt[layer] *
                                  naxes_MASTER[0];
                }
                else // non real time
                {

                    xpos[layer] += vxpix[layer] * atmturbconf.TimeStep;
                    ypos[layer] += vypix[layer] * atmturbconf.TimeStep;
                }



                long xref = (long)(xpos[layer]);
                long yref = (long)(ypos[layer]);

                while(xpos[layer] < 0)
                {
                    xpos[layer] += 1.0 * naxes_MASTER[0];
                    xposfcnt[layer]++;
                    xref = (long)(xpos[layer]);
                }
                while(xpos[layer] > 1.0 * naxes_MASTER[0])
                {
                    xpos[layer] -= 1.0 * naxes_MASTER[0];
                    xposfcnt[layer]--;
                    xref = (long)(xpos[layer]);
                }

                while(ypos[layer] < 0)
                {
                    ypos[layer] += 1.0 * naxes_MASTER[1];
                    yposfcnt[layer]++;
                    yref = (long)(ypos[layer]);
                }
                while(ypos[layer] > 1.0 * naxes_MASTER[1])
                {
                    ypos[layer] -= 1.0 * naxes_MASTER[1];
                    yposfcnt[layer]--;
                    yref = (long)(ypos[layer]);
                }


                if(xref == naxes_MASTER[0])
                {
                    xref = 0;
                }
                if(yref == naxes_MASTER[1])
                {
                    yref = 0;
                }

                long iimax = naxes_MASTER[0] - xref;
                long jjmax = naxes_MASTER[1] - yref;
                if(iimax > naxes[0])
                {
                    iimax = naxes[0];
                }
                if(jjmax > naxes[1])
                {
                    jjmax = naxes[1];
                }
                //xrefm = xref - naxes_MASTER[0];
                //yrefm = yref - naxes_MASTER[1];


                /* make wavefront */
                if(WFprecision == 0) // floating point precision
                {
                    if(BICUBIC == 0)
                    {
                        // bilinear interpolation
                        for(uint32_t ii = 0; ii < naxes[0]; ii++)
                            for(uint32_t jj = 0; jj < naxes[1]; jj++)
                            {
                                double iimf = fmod((xpos[layer] + ii), 1.0 * naxes_MASTER[0]);
                                double jjmf = fmod((ypos[layer] + jj), 1.0 * naxes_MASTER[1]);
                                long iim = (long)(iimf);
                                long jjm = (long)(jjmf);
                                double iifrac = iimf - iim;
                                double jjfrac = jjmf - jjm;
                                long iim1 = iim + 1;
                                long jjm1 = jjm + 1;

                                if(iim == atmturbconf.MasterSize)
                                {
                                    iim = 0;
                                }
                                if(jjm == atmturbconf.MasterSize)
                                {
                                    jjm = 0;
                                }
                                if(iim1 > atmturbconf.MasterSize - 1)
                                {
                                    iim1 -= atmturbconf.MasterSize;
                                }
                                if(jjm1 > atmturbconf.MasterSize - 1)
                                {
                                    jjm1 -= atmturbconf.MasterSize;
                                }

                                double value = (1.0 - iifrac) * (1.0 - jjfrac) *
                                               data.image[ID_TML[layer]].array.F[jjm
                                                       * naxes_MASTER[0] + iim];
                                value += (1.0 - iifrac) * (jjfrac) * data.image[ID_TML[layer]].array.F[jjm1 *
                                         naxes_MASTER[0] + iim];
                                value += (iifrac) * (jjfrac) * data.image[ID_TML[layer]].array.F[jjm1 *
                                         naxes_MASTER[0] + iim1];
                                value += (iifrac) * (1.0 - jjfrac) * data.image[ID_TML[layer]].array.F[jjm *
                                         naxes_MASTER[0] + iim1];

                                data.image[ID_array1].array.F[jj * naxes[0] + ii] += value;
                                if(atmturbconf.flag_WFampl == 1)
                                {
                                    float re = data.image[ID_array2].array.CF[jj * naxes[0] + ii].re;
                                    float im = data.image[ID_array2].array.CF[jj * naxes[0] + ii].im;

                                    data.image[ID_array2].array.CF[jj * naxes[0] + ii].re =
                                        re * cos(value) - im * sin(value);

                                    data.image[ID_array2].array.CF[jj * naxes[0] + ii].im =
                                        re * sin(value) + im * cos(value);
                                }
                            }

                        //fpxypos = fopen("xypos.log", "a");
                        //fprintf(fpxypos, "%5ld %4ld    %10.8f %10.8f      %5ld %10.8f %5ld %10.8f    %10.8f %10.8f  %.18g        %.18g   %.18g   %.18g   %.18g\n", vindex, layer, xpos[layer], ypos[layer], iim, iifrac, jjm, jjfrac, 1.0*iim+iifrac, 1.0*jjm+jjfrac, value, data.image[ID_TML[layer]].array.F[jjm*naxes_MASTER[0]+iim], data.image[ID_TML[layer]].array.F[jjm1*naxes_MASTER[0]+iim], data.image[ID_TML[layer]].array.F[jjm1*naxes_MASTER[0]+iim1], data.image[ID_TML[layer]].array.F[jjm*naxes_MASTER[0]+iim1]);
                        //fclose(fpxypos);
                    }
                    else
                    {
                        // bicubic interpolation
                        for(uint32_t ii = 0; ii < naxes[0]; ii++)
                            for(uint32_t jj = 0; jj < naxes[1]; jj++)
                            {
                                double a00, a01, a02, a03;
                                double a10, a11, a12, a13;
                                double a20, a21, a22, a23;
                                double a30, a31, a32, a33;
                                double p00, p01, p02, p03;
                                double p10, p11, p12, p13;
                                double p20, p21, p22, p23;
                                double p30, p31, p32, p33;

                                double iimf = fmod((xpos[layer] + ii), 1.0 * naxes_MASTER[0]);
                                double jjmf = fmod((ypos[layer] + jj), 1.0 * naxes_MASTER[1]);

                                long iim = (long)(iimf);
                                long jjm = (long)(jjmf);

                                double x = iimf - iim;
                                double y = jjmf - jjm;


                                long iim0 = iim - 1;
                                long iim1 = iim;
                                long iim2 = iim + 1;
                                long iim3 = iim + 2;
                                if(iim1 > atmturbconf.MasterSize - 1)
                                {
                                    iim1 -= atmturbconf.MasterSize;
                                }
                                if(iim2 > atmturbconf.MasterSize - 1)
                                {
                                    iim2 -= atmturbconf.MasterSize;
                                }
                                if(iim3 > atmturbconf.MasterSize - 1)
                                {
                                    iim3 -= atmturbconf.MasterSize;
                                }
                                if(iim0 < 0)
                                {
                                    iim0 += atmturbconf.MasterSize;
                                }

                                long jjm0 = jjm - 1;
                                long jjm1 = jjm;
                                long jjm2 = jjm + 1;
                                long jjm3 = jjm + 2;
                                if(jjm1 > atmturbconf.MasterSize - 1)
                                {
                                    iim1 -= atmturbconf.MasterSize;
                                }
                                if(jjm2 > atmturbconf.MasterSize - 1)
                                {
                                    jjm2 -= atmturbconf.MasterSize;
                                }
                                if(jjm3 > atmturbconf.MasterSize - 1)
                                {
                                    jjm3 -= atmturbconf.MasterSize;
                                }
                                if(jjm0 < 0)
                                {
                                    jjm0 += atmturbconf.MasterSize;
                                }

                                /*assert(iim0>=0);
                                assert(iim1>=0);
                                assert(iim2>=0);
                                assert(iim3>=0);
                                assert(iim0<MasterSize);
                                assert(iim1<MasterSize);
                                assert(iim2<MasterSize);
                                assert(iim3<MasterSize);
                                assert(jjm0>=0);
                                assert(jjm1>=0);
                                assert(jjm2>=0);
                                assert(jjm3>=0);
                                assert(jjm0<MasterSize);
                                assert(jjm1<MasterSize);
                                assert(jjm2<MasterSize);
                                assert(jjm3<MasterSize);
                                */

                                p00 = data.image[ID_TML[layer]].array.F[jjm0 * naxes_MASTER[0] + iim0];
                                p01 = data.image[ID_TML[layer]].array.F[jjm1 * naxes_MASTER[0] + iim0];
                                p02 = data.image[ID_TML[layer]].array.F[jjm2 * naxes_MASTER[0] + iim0];
                                p03 = data.image[ID_TML[layer]].array.F[jjm3 * naxes_MASTER[0] + iim0];

                                p10 = data.image[ID_TML[layer]].array.F[jjm0 * naxes_MASTER[0] + iim1];
                                p11 = data.image[ID_TML[layer]].array.F[jjm1 * naxes_MASTER[0] + iim1];
                                p12 = data.image[ID_TML[layer]].array.F[jjm2 * naxes_MASTER[0] + iim1];
                                p13 = data.image[ID_TML[layer]].array.F[jjm3 * naxes_MASTER[0] + iim1];

                                p20 = data.image[ID_TML[layer]].array.F[jjm0 * naxes_MASTER[0] + iim2];
                                p21 = data.image[ID_TML[layer]].array.F[jjm1 * naxes_MASTER[0] + iim2];
                                p22 = data.image[ID_TML[layer]].array.F[jjm2 * naxes_MASTER[0] + iim2];
                                p23 = data.image[ID_TML[layer]].array.F[jjm3 * naxes_MASTER[0] + iim2];

                                p30 = data.image[ID_TML[layer]].array.F[jjm0 * naxes_MASTER[0] + iim3];
                                p31 = data.image[ID_TML[layer]].array.F[jjm1 * naxes_MASTER[0] + iim3];
                                p32 = data.image[ID_TML[layer]].array.F[jjm2 * naxes_MASTER[0] + iim3];
                                p33 = data.image[ID_TML[layer]].array.F[jjm3 * naxes_MASTER[0] + iim3];


                                a00 = p11;
                                a01 = -.5 * p10 + .5 * p12;
                                a02 = p10 - 2.5 * p11 + 2 * p12 - .5 * p13;
                                a03 = -.5 * p10 + 1.5 * p11 - 1.5 * p12 + .5 * p13;
                                a10 = -.5 * p01 + .5 * p21;
                                a11 = .25 * p00 - .25 * p02 - .25 * p20 + .25 * p22;
                                a12 = -.5 * p00 + 1.25 * p01 - p02 + .25 * p03 + .5 * p20 - 1.25 * p21 + p22 -
                                      .25 * p23;
                                a13 = .25 * p00 - .75 * p01 + .75 * p02 - .25 * p03 - .25 * p20 + .75 * p21 -
                                      .75 * p22 + .25 * p23;
                                a20 = p01 - 2.5 * p11 + 2 * p21 - .5 * p31;
                                a21 = -.5 * p00 + .5 * p02 + 1.25 * p10 - 1.25 * p12 - p20 + p22 + .25 * p30 -
                                      .25 * p32;
                                a22 = p00 - 2.5 * p01 + 2 * p02 - .5 * p03 - 2.5 * p10 + 6.25 * p11 - 5 * p12 +
                                      1.25 * p13 + 2 * p20 - 5 * p21 + 4 * p22 - p23 - .5 * p30 + 1.25 * p31 - p32 +
                                      .25 * p33;
                                a23 = -.5 * p00 + 1.5 * p01 - 1.5 * p02 + .5 * p03 + 1.25 * p10 - 3.75 * p11 +
                                      3.75 * p12 - 1.25 * p13 - p20 + 3 * p21 - 3 * p22 + p23 + .25 * p30 - .75 * p31
                                      + .75 * p32 - .25 * p33;
                                a30 = -.5 * p01 + 1.5 * p11 - 1.5 * p21 + .5 * p31;
                                a31 = .25 * p00 - .25 * p02 - .75 * p10 + .75 * p12 + .75 * p20 - .75 * p22 -
                                      .25 * p30 + .25 * p32;
                                a32 = -.5 * p00 + 1.25 * p01 - p02 + .25 * p03 + 1.5 * p10 - 3.75 * p11 + 3 *
                                      p12 - .75 * p13 - 1.5 * p20 + 3.75 * p21 - 3 * p22 + .75 * p23 + .5 * p30 - 1.25
                                      * p31 + p32 - .25 * p33;
                                a33 = .25 * p00 - .75 * p01 + .75 * p02 - .25 * p03 - .75 * p10 + 2.25 * p11 -
                                      2.25 * p12 + .75 * p13 + .75 * p20 - 2.25 * p21 + 2.25 * p22 - .75 * p23 - .25 *
                                      p30 + .75 * p31 - .75 * p32 + .25 * p33;

                                double x2 = x * x;
                                double x3 = x2 * x;
                                double y2 = y * y;
                                double y3 = y2 * y;

                                double value = (a00 + a01 * y + a02 * y2 + a03 * y3)
                                               + (a10 + a11 * y + a12 * y2 + a13 * y3) * x
                                               + (a20 + a21 * y + a22 * y2 + a23 * y3) * x2
                                               + (a30 + a31 * y + a32 * y2 + a33 * y3) * x3;


                                data.image[ID_array1].array.F[jj * naxes[0] + ii] += value;
                                if(atmturbconf.flag_WFampl == 1)
                                {
                                    float re = data.image[ID_array2].array.CF[jj * naxes[0] + ii].re;
                                    float im = data.image[ID_array2].array.CF[jj * naxes[0] + ii].im;
                                    data.image[ID_array2].array.CF[jj * naxes[0] + ii].re =
                                        re * cos(value) - im * sin(value);
                                    data.image[ID_array2].array.CF[jj * naxes[0] + ii].im =
                                        re * sin(value) + im * cos(value);
                                }


                                //fpxypos = fopen("xypos3.log", "a");
                                //fprintf(fpxypos, "%5ld %4ld    %10.8f %10.8f      %5ld %10.8f %5ld %10.8f    %10.8f %10.8f  %.18g        %.18g   %.18g   %.18g   %.18g\n", vindex, layer, xpos[layer], ypos[layer], iim, x, jjm, y, 1.0*iim+x, 1.0*jjm+y, value, data.image[ID_TML[layer]].array.F[jjm*naxes_MASTER[0]+iim], data.image[ID_TML[layer]].array.F[jjm1*naxes_MASTER[0]+iim], data.image[ID_TML[layer]].array.F[jjm1*naxes_MASTER[0]+iim1], data.image[ID_TML[layer]].array.F[jjm*naxes_MASTER[0]+iim1]);
                                //fclose(fpxypos);
                            }
                    }

                }
                else // double precision
                {
                    if(BICUBIC == 0)
                    {
                        // bilinear interpolation
                        for(uint32_t ii = 0; ii < naxes[0]; ii++)
                            for(uint32_t jj = 0; jj < naxes[1]; jj++)
                            {
                                double iimf = fmod((xpos[layer] + ii), 1.0 * naxes_MASTER[0]);
                                double jjmf = fmod((ypos[layer] + jj), 1.0 * naxes_MASTER[1]);
                                long iim = (long)(iimf);
                                long jjm = (long)(jjmf);
                                double iifrac = iimf - iim;
                                double jjfrac = jjmf - jjm;
                                long iim1 = iim + 1;
                                long jjm1 = jjm + 1;
                                if(iim == atmturbconf.MasterSize)
                                {
                                    iim = 0;
                                }
                                if(jjm == atmturbconf.MasterSize)
                                {
                                    jjm = 0;
                                }
                                if(iim1 > atmturbconf.MasterSize - 1)
                                {
                                    iim1 -= atmturbconf.MasterSize;
                                }
                                if(jjm1 > atmturbconf.MasterSize - 1)
                                {
                                    jjm1 -= atmturbconf.MasterSize;
                                }

                                double value = (1.0 - iifrac) * (1.0 - jjfrac) *
                                               data.image[ID_TML[layer]].array.D[jjm
                                                       * naxes_MASTER[0] + iim];
                                value += (1.0 - iifrac) * (jjfrac) * data.image[ID_TML[layer]].array.D[jjm1 *
                                         naxes_MASTER[0] + iim];
                                value += (iifrac) * (jjfrac) * data.image[ID_TML[layer]].array.D[jjm1 *
                                         naxes_MASTER[0] + iim1];
                                value += (iifrac) * (1.0 - jjfrac) * data.image[ID_TML[layer]].array.D[jjm *
                                         naxes_MASTER[0] + iim1];

                                data.image[ID_array1].array.D[jj * naxes[0] + ii] += value;
                                if(atmturbconf.flag_WFampl == 1)
                                {
                                    double re = data.image[ID_array2].array.CD[jj * naxes[0] + ii].re;
                                    double im = data.image[ID_array2].array.CD[jj * naxes[0] + ii].im;
                                    data.image[ID_array2].array.CD[jj * naxes[0] + ii].re =
                                        re * cos(value) - im * sin(value);
                                    data.image[ID_array2].array.CD[jj * naxes[0] + ii].im =
                                        re * sin(value) + im * cos(value);
                                }
                            }
                    }
                    else
                    {
                        // bicubic interpolation
                        for(uint32_t ii = 0; ii < naxes[0]; ii++)
                            for(uint32_t jj = 0; jj < naxes[1]; jj++)
                            {
                                double a00, a01, a02, a03;
                                double a10, a11, a12, a13;
                                double a20, a21, a22, a23;
                                double a30, a31, a32, a33;
                                double p00, p01, p02, p03;
                                double p10, p11, p12, p13;
                                double p20, p21, p22, p23;
                                double p30, p31, p32, p33;


                                double iimf = fmod((xpos[layer] + ii), 1.0 * naxes_MASTER[0]);
                                double jjmf = fmod((ypos[layer] + jj), 1.0 * naxes_MASTER[1]);

                                long iim = (long)(iimf);
                                long jjm = (long)(jjmf);

                                double x = iimf - iim;
                                double y = jjmf - jjm;


                                long iim0 = iim - 1;
                                long iim1 = iim;
                                long iim2 = iim + 1;
                                long iim3 = iim + 2;
                                if(iim1 > atmturbconf.MasterSize - 1)
                                {
                                    iim1 -= atmturbconf.MasterSize;
                                }
                                if(iim2 > atmturbconf.MasterSize - 1)
                                {
                                    iim2 -= atmturbconf.MasterSize;
                                }
                                if(iim3 > atmturbconf.MasterSize - 1)
                                {
                                    iim3 -= atmturbconf.MasterSize;
                                }
                                if(iim0 < 0)
                                {
                                    iim0 += atmturbconf.MasterSize;
                                }

                                long jjm0 = jjm - 1;
                                long jjm1 = jjm;
                                long jjm2 = jjm + 1;
                                long jjm3 = jjm + 2;
                                if(jjm1 > atmturbconf.MasterSize - 1)
                                {
                                    iim1 -= atmturbconf.MasterSize;
                                }
                                if(jjm2 > atmturbconf.MasterSize - 1)
                                {
                                    jjm2 -= atmturbconf.MasterSize;
                                }
                                if(jjm3 > atmturbconf.MasterSize - 1)
                                {
                                    jjm3 -= atmturbconf.MasterSize;
                                }
                                if(jjm0 < 0)
                                {
                                    jjm0 += atmturbconf.MasterSize;
                                }

                                /*assert(iim0>=0);
                                assert(iim1>=0);
                                assert(iim2>=0);
                                assert(iim3>=0);
                                assert(iim0<MasterSize);
                                assert(iim1<MasterSize);
                                assert(iim2<MasterSize);
                                assert(iim3<MasterSize);
                                assert(jjm0>=0);
                                assert(jjm1>=0);
                                assert(jjm2>=0);
                                assert(jjm3>=0);
                                assert(jjm0<MasterSize);
                                assert(jjm1<MasterSize);
                                assert(jjm2<MasterSize);
                                assert(jjm3<MasterSize);
                                */
                                p00 = data.image[ID_TML[layer]].array.D[jjm0 * naxes_MASTER[0] + iim0];
                                p01 = data.image[ID_TML[layer]].array.D[jjm1 * naxes_MASTER[0] + iim0];
                                p02 = data.image[ID_TML[layer]].array.D[jjm2 * naxes_MASTER[0] + iim0];
                                p03 = data.image[ID_TML[layer]].array.D[jjm3 * naxes_MASTER[0] + iim0];

                                p10 = data.image[ID_TML[layer]].array.D[jjm0 * naxes_MASTER[0] + iim1];
                                p11 = data.image[ID_TML[layer]].array.D[jjm1 * naxes_MASTER[0] + iim1];
                                p12 = data.image[ID_TML[layer]].array.D[jjm2 * naxes_MASTER[0] + iim1];
                                p13 = data.image[ID_TML[layer]].array.D[jjm3 * naxes_MASTER[0] + iim1];

                                p20 = data.image[ID_TML[layer]].array.D[jjm0 * naxes_MASTER[0] + iim2];
                                p21 = data.image[ID_TML[layer]].array.D[jjm1 * naxes_MASTER[0] + iim2];
                                p22 = data.image[ID_TML[layer]].array.D[jjm2 * naxes_MASTER[0] + iim2];
                                p23 = data.image[ID_TML[layer]].array.D[jjm3 * naxes_MASTER[0] + iim2];

                                p30 = data.image[ID_TML[layer]].array.D[jjm0 * naxes_MASTER[0] + iim3];
                                p31 = data.image[ID_TML[layer]].array.D[jjm1 * naxes_MASTER[0] + iim3];
                                p32 = data.image[ID_TML[layer]].array.D[jjm2 * naxes_MASTER[0] + iim3];
                                p33 = data.image[ID_TML[layer]].array.D[jjm3 * naxes_MASTER[0] + iim3];


                                a00 = p11;
                                a01 = -.5 * p10 + .5 * p12;
                                a02 = p10 - 2.5 * p11 + 2 * p12 - .5 * p13;
                                a03 = -.5 * p10 + 1.5 * p11 - 1.5 * p12 + .5 * p13;
                                a10 = -.5 * p01 + .5 * p21;
                                a11 = .25 * p00 - .25 * p02 - .25 * p20 + .25 * p22;
                                a12 = -.5 * p00 + 1.25 * p01 - p02 + .25 * p03 + .5 * p20 - 1.25 * p21 + p22 -
                                      .25 * p23;
                                a13 = .25 * p00 - .75 * p01 + .75 * p02 - .25 * p03 - .25 * p20 + .75 * p21 -
                                      .75 * p22 + .25 * p23;
                                a20 = p01 - 2.5 * p11 + 2 * p21 - .5 * p31;
                                a21 = -.5 * p00 + .5 * p02 + 1.25 * p10 - 1.25 * p12 - p20 + p22 + .25 * p30 -
                                      .25 * p32;
                                a22 = p00 - 2.5 * p01 + 2 * p02 - .5 * p03 - 2.5 * p10 + 6.25 * p11 - 5 * p12 +
                                      1.25 * p13 + 2 * p20 - 5 * p21 + 4 * p22 - p23 - .5 * p30 + 1.25 * p31 - p32 +
                                      .25 * p33;
                                a23 = -.5 * p00 + 1.5 * p01 - 1.5 * p02 + .5 * p03 + 1.25 * p10 - 3.75 * p11 +
                                      3.75 * p12 - 1.25 * p13 - p20 + 3 * p21 - 3 * p22 + p23 + .25 * p30 - .75 * p31
                                      + .75 * p32 - .25 * p33;
                                a30 = -.5 * p01 + 1.5 * p11 - 1.5 * p21 + .5 * p31;
                                a31 = .25 * p00 - .25 * p02 - .75 * p10 + .75 * p12 + .75 * p20 - .75 * p22 -
                                      .25 * p30 + .25 * p32;
                                a32 = -.5 * p00 + 1.25 * p01 - p02 + .25 * p03 + 1.5 * p10 - 3.75 * p11 + 3 *
                                      p12 - .75 * p13 - 1.5 * p20 + 3.75 * p21 - 3 * p22 + .75 * p23 + .5 * p30 - 1.25
                                      * p31 + p32 - .25 * p33;
                                a33 = .25 * p00 - .75 * p01 + .75 * p02 - .25 * p03 - .75 * p10 + 2.25 * p11 -
                                      2.25 * p12 + .75 * p13 + .75 * p20 - 2.25 * p21 + 2.25 * p22 - .75 * p23 - .25 *
                                      p30 + .75 * p31 - .75 * p32 + .25 * p33;

                                double x2 = x * x;
                                double x3 = x2 * x;
                                double y2 = y * y;
                                double y3 = y2 * y;

                                double value = (a00 + a01 * y + a02 * y2 + a03 * y3) + (a10 + a11 * y + a12 * y2
                                               + a13
                                               * y3) * x + (a20 + a21 * y + a22 * y2 + a23 * y3) * x2 +
                                               (a30 + a31 * y + a32 * y2 + a33 * y3) * x3;


                                data.image[ID_array1].array.D[jj * naxes[0] + ii] += value;
                                if(atmturbconf.flag_WFampl == 1)
                                {
                                    double re = data.image[ID_array2].array.CD[jj * naxes[0] + ii].re;
                                    double im = data.image[ID_array2].array.CD[jj * naxes[0] + ii].im;
                                    data.image[ID_array2].array.CD[jj * naxes[0] + ii].re =
                                        re * cos(value) - im * sin(value);
                                    data.image[ID_array2].array.CD[jj * naxes[0] + ii].im =
                                        re * sin(value) + im * cos(value);
                                }
                            }
                    }
                }



                /* make swavefront */
                if(atmturbconf.flag_SWF_make == 1)
                {
                    if(WFprecision == 0)
                    {
                        if(BICUBIC == 0)
                        {
                            for(uint32_t ii = 0; ii < naxes[0]; ii++)
                                for(uint32_t jj = 0; jj < naxes[1]; jj++)
                                {
                                    double iimf = fmod((xpos[layer] + ii), 1.0 * naxes_MASTER[0]);
                                    double jjmf = fmod((ypos[layer] + jj), 1.0 * naxes_MASTER[1]);
                                    long iim = (long)(iimf);
                                    long jjm = (long)(jjmf);
                                    double iifrac = iimf - iim;
                                    double jjfrac = jjmf - jjm;
                                    long iim1 = iim + 1;
                                    long jjm1 = jjm + 1;
                                    if(iim == atmturbconf.MasterSize)
                                    {
                                        iim = 0;
                                    }
                                    if(jjm == atmturbconf.MasterSize)
                                    {
                                        jjm = 0;
                                    }
                                    if(iim1 > atmturbconf.MasterSize - 1)
                                    {
                                        iim1 -= atmturbconf.MasterSize;
                                    }
                                    if(jjm1 > atmturbconf.MasterSize - 1)
                                    {
                                        jjm1 -= atmturbconf.MasterSize;
                                    }

                                    double value = (1.0 - iifrac) * (1.0 - jjfrac) *
                                                   data.image[ID_TML[layer]].array.F[jjm
                                                           * naxes_MASTER[0] + iim];
                                    value += (1.0 - iifrac) * (jjfrac) * data.image[ID_TML[layer]].array.F[jjm1 *
                                             naxes_MASTER[0] + iim];
                                    value += (iifrac) * (jjfrac) * data.image[ID_TML[layer]].array.F[jjm1 *
                                             naxes_MASTER[0] + iim1];
                                    value += (iifrac) * (1.0 - jjfrac) * data.image[ID_TML[layer]].array.F[jjm *
                                             naxes_MASTER[0] + iim1];

                                    value *= Scoeff;  // multiplicative coeff to go from ref lambda to science lambda

                                    data.image[ID_sarray1].array.F[jj * naxes[0] + ii] += value;

                                    if(atmturbconf.flag_WFampl == 1)
                                    {
                                        float re = data.image[ID_sarray2].array.CF[jj * naxes[0] + ii].re;
                                        float im = data.image[ID_sarray2].array.CF[jj * naxes[0] + ii].im;
                                        data.image[ID_sarray2].array.CF[jj * naxes[0] + ii].re =
                                            re * cos(value) - im * sin(value);
                                        data.image[ID_sarray2].array.CF[jj * naxes[0] + ii].im =
                                            re * sin(value) + im * cos(value);
                                    }
                                }
                        }
                        else
                        {
                            // bicubic interpolation
                            for(uint32_t ii = 0; ii < naxes[0]; ii++)
                                for(uint32_t jj = 0; jj < naxes[1]; jj++)
                                {
                                    double a00, a01, a02, a03;
                                    double a10, a11, a12, a13;
                                    double a20, a21, a22, a23;
                                    double a30, a31, a32, a33;
                                    double p00, p01, p02, p03;
                                    double p10, p11, p12, p13;
                                    double p20, p21, p22, p23;
                                    double p30, p31, p32, p33;


                                    double iimf = fmod((xpos[layer] + ii), 1.0 * naxes_MASTER[0]);
                                    double jjmf = fmod((ypos[layer] + jj), 1.0 * naxes_MASTER[1]);

                                    long iim = (long)(iimf);
                                    long jjm = (long)(jjmf);

                                    double x = iimf - iim;
                                    double y = jjmf - jjm;


                                    long iim0 = iim - 1;
                                    long iim1 = iim;
                                    long iim2 = iim + 1;
                                    long iim3 = iim + 2;
                                    if(iim1 > atmturbconf.MasterSize - 1)
                                    {
                                        iim1 -= atmturbconf.MasterSize;
                                    }
                                    if(iim2 > atmturbconf.MasterSize - 1)
                                    {
                                        iim2 -= atmturbconf.MasterSize;
                                    }
                                    if(iim3 > atmturbconf.MasterSize - 1)
                                    {
                                        iim3 -= atmturbconf.MasterSize;
                                    }
                                    if(iim0 < 0)
                                    {
                                        iim0 += atmturbconf.MasterSize;
                                    }

                                    long jjm0 = jjm - 1;
                                    long jjm1 = jjm;
                                    long jjm2 = jjm + 1;
                                    long jjm3 = jjm + 2;
                                    if(jjm1 > atmturbconf.MasterSize - 1)
                                    {
                                        iim1 -= atmturbconf.MasterSize;
                                    }
                                    if(jjm2 > atmturbconf.MasterSize - 1)
                                    {
                                        jjm2 -= atmturbconf.MasterSize;
                                    }
                                    if(jjm3 > atmturbconf.MasterSize - 1)
                                    {
                                        jjm3 -= atmturbconf.MasterSize;
                                    }
                                    if(jjm0 < 0)
                                    {
                                        jjm0 += atmturbconf.MasterSize;
                                    }

                                    /*assert(iim0>=0);
                                    assert(iim1>=0);
                                    assert(iim2>=0);
                                    assert(iim3>=0);
                                    assert(iim0<MasterSize);
                                    assert(iim1<MasterSize);
                                    assert(iim2<MasterSize);
                                    assert(iim3<MasterSize);
                                    assert(jjm0>=0);
                                    assert(jjm1>=0);
                                    assert(jjm2>=0);
                                    assert(jjm3>=0);
                                    assert(jjm0<MasterSize);
                                    assert(jjm1<MasterSize);
                                    assert(jjm2<MasterSize);
                                    assert(jjm3<MasterSize);
                                    */
                                    p00 = data.image[ID_TML[layer]].array.F[jjm0 * naxes_MASTER[0] + iim0];
                                    p01 = data.image[ID_TML[layer]].array.F[jjm1 * naxes_MASTER[0] + iim0];
                                    p02 = data.image[ID_TML[layer]].array.F[jjm2 * naxes_MASTER[0] + iim0];
                                    p03 = data.image[ID_TML[layer]].array.F[jjm3 * naxes_MASTER[0] + iim0];

                                    p10 = data.image[ID_TML[layer]].array.F[jjm0 * naxes_MASTER[0] + iim1];
                                    p11 = data.image[ID_TML[layer]].array.F[jjm1 * naxes_MASTER[0] + iim1];
                                    p12 = data.image[ID_TML[layer]].array.F[jjm2 * naxes_MASTER[0] + iim1];
                                    p13 = data.image[ID_TML[layer]].array.F[jjm3 * naxes_MASTER[0] + iim1];

                                    p20 = data.image[ID_TML[layer]].array.F[jjm0 * naxes_MASTER[0] + iim2];
                                    p21 = data.image[ID_TML[layer]].array.F[jjm1 * naxes_MASTER[0] + iim2];
                                    p22 = data.image[ID_TML[layer]].array.F[jjm2 * naxes_MASTER[0] + iim2];
                                    p23 = data.image[ID_TML[layer]].array.F[jjm3 * naxes_MASTER[0] + iim2];

                                    p30 = data.image[ID_TML[layer]].array.F[jjm0 * naxes_MASTER[0] + iim3];
                                    p31 = data.image[ID_TML[layer]].array.F[jjm1 * naxes_MASTER[0] + iim3];
                                    p32 = data.image[ID_TML[layer]].array.F[jjm2 * naxes_MASTER[0] + iim3];
                                    p33 = data.image[ID_TML[layer]].array.F[jjm3 * naxes_MASTER[0] + iim3];


                                    a00 = p11;
                                    a01 = -.5 * p10 + .5 * p12;
                                    a02 = p10 - 2.5 * p11 + 2 * p12 - .5 * p13;
                                    a03 = -.5 * p10 + 1.5 * p11 - 1.5 * p12 + .5 * p13;
                                    a10 = -.5 * p01 + .5 * p21;
                                    a11 = .25 * p00 - .25 * p02 - .25 * p20 + .25 * p22;
                                    a12 = -.5 * p00 + 1.25 * p01 - p02 + .25 * p03 + .5 * p20 - 1.25 * p21 + p22 -
                                          .25 * p23;
                                    a13 = .25 * p00 - .75 * p01 + .75 * p02 - .25 * p03 - .25 * p20 + .75 * p21 -
                                          .75 * p22 + .25 * p23;
                                    a20 = p01 - 2.5 * p11 + 2 * p21 - .5 * p31;
                                    a21 = -.5 * p00 + .5 * p02 + 1.25 * p10 - 1.25 * p12 - p20 + p22 + .25 * p30 -
                                          .25 * p32;
                                    a22 = p00 - 2.5 * p01 + 2 * p02 - .5 * p03 - 2.5 * p10 + 6.25 * p11 - 5 * p12 +
                                          1.25 * p13 + 2 * p20 - 5 * p21 + 4 * p22 - p23 - .5 * p30 + 1.25 * p31 - p32 +
                                          .25 * p33;
                                    a23 = -.5 * p00 + 1.5 * p01 - 1.5 * p02 + .5 * p03 + 1.25 * p10 - 3.75 * p11 +
                                          3.75 * p12 - 1.25 * p13 - p20 + 3 * p21 - 3 * p22 + p23 + .25 * p30 - .75 * p31
                                          + .75 * p32 - .25 * p33;
                                    a30 = -.5 * p01 + 1.5 * p11 - 1.5 * p21 + .5 * p31;
                                    a31 = .25 * p00 - .25 * p02 - .75 * p10 + .75 * p12 + .75 * p20 - .75 * p22 -
                                          .25 * p30 + .25 * p32;
                                    a32 = -.5 * p00 + 1.25 * p01 - p02 + .25 * p03 + 1.5 * p10 - 3.75 * p11 + 3 *
                                          p12 - .75 * p13 - 1.5 * p20 + 3.75 * p21 - 3 * p22 + .75 * p23 + .5 * p30 - 1.25
                                          * p31 + p32 - .25 * p33;
                                    a33 = .25 * p00 - .75 * p01 + .75 * p02 - .25 * p03 - .75 * p10 + 2.25 * p11 -
                                          2.25 * p12 + .75 * p13 + .75 * p20 - 2.25 * p21 + 2.25 * p22 - .75 * p23 - .25 *
                                          p30 + .75 * p31 - .75 * p32 + .25 * p33;

                                    double x2 = x * x;
                                    double x3 = x2 * x;
                                    double y2 = y * y;
                                    double y3 = y2 * y;

                                    double value = (a00 + a01 * y + a02 * y2 + a03 * y3) + (a10 + a11 * y + a12 * y2
                                                   + a13
                                                   * y3) * x + (a20 + a21 * y + a22 * y2 + a23 * y3) * x2 +
                                                   (a30 + a31 * y + a32 * y2 + a33 * y3) * x3;


                                    value *= Scoeff;  // multiplicative coeff to go from ref lambda to science lambda

                                    data.image[ID_sarray1].array.F[jj * naxes[0] + ii] += value;

                                    if(atmturbconf.flag_WFampl == 1)
                                    {
                                        float re = data.image[ID_sarray2].array.CF[jj * naxes[0] + ii].re;
                                        float im = data.image[ID_sarray2].array.CF[jj * naxes[0] + ii].im;
                                        data.image[ID_sarray2].array.CF[jj * naxes[0] + ii].re =
                                            re * cos(value) - im * sin(value);
                                        data.image[ID_sarray2].array.CF[jj * naxes[0] + ii].im =
                                            re * sin(value) + im * cos(value);
                                    }
                                }
                        }
                    }
                    else // double precision
                    {
                        if(BICUBIC == 0)
                        {
                            for(uint32_t ii = 0; ii < naxes[0]; ii++)
                                for(uint32_t jj = 0; jj < naxes[1]; jj++)
                                {
                                    double iimf = fmod((xpos[layer] + ii), 1.0 * naxes_MASTER[0]);
                                    double jjmf = fmod((ypos[layer] + jj), 1.0 * naxes_MASTER[1]);
                                    long iim = (long)(iimf);
                                    long jjm = (long)(jjmf);
                                    double iifrac = iimf - iim;
                                    double jjfrac = jjmf - jjm;
                                    long iim1 = iim + 1;
                                    long jjm1 = jjm + 1;
                                    if(iim == atmturbconf.MasterSize)
                                    {
                                        iim = 0;
                                    }
                                    if(jjm == atmturbconf.MasterSize)
                                    {
                                        jjm = 0;
                                    }
                                    if(iim1 > atmturbconf.MasterSize - 1)
                                    {
                                        iim1 -= atmturbconf.MasterSize;
                                    }
                                    if(jjm1 > atmturbconf.MasterSize - 1)
                                    {
                                        jjm1 -= atmturbconf.MasterSize;
                                    }

                                    double value = (1.0 - iifrac) * (1.0 - jjfrac) *
                                                   data.image[ID_TML[layer]].array.D[jjm
                                                           * naxes_MASTER[0] + iim];
                                    value += (1.0 - iifrac) * (jjfrac) * data.image[ID_TML[layer]].array.D[jjm1 *
                                             naxes_MASTER[0] + iim];
                                    value += (iifrac) * (jjfrac) * data.image[ID_TML[layer]].array.D[jjm1 *
                                             naxes_MASTER[0] + iim1];
                                    value += (iifrac) * (1.0 - jjfrac) * data.image[ID_TML[layer]].array.D[jjm *
                                             naxes_MASTER[0] + iim1];

                                    value *= Scoeff;  // multiplicative coeff to go from ref lambda to science lambda

                                    data.image[ID_sarray1].array.D[jj * naxes[0] + ii] += value;

                                    if(atmturbconf.flag_WFampl == 1)
                                    {
                                        double re = data.image[ID_sarray2].array.CD[jj * naxes[0] + ii].re;
                                        double im = data.image[ID_sarray2].array.CD[jj * naxes[0] + ii].im;
                                        data.image[ID_sarray2].array.CD[jj * naxes[0] + ii].re =
                                            re * cos(value) - im * sin(value);
                                        data.image[ID_sarray2].array.CD[jj * naxes[0] + ii].im =
                                            re * sin(value) + im * cos(value);
                                    }
                                }
                        }
                        else
                        {
                            // bicubic interpolation
                            for(uint32_t ii = 0; ii < naxes[0]; ii++)
                                for(uint32_t jj = 0; jj < naxes[1]; jj++)
                                {
                                    double a00, a01, a02, a03;
                                    double a10, a11, a12, a13;
                                    double a20, a21, a22, a23;
                                    double a30, a31, a32, a33;
                                    double p00, p01, p02, p03;
                                    double p10, p11, p12, p13;
                                    double p20, p21, p22, p23;
                                    double p30, p31, p32, p33;


                                    double iimf = fmod((xpos[layer] + ii), 1.0 * naxes_MASTER[0]);
                                    double jjmf = fmod((ypos[layer] + jj), 1.0 * naxes_MASTER[1]);

                                    long iim = (long)(iimf);
                                    long jjm = (long)(jjmf);

                                    double x = iimf - iim;
                                    double y = jjmf - jjm;


                                    long iim0 = iim - 1;
                                    long iim1 = iim;
                                    long iim2 = iim + 1;
                                    long iim3 = iim + 2;
                                    if(iim1 > atmturbconf.MasterSize - 1)
                                    {
                                        iim1 -= atmturbconf.MasterSize;
                                    }
                                    if(iim2 > atmturbconf.MasterSize - 1)
                                    {
                                        iim2 -= atmturbconf.MasterSize;
                                    }
                                    if(iim3 > atmturbconf.MasterSize - 1)
                                    {
                                        iim3 -= atmturbconf.MasterSize;
                                    }
                                    if(iim0 < 0)
                                    {
                                        iim0 += atmturbconf.MasterSize;
                                    }

                                    long jjm0 = jjm - 1;
                                    long jjm1 = jjm;
                                    long jjm2 = jjm + 1;
                                    long jjm3 = jjm + 2;
                                    if(jjm1 > atmturbconf.MasterSize - 1)
                                    {
                                        iim1 -= atmturbconf.MasterSize;
                                    }
                                    if(jjm2 > atmturbconf.MasterSize - 1)
                                    {
                                        jjm2 -= atmturbconf.MasterSize;
                                    }
                                    if(jjm3 > atmturbconf.MasterSize - 1)
                                    {
                                        jjm3 -= atmturbconf.MasterSize;
                                    }
                                    if(jjm0 < 0)
                                    {
                                        jjm0 += atmturbconf.MasterSize;
                                    }

                                    /*assert(iim0>=0);
                                    assert(iim1>=0);
                                    assert(iim2>=0);
                                    assert(iim3>=0);
                                    assert(iim0<MasterSize);
                                    assert(iim1<MasterSize);
                                    assert(iim2<MasterSize);
                                    assert(iim3<MasterSize);
                                    assert(jjm0>=0);
                                    assert(jjm1>=0);
                                    assert(jjm2>=0);
                                    assert(jjm3>=0);
                                    assert(jjm0<MasterSize);
                                    assert(jjm1<MasterSize);
                                    assert(jjm2<MasterSize);
                                    assert(jjm3<MasterSize);
                                    */


                                    p00 = data.image[ID_TML[layer]].array.D[jjm0 * naxes_MASTER[0] + iim0];
                                    p01 = data.image[ID_TML[layer]].array.D[jjm1 * naxes_MASTER[0] + iim0];
                                    p02 = data.image[ID_TML[layer]].array.D[jjm2 * naxes_MASTER[0] + iim0];
                                    p03 = data.image[ID_TML[layer]].array.D[jjm3 * naxes_MASTER[0] + iim0];

                                    p10 = data.image[ID_TML[layer]].array.D[jjm0 * naxes_MASTER[0] + iim1];
                                    p11 = data.image[ID_TML[layer]].array.D[jjm1 * naxes_MASTER[0] + iim1];
                                    p12 = data.image[ID_TML[layer]].array.D[jjm2 * naxes_MASTER[0] + iim1];
                                    p13 = data.image[ID_TML[layer]].array.D[jjm3 * naxes_MASTER[0] + iim1];

                                    p20 = data.image[ID_TML[layer]].array.D[jjm0 * naxes_MASTER[0] + iim2];
                                    p21 = data.image[ID_TML[layer]].array.D[jjm1 * naxes_MASTER[0] + iim2];
                                    p22 = data.image[ID_TML[layer]].array.D[jjm2 * naxes_MASTER[0] + iim2];
                                    p23 = data.image[ID_TML[layer]].array.D[jjm3 * naxes_MASTER[0] + iim2];

                                    p30 = data.image[ID_TML[layer]].array.D[jjm0 * naxes_MASTER[0] + iim3];
                                    p31 = data.image[ID_TML[layer]].array.D[jjm1 * naxes_MASTER[0] + iim3];
                                    p32 = data.image[ID_TML[layer]].array.D[jjm2 * naxes_MASTER[0] + iim3];
                                    p33 = data.image[ID_TML[layer]].array.D[jjm3 * naxes_MASTER[0] + iim3];


                                    a00 = p11;
                                    a01 = -.5 * p10 + .5 * p12;
                                    a02 = p10 - 2.5 * p11 + 2 * p12 - .5 * p13;
                                    a03 = -.5 * p10 + 1.5 * p11 - 1.5 * p12 + .5 * p13;
                                    a10 = -.5 * p01 + .5 * p21;
                                    a11 = .25 * p00 - .25 * p02 - .25 * p20 + .25 * p22;
                                    a12 = -.5 * p00 + 1.25 * p01 - p02 + .25 * p03 + .5 * p20 - 1.25 * p21 + p22 -
                                          .25 * p23;
                                    a13 = .25 * p00 - .75 * p01 + .75 * p02 - .25 * p03 - .25 * p20 + .75 * p21 -
                                          .75 * p22 + .25 * p23;
                                    a20 = p01 - 2.5 * p11 + 2 * p21 - .5 * p31;
                                    a21 = -.5 * p00 + .5 * p02 + 1.25 * p10 - 1.25 * p12 - p20 + p22 + .25 * p30 -
                                          .25 * p32;
                                    a22 = p00 - 2.5 * p01 + 2 * p02 - .5 * p03 - 2.5 * p10 + 6.25 * p11 - 5 * p12 +
                                          1.25 * p13 + 2 * p20 - 5 * p21 + 4 * p22 - p23 - .5 * p30 + 1.25 * p31 - p32 +
                                          .25 * p33;
                                    a23 = -.5 * p00 + 1.5 * p01 - 1.5 * p02 + .5 * p03 + 1.25 * p10 - 3.75 * p11 +
                                          3.75 * p12 - 1.25 * p13 - p20 + 3 * p21 - 3 * p22 + p23 + .25 * p30 - .75 * p31
                                          + .75 * p32 - .25 * p33;
                                    a30 = -.5 * p01 + 1.5 * p11 - 1.5 * p21 + .5 * p31;
                                    a31 = .25 * p00 - .25 * p02 - .75 * p10 + .75 * p12 + .75 * p20 - .75 * p22 -
                                          .25 * p30 + .25 * p32;
                                    a32 = -.5 * p00 + 1.25 * p01 - p02 + .25 * p03 + 1.5 * p10 - 3.75 * p11 + 3 *
                                          p12 - .75 * p13 - 1.5 * p20 + 3.75 * p21 - 3 * p22 + .75 * p23 + .5 * p30 - 1.25
                                          * p31 + p32 - .25 * p33;
                                    a33 = .25 * p00 - .75 * p01 + .75 * p02 - .25 * p03 - .75 * p10 + 2.25 * p11 -
                                          2.25 * p12 + .75 * p13 + .75 * p20 - 2.25 * p21 + 2.25 * p22 - .75 * p23 - .25 *
                                          p30 + .75 * p31 - .75 * p32 + .25 * p33;

                                    double x2 = x * x;
                                    double x3 = x2 * x;
                                    double y2 = y * y;
                                    double y3 = y2 * y;

                                    double value = (a00 + a01 * y + a02 * y2 + a03 * y3) + (a10 + a11 * y + a12 * y2
                                                   + a13
                                                   * y3) * x + (a20 + a21 * y + a22 * y2 + a23 * y3) * x2 +
                                                   (a30 + a31 * y + a32 * y2 + a33 * y3) * x3;


                                    value *= Scoeff;  // multiplicative coeff to go from ref lambda to science lambda

                                    data.image[ID_sarray1].array.D[jj * naxes[0] + ii] += value;

                                    if(atmturbconf.flag_WFampl == 1)
                                    {
                                        double re = data.image[ID_sarray2].array.CD[jj * naxes[0] + ii].re;
                                        double im = data.image[ID_sarray2].array.CD[jj * naxes[0] + ii].im;
                                        data.image[ID_sarray2].array.CD[jj * naxes[0] + ii].re =
                                            re * cos(value) - im * sin(value);
                                        data.image[ID_sarray2].array.CD[jj * naxes[0] + ii].im =
                                            re * sin(value) + im * cos(value);
                                    }
                                }
                        }
                    }
                }
            }




            // REFERENCE LAMBDA
            if(WFprecision == 0)
            {
                for(uint32_t ii = 0; ii < naxesout[0]; ii++)
                    for(uint32_t jj = 0; jj < naxesout[1]; jj++)
                    {
                        long ii1 = ii + (naxes[0] - naxesout[0]) / 2;
                        long jj1 = jj + (naxes[1] - naxesout[1]) / 2;
                        data.image[IDout_array_pha].array.F[frame * naxesout[0]*naxesout[1] + jj *
                                                            naxesout[0] + ii] = data.image[ID_array1].array.F[jj1 * naxes[0] + ii1];
                    }
            }
            else
            {
                for(uint32_t ii = 0; ii < naxesout[0]; ii++)
                    for(uint32_t jj = 0; jj < naxesout[1]; jj++)
                    {
                        long ii1 = ii + (naxes[0] - naxesout[0]) / 2;
                        long jj1 = jj + (naxes[1] - naxesout[1]) / 2;
                        data.image[IDout_array_pha].array.D[frame * naxesout[0]*naxesout[1] + jj *
                                                            naxesout[0] + ii] = data.image[ID_array1].array.D[jj1 * naxes[0] + ii1];
                    }
            }

            if(WFprecision == 0)
            {
                if(atmturbconf.flag_WFampl == 1)
                {
                    for(uint32_t ii = 0; ii < naxesout[0]; ii++)
                        for(uint32_t jj = 0; jj < naxesout[1]; jj++)
                        {
                            long ii1 = ii + (naxes[0] - naxesout[0]) / 2;
                            long jj1 = jj + (naxes[1] - naxesout[1]) / 2;
                            array[frame * naxesout[0]*naxesout[1] + jj * naxesout[0] + ii].re =
                                data.image[ID_array2].array.CF[jj1 * naxes[0] + ii1].re;
                            array[frame * naxesout[0]*naxesout[1] + jj * naxesout[0] + ii].im =
                                data.image[ID_array2].array.CF[jj1 * naxes[0] + ii1].im;

                            float re = array[frame * naxesout[0] * naxesout[1] + jj * naxesout[0] + ii].re;
                            float im = array[frame * naxesout[0] * naxesout[1] + jj * naxesout[0] + ii].im;
                            data.image[IDout_array_amp].array.F[frame * naxesout[0]*naxesout[1] + jj *
                                                                naxesout[0] + ii] = sqrt(re * re + im * im);
                            float pha = atan2(im, re);
                            data.image[IDout_array_pha].array.F[frame * naxesout[0]*naxesout[1] + jj *
                                                                naxesout[0] + ii] = pha + 2.0 * M_PI * ((long)(
                                                                        data.image[IDout_array_pha].array.F[frame * naxesout[0] * naxesout[1] + jj *
                                                                                naxesout[0] + ii] / 2.0 / M_PI + 1000.5) - 1000.0);
                        }
                }
            }
            else
            {
                for(uint32_t ii = 0; ii < naxesout[0]; ii++)
                    for(uint32_t jj = 0; jj < naxesout[1]; jj++)
                    {
                        long ii1 = ii + (naxes[0] - naxesout[0]) / 2;
                        long jj1 = jj + (naxes[1] - naxesout[1]) / 2;
                        array_double[frame * naxesout[0]*naxesout[1] + jj * naxesout[0] + ii].re =
                            data.image[ID_array2].array.CD[jj1 * naxes[0] + ii1].re;
                        array_double[frame * naxesout[0]*naxesout[1] + jj * naxesout[0] + ii].im =
                            data.image[ID_array2].array.CD[jj1 * naxes[0] + ii1].im;

                        double re = array_double[frame * naxesout[0] * naxesout[1] + jj * naxesout[0] +
                                                       ii].re;
                        double im = array_double[frame * naxesout[0] * naxesout[1] + jj * naxesout[0] +
                                                       ii].im;
                        data.image[IDout_array_amp].array.D[frame * naxesout[0]*naxesout[1] + jj *
                                                            naxesout[0] + ii] = sqrt(re * re + im * im);
                        double pha = atan2(im, re);
                        data.image[IDout_array_pha].array.D[frame * naxesout[0]*naxesout[1] + jj *
                                                            naxesout[0] + ii] = pha + 2.0 * M_PI * ((long)(
                                                                    data.image[IDout_array_pha].array.D[frame * naxesout[0] * naxesout[1] + jj *
                                                                            naxesout[0] + ii] / 2.0 / M_PI + 1000.5) - 1000.0);
                    }
            }





            // WRITE CURRENT WF TO SHARED MEMORY

            if(atmturbconf.flag_SWF_SHMoutput == 1)
            {
                if(atmturbconf.flag_WFampl == 0)
                {
                    if(WFprecision == 0)
                    {
                        for(uint64_t ii = 0; ii < naxesout[0]*naxesout[1]; ii++)
                        {
                            data.image[IDshmpha].array.F[ii] = data.image[ID_array1].array.F[frame *
                                                               naxesout[0] * naxesout[1] + ii];
                        }
                    }
                    else
                    {
                        for(uint64_t ii = 0; ii < naxesout[0]*naxesout[1]; ii++)
                        {
                            data.image[IDshmpha].array.D[ii] = data.image[ID_array1].array.D[frame *
                                                               naxesout[0] * naxesout[1] + ii];
                        }
                    }
                }
                else
                {
                    if(WFprecision == 0)
                    {
                        for(uint64_t ii = 0; ii < naxesout[0]*naxesout[1]; ii++)
                        {
                            data.image[IDshmpha].array.F[ii] = data.image[IDout_array_pha].array.F[frame *
                                                               naxesout[0] * naxesout[1] + ii];
                            data.image[IDshmamp].array.F[ii] = data.image[IDout_array_amp].array.F[frame *
                                                               naxesout[0] * naxesout[1] + ii];
                        }
                    }
                    else
                    {
                        for(uint64_t ii = 0; ii < naxesout[0]*naxesout[1]; ii++)
                        {
                            data.image[IDshmpha].array.D[ii] = data.image[IDout_array_pha].array.D[frame *
                                                               naxesout[0] * naxesout[1] + ii];
                            data.image[IDshmamp].array.D[ii] = data.image[IDout_array_amp].array.D[frame *
                                                               naxesout[0] * naxesout[1] + ii];
                        }
                    }
                }
            }



            // SCIENCE LAMBDA
            if(atmturbconf.flag_SWF_make == 1)
            {
                if(WFprecision == 0)
                {
                    for(uint32_t ii = 0; ii < naxesout[0]; ii++)
                        for(uint32_t jj = 0; jj < naxesout[1]; jj++)
                        {
                            long ii1 = ii + (naxes[0] - naxesout[0]) / 2;
                            long jj1 = jj + (naxes[1] - naxesout[1]) / 2;
                            data.image[IDout_sarray_pha].array.F[frame * naxesout[0]*naxesout[1] + jj *
                                                                 naxesout[0] + ii] = data.image[ID_sarray1].array.F[jj1 * naxes[0] + ii1];
                        }
                }
                else
                {
                    for(uint32_t ii = 0; ii < naxesout[0]; ii++)
                        for(uint32_t jj = 0; jj < naxesout[1]; jj++)
                        {
                            long ii1 = ii + (naxes[0] - naxesout[0]) / 2;
                            long jj1 = jj + (naxes[1] - naxesout[1]) / 2;
                            data.image[IDout_sarray_pha].array.D[frame * naxesout[0]*naxesout[1] + jj *
                                                                 naxesout[0] + ii] = data.image[ID_sarray1].array.D[jj1 * naxes[0] + ii1];
                        }
                }


                if(atmturbconf.flag_WFampl == 1)
                {
                    if(WFprecision == 0)
                    {
                        for(uint64_t ii2 = 0; ii2 < xsizepeakpha * ysizepeakpha; ii2++)
                        {
                            data.image[IDpeakpha_re_bin].array.F[ii2] = 0.0;
                            data.image[IDpeakpha_im_bin].array.F[ii2] = 0.0;
                            data.image[IDpeakpha_bin].array.F[ii2] = 0.0;
                            data.image[IDpeakpha_bin_ch].array.F[ii2] = 0.0;
                        }
                    }
                    else
                    {
                        for(uint64_t ii2 = 0; ii2 < xsizepeakpha * ysizepeakpha; ii2++)
                        {
                            data.image[IDpeakpha_re_bin].array.D[ii2] = 0.0;
                            data.image[IDpeakpha_im_bin].array.D[ii2] = 0.0;
                            data.image[IDpeakpha_bin].array.D[ii2] = 0.0;
                            data.image[IDpeakpha_bin_ch].array.D[ii2] = 0.0;
                        }
                    }

                    //peakpha_re = 0.0;
                    //peakpha_im = 0.0;

                    if(WFprecision == 0)
                    {
                        for(uint32_t ii = 0; ii < naxesout[0]; ii++)
                            for(uint32_t jj = 0; jj < naxesout[1]; jj++)
                            {
                                long ii1 = ii + (naxes[0] - naxesout[0]) / 2;
                                long jj1 = jj + (naxes[1] - naxesout[1]) / 2;
                                sarray[frame * naxesout[0]*naxesout[1] + jj * naxesout[0] + ii].re =
                                    data.image[ID_sarray2].array.CF[jj1 * naxes[0] + ii1].re;
                                sarray[frame * naxesout[0]*naxesout[1] + jj * naxesout[0] + ii].im =
                                    data.image[ID_sarray2].array.CF[jj1 * naxes[0] + ii1].im;
                                float re = sarray[frame * naxesout[0] * naxesout[1] + jj * naxesout[0] + ii].re;
                                float im = sarray[frame * naxesout[0] * naxesout[1] + jj * naxesout[0] + ii].im;
                                data.image[IDout_sarray_amp].array.F[frame * naxesout[0]*naxesout[1] + jj *
                                                                     naxesout[0] + ii] = sqrt(re * re + im * im);
                                float pha = atan2(im, re);
                                data.image[IDpeakpha_re].array.F[jj * naxesout[0] + ii] = cos(
                                            data.image[IDout_sarray_pha].array.F[frame * naxesout[0] * naxesout[1] + jj *
                                                    naxesout[0] + ii] - pha);
                                data.image[IDpeakpha_im].array.F[jj * naxesout[0] + ii] = sin(
                                            data.image[IDout_sarray_pha].array.F[frame * naxesout[0] * naxesout[1] + jj *
                                                    naxesout[0] + ii] - pha);
                                data.image[IDout_sarray_pha].array.F[frame * naxesout[0]*naxesout[1] + jj *
                                                                     naxesout[0] + ii] = pha;
                                long ii2 = (long)(1.0 * ii / naxesout[0] * xsizepeakpha);
                                long jj2 = (long)(1.0 * jj / naxesout[1] * ysizepeakpha);

                                if((ii2 < xsizepeakpha) && (jj2 < ysizepeakpha))
                                {
                                    data.image[IDpeakpha_re_bin].array.F[jj2 * xsizepeakpha + ii2] += cos(
                                                data.image[ID_sarray1].array.F[jj1 * naxes[0] + ii1] - pha);
                                    data.image[IDpeakpha_im_bin].array.F[jj2 * xsizepeakpha + ii2] += sin(
                                                data.image[ID_sarray1].array.F[jj1 * naxes[0] + ii1] - pha);
                                }
                            }
                    }
                    else
                    {
                        for(uint32_t ii = 0; ii < naxesout[0]; ii++)
                            for(uint32_t jj = 0; jj < naxesout[1]; jj++)
                            {
                                long ii1 = ii + (naxes[0] - naxesout[0]) / 2;
                                long jj1 = jj + (naxes[1] - naxesout[1]) / 2;
                                sarray_double[frame * naxesout[0]*naxesout[1] + jj * naxesout[0] + ii].re =
                                    data.image[ID_sarray2].array.CD[jj1 * naxes[0] + ii1].re;
                                sarray_double[frame * naxesout[0]*naxesout[1] + jj * naxesout[0] + ii].im =
                                    data.image[ID_sarray2].array.CD[jj1 * naxes[0] + ii1].im;
                                double re = sarray_double[frame * naxesout[0] * naxesout[1] + jj * naxesout[0] +
                                                                ii].re;
                                double im = sarray_double[frame * naxesout[0] * naxesout[1] + jj * naxesout[0] +
                                                                ii].im;
                                data.image[IDout_sarray_amp].array.D[frame * naxesout[0]*naxesout[1] + jj *
                                                                     naxesout[0] + ii] = sqrt(re * re + im * im);
                                double pha = atan2(im, re);
                                data.image[IDpeakpha_re].array.D[jj * naxesout[0] + ii] = cos(
                                            data.image[IDout_sarray_pha].array.D[frame * naxesout[0] * naxesout[1] + jj *
                                                    naxesout[0] + ii] - pha);
                                data.image[IDpeakpha_im].array.D[jj * naxesout[0] + ii] = sin(
                                            data.image[IDout_sarray_pha].array.D[frame * naxesout[0] * naxesout[1] + jj *
                                                    naxesout[0] + ii] - pha);
                                data.image[IDout_sarray_pha].array.D[frame * naxesout[0]*naxesout[1] + jj *
                                                                     naxesout[0] + ii] = pha;
                                long ii2 = (long)(1.0 * ii / naxesout[0] * xsizepeakpha);
                                long jj2 = (long)(1.0 * jj / naxesout[1] * ysizepeakpha);

                                if((ii2 < xsizepeakpha) && (jj2 < ysizepeakpha))
                                {
                                    data.image[IDpeakpha_re_bin].array.D[jj2 * xsizepeakpha + ii2] += cos(
                                                data.image[ID_sarray1].array.D[jj1 * naxes[0] + ii1] - pha);
                                    data.image[IDpeakpha_im_bin].array.D[jj2 * xsizepeakpha + ii2] += sin(
                                                data.image[ID_sarray1].array.D[jj1 * naxes[0] + ii1] - pha);
                                }
                            }
                    }


                    //peakpha = atan2(peakpha_im, peakpha_re);
                    //printf("peak pha = %lf\n", peakpha/2.0/M_PI);


                    //peakpha = 0.0;
                    if(WFprecision == 0)
                    {
                        for(uint64_t ii2 = 0; ii2 < xsizepeakpha * ysizepeakpha; ii2++)
                        {
                            data.image[IDpeakpha_bin].array.F[ii2] = atan2(
                                        data.image[IDpeakpha_im_bin].array.F[ii2],
                                        data.image[IDpeakpha_re_bin].array.F[ii2]);
                            //	while(data.image[IDpeakpha_bin].array.F[ii2]<0.0)
                            //	data.image[IDpeakpha_bin].array.F[ii2] += 2.0*M_PI;
                        }
                    }
                    else
                    {
                        for(uint64_t ii2 = 0; ii2 < xsizepeakpha * ysizepeakpha; ii2++)
                        {
                            data.image[IDpeakpha_bin].array.D[ii2] = atan2(
                                        data.image[IDpeakpha_im_bin].array.D[ii2],
                                        data.image[IDpeakpha_re_bin].array.D[ii2]);
                        }
                    }


                    long chcnt = 1;
                    long chcnt0cnt = 0;
                    long chiter = 0;
                    double plim = 1.0;
                    while((plim > 0.51) && (chiter < 1000) && (chcnt0cnt < 5))
                    {
                        chiter++;
                        chcnt = 0;
                        if(WFprecision == 0)
                        {
                            for(uint64_t ii2 = 0; ii2 < xsizepeakpha * ysizepeakpha; ii2++)
                            {
                                data.image[IDpeakpha_bin_ch].array.F[ii2] = 0.0;
                            }
                        }
                        else
                        {
                            for(uint64_t ii2 = 0; ii2 < xsizepeakpha * ysizepeakpha; ii2++)
                            {
                                data.image[IDpeakpha_bin_ch].array.D[ii2] = 0.0;
                            }
                        }

                        if(WFprecision == 0)
                        {
                            for(uint32_t ii2 = 0; ii2 < xsizepeakpha - 1; ii2++)
                                for(uint32_t jj2 = 0; jj2 < ysizepeakpha; jj2++)
                                {
                                    uint64_t index1 = jj2 * xsizepeakpha + ii2;
                                    uint64_t index2 = jj2 * xsizepeakpha + ii2 + 1;
                                    float pv1 = data.image[IDpeakpha_bin].array.F[index1];
                                    float pv2 = data.image[IDpeakpha_bin].array.F[index2];
                                    if(pv2 > pv1 + M_PI)
                                    {
                                        data.image[IDpeakpha_bin_ch].array.F[index1] += 0.2;
                                        data.image[IDpeakpha_bin_ch].array.F[index2] -= 0.2;
                                    }
                                    if(pv1 > pv2 + M_PI)
                                    {
                                        data.image[IDpeakpha_bin_ch].array.F[index1] -= 0.2;
                                        data.image[IDpeakpha_bin_ch].array.F[index2] += 0.2;
                                    }
                                }

                            for(uint32_t ii2 = 0; ii2 < xsizepeakpha; ii2++)
                                for(uint32_t jj2 = 0; jj2 < ysizepeakpha - 1; jj2++)
                                {
                                    uint64_t index1 = jj2 * xsizepeakpha + ii2;
                                    uint64_t index2 = (jj2 + 1) * xsizepeakpha + ii2;
                                    float pv1 = data.image[IDpeakpha_bin].array.F[index1];
                                    float pv2 = data.image[IDpeakpha_bin].array.F[index2];
                                    if(pv2 > pv1 + M_PI)
                                    {
                                        data.image[IDpeakpha_bin_ch].array.F[index1] += 0.2;
                                        data.image[IDpeakpha_bin_ch].array.F[index2] -= 0.2;
                                    }
                                    if(pv1 > pv2 + M_PI)
                                    {
                                        data.image[IDpeakpha_bin_ch].array.F[index1] -= 0.2;
                                        data.image[IDpeakpha_bin_ch].array.F[index2] += 0.2;
                                    }
                                }


                            for(uint32_t ii2 = 0; ii2 < xsizepeakpha - 1; ii2++)
                                for(uint32_t jj2 = 0; jj2 < ysizepeakpha - 1; jj2++)
                                {
                                    uint64_t index1 = jj2 * xsizepeakpha + ii2;
                                    uint64_t index2 = (jj2 + 1) * xsizepeakpha + ii2 + 1;
                                    float pv1 = data.image[IDpeakpha_bin].array.F[index1];
                                    float pv2 = data.image[IDpeakpha_bin].array.F[index2];
                                    if(pv2 > pv1 + M_PI)
                                    {
                                        data.image[IDpeakpha_bin_ch].array.F[index1] += 0.2;
                                        data.image[IDpeakpha_bin_ch].array.F[index2] -= 0.2;
                                    }
                                    if(pv1 > pv2 + M_PI)
                                    {
                                        data.image[IDpeakpha_bin_ch].array.F[index1] -= 0.2;
                                        data.image[IDpeakpha_bin_ch].array.F[index2] += 0.2;
                                    }
                                }


                            for(uint32_t ii2 = 0; ii2 < xsizepeakpha - 1; ii2++)
                                for(uint32_t jj2 = 0; jj2 < ysizepeakpha - 1; jj2++)
                                {
                                    uint64_t index1 = (jj2 + 1) * xsizepeakpha + ii2;
                                    uint64_t index2 = jj2 * xsizepeakpha + ii2 + 1;
                                    float pv1 = data.image[IDpeakpha_bin].array.F[index1];
                                    float pv2 = data.image[IDpeakpha_bin].array.F[index2];
                                    if(pv2 > pv1 + M_PI)
                                    {
                                        data.image[IDpeakpha_bin_ch].array.F[index1] += 0.2;
                                        data.image[IDpeakpha_bin_ch].array.F[index2] -= 0.2;
                                    }
                                    if(pv1 > pv2 + M_PI)
                                    {
                                        data.image[IDpeakpha_bin_ch].array.F[index1] -= 0.2;
                                        data.image[IDpeakpha_bin_ch].array.F[index2] += 0.2;
                                    }
                                }



                            for(uint32_t ii2 = 0; ii2 < xsizepeakpha - 2; ii2++)
                                for(uint32_t jj2 = 0; jj2 < ysizepeakpha; jj2++)
                                {
                                    uint64_t index1 = jj2 * xsizepeakpha + ii2;
                                    uint64_t index2 = jj2 * xsizepeakpha + ii2 + 2;
                                    float pv1 = data.image[IDpeakpha_bin].array.F[index1];
                                    float pv2 = data.image[IDpeakpha_bin].array.F[index2];
                                    if(pv2 > pv1 + M_PI)
                                    {
                                        data.image[IDpeakpha_bin_ch].array.F[index1] += pcoeff2;
                                        data.image[IDpeakpha_bin_ch].array.F[index2] -= pcoeff2;
                                    }
                                    if(pv1 > pv2 + M_PI)
                                    {
                                        data.image[IDpeakpha_bin_ch].array.F[index1] -= pcoeff2;
                                        data.image[IDpeakpha_bin_ch].array.F[index2] += pcoeff2;
                                    }
                                }


                            for(uint32_t ii2 = 0; ii2 < xsizepeakpha; ii2++)
                                for(uint32_t jj2 = 0; jj2 < ysizepeakpha - 2; jj2++)
                                {
                                    uint64_t index1 = (jj2 + 2) * xsizepeakpha + ii2;
                                    uint64_t index2 = jj2 * xsizepeakpha + ii2;
                                    float pv1 = data.image[IDpeakpha_bin].array.F[index1];
                                    float pv2 = data.image[IDpeakpha_bin].array.F[index2];
                                    if(pv2 > pv1 + M_PI)
                                    {
                                        data.image[IDpeakpha_bin_ch].array.F[index1] += pcoeff2;
                                        data.image[IDpeakpha_bin_ch].array.F[index2] -= pcoeff2;
                                    }
                                    if(pv1 > pv2 + M_PI)
                                    {
                                        data.image[IDpeakpha_bin_ch].array.F[index1] -= pcoeff2;
                                        data.image[IDpeakpha_bin_ch].array.F[index2] += pcoeff2;
                                    }
                                }

                            for(uint32_t ii2 = 0; ii2 < xsizepeakpha - 1; ii2++)
                                for(uint32_t jj2 = 0; jj2 < ysizepeakpha - 2; jj2++)
                                {
                                    uint64_t index1 = (jj2 + 2) * xsizepeakpha + ii2 + 1;
                                    uint64_t index2 = jj2 * xsizepeakpha + ii2;
                                    float pv1 = data.image[IDpeakpha_bin].array.F[index1];
                                    float pv2 = data.image[IDpeakpha_bin].array.F[index2];
                                    if(pv2 > pv1 + M_PI)
                                    {
                                        data.image[IDpeakpha_bin_ch].array.F[index1] += pcoeff2;
                                        data.image[IDpeakpha_bin_ch].array.F[index2] -= pcoeff2;
                                    }
                                    if(pv1 > pv2 + M_PI)
                                    {
                                        data.image[IDpeakpha_bin_ch].array.F[index1] -= pcoeff2;
                                        data.image[IDpeakpha_bin_ch].array.F[index2] += pcoeff2;
                                    }
                                }

                            for(uint32_t ii2 = 0; ii2 < xsizepeakpha - 1; ii2++)
                                for(uint32_t jj2 = 0; jj2 < ysizepeakpha - 2; jj2++)
                                {
                                    uint64_t index1 = (jj2 + 2) * xsizepeakpha + ii2;
                                    uint64_t index2 = jj2 * xsizepeakpha + ii2 + 1;
                                    float pv1 = data.image[IDpeakpha_bin].array.F[index1];
                                    float pv2 = data.image[IDpeakpha_bin].array.F[index2];
                                    if(pv2 > pv1 + M_PI)
                                    {
                                        data.image[IDpeakpha_bin_ch].array.F[index1] += pcoeff2;
                                        data.image[IDpeakpha_bin_ch].array.F[index2] -= pcoeff2;
                                    }
                                    if(pv1 > pv2 + M_PI)
                                    {
                                        data.image[IDpeakpha_bin_ch].array.F[index1] -= pcoeff2;
                                        data.image[IDpeakpha_bin_ch].array.F[index2] += pcoeff2;
                                    }
                                }

                            for(uint32_t ii2 = 0; ii2 < xsizepeakpha - 2; ii2++)
                                for(uint32_t jj2 = 0; jj2 < ysizepeakpha - 2; jj2++)
                                {
                                    uint64_t index1 = (jj2 + 2) * xsizepeakpha + ii2 + 2;
                                    uint64_t index2 = jj2 * xsizepeakpha + ii2;
                                    float pv1 = data.image[IDpeakpha_bin].array.F[index1];
                                    float pv2 = data.image[IDpeakpha_bin].array.F[index2];
                                    if(pv2 > pv1 + M_PI)
                                    {
                                        data.image[IDpeakpha_bin_ch].array.F[index1] += pcoeff2;
                                        data.image[IDpeakpha_bin_ch].array.F[index2] -= pcoeff2;
                                    }
                                    if(pv1 > pv2 + M_PI)
                                    {
                                        data.image[IDpeakpha_bin_ch].array.F[index1] -= pcoeff2;
                                        data.image[IDpeakpha_bin_ch].array.F[index2] += pcoeff2;
                                    }
                                }

                            for(uint32_t ii2 = 0; ii2 < xsizepeakpha - 2; ii2++)
                                for(uint32_t jj2 = 0; jj2 < ysizepeakpha - 2; jj2++)
                                {
                                    uint64_t index1 = (jj2 + 2) * xsizepeakpha + ii2;
                                    uint64_t index2 = jj2 * xsizepeakpha + ii2 + 2;
                                    float pv1 = data.image[IDpeakpha_bin].array.F[index1];
                                    float pv2 = data.image[IDpeakpha_bin].array.F[index2];
                                    if(pv2 > pv1 + M_PI)
                                    {
                                        data.image[IDpeakpha_bin_ch].array.F[index1] += pcoeff2;
                                        data.image[IDpeakpha_bin_ch].array.F[index2] -= pcoeff2;
                                    }
                                    if(pv1 > pv2 + M_PI)
                                    {
                                        data.image[IDpeakpha_bin_ch].array.F[index1] -= pcoeff2;
                                        data.image[IDpeakpha_bin_ch].array.F[index2] += pcoeff2;
                                    }
                                }

                            for(uint32_t ii2 = 0; ii2 < xsizepeakpha - 2; ii2++)
                                for(uint32_t jj2 = 0; jj2 < ysizepeakpha - 1; jj2++)
                                {
                                    uint64_t index1 = (jj2 + 1) * xsizepeakpha + ii2;
                                    uint64_t index2 = jj2 * xsizepeakpha + ii2 + 2;
                                    float pv1 = data.image[IDpeakpha_bin].array.F[index1];
                                    float pv2 = data.image[IDpeakpha_bin].array.F[index2];
                                    if(pv2 > pv1 + M_PI)
                                    {
                                        data.image[IDpeakpha_bin_ch].array.F[index1] += pcoeff2;
                                        data.image[IDpeakpha_bin_ch].array.F[index2] -= pcoeff2;
                                    }
                                    if(pv1 > pv2 + M_PI)
                                    {
                                        data.image[IDpeakpha_bin_ch].array.F[index1] -= pcoeff2;
                                        data.image[IDpeakpha_bin_ch].array.F[index2] += pcoeff2;
                                    }
                                }

                            for(uint32_t ii2 = 0; ii2 < xsizepeakpha - 2; ii2++)
                                for(uint32_t jj2 = 0; jj2 < ysizepeakpha - 1; jj2++)
                                {
                                    uint64_t index1 = (jj2 + 1) * xsizepeakpha + ii2 + 2;
                                    uint64_t index2 = jj2 * xsizepeakpha + ii2;
                                    float pv1 = data.image[IDpeakpha_bin].array.F[index1];
                                    float pv2 = data.image[IDpeakpha_bin].array.F[index2];
                                    if(pv2 > pv1 + M_PI)
                                    {
                                        data.image[IDpeakpha_bin_ch].array.F[index1] += pcoeff2;
                                        data.image[IDpeakpha_bin_ch].array.F[index2] -= pcoeff2;
                                    }
                                    if(pv1 > pv2 + M_PI)
                                    {
                                        data.image[IDpeakpha_bin_ch].array.F[index1] -= pcoeff2;
                                        data.image[IDpeakpha_bin_ch].array.F[index2] += pcoeff2;
                                    }
                                }



                            plim = 0.0;
                            for(uint64_t ii2 = 0; ii2 < xsizepeakpha * ysizepeakpha; ii2++)
                                if(fabs(data.image[IDpeakpha_bin_ch].array.F[ii2]) > plim)
                                {
                                    plim = fabs(data.image[IDpeakpha_bin_ch].array.F[ii2]);
                                }
                            plim -= 0.001;

                            if(plim < 0.5)
                            {
                                plim = 0.5;
                            }

                            //	plim = 2.39;

                            //                        save_fits("peakpha_bin", "peakpha_bin.fits");

                            for(uint32_t ii2 = 0; ii2 < xsizepeakpha; ii2++)
                                for(uint32_t jj2 = 0; jj2 < ysizepeakpha; jj2++)
                                {
                                    if(data.image[IDpeakpha_bin_ch].array.F[jj2 * xsizepeakpha + ii2] > plim)
                                    {
                                        if((ii2 > 1) && (jj2 > 1) && (ii2 < xsizepeakpha - 2)
                                                && (jj2 < ysizepeakpha - 2))
                                        {
                                            chcnt ++;
                                        }
                                        data.image[IDpeakpha_bin].array.F[jj2 * xsizepeakpha + ii2] += 2.0 * M_PI;
                                    }
                                    if(data.image[IDpeakpha_bin_ch].array.F[jj2 * xsizepeakpha + ii2] < -plim)
                                    {
                                        if((ii2 > 1) && (jj2 > 1) && (ii2 < xsizepeakpha - 2)
                                                && (jj2 < ysizepeakpha - 2))
                                        {
                                            chcnt ++;
                                        }
                                        data.image[IDpeakpha_bin].array.F[jj2 * xsizepeakpha + ii2] -= 2.0 * M_PI;
                                    }

                                }
                        }
                        else
                        {
                            for(uint32_t ii2 = 0; ii2 < xsizepeakpha - 1; ii2++)
                                for(uint32_t jj2 = 0; jj2 < ysizepeakpha; jj2++)
                                {
                                    uint64_t index1 = jj2 * xsizepeakpha + ii2;
                                    uint64_t index2 = jj2 * xsizepeakpha + ii2 + 1;
                                    double pv1 = data.image[IDpeakpha_bin].array.D[index1];
                                    double pv2 = data.image[IDpeakpha_bin].array.D[index2];
                                    if(pv2 > pv1 + M_PI)
                                    {
                                        data.image[IDpeakpha_bin_ch].array.D[index1] += 0.2;
                                        data.image[IDpeakpha_bin_ch].array.D[index2] -= 0.2;
                                    }
                                    if(pv1 > pv2 + M_PI)
                                    {
                                        data.image[IDpeakpha_bin_ch].array.D[index1] -= 0.2;
                                        data.image[IDpeakpha_bin_ch].array.D[index2] += 0.2;
                                    }
                                }

                            for(uint32_t ii2 = 0; ii2 < xsizepeakpha; ii2++)
                                for(uint32_t jj2 = 0; jj2 < ysizepeakpha - 1; jj2++)
                                {
                                    uint64_t index1 = jj2 * xsizepeakpha + ii2;
                                    uint64_t index2 = (jj2 + 1) * xsizepeakpha + ii2;
                                    double pv1 = data.image[IDpeakpha_bin].array.D[index1];
                                    double pv2 = data.image[IDpeakpha_bin].array.D[index2];
                                    if(pv2 > pv1 + M_PI)
                                    {
                                        data.image[IDpeakpha_bin_ch].array.D[index1] += 0.2;
                                        data.image[IDpeakpha_bin_ch].array.D[index2] -= 0.2;
                                    }
                                    if(pv1 > pv2 + M_PI)
                                    {
                                        data.image[IDpeakpha_bin_ch].array.D[index1] -= 0.2;
                                        data.image[IDpeakpha_bin_ch].array.D[index2] += 0.2;
                                    }
                                }


                            for(uint32_t ii2 = 0; ii2 < xsizepeakpha - 1; ii2++)
                                for(uint32_t jj2 = 0; jj2 < ysizepeakpha - 1; jj2++)
                                {
                                    uint64_t index1 = jj2 * xsizepeakpha + ii2;
                                    uint64_t index2 = (jj2 + 1) * xsizepeakpha + ii2 + 1;
                                    double pv1 = data.image[IDpeakpha_bin].array.D[index1];
                                    double pv2 = data.image[IDpeakpha_bin].array.D[index2];
                                    if(pv2 > pv1 + M_PI)
                                    {
                                        data.image[IDpeakpha_bin_ch].array.D[index1] += 0.2;
                                        data.image[IDpeakpha_bin_ch].array.D[index2] -= 0.2;
                                    }
                                    if(pv1 > pv2 + M_PI)
                                    {
                                        data.image[IDpeakpha_bin_ch].array.D[index1] -= 0.2;
                                        data.image[IDpeakpha_bin_ch].array.D[index2] += 0.2;
                                    }
                                }


                            for(uint32_t ii2 = 0; ii2 < xsizepeakpha - 1; ii2++)
                                for(uint32_t jj2 = 0; jj2 < ysizepeakpha - 1; jj2++)
                                {
                                    uint64_t index1 = (jj2 + 1) * xsizepeakpha + ii2;
                                    uint64_t index2 = jj2 * xsizepeakpha + ii2 + 1;
                                    double pv1 = data.image[IDpeakpha_bin].array.D[index1];
                                    double pv2 = data.image[IDpeakpha_bin].array.D[index2];
                                    if(pv2 > pv1 + M_PI)
                                    {
                                        data.image[IDpeakpha_bin_ch].array.D[index1] += 0.2;
                                        data.image[IDpeakpha_bin_ch].array.D[index2] -= 0.2;
                                    }
                                    if(pv1 > pv2 + M_PI)
                                    {
                                        data.image[IDpeakpha_bin_ch].array.D[index1] -= 0.2;
                                        data.image[IDpeakpha_bin_ch].array.D[index2] += 0.2;
                                    }
                                }


                            for(uint32_t ii2 = 0; ii2 < xsizepeakpha - 2; ii2++)
                                for(uint32_t jj2 = 0; jj2 < ysizepeakpha; jj2++)
                                {
                                    uint64_t index1 = jj2 * xsizepeakpha + ii2;
                                    uint64_t index2 = jj2 * xsizepeakpha + ii2 + 2;
                                    double pv1 = data.image[IDpeakpha_bin].array.D[index1];
                                    double pv2 = data.image[IDpeakpha_bin].array.D[index2];
                                    if(pv2 > pv1 + M_PI)
                                    {
                                        data.image[IDpeakpha_bin_ch].array.D[index1] += pcoeff2;
                                        data.image[IDpeakpha_bin_ch].array.D[index2] -= pcoeff2;
                                    }
                                    if(pv1 > pv2 + M_PI)
                                    {
                                        data.image[IDpeakpha_bin_ch].array.D[index1] -= pcoeff2;
                                        data.image[IDpeakpha_bin_ch].array.D[index2] += pcoeff2;
                                    }
                                }


                            for(uint32_t ii2 = 0; ii2 < xsizepeakpha; ii2++)
                                for(uint32_t jj2 = 0; jj2 < ysizepeakpha - 2; jj2++)
                                {
                                    uint64_t index1 = (jj2 + 2) * xsizepeakpha + ii2;
                                    uint64_t index2 = jj2 * xsizepeakpha + ii2;
                                    double pv1 = data.image[IDpeakpha_bin].array.D[index1];
                                    double pv2 = data.image[IDpeakpha_bin].array.D[index2];
                                    if(pv2 > pv1 + M_PI)
                                    {
                                        data.image[IDpeakpha_bin_ch].array.D[index1] += pcoeff2;
                                        data.image[IDpeakpha_bin_ch].array.D[index2] -= pcoeff2;
                                    }
                                    if(pv1 > pv2 + M_PI)
                                    {
                                        data.image[IDpeakpha_bin_ch].array.D[index1] -= pcoeff2;
                                        data.image[IDpeakpha_bin_ch].array.D[index2] += pcoeff2;
                                    }
                                }

                            for(uint32_t ii2 = 0; ii2 < xsizepeakpha - 1; ii2++)
                                for(uint32_t jj2 = 0; jj2 < ysizepeakpha - 2; jj2++)
                                {
                                    uint64_t index1 = (jj2 + 2) * xsizepeakpha + ii2 + 1;
                                    uint64_t index2 = jj2 * xsizepeakpha + ii2;
                                    double pv1 = data.image[IDpeakpha_bin].array.D[index1];
                                    double pv2 = data.image[IDpeakpha_bin].array.D[index2];
                                    if(pv2 > pv1 + M_PI)
                                    {
                                        data.image[IDpeakpha_bin_ch].array.D[index1] += pcoeff2;
                                        data.image[IDpeakpha_bin_ch].array.D[index2] -= pcoeff2;
                                    }
                                    if(pv1 > pv2 + M_PI)
                                    {
                                        data.image[IDpeakpha_bin_ch].array.D[index1] -= pcoeff2;
                                        data.image[IDpeakpha_bin_ch].array.D[index2] += pcoeff2;
                                    }
                                }

                            for(uint32_t ii2 = 0; ii2 < xsizepeakpha - 1; ii2++)
                                for(uint32_t jj2 = 0; jj2 < ysizepeakpha - 2; jj2++)
                                {
                                    uint64_t index1 = (jj2 + 2) * xsizepeakpha + ii2;
                                    uint64_t index2 = jj2 * xsizepeakpha + ii2 + 1;
                                    double pv1 = data.image[IDpeakpha_bin].array.D[index1];
                                    double pv2 = data.image[IDpeakpha_bin].array.D[index2];
                                    if(pv2 > pv1 + M_PI)
                                    {
                                        data.image[IDpeakpha_bin_ch].array.D[index1] += pcoeff2;
                                        data.image[IDpeakpha_bin_ch].array.D[index2] -= pcoeff2;
                                    }
                                    if(pv1 > pv2 + M_PI)
                                    {
                                        data.image[IDpeakpha_bin_ch].array.D[index1] -= pcoeff2;
                                        data.image[IDpeakpha_bin_ch].array.D[index2] += pcoeff2;
                                    }
                                }

                            for(uint32_t ii2 = 0; ii2 < xsizepeakpha - 2; ii2++)
                                for(uint32_t jj2 = 0; jj2 < ysizepeakpha - 2; jj2++)
                                {
                                    uint64_t index1 = (jj2 + 2) * xsizepeakpha + ii2 + 2;
                                    uint64_t index2 = jj2 * xsizepeakpha + ii2;
                                    double pv1 = data.image[IDpeakpha_bin].array.D[index1];
                                    double pv2 = data.image[IDpeakpha_bin].array.D[index2];
                                    if(pv2 > pv1 + M_PI)
                                    {
                                        data.image[IDpeakpha_bin_ch].array.D[index1] += pcoeff2;
                                        data.image[IDpeakpha_bin_ch].array.D[index2] -= pcoeff2;
                                    }
                                    if(pv1 > pv2 + M_PI)
                                    {
                                        data.image[IDpeakpha_bin_ch].array.D[index1] -= pcoeff2;
                                        data.image[IDpeakpha_bin_ch].array.D[index2] += pcoeff2;
                                    }
                                }

                            for(uint32_t ii2 = 0; ii2 < xsizepeakpha - 2; ii2++)
                                for(uint32_t jj2 = 0; jj2 < ysizepeakpha - 2; jj2++)
                                {
                                    uint64_t index1 = (jj2 + 2) * xsizepeakpha + ii2;
                                    uint64_t index2 = jj2 * xsizepeakpha + ii2 + 2;
                                    double pv1 = data.image[IDpeakpha_bin].array.D[index1];
                                    double pv2 = data.image[IDpeakpha_bin].array.D[index2];
                                    if(pv2 > pv1 + M_PI)
                                    {
                                        data.image[IDpeakpha_bin_ch].array.D[index1] += pcoeff2;
                                        data.image[IDpeakpha_bin_ch].array.D[index2] -= pcoeff2;
                                    }
                                    if(pv1 > pv2 + M_PI)
                                    {
                                        data.image[IDpeakpha_bin_ch].array.D[index1] -= pcoeff2;
                                        data.image[IDpeakpha_bin_ch].array.D[index2] += pcoeff2;
                                    }
                                }

                            for(uint32_t ii2 = 0; ii2 < xsizepeakpha - 2; ii2++)
                                for(uint32_t jj2 = 0; jj2 < ysizepeakpha - 1; jj2++)
                                {
                                    uint64_t index1 = (jj2 + 1) * xsizepeakpha + ii2;
                                    uint64_t index2 = jj2 * xsizepeakpha + ii2 + 2;
                                    double pv1 = data.image[IDpeakpha_bin].array.D[index1];
                                    double pv2 = data.image[IDpeakpha_bin].array.D[index2];
                                    if(pv2 > pv1 + M_PI)
                                    {
                                        data.image[IDpeakpha_bin_ch].array.D[index1] += pcoeff2;
                                        data.image[IDpeakpha_bin_ch].array.D[index2] -= pcoeff2;
                                    }
                                    if(pv1 > pv2 + M_PI)
                                    {
                                        data.image[IDpeakpha_bin_ch].array.D[index1] -= pcoeff2;
                                        data.image[IDpeakpha_bin_ch].array.D[index2] += pcoeff2;
                                    }
                                }

                            for(uint32_t ii2 = 0; ii2 < xsizepeakpha - 2; ii2++)
                                for(uint32_t jj2 = 0; jj2 < ysizepeakpha - 1; jj2++)
                                {
                                    uint64_t index1 = (jj2 + 1) * xsizepeakpha + ii2 + 2;
                                    uint64_t index2 = jj2 * xsizepeakpha + ii2;
                                    double pv1 = data.image[IDpeakpha_bin].array.D[index1];
                                    double pv2 = data.image[IDpeakpha_bin].array.D[index2];
                                    if(pv2 > pv1 + M_PI)
                                    {
                                        data.image[IDpeakpha_bin_ch].array.D[index1] += pcoeff2;
                                        data.image[IDpeakpha_bin_ch].array.D[index2] -= pcoeff2;
                                    }
                                    if(pv1 > pv2 + M_PI)
                                    {
                                        data.image[IDpeakpha_bin_ch].array.D[index1] -= pcoeff2;
                                        data.image[IDpeakpha_bin_ch].array.D[index2] += pcoeff2;
                                    }
                                }



                            plim = 0.0;
                            for(uint64_t ii2 = 0; ii2 < xsizepeakpha * ysizepeakpha; ii2++)
                                if(fabs(data.image[IDpeakpha_bin_ch].array.D[ii2]) > plim)
                                {
                                    plim = fabs(data.image[IDpeakpha_bin_ch].array.D[ii2]);
                                }
                            plim -= 0.001;

                            if(plim < 0.5)
                            {
                                plim = 0.5;
                            }

                            //	plim = 2.39;

                            //                        save_fits("peakpha_bin", "peakpha_bin.fits");

                            for(uint32_t ii2 = 0; ii2 < xsizepeakpha; ii2++)
                                for(uint32_t jj2 = 0; jj2 < ysizepeakpha; jj2++)
                                {
                                    if(data.image[IDpeakpha_bin_ch].array.D[jj2 * xsizepeakpha + ii2] > plim)
                                    {
                                        if((ii2 > 1) && (jj2 > 1) && (ii2 < xsizepeakpha - 2)
                                                && (jj2 < ysizepeakpha - 2))
                                        {
                                            chcnt ++;
                                        }
                                        data.image[IDpeakpha_bin].array.D[jj2 * xsizepeakpha + ii2] += 2.0 * M_PI;
                                    }
                                    if(data.image[IDpeakpha_bin_ch].array.D[jj2 * xsizepeakpha + ii2] < -plim)
                                    {
                                        if((ii2 > 1) && (jj2 > 1) && (ii2 < xsizepeakpha - 2)
                                                && (jj2 < ysizepeakpha - 2))
                                        {
                                            chcnt ++;
                                        }
                                        data.image[IDpeakpha_bin].array.D[jj2 * xsizepeakpha + ii2] -= 2.0 * M_PI;
                                    }

                                }


                        }
                        //                    printf("chiter = %ld    [%ld]  %f\n", chiter, chcnt, plim);
                        //                      save_fits("peakpha_bin_ch", "peakpha_bin_ch.fits");

                        if(chcnt == 0)
                        {
                            chcnt0cnt++;
                        }
                        else
                        {
                            chcnt0cnt = 0;
                        }
                    }


                    /*          list_image_ID();
                              save_fits("peakpha_bin", "peakpha2_bin.fits");
                              save_fits("peakphare_bin", "peakpha_re_bin.fits");
                              save_fits("peakphaim_bin", "peakpha_im_bin.fits");
                    */

                    if(WFprecision == 0)
                    {
                        for(uint32_t ii = 0; ii < naxesout[0]; ii++)
                            for(uint32_t jj = 0; jj < naxesout[1]; jj++)
                            {
                                double peakpha = 0.0;
                                long ii2 = (long)(1.0 * ii / naxesout[0] * xsizepeakpha);
                                long jj2 = (long)(1.0 * jj / naxesout[1] * ysizepeakpha);

                                if((ii2 < xsizepeakpha) && (jj2 < ysizepeakpha))
                                {
                                    peakpha = data.image[IDpeakpha_bin].array.F[jj2 * xsizepeakpha + ii2];
                                }

                                long ii1 = ii + (naxes[0] - naxesout[0]) / 2;
                                long jj1 = jj + (naxes[1] - naxesout[1]) / 2;

                                float pha = data.image[IDout_sarray_pha].array.F[frame * naxesout[0] *
                                            naxesout[1] +
                                            jj * naxesout[0] + ii];
                                data.image[IDout_sarray_pha].array.F[frame * naxesout[0]*naxesout[1] + jj *
                                                                     naxesout[0] + ii] = pha + 2.0 * M_PI * ((long)((
                                                                             data.image[ID_sarray1].array.F[jj1 * naxes[0] + ii1] - pha - peakpha) / 2.0 /
                                                                             M_PI + 1000.5) - 1000.0);
                            }
                    }
                    else
                    {
                        for(uint32_t ii = 0; ii < naxesout[0]; ii++)
                            for(uint32_t jj = 0; jj < naxesout[1]; jj++)
                            {
                                double peakpha = 0.0;
                                long ii2 = (long)(1.0 * ii / naxesout[0] * xsizepeakpha);
                                long jj2 = (long)(1.0 * jj / naxesout[1] * ysizepeakpha);

                                if((ii2 < xsizepeakpha) && (jj2 < ysizepeakpha))
                                {
                                    peakpha = data.image[IDpeakpha_bin].array.D[jj2 * xsizepeakpha + ii2];
                                }

                                long ii1 = ii + (naxes[0] - naxesout[0]) / 2;
                                long jj1 = jj + (naxes[1] - naxesout[1]) / 2;

                                double pha = data.image[IDout_sarray_pha].array.D[frame * naxesout[0] *
                                             naxesout[1] +
                                             jj * naxesout[0] + ii];
                                data.image[IDout_sarray_pha].array.D[frame * naxesout[0]*naxesout[1] + jj *
                                                                     naxesout[0] + ii] = pha + 2.0 * M_PI * ((long)((
                                                                             data.image[ID_sarray1].array.D[jj1 * naxes[0] + ii1] - pha - peakpha) / 2.0 /
                                                                             M_PI + 1000.5) - 1000.0);
                            }
                    }
                }
            }



            /*
                        if(make_cwavefront==1)
                        {
                            for(ii=0; ii<naxesout[0]; ii++)
                                for(jj=0; jj<naxesout[1]; jj++)
                                {
                                    ii1 = ii+(naxes[0]-naxesout[0])/2;
                                    jj1 = jj+(naxes[1]-naxesout[1])/2;
                                    carray[frame*naxesout[0]*naxesout[1]+jj*naxesout[0]+ii].re = data.image[ID_carray2].array.CF[jj1*naxes[0]+ii1].re;
                                    carray[frame*naxesout[0]*naxesout[1]+jj*naxesout[0]+ii].im = data.image[ID_carray2].array.CF[jj1*naxes[0]+ii1].im;
                                }
                        }

            */


            // WRITE CURRENT WF TO SHARED MEMORY
            if((atmturbconf.flag_SHMoutput == 1) && (atmturbconf.flag_SWF_make == 1))
            {
                double coeff = 1.0;

                switch(atmturbconf.flag_SWF_SHMouputM)
                {
                case 1 :
                    coeff = SLAMBDA / 2.0 / M_PI;
                    break;
                case 2 :
                    coeff = SLAMBDA / 2.0 / M_PI * 1e6;
                    break;
                default :
                    coeff = 1.0;
                    break;
                }



                if(atmturbconf.flag_WFampl == 0)
                {
                    data.image[IDshmspha].md[0].write = 1;
                    data.image[IDshmspha].kw[0].value.numf = tnowdouble;
                    if(WFprecision == 0)
                    {
                        for(uint64_t ii = 0; ii < naxesout[0]*naxesout[1]; ii++)
                        {
                            data.image[IDshmspha].array.F[ii] = data.image[ID_sarray1].array.F[frame *
                                                                naxesout[0] * naxesout[1] + ii] * coeff;
                        }
                    }
                    else
                    {
                        for(uint64_t ii = 0; ii < naxesout[0]*naxesout[1]; ii++)
                        {
                            data.image[IDshmspha].array.D[ii] = data.image[ID_sarray1].array.D[frame *
                                                                naxesout[0] * naxesout[1] + ii] * coeff;
                        }
                    }
                    data.image[IDshmspha].md[0].cnt0++;
                    data.image[IDshmspha].md[0].write = 0;
                }
                else
                {
                    data.image[IDshmspha].md[0].write = 1;
                    data.image[IDshmsamp].md[0].write = 1;
                    data.image[IDshmspha].kw[0].value.numf = tnowdouble;
                    data.image[IDshmsamp].kw[0].value.numf = tnowdouble;
                    if(WFprecision == 0)
                    {
                        for(uint64_t ii = 0; ii < naxesout[0]*naxesout[1]; ii++)
                        {
                            data.image[IDshmspha].array.F[ii] = data.image[IDout_sarray_pha].array.F[frame *
                                                                naxesout[0] * naxesout[1] + ii] * coeff;
                            data.image[IDshmsamp].array.F[ii] = data.image[IDout_sarray_amp].array.F[frame *
                                                                naxesout[0] * naxesout[1] + ii];
                        }
                    }
                    else
                    {
                        for(uint64_t ii = 0; ii < naxesout[0]*naxesout[1]; ii++)
                        {
                            data.image[IDshmspha].array.D[ii] = data.image[IDout_sarray_pha].array.D[frame *
                                                                naxesout[0] * naxesout[1] + ii] * coeff;
                            data.image[IDshmsamp].array.D[ii] = data.image[IDout_sarray_amp].array.D[frame *
                                                                naxesout[0] * naxesout[1] + ii];
                        }
                    }
                    data.image[IDshmspha].md[0].cnt0++;
                    data.image[IDshmsamp].md[0].cnt0++;
                    data.image[IDshmspha].md[0].write = 0;
                    data.image[IDshmsamp].md[0].write = 0;
                }

            }
        }


        if(atmturbconf.flag_WFoutput == 1) // WRITE REFERENCE LAMBDA
        {
            char fnameoutpha[STRINGMAXLEN_FILENAME];
            WRITE_FILENAME(fnameoutpha, "%s%08ld.%09ld.pha.fits", atmturbconf.WFfileprefix,
                           cubeindex, (long)(1.0e12 * SLAMBDA + 0.5));

            if(WFprecision == 0)
            {
                save_fl_fits("outarraypha", fnameoutpha);
            }
            else
            {
                save_db_fits("outarraypha", fnameoutpha);
            }


            if(atmturbconf.flag_WFampl == 1)
            {
                char fnameoutamp[STRINGMAXLEN_FILENAME];
                WRITE_FILENAME(fnameoutamp, "%s%08ld.%09ld.amp.fits", atmturbconf.WFfileprefix,
                               cubeindex, (long)(1.0e12 * SLAMBDA + 0.5));

                if(WFprecision == 0)
                {
                    save_fl_fits("outarrayamp", fnameoutamp);
                }
                else
                {
                    save_db_fits("outarrayamp", fnameoutamp);
                }
            }
        }
        else if(atmturbconf.flag_WFoutput == 0)
        {
            // CREATE EMPTY FILES
            char fnameoutpha[STRINGMAXLEN_FILENAME];
            WRITE_FILENAME(fnameoutpha, "%s%08ld.%09ld.pha.fits", atmturbconf.WFfileprefix,
                           cubeindex, (long)(1.0e12 * SLAMBDA + 0.5));
            EXECUTE_SYSTEM_COMMAND("touch %s", fnameoutpha);

            if(atmturbconf.flag_WFampl == 1)
            {
                char fnameoutamp[STRINGMAXLEN_FILENAME];
                WRITE_FILENAME(fnameoutamp, "%s%08ld.%09ld.amp.fits", atmturbconf.WFfileprefix,
                               cubeindex, (long)(1.0e12 * SLAMBDA + 0.5));
                EXECUTE_SYSTEM_COMMAND("touch %s", fnameoutamp);
            }
        }


        if((atmturbconf.flag_SWF_make == 1)
                && (atmturbconf.flag_SWF_Wite2Disk == 1)) // WRITE SCIENCE LAMBDA
        {
            printf("WRITING WAVEFRONT FILE ...");
            fflush(stdout);

            char fnameoutpha[STRINGMAXLEN_FILENAME];
            WRITE_FILENAME(fnameoutpha, "%s%08ld.%09ld.pha.fits", atmturbconf.SWFfileprefix,
                           cubeindex, (long)(1.0e12 * SLAMBDA + 0.5));

            if(WFprecision == 0)
            {
                save_fl_fits("outsarraypha", fnameoutpha);
            }
            else
            {
                save_db_fits("outsarraypha", fnameoutpha);
            }

            printf(" - ");
            fflush(stdout);

            if(atmturbconf.flag_WFampl == 1)
            {
                char fnameoutamp[STRINGMAXLEN_FILENAME];
                WRITE_FILENAME(fnameoutamp, "%s%08ld.%09ld.amp.fits", atmturbconf.SWFfileprefix,
                               cubeindex, (long)(1.0e12 * SLAMBDA + 0.5));

                if(WFprecision == 0)
                {
                    save_fl_fits("outsarrayamp", fnameoutamp);
                }
                else
                {
                    save_db_fits("outsarrayamp", fnameoutamp);
                }
            }

            printf("\n");
            fflush(stdout);
        }

    }

    delete_image_ID("array1", DELETE_IMAGE_ERRMODE_WARNING);
    if(WFprecision == 0)
    {
        free(array);
    }
    else
    {
        free(array_double);
    }

    if(atmturbconf.flag_SWF_make == 1)
    {
        delete_image_ID("sarray1", DELETE_IMAGE_ERRMODE_WARNING);
        if(WFprecision == 0)
        {
            free(sarray);
        }
        else
        {
            free(sarray_double);
        }
    }
    /*  if(make_cwavefront==1)
      {
          delete_image_ID("carray1");
          free(carray);
      }*/

    free(SLAYER_ALT);
    free(super_layer_index);
    free(xpos);
    free(ypos);
    free(xpos0);
    free(ypos0);
    free(xposfcnt);
    free(yposfcnt);
    free(vxpix);
    free(vypix);


    free(LAYER_ALT);
    free(LAYER_CN2);
    free(LAYER_SPD);
    free(LAYER_DIR);
    free(LAYER_OUTERSCALE);
    free(LAYER_INNERSCALE);
    free(LAYER_SIGMAWSPEED);
    free(LAYER_LWIND);
    free(ID_TM);
    free(ID_TML);

    free(naxes);
    free(naxesout);

    return RETURN_SUCCESS;
}







