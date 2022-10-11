/**
 * @file    ReadConf.c
 * @brief   Read configuration file from disk
 *
 *
 */

#include <stdint.h>
#include <string.h>

#include "CommandLineInterface/CLIcore.h"

#include "COREMOD_iofits/COREMOD_iofits.h" // loadfits

#include "COREMOD_tools/fileutils.h" // read_config_parameter

#include "AtmosphericTurbulence_conf.h"

errno_t AtmosphericTurbulence_ReadConf(const char *restrict fnameconf,
                                       ATMTURBCONF *atmturbconf)
{
    char KEYWORD[200];
    char CONTENT[200];

    // ------------ TURBULENCE AND ATMOSPHERE PARAMETERS -------------------------------------------

    atmturbconf->lambda =
        read_config_parameter_float(fnameconf, "TURBULENCE_REF_WAVEL") *
        0.000001;
    atmturbconf->seeing =
        read_config_parameter_float(fnameconf, "TURBULENCE_SEEING");
    read_config_parameter(fnameconf,
                          "TURBULENCE_PROF_FILE",
                          atmturbconf->turbulenceprof_fname);
    atmturbconf->zenithangle =
        read_config_parameter_float(fnameconf, "ZENITH_ANGLE");
    atmturbconf->sourceXpos =
        read_config_parameter_float(fnameconf, "SOURCE_XPOS");
    atmturbconf->sourceYpos =
        read_config_parameter_float(fnameconf, "SOURCE_YPOS");

    // ------------ LOW_FIDELITY OUTPUT AT REF LAMBDA ------------------------------------------------

    if(read_config_parameter_exists(fnameconf, "WFOUTPUT") == 1)
    {
        atmturbconf->flag_WFoutput =
            read_config_parameter_int(fnameconf, "WFOUTPUT");
    }
    read_config_parameter(fnameconf,
                          "WF_FILE_PREFIX",
                          atmturbconf->WFfileprefix);
    if(read_config_parameter_exists(fnameconf, "SHM_OUTPUT") == 1)
    {
        atmturbconf->flag_SHMoutput =
            read_config_parameter_int(fnameconf, "SHM_OUTPUT");
    }

    // ------------- HIGH FIDELITY OUTPUT WAVELENGTH ---------------------------------

    atmturbconf->flag_SWF_make =
        read_config_parameter_int(fnameconf, "MAKE_SWAVEFRONT");
    atmturbconf->flag_SWF_Wite2Disk =
        read_config_parameter_int(fnameconf, "SWF_WRITE2DISK");
    read_config_parameter(fnameconf,
                          "SWF_FILE_PREFIX",
                          atmturbconf->SWFfileprefix);
    if(read_config_parameter_exists(fnameconf, "SHM_SOUTPUT") == 1)
    {
        atmturbconf->flag_SWF_SHMoutput =
            read_config_parameter_int(fnameconf, "SHM_SOUTPUT");
    }
    read_config_parameter(fnameconf, "SHM_SPREFIX", atmturbconf->SWFSHMprefix);
    atmturbconf->flag_SWF_SHMouputM =
        read_config_parameter_int(fnameconf, "SHM_SOUTPUTM");

    // ------------ OUTPUT PARAMETERS --------------------------------------------

    atmturbconf->WFsize = read_config_parameter_int(fnameconf, "WFsize");
    atmturbconf->PupilScale =
        read_config_parameter_float(fnameconf, "PUPIL_SCALE");

    // ------------- TIMING ---------------------------------------

    atmturbconf->flag_RealTime =
        read_config_parameter_int(fnameconf, "REALTIME");
    atmturbconf->RealTimeFactor =
        read_config_parameter_float(fnameconf, "REALTIMEFACTOR");
    atmturbconf->TimeStep =
        read_config_parameter_float(fnameconf, "WFTIME_STEP");
    atmturbconf->TimeSpanCube =
        read_config_parameter_float(fnameconf, "TIME_CUBESPAN");
    atmturbconf->NBTimeCube =
        read_config_parameter_float(fnameconf, "TIME_NBCUBE");
    atmturbconf->TimeDelayus =
        read_config_parameter_int(fnameconf, "SIMTDELAY");
    atmturbconf->flag_WaitSem =
        read_config_parameter_int(fnameconf, "WAITFORSEM");
    read_config_parameter(fnameconf, "WAITSEMIMNAME", atmturbconf->WaitSemName);

    // ------------ COMPUTATION PARAMETERS, MODES --------------------------------

    if(read_config_parameter_exists(fnameconf, "SKIP_EXISTING") == 1)
    {
        atmturbconf->flag_SkipExisting =
            read_config_parameter_int(fnameconf, "SKIP_EXISTING");
    }
    else
    {
        atmturbconf->flag_SkipExisting = 0;
    }

    atmturbconf->WFrawSize =
        read_config_parameter_int(fnameconf, "WF_RAW_SIZE");
    atmturbconf->MasterSize =
        read_config_parameter_int(fnameconf, "MASTER_SIZE");

    //  ------------ WAVEFRONT AMPLITUDE -------------------------------------------

    atmturbconf->flag_WFampl =
        read_config_parameter_int(fnameconf, "WAVEFRONT_AMPLITUDE");

    if(atmturbconf->flag_WFampl == 1)
    {
        atmturbconf->flag_FresnelProp =
            read_config_parameter_int(fnameconf, "FRESNEL_PROPAGATION");
    }
    else
    {
        atmturbconf->flag_FresnelProp = 0;
    }
    atmturbconf->FresnelPropBin =
        read_config_parameter_float(fnameconf, "FRESNEL_PROPAGATION_BIN");

    // ------------ POSTPROCESSING --------------

    strcpy(KEYWORD, "PUPIL_AMPL_FILE");
    if(read_config_parameter_exists(fnameconf, "PUPIL_AMPL_FILE") == 1)
    {
        read_config_parameter(fnameconf, "PUPIL_AMPL_FILE", CONTENT);
        load_fits(CONTENT, "ST_pa", LOADFITS_ERRMODE_WARNING, NULL);
    }

    strcpy(KEYWORD, "PUPIL_PHA_FILE");
    if(read_config_parameter_exists(fnameconf, "PUPIL_PHA_FILE") == 1)
    {
        read_config_parameter(fnameconf, "PUPIL_PHA_FILE", CONTENT);
        load_fits(CONTENT, "ST_pp", LOADFITS_ERRMODE_WARNING, NULL);
    }

    return RETURN_SUCCESS;
}
