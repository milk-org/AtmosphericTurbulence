/**
 * @file    ReadConf.c
 * @brief   Read configuration file from disk
 *
 *
 */


#include "CommandLineInterface/CLIcore.h"


errno_t AtmosphericTurbulence_ReadConf()
{
    char KEYWORD[200];
    char CONTENT[200];
    
    
   // ------------ TURBULENCE AND ATMOSPHERE PARAMETERS -------------------------------------------
    
    strcpy(KEYWORD,"TURBULENCE_REF_WAVEL");
    read_config_parameter(CONFFILE, KEYWORD, CONTENT);
    CONF_LAMBDA = atof(CONTENT)*0.000001;
    
    strcpy(KEYWORD,"TURBULENCE_SEEING");
    read_config_parameter(CONFFILE,KEYWORD,CONTENT);
    CONF_SEEING = atof(CONTENT);

    strcpy(KEYWORD,"TURBULENCE_PROF_FILE");
    read_config_parameter(CONFFILE, KEYWORD, CONF_TURBULENCE_PROF_FILE);
 
    strcpy(KEYWORD,"ZENITH_ANGLE");
    read_config_parameter(CONFFILE,KEYWORD,CONTENT);
    CONF_ZANGLE = atof(CONTENT);

    strcpy(KEYWORD,"SOURCE_XPOS");
    read_config_parameter(CONFFILE,KEYWORD,CONTENT);
    CONF_SOURCE_Xpos = atof(CONTENT);

    strcpy(KEYWORD,"SOURCE_YPOS");
    read_config_parameter(CONFFILE,KEYWORD,CONTENT);
    CONF_SOURCE_Ypos = atof(CONTENT);

 
    
    // ------------ LOW_FIDELITY OUTPUT AT REF LAMBDA ------------------------------------------------
    
    strcpy(KEYWORD,"WFOUTPUT");
    if(read_config_parameter_exists(CONFFILE, KEYWORD)==1)
    {
        read_config_parameter(CONFFILE,KEYWORD, CONTENT);
        CONF_WFOUTPUT = atoi(CONTENT);
    }
    
    strcpy(KEYWORD,"WF_FILE_PREFIX");
    read_config_parameter(CONFFILE,KEYWORD, CONF_WF_FILE_PREFIX);

    strcpy(KEYWORD,"SHM_OUTPUT");
    if(read_config_parameter_exists(CONFFILE, KEYWORD)==1)
    {
        read_config_parameter(CONFFILE, KEYWORD, CONTENT);
        CONF_SHM_OUTPUT = atoi(CONTENT);
    }

    
    // ------------- HIGH FIDELITY OUTPUT WAVELENGTH ---------------------------------
    
    strcpy(KEYWORD,"MAKE_SWAVEFRONT");
    read_config_parameter(CONFFILE, KEYWORD, CONTENT);
    CONF_MAKE_SWAVEFRONT = atoi(CONTENT);

//    strcpy(KEYWORD,"SLAMBDA");
 //   read_config_parameter(CONFFILE, KEYWORD, CONTENT);
 //   CONF_SLAMBDA = atof(CONTENT)*0.000001;

    strcpy(KEYWORD,"SWF_WRITE2DISK");
    read_config_parameter(CONFFILE, KEYWORD, CONTENT);
    CONF_SWF_WRITE2DISK = atoi(CONTENT);

    strcpy(KEYWORD,"SWF_FILE_PREFIX");
    read_config_parameter(CONFFILE, KEYWORD, CONF_SWF_FILE_PREFIX);

    strcpy(KEYWORD,"SHM_SOUTPUT");
    if(read_config_parameter_exists(CONFFILE, KEYWORD)==1)
    {
        read_config_parameter(CONFFILE, KEYWORD, CONTENT);
        CONF_SHM_SOUTPUT = atoi(CONTENT);
    }
    strcpy(KEYWORD,"SHM_SPREFIX");
    read_config_parameter(CONFFILE,KEYWORD, CONF_SHM_SPREFIX);

    strcpy(KEYWORD,"SHM_SOUTPUTM");
    read_config_parameter(CONFFILE, KEYWORD, CONTENT);
    CONF_SHM_SOUTPUTM = atoi(CONTENT);


    // ------------ OUTPUT PARAMETERS --------------------------------------------
    
    strcpy(KEYWORD,"PUPIL_SCALE");
    read_config_parameter(CONFFILE, KEYWORD, CONTENT);
    CONF_PUPIL_SCALE = atof(CONTENT);
    
    strcpy(KEYWORD,"WFsize");
    read_config_parameter(CONFFILE, KEYWORD, CONTENT);
    CONF_WFsize = atol(CONTENT);

 
    
    // ------------- TIMING ---------------------------------------
    
    strcpy(KEYWORD,"REALTIME");
    read_config_parameter(CONFFILE, KEYWORD, CONTENT);
    CONF_ATMWF_REALTIME = atoi(CONTENT);
    
    strcpy(KEYWORD,"REALTIMEFACTOR");
    read_config_parameter(CONFFILE, KEYWORD, CONTENT);
    CONF_ATMWF_REALTIMEFACTOR = atof(CONTENT);

    strcpy(KEYWORD, "WFTIME_STEP");
    read_config_parameter(CONFFILE, KEYWORD, CONTENT);
    CONF_WFTIME_STEP = atof(CONTENT);

    strcpy(KEYWORD, "TIME_SPAN");
    read_config_parameter(CONFFILE, KEYWORD, CONTENT);
    CONF_TIME_SPAN = atof(CONTENT);

    strcpy(KEYWORD, "NB_TSPAN");
    read_config_parameter(CONFFILE, KEYWORD, CONTENT);
    CONF_NB_TSPAN = atol(CONTENT);

    strcpy(KEYWORD, "SIMTDELAY");
    read_config_parameter(CONFFILE, KEYWORD, CONTENT);
    CONF_SIMTDELAY = atol(CONTENT);

    
    strcpy(KEYWORD, "WAITFORSEM");
    read_config_parameter(CONFFILE, KEYWORD, CONTENT);
    CONF_WAITFORSEM = atoi(CONTENT);
    
    strcpy(KEYWORD, "WAITSEMIMNAME");
    read_config_parameter(CONFFILE, KEYWORD, CONF_WAITSEMIMNAME);
    
    
    
    // ------------ COMPUTATION PARAMETERS, MODES --------------------------------
    
    
    strcpy(KEYWORD,"SKIP_EXISTING");
    if(read_config_parameter_exists(CONFFILE,KEYWORD)==1)
    {
        read_config_parameter(CONFFILE,KEYWORD,CONTENT);
        CONF_SKIP_EXISTING = 1;
    }
    else
        CONF_SKIP_EXISTING = 0;

    strcpy(KEYWORD,"WF_RAW_SIZE");
    read_config_parameter(CONFFILE, KEYWORD, CONTENT);
    CONF_WF_RAW_SIZE = atol(CONTENT);
 
    strcpy(KEYWORD,"MASTER_SIZE");
    read_config_parameter(CONFFILE, KEYWORD,CONTENT);
    CONF_MASTER_SIZE = atol(CONTENT);



    //  ------------ WAVEFRONT AMPLITUDE -------------------------------------------

    
    strcpy(KEYWORD,"WAVEFRONT_AMPLITUDE");
    read_config_parameter(CONFFILE, KEYWORD, CONTENT);
    CONF_WAVEFRONT_AMPLITUDE = atoi(CONTENT);

	if(CONF_WAVEFRONT_AMPLITUDE == 1)
	{
		strcpy(KEYWORD,"FRESNEL_PROPAGATION");
		read_config_parameter(CONFFILE, KEYWORD, CONTENT);
		CONF_FRESNEL_PROPAGATION = atoi(CONTENT);
	}
	else
		CONF_FRESNEL_PROPAGATION = 0;
		
    strcpy(KEYWORD,"FRESNEL_PROPAGATION_BIN");
    read_config_parameter(CONFFILE, KEYWORD, CONTENT);
    CONF_FRESNEL_PROPAGATION_BIN = atof(CONTENT);

    // ------------ POSTPROCESSING --------------
 

    strcpy(KEYWORD,"PUPIL_AMPL_FILE");
    if(read_config_parameter_exists(CONFFILE,KEYWORD)==1)
    {
        read_config_parameter(CONFFILE, KEYWORD, CONTENT);
        load_fits(CONTENT, "ST_pa", 1);
    }

    strcpy(KEYWORD,"PUPIL_PHA_FILE");
    if(read_config_parameter_exists(CONFFILE,KEYWORD)==1)
    {
        read_config_parameter(CONFFILE,KEYWORD,CONTENT);
        load_fits(CONTENT, "ST_pp", 1);
    }


    
    return RETURN_SUCCESS;
}

