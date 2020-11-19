/**
 * @file    AtmosphericTurbulence_conf.h
 *
 *
 */

#ifndef _ATMOSPHERICTURBULENCE_CONF_H
#define _ATMOSPHERICTURBULENCE_CONF_H





typedef struct {
	
	// ------------ TURBULENCE AND ATMOSPHERE PARAMETERS ---------------
	
	// reference lambda [m]
	float lambda; // CONF_LAMBDA
	
	// atmospheric seeing [arcsec] at lambda 
	float seeing; // CONF_SEEING
	
	// turbulence profile file name
	char turbulenceprof_fname[STRINGMAXLEN_FILENAME]; // CONF_TURBULENCE_PROF_FILE

	// zenith angle [rad]
	float zenithangle; // CONF_ZANGLE
	
	float sourceXpos; // CONF_SOURCE_Xpos
	float sourceYpos; // CONF_SOURCE_Ypos
	
	
	// ------------ LOW_FIDELITY OUTPUT AT REF LAMBDA ------------------
	
	int  flag_WFoutput; // CONF_WFOUTPUT 1
	char WFfileprefix[STRINGMAXLEN_FILENAME]; // CONF_WF_FILE_PREFIX
	int  flag_SHMoutput; // CONF_SHM_OUTPUT 0
	
	
	// ------------- HIGH FIDELITY OUTPUT WAVELENGTH -------------------
	
	int  flag_SWF_make; // CONF_MAKE_SWAVEFRONT 0
	int  flag_SWF_Wite2Disk; // CONF_SWF_WRITE2DISK 0
	char SWFfileprefix[STRINGMAXLEN_FILENAME]; // CONF_SWF_FILE_PREFIX
	
	int  flag_SWF_SHMoutput; // CONF_SHM_SOUTPUT
	char SWFSHMprefix[STRINGMAXLEN_FILENAME]; // CONF_SHM_SPREFIX
	int  flag_SWF_SHMouputM; // CONF_SHM_SOUTPUTM 0 output in [meter]
	
	// ------------ OUTPUT PARAMETERS ----------------------------------
	uint32_t WFsize;
	float    PupilScale;  // CONF_PUPIL_SCALE
	
	// ------------- TIMING --------------------------------------------
	int   flag_RealTime;   // CONF_ATMWF_REALTIME
	float RealTimeFactor;  // CONF_ATMWF_REALTIMEFACTOR 1.0
	float TimeStep;        // CONF_WFTIME_STEP
	float  TimeSpanCube;       // CONF_TIME_SPAN
	long  NBTimeCube;          // CONF_TIME_SPAN
	long  TimeDelayus;     // CONF_SIMTDELAY 0
	int   flag_WaitSem;    // CONF_WAITFORSEM
	char  WaitSemName[STRINGMAXLEN_IMGNAME];    // CONF_WAITSEMIMNAME
	
	// ------------ COMPUTATION PARAMETERS, MODES ----------------------
	int      flag_SkipExisting;  // CONF_SKIP_EXISTING
	uint32_t WFrawSize;          // CONF_WF_RAW_SIZE
	uint32_t MasterSize;         // CONF_MASTER_SIZE
	
	// ------------ WAVEFRONT AMPLITUDE --------------------------------
	int     flag_WFampl;          // CONF_WAVEFRONT_AMPLITUDE
	int     flag_FresnelProp;    // CONF_FRESNEL_PROPAGATION
	float   FresnelPropBin;      // CONF_FRESNEL_PROPAGATION_BIN
	
	
	
  
} ATMTURBCONF;



#endif
