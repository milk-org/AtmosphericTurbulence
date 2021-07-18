set(LIBNAME "milkAtmosphericTurbulence")

message("")
message(" SRCNAME = ${SRCNAME} -> LIBNAME = ${LIBNAME}")

# list source files (.c) other than modulename.c
set(SOURCEFILES
	makeHV_CN2prof.c
	make_master_turbulence_screen.c
	make_master_turbulence_screen_pow.c
	make_AtmosphericTurbulence_vonKarmanWind.c
	AtmosphericTurbulence_prediction.c
	make_wavefront_sequence.c
	ReadConf.c
	Z_Air.c
	)

# list include files (.h) that should be installed on system
set(INCLUDEFILES
	)

# list scripts that should be installed on system
set(SCRIPTFILES
	scripts/milk-atmturb-runturb
	scripts/milk-atmturb-mkconf
)

set(LINKLIBS
	CLIcore
)


# DEFAULT SETTINGS
# Do not change unless needed
# =====================================================================

# SRCNAME is current directory (last subdir)
set(NAME0 "${CMAKE_CURRENT_SOURCE_DIR}")
string(REPLACE "/" " " NAME1 ${NAME0})
string(REGEX MATCH "[a-zA-Z0-9_]+$" SRCNAME "${NAME0}")
message(" SRCNAME = ${SRCNAME}")

project(lib_${LIBNAME}_project)

# Library can be compiled from multiple source files
# Convention: the main souce file is named <libname>.c
#
add_library(${LIBNAME} SHARED ${SRCNAME}.c ${SOURCEFILES})

target_include_directories(${LIBNAME} PUBLIC ${PROJECT_SOURCE_DIR}/.. ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(${LIBNAME} PRIVATE ${LINKLIBS})

install(TARGETS ${LIBNAME} DESTINATION lib)
install(FILES ${SRCNAME}.h ${INCLUDEFILES} DESTINATION include/${SRCNAME})
install(FILES ${INCLUDEFILES} DESTINATION include/${SRCNAME})
install(PROGRAMS ${SCRIPTFILES} DESTINATION bin)

