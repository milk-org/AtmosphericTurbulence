# list source files (.c) other than modulename.c
set(SOURCEFILES
	make_master_turbulence_screen.c
	make_master_turbulence_screen_pow.c)

# list include files (.h) that should be installed on system
set(INCLUDEFILES
	)

# list scripts that should be installed on system
set(SCRIPTFILES
	scripts/milk-semtestspeed)




# DEFAULT SETTINGS 
# Do not change unless needed
# =====================================================================

# SRCNAME is current directory (last subdir)
set(NAME0 "${CMAKE_CURRENT_SOURCE_DIR}")
string(REPLACE "/" " " NAME1 ${NAME0})
string(REGEX MATCH "[a-zA-Z0-9_]+$" SRCNAME "${NAME0}")
message(" SRCNAME = ${SRCNAME}")

# LIBNAME constructed by removing spaces
string(REPLACE "_" "" LIBNAME "${PROJECT_NAME}${SRCNAME}")
message(" LIBNAME = ${LIBNAME}")



project(lib_${LIBNAME}_project)

include_directories ("${PROJECT_SOURCE_DIR}/src")
include_directories ("${PROJECT_SOURCE_DIR}/..")


# Library can be compiled from multiple source files
# Convention: the main souce file is named <libname>.c
# 
add_library(${LIBNAME} SHARED ${SRCNAME}.c ${SOURCEFILES})

install(TARGETS ${LIBNAME} DESTINATION lib)
install(FILES ${SRCNAME}.h ${INCLUDEFILES} DESTINATION include/${SRCNAME})

install(PROGRAMS
	${SCRIPTFILES}
	DESTINATION bin)

