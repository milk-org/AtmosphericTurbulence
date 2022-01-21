/**
 * @file    ReadConf.h
 * @brief   Read configuration file from disk
 *
 *
 */

#ifndef _ATMOSPHERICTURBULENCE_READCONF_H
#define _ATMOSPHERICTURBULENCE_READCONF_H

errno_t AtmosphericTurbulence_ReadConf(const char *restrict fnameconf,
                                       ATMTURBCONF *atmturbconf);

#endif
