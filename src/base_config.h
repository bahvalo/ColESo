//======================================================================================================================
// Lowest-level header with build configuration settings, MUST BE included everywhere!
//======================================================================================================================
#pragma once
#ifndef BASE_CONFIG_HEADER
#define BASE_CONFIG_HEADER

// Enable this to compile ColESo with dd and qd precision routines
// (when building using a makefile, this flag is set automatically)
#define EXTRAPRECISION_COLESO

//======================================================================================================================
//
//======================================================================================================================

//C++11 specific things
#if __cplusplus >= 201103
    #define OVERRIDE override
#else
    #define OVERRIDE
#endif

// Compiler warnings settings - НИКАКИХ НОВЫХ ВОРНИНГОВ СЮДА НЕ ПИХАТЬ!
#pragma warning (disable: 4127) // отключить ругань майкрософтовского компеллера с W4 на константы в условном операторе
#pragma warning (disable: 6993) // ворнинг аналайзера об игнорировании ОпенМП
#pragma warning (disable: 4068) // ворнинг аналайзера об игнорировании незнакомой прагмы ivdep
#pragma warning (disable: 6211) // ворнинг на каждый new, что надо ловить эксепшн. пока отключено. 

// Compiler-dependent defines, if not set at compiler call
#ifndef _CRT_SECURE_NO_DEPRECATE
    #define _CRT_SECURE_NO_DEPRECATE 1
#endif
#ifndef _CRT_SECURE_NO_WARNINGS
    #define _CRT_SECURE_NO_WARNINGS  1
#endif
#ifndef __USE_MINGW_ANSI_STDIO
    #define __USE_MINGW_ANSI_STDIO // for %llu in fprintf 
#endif

// Default system includes
#include <stdio.h>
#include <stdlib.h>

typedef double NativeDouble;  // remember system type for the case of later redefinition

#endif    
