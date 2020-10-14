//======================================================================================================================
// File with personal configuration settings, MUST BE included everywhere
//======================================================================================================================

#ifndef PERSONAL_HEADER
#define PERSONAL_HEADER

// Disable some warnings
#define _CRT_SECURE_NO_DEPRECATE 1
#define _CRT_SECURE_NO_WARNINGS  1
#pragma warning (disable: 4127) // Microsoft's compiler's warning about if(1) and if(0)
#pragma warning (disable: 6211) // No try/catch block for new[] operators
#pragma warning (disable: 4068) // Microsoft's compiler's warning about unknown pragma 
#pragma GCC diagnostic ignored "-Wformat-security" // turn off "format not a string literal and no format arguments" gcc warning

// Default system includes
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <string.h>
#include <vector>

// Print, exit and memory functions. The user can substitute his own functions here
#define crash(...) exit(printf(__VA_ARGS__))
std::string stringprintf(const char* fmt, ...);
#define ASSERT(X,...) if(!(X)) crash( (stringprintf("Assert (%s) failed: ", #X) +  stringprintf(__VA_ARGS__)).c_str() );
#define pprintf(...) printf(__VA_ARGS__)
#define pprintf0(...) printf(__VA_ARGS__)
#define printf0(...) printf(__VA_ARGS__)
template <typename T> inline T* GimmeMem(size_t N, const char* = 0){ return new T[N]; }
template <typename T> inline void FreeMem(T*& p) { delete[] p; p = 0; }

// Miscellaneous
#define OVERRIDE // #define OVERRIDE override if C++11 standart is supported
#define NativeDouble double // Do not change this!

// Uncomment this for use double-double or quad-double precision (if you have QD library)
//#define EXTRAPRECISION_COLESO

#endif
