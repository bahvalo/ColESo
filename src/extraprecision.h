//======================================================================================================================
// Additional macro definition for compatibility with DD/QD computations
// This is a casual header, it does no damage
//======================================================================================================================
#pragma once
#ifndef EXTRAPRECISION_HEADER
#define EXTRAPRECISION_HEADER

#include "base_config.h"
#include "base_const.h" // extraprecision.h contains specifications of the templates declared in base_const.h
#include <stdarg.h>
#include <string.h>
using namespace std;

// Check that QD library is available
#ifdef NOLIB_QD
    #error Quad double library is required but not linked
#endif

#ifndef double
    #pragma warning (disable: 4800) // 'int' : forcing value to bool 'true' or 'false' (performance warning)
    #include "qd/dd_real.h"
    #include "qd/qd_real.h"
#endif


// Nanotechnology
int IsNaN(dd_real x); 
int IsNaN(qd_real x); 

//======================================================================================================================
// Formatted output with correction of floating-point data
//======================================================================================================================

// Formatted output to std::string without the correction
string _stringprintf(const char* fmt, ...);

//----------------------------------------------------------------------------------------------------------------------
// Formatted output to std::string with correction of floating-point data
//----------------------------------------------------------------------------------------------------------------------
template<typename fpv, int numpercent>
inline string vstringprintf(const char* fmt, va_list &ap) {
    char _fmt[1000];
    string S;
    
    int on = 0;
    const char* F = fmt;
    const char* FB = NULL;
    while(*F) {
        if(on==0) {
            if(*F=='%') { FB = F; on=1; }
            else S += *F;
        }
        else {
            switch(*F) {
            case '%': on = 0; for(int j=0; j<numpercent; j++) S += "%"; break;
            case 'f': case 'F': case 'e': case 'E': 
            case 'g': case 'G': case 'a': case 'A': {
                fpv V = va_arg(ap, fpv);
                strncpy(_fmt, FB, F-FB+1); _fmt[F-FB+1]=0;
                S += _stringprintf(_fmt, static_cast<NativeDouble>(V));
                on = 0;
                break; }
            case 'i': case 'I': case 'd': case 'D': case 'u':
            case 'U': case 'o': case 'x': case 'X': {
                int V = va_arg(ap, int);
                strncpy(_fmt, FB, F-FB+1); _fmt[F-FB+1]=0;
                S += _stringprintf(_fmt, V);
                on = 0;
                break; }
            case 'c': case 'C': {
                int V = va_arg(ap, int);
                strncpy(_fmt, FB, F-FB+1); _fmt[F-FB+1]=0;
                S += _stringprintf(_fmt, char(V));
                on = 0;
                break; }
            case 's': case 'S': case 'p': {
                char* V = va_arg(ap, char*);
                strncpy(_fmt, FB, F-FB+1); _fmt[F-FB+1]=0;
                S += _stringprintf(_fmt, V);
                on = 0;
                break; }
            default: ;
            }
        }
        F++;
    }
    return S; 
}

template<typename fpv>
inline string _vstringprintf(const char* fmt, ...) {
    va_list ap;
    va_start(ap,fmt);
    string S = vstringprintf<fpv, 1>(fmt, ap);
    va_end(ap);
    return S;
}

inline int fprintf_safe(FILE* pF, const char *fmt,...){ // Write to log file
    va_list ap;
    #ifdef FIX_ARGS
        FIX_ARGS;
    #endif
    va_start(ap,fmt);
    int r=vfprintf(pF,fmt,ap);
    va_end(ap);
    return(r);
}
//======================================================================================================================


//======================================================================================================================
// Some arithmetic subroutines for extra precision mode
//======================================================================================================================

template<> inline dd_real GetPiNumber() { return dd_real::_pi; }
template<> inline dd_real GetPiNumber2() { return dd_real::_2pi; }
template<> inline qd_real GetPiNumber() { return qd_real::_pi; }
template<> inline qd_real GetPiNumber2() { return qd_real::_2pi; }
template<> inline dd_real GetLn2() { return dd_real::_log2; }
template<> inline qd_real GetLn2() { return qd_real::_log2; }

template<> inline dd_real get_eps() { return dd_real::_eps; } // smallest such that 1.0+eps != 1.0
template<> inline qd_real get_eps() { return qd_real::_eps; } // smallest such that 1.0+eps != 1.0

#endif
