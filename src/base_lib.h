//======================================================================================================================
//Basic macro definitions 
//======================================================================================================================
#pragma once
#ifndef BASE_INFRA_HEADER
#define BASE_INFRA_HEADER

#include "base_config.h"
#include "base_macro.h"
#include <vector>
#include <string>
using namespace std;

// --------------------------------------------------------------------------------------------------------------------
//   Printing & log interfaces
// --------------------------------------------------------------------------------------------------------------------

#define pprintf printf
#define printf0 printf
#define pprintf0 printf

// Formatted output to std::string conversion
std::string correctString(const char *s);
std::string stringprintf(const char* fmt, ...);

// Abnormal termination with writing error message to log, stdout and stderr
#define crash(...) exit(printf(__VA_ARGS__))
#define ASSERT(X,...) if(!(X)) crash("Assert (%s) failed: ", #X);
#define SAFE_ASSERT ASSERT

// Allocation wrapper - arrays
template <typename T> inline T* GimmeMem(size_t N, const char* = 0){ return new T[N]; }
template <typename T> inline void FreeMem(T*& p) { delete[] p; p = 0; }
template <typename T> inline T* GimmeMemSingle(const char* = 0){ return new T; }
template <typename T> inline void FreeMemSingle(T*& p) { delete p; p = 0; }


// --------------------------------------------------------------------------------------------------------------------
// NaN-otechnology
// --------------------------------------------------------------------------------------------------------------------
inline int IsNaN(NativeDouble x) {
    unsigned char* X = (unsigned char*) &x;
    if(((X[6] | 0x0F) & X[7] & 0x7F) == 0x7F) return 1;
    return !(((X[6]^0xF0)&0xF7) || ((X[7]^0x7F)&0x7F));
}
template<typename T>
inline void MakeNaN(T& x) {
    unsigned char* X = (unsigned char*) &x;
    for(unsigned int i=0; i<sizeof(T); i++) X[i]=0xFF;
}

//----------------------------------------------------------------------------------------------------------------------
// Глобальные инлайнеры
//----------------------------------------------------------------------------------------------------------------------

template<typename T> inline void SWAP(T &a, T &b){ T swapbuf = a; a = b; b = swapbuf;}

//to silence unused warnings - only for virt functions with given interface!
inline void UnuseIt(double x){x = 0.0;}
inline void UnuseIt(void* x){x = 0;}
inline void UnuseIt(int x){x = 0;}

#endif
