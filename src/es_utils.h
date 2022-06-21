//======================================================================================================================
// EXACT SOLUTION MODULE
//======================================================================================================================
#pragma once
#ifndef ES_UTILS_HEADER
#define ES_UTILS_HEADER

#include "base_containers.h"

enum tGItype {
    GI_LEGENDRE = 0,  // Gauss-Legendre quadratures for int_0^1 f(x) dx
    GI_JACOBI1  = 1   // Gauss-Jacobi quadratures for int_0^1 f(x)/sqrt(x) dx
};

template<typename fpv>
void GaussPointsInit(int NumPoints, fpv* Nodes, fpv* Weights, tGItype gi_type);

// calculating a^n where 'a' is a double and 'n' is a non-negative const int
template<int n> inline double pow_to_const_int(double a) { return pow_to_const_int<n-1>(a) * a; }
template<> inline double pow_to_const_int<0>(double) { return 1.0; }

template<typename fpv>
struct tGaussIntegrator {
    int GR;
    fpv *GN, *GC;

    tGaussIntegrator() { GN=GC=NULL; GR=4*sizeof(fpv); }
    tGaussIntegrator(int _GR) { GN=GC=NULL; GR=0; Init(_GR); }
    tGaussIntegrator(const tGaussIntegrator<fpv>& gi) { GN=GC=NULL; *this = gi; }
    inline void Clear() { GN=GC=NULL; GR=0; } // not a destructor - just nullify pointers and data
    inline void Free() { if(GN!=NULL) delete[] GN; GN=GC=NULL; } // do not touch GR
    ~tGaussIntegrator() { Free(); }

    inline void Init(int _GR = 0, tGItype GItype = GI_LEGENDRE) { // выделение памяти и определение узлов и коэффициентов квадратуры
        Free();
        if(_GR==0) _GR = GR;
        if(_GR <= 0 || _GR > 0x7FFF) return;
        GR = _GR;
        GN = new fpv[GR*2];
        GC = GN + GR;
        GaussPointsInit<fpv>(GR, GN, GC, GItype);
    }
    inline tGaussIntegrator<fpv>& operator=(const tGaussIntegrator<fpv> &gi){
        if(this == &gi) return *this;
        Free();
        GR = gi.GR;
        if(gi.GN!=NULL) { // if the source array was allocated, allocate the destination
            GN = new fpv[GR*2];
            GC = GN + GR;
            for(int i=0; i<GR*2; ++i) GN[i]=gi.GN[i];
        }
        return *this;
    }
};

// Calculation using compound Gauss formula of the integral
// int exp(-alpha^2 (x-x0)^2/2) (x-x0)^k h(x) x^gamma dx, x = xmin..xmax,
// where k=0,1,... (not big), alpha>0, gamma=0 (mode=0) or -1/2 (mode=1),
// h -- analytic function on [0..xmax] with convergence raduis >=r at each point, q = 1/r.
// If mode==1, xmin should be zero
template<typename fpv>
struct tCompoundGaussIntegrator {
    tGaussIntegrator<fpv> GLI, GJI; // Gauss-Legendre and Gauss-Jacobi integrators
    double reducer; // Workaround: multiply number of quadrature nodes by this value (set < 1 to reduce computation time)
    tCompoundGaussIntegrator() { reducer = 1.0; }
    inline void Clear() { GLI.Clear(); GJI.Clear(); } // not a destructor - just nullify pointers and data
    inline void Free() { GLI.Free(); GJI.Free(); } // the same as destructor
    inline void Init(int GR = 0) {
        if(GR<=0) GR = sizeof(fpv)*4;
        GLI.Init(GR, GI_LEGENDRE);
        GJI.Init(GR, GI_JACOBI1);
    }
    // Automatic integration subroutine
    template<unsigned int N>
    tFixBlock<fpv, N> Integrate(fpv xmin, fpv xmax, fpv alpha, fpv x0, int k, int mode,
                                fpv M, fpv q, tFixBlock<fpv ,N> (*fun_ptr)(fpv, void*), void*) const;
};

#endif
