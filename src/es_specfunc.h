//======================================================================================================================
// SPECIAL FUNCTIONS
//======================================================================================================================
#pragma once
#ifndef ES_SPECFUNC_HEADER
#define ES_SPECFUNC_HEADER

#include "base_config.h"

// Bessel functions of the first and second kind
template<typename fpv> fpv BesselJ0(fpv arg);
template<typename fpv> fpv BesselJ1(fpv arg);
template<typename fpv> fpv BesselN0(fpv arg);
template<typename fpv> fpv BesselN1(fpv arg);

// Bessel function of the first kind and its derivatives
NativeDouble BesselJ(int l, NativeDouble x, NativeDouble* dbdx = NULL, NativeDouble* d2bdx2 = NULL,
                     NativeDouble* d2bdx3 = NULL, NativeDouble* d2bdx4 = NULL, NativeDouble* d2bdx5 = NULL);

// k-th zero of J'_n(x), where k=RadialMode, n=AngularMode
NativeDouble BesselPrimeZero(int AngularMode, int RadialMode, int log=0);

// Calculation of Bessel functions J_n(x) for n=0,...,nmax-1
// ATTENTION! Output array J must be of size S=nmax+MIN(nmax,65)
// Buffer array 'buf' must be of size S; if buf==NULL, we assume that J is of size 2*S
void BesselFunctionsJ(NativeDouble x, int nmax, NativeDouble *J, NativeDouble* buf = NULL);

// Calculation of Bessel functions Y_n(x) for n=0,...,nmax-1
// Output array Y must be of size nmax
// Returns n+1 where n is the maximal index such that Y_n(x) is successfully calculated
int BesselFunctionsY(NativeDouble x, int nmax, NativeDouble *Y);

// Bessel function of complex argument
template<typename real>
void BesselJComplex(int l, real x, real y, real* ReBesselJ, real* ImBesselJ, real* ReDBesselJ = NULL, real* ImDBesselJ = NULL);

// Radial part of eigenfunctions of the Laplace operator in a circle exterior (r > y)
// with the Neumann conditions at r = rc
// These functions are given by  u_{nu,k}(r) = phi(nu, k*r, k*y), r >= y,
// phi(nu,x,y) = (- N'_nu(y) J_nu(x) + J'_nu(y) N_nu(x)) / sqrt((J'_nu(y))^2 + (N'_nu(y))^2), x >= y.
// This subroutine calculates the values of phi(nu,x,y) and dphi/dx(nu,x,y) for given x,y and each nu=0,...,Nmax-1
void BesselOuterFunction(NativeDouble x, NativeDouble y, int Nmax, NativeDouble* phi, NativeDouble* dphi);

// Modified Bessel functions for n = 0, ..., n1-1
template<typename fpv> void BesselI(fpv x, int n1, fpv eps, int s, fpv *t);

void PrintBesselDerivatives(int l);

NativeDouble cephes_erf(NativeDouble x); // erf

template<typename fpv> void Fresnl(fpv x, fpv& S, fpv& C); // Fresnel integrals

NativeDouble expn( int n, NativeDouble x ); // int_1^\infty exp(-x*t) / t^n dt


#endif
