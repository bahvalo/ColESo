//======================================================================================================================
// Geometry library level 0 - basic geometric subroutines
//======================================================================================================================
#pragma once
#ifndef GEOM_L0_HEADER
#define GEOM_L0_HEADER
#include "base_macro.h"
#include "base_const.h"
#include <math.h>

// Vector product
template<typename fpv>
inline void GL_Vector_Product(const fpv *A, const fpv *B, fpv *CC){
    CC[0] = A[1]*B[2] - A[2]*B[1];
    CC[1] = A[2]*B[0] - A[0]*B[2];
    CC[2] = A[0]*B[1] - A[1]*B[0];
}

// Вычисление биномиальных коэффициентов
inline int Factorial(int n) {int i, f = 1; for(i=2; i<=n; i++) f*=i; return f;}
inline int BinomialCoeff0(int n, int k){return Factorial(n) / Factorial(k) / Factorial(n-k);}

template<typename fpv>
inline void RoundToCentre(fpv &R, fpv limit) {
    if(limit > 0.5 * huge) return;
    R /= limit;
    R -= floor(R); // Взяли дробную часть, в промежутке (0,1)
    if(R>0.5) R -= 1.0; // Перевели в промежуток (-0.5,0.5)
    R *= limit;
}
// угол между осью х и направлением на точку (х,у)
template <class fpv>
fpv GetAngle(fpv x, fpv y) {
    fpv absx = fabs(x), absy = fabs(y);
    // If x=y=0, angle is undefined. This subroutine returns 0
    if(absx < get_min_value<fpv>() && absy < get_min_value<fpv>()) {
        return fpv(0.0);
    }
    if(absx >= absy) {
        fpv Phi = atan(y/x);
        if(x < 0.0) return Phi + GetPiNumber<fpv>();
        else {
            if(Phi<0.0) return Phi + GetPiNumber2<fpv>();
            else return Phi;
        }
    }
    else {
        fpv Phi = -atan(x/y);
        if(y < 0.0) return Phi + 1.5*GetPiNumber<fpv>(); // 1.5 представляется в памяти точно
        else return Phi + 0.5*GetPiNumber<fpv>(); // 0.5 представляется в памяти точно
    }
}

inline double GetAngle(double x, double y, double NullAngle) {
    if(fabs(x) < get_min_value<double>() && fabs(y) < get_min_value<double>()) return NullAngle;
    double Phi = GetAngle(x, y);
    if(fabs(Phi - NullAngle) < GetPiNumber<double>()) return Phi;
    Phi -= NullAngle;
    RoundToCentre(Phi, Pi2);
    return Phi + NullAngle;
}
inline double GetAngle(double x, double y, double* r0) { return GetAngle(x-r0[0], y-r0[1]); }
inline double GetAngle(double x, double y, double* r0, double NullAngle) { return GetAngle(x-r0[0], y-r0[1], NullAngle); }

// Rotation of vectors
// 2D - rotation of vec on angle phi. And version with given sin, cos of phi already calculated
template<typename fpv>
inline void RotateVector2D(const fpv* vec, fpv sinphi, fpv cosphi, fpv* out){
    fpv v[2]={vec[0], vec[1]};
    out[0] = v[0]*cosphi - v[1]*sinphi;
    out[1] = v[0]*sinphi + v[1]*cosphi;
}
template<typename fpv>
inline void RotateVector2D(const fpv* vec, fpv phi, fpv* out){ RotateVector2D(vec,sin(phi),cos(phi),out); }
template<typename fpv>
inline void RotateVector2D(fpv* vec, fpv phi) {                RotateVector2D(vec, phi, vec); }
template<typename fpv>
inline void RotateVector2D(fpv* vec, fpv sinphi, fpv cosphi) { RotateVector2D(vec, sinphi, cosphi, vec); }

// 3D - rotation of vector vec around vector e on angle phi. And version with given sin, cos of phi already calculated
template<typename fpv>
inline void RotateVector(const fpv* vec, const fpv* e, fpv sinphi, fpv cosphi, fpv* out) {
    fpv er = VDOT(e, vec), evr[3];
    GL_Vector_Product(e, vec, evr);
    out[0] = vec[0]*cosphi + e[0]*er*(1.0-cosphi) + evr[0]*sinphi;
    out[1] = vec[1]*cosphi + e[1]*er*(1.0-cosphi) + evr[1]*sinphi;
    out[2] = vec[2]*cosphi + e[2]*er*(1.0-cosphi) + evr[2]*sinphi;
}
inline void RotateVector(const double* vec, const double* e, double phi, double* out){
    RotateVector(vec, e, sin(phi), cos(phi), out);
}

#endif
