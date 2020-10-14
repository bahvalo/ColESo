//======================================================================================================================
// Geometry library level 0 - basic geometric subroutnies
//======================================================================================================================
#pragma once
#ifndef GEOM_L0_HEADER
#define GEOM_L0_HEADER
#include "const.h"
#include <math.h>
//----------------------------------------------------------------------------------------------------------------------
// Inline functions  (operations are for vector of 3 components)
//----------------------------------------------------------------------------------------------------------------------

// Vector absolute value (L2 norm)
inline double GL_Absolute(const double *V){ return sqrt(VDOT(V,V));}

// Dot product 
inline double GL_Scalar_Product(const double *V1, const double *V2){ return VDOT(V1,V2); }

// Vector product
template<typename fpv>
inline void GL_Vector_Product(const fpv *A, const fpv *B, fpv *CC){
    CC[0] = A[1]*B[2] - A[2]*B[1];
    CC[1] = A[2]*B[0] - A[0]*B[2];
    CC[2] = A[0]*B[1] - A[1]*B[0];
}
 
// Absolute value of vector product 
inline double GL_Absolute_Vector_Product(double* V1, double* V2) {
    double V3[3];
    GL_Vector_Product(V1, V2, V3);
    return GL_Absolute(V3);
}

// Triple product
inline double GL_ScalarTriple(const double *P1, const double *P2, const double *P3){ 
    return P1[0]*(P2[1]*P3[2]-P3[1]*P2[2]) + P1[1]*(P2[2]*P3[0]-P2[0]*P3[2]) + P1[2]*(P2[0]*P3[1]-P2[1]*P3[0]);
}

// Cosine of angle between 2 vectors 
inline double GL_Cos_2Vectors(const double *V1, const double *V2){ 
    return VDOT(V1,V2)/sqrt(VDOT(V1,V1)*VDOT(V2,V2));
}

// Distance between 2 points
inline double GL_Dist3D(const double *v1, const double *v2){ // 3D version
    return sqrt((v1[0]-v2[0])*(v1[0]-v2[0])+(v1[1]-v2[1])*(v1[1]-v2[1])+(v1[2]-v2[2])*(v1[2]-v2[2]));
}
inline double GL_Dist2D(const double *v1, const double *v2){  // 2D version 
    return sqrt((v1[0]-v2[0])*(v1[0]-v2[0])+(v1[1]-v2[1])*(v1[1]-v2[1]));
}

// Area of triangle 
inline double GL_S_Triangle(const double* P1, const double* P2, const double* P3){ //3D version 
    double V[3] = { P1[1]*(P2[2]-P3[2]) - P2[1]*(P1[2]-P3[2]) + P3[1]*(P1[2]-P2[2]),
                    P1[2]*(P2[0]-P3[0]) - P2[2]*(P1[0]-P3[0]) + P3[2]*(P1[0]-P2[0]),
                    P1[0]*(P2[1]-P3[1]) - P2[0]*(P1[1]-P3[1]) + P3[0]*(P1[1]-P2[1])}; 
    return sqrt(VDOT(V,V)*0.25);
}
inline double GL_S_Triangle_2D(const double* P1, const double* P2, const double* P3){ //2D version
    return 0.5*fabs((P3[1]-P1[1])*(P2[0]-P1[0]) - (P3[0]-P1[0])*(P2[1]-P1[1])); 
}

// Volume of tetrahedron 
inline double GL_V_Tetrahedron_Signed(const double *P1, const double *P2, const double *P3, const double *P4){ 
    double V2[3] = {P2[0]-P1[0], P2[1]-P1[1], P2[2]-P1[2]};
    double V3[3] = {P3[0]-P1[0], P3[1]-P1[1], P3[2]-P1[2]};
    double V4[3] = {P4[0]-P1[0], P4[1]-P1[1], P4[2]-P1[2]};
    return -GL_ScalarTriple(V2, V3, V4)*C1_6;
}
inline double GL_V_Tetrahedron(const double *P1, const double *P2, const double *P3, const double *P4){
    return (fabs(GL_V_Tetrahedron_Signed(P1,P2,P3,P4)));
}

// Gets plane (coefficients A,B,C,D) by 3 given points 
inline void GL_Plane_3Points(const double *P1, const double *P2, const double *P3,
                      double &A, double &B, double &C, double &D){
    A = P1[1]*(P2[2]-P3[2]) - P2[1]*(P1[2]-P3[2]) + P3[1]*(P1[2]-P2[2]);
    B = P1[2]*(P2[0]-P3[0]) - P2[2]*(P1[0]-P3[0]) + P3[2]*(P1[0]-P2[0]);
    C = P1[0]*(P2[1]-P3[1]) - P2[0]*(P1[1]-P3[1]) + P3[0]*(P1[1]-P2[1]);
    D =-P1[0]*(P2[1]*P3[2]-P3[1]*P2[2]) + P2[0]*(P1[1]*P3[2]-P3[1]*P1[2])-P3[0]*(P1[1]*P2[2]-P2[1]*P1[2]);
}

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

void RotateSymmetricTensor2D(const double* vec, double phi, double* out);
void RotateSymmetricTensor3D(const double* vec, double phi, double* out);

// Normalize vector and put the magnitude in vec[3] position
inline void NormalizeVector3D(double* vec){ // Convert {xyz} to {xyzS}
    if( (vec[3]=sqrt(VDOT(vec,vec))) < tiny) return;
    for(int icoor=0; icoor<3; icoor++) vec[icoor] /= vec[3];
}

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
    double Phi = GetAngle(x, y) - NullAngle;
    RoundToCentre(Phi, Pi2);
    return Phi + NullAngle;
}
inline double GetAngle(double x, double y, double* r0) { return GetAngle(x-r0[0], y-r0[1]); }
inline double GetAngle(double x, double y, double* r0, double NullAngle) { return GetAngle(x-r0[0], y-r0[1], NullAngle); }

//----------------------------------------------------------------------------------------------------------------------
// Heavy  functions 
//----------------------------------------------------------------------------------------------------------------------

// Finds intersection point between line and plane 
void GL_Point_Line_Plane(const double *P1, const double *P2, double A, double B, double C, double D, double *P);

// Finds circumcenter CC of triangle (R1,R2,R3). 
// TriangleCentreError - optional, returns numerical error 
// InsideOnly - optional, 1 - returns point that belong to triangle (point on the edge if obtuse triangle)
//                        0 - returns actual point, wherever 
// return value - 1 for obtuse triangle in case InsideOnly=1, 0 otherwise
int GL_TriangleCircumcentre(const double* R1, const double* R2, const double* R3, double* CC, 
                            double* TriangleCentreError=NULL, int InsideOnly=1);

// Finds circumcenter C of tetrahedron (P1,P2,P3,P4)
// InsideOnly - optional, 1 - returns point that belong to tetrahedron (point on the edge if obtuse tetrahedron)
//                        0 - returns actual point, wherever it is
// return value - 1 for obtuse tetrahedron in case of InsideOnly=1, 0 otherwise
int GL_TetrahedronCircumcentre(const double* P1, const double* P2, const double* P3, const double* P4, 
                               double* C, int InsideOnly=1);

// Checks if point (x,y) is outside of circular segment (A,B,R), returns 1 if ourside. 
int GL_PointOutOfCircularSegment(double x, double y, double* A, double* B, double R); //2D only

// Finds projection (xout,yout) of point (x,y) on segment (c1,c2)
double GL_ProjectionPointSegment(double x, double y, const double* c1, const double* c2, double& xout, double& yout);//2D only

// Finds projection point (mp) of point m on the segment(edge) [c1,c2].
// Return distance from mp to m if mp lies on segment and c1 or c2 overwise.
double GL_ProjectionPointSegment(const double *m, const double* c1, const double* c2, double* mp); //3D 
 
#endif
