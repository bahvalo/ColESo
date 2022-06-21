// *********************************************************************************************************************
// *****                                                                                                           *****
// *****                                  Collection of exact solutions (ColESo)                                   *****
// *****                    Incompressible irrotational flows around a cylinder and a sphere                       *****
// *****                                                                                                           *****
// *********************************************************************************************************************

#include "coleso.h"
#include "base_parser.h"
#include "geom_primitive.h"

//======================================================================================================================
void s_CurlFreeCylinder::ReadParams(tFileBuffer& FB) {
    tParamManager PM;
    PM.Request(X, "X");
    PM.Request(Y, "Y");
    PM.Request(R, "R");
    PM.Request(Gamma, "Gamma");  // circulation
    PM.Request(SoundVel, "SoundVel");
    PM.Request(FlowVel, "FlowVel");
    PM.Request(gam, "gam");      // specific ratio
    PM.ReadParamsFromBuffer(FB);
}
//======================================================================================================================

//======================================================================================================================
void s_CurlFreeCylinder::Init() {
    if(R<=0.0 || SoundVel<=0.0 || FlowVel<=0.0 || gam<=1.0) crash("s_CurlFreeCylinder::Init: wrong parameters");
}
//======================================================================================================================

//======================================================================================================================
void s_CurlFreeCylinder::PointValue(double, const double* coor, double* V) const {
    double x = coor[0] - X;
    double y = coor[1] - Y;
    if(x*x + y*y < tiny) {
        V[Var_U] = V[Var_V] = 0.0;
    }
    else {
        double phi = GetAngle(x, y);
        double r = R*R/(x*x + y*y);
        V[Var_U] = FlowVel * (1.0 - r*cos(2.*phi));
        V[Var_V] = - FlowVel * r * sin(2.*phi);
        double m = Gamma / (Pi2 * (x*x + y*y));
        V[Var_U] -= m * y;
        V[Var_V] += m * x;
    }
    V[Var_W] = 0.0;
    V[Var_R] = 1.0 + 0.5*(SQR(FlowVel)-SQR(V[Var_U])-SQR(V[Var_V])) / (SoundVel*SoundVel);
    V[Var_P] = SoundVel*SoundVel/gam + 0.5*(SQR(FlowVel)-SQR(V[Var_U])-SQR(V[Var_V]));
}
//======================================================================================================================

//======================================================================================================================
void s_PotentialSphere::ReadParams(tFileBuffer& FB) {
    tParamManager PM;
    PM.Request(X, "X");
    PM.Request(Y, "Y");
    PM.Request(Z, "Z");
    PM.Request(R, "R");
    PM.Request(SoundVel, "SoundVel");
    PM.Request(FlowVel, "FlowVel");
    PM.Request(gam, "gam");
    PM.ReadParamsFromBuffer(FB);
}
//======================================================================================================================

//======================================================================================================================
void s_PotentialSphere::Init() {
    if(R<=0.0 || SoundVel<=0.0 || FlowVel<=0.0 || gam<=1.0) crash("s_PotentialSphere::Init: wrong parameters");
}
//======================================================================================================================

//======================================================================================================================
void s_PotentialSphere::PointValue(double, const double* coor, double* V) const {
    double x = coor[0] - X;
    double y = coor[1] - Y;
    double z = coor[2] - Z;
    double r = sqrt(x*x + y*y + z*z);
    if(r < tiny) {
        V[Var_U] = V[Var_V] = V[Var_W] = 0.0;
    }
    else {
        double r5 = r*r*r*r*r;
        V[Var_U] = FlowVel * (1.0 - 0.5*R*R*R/r5*(3.0*x*x-r*r));
        V[Var_V] = - FlowVel * 1.5 * R*R*R*x*y/r5;
        V[Var_W] = - FlowVel * 1.5 * R*R*R*x*z/r5;
    }
    double dKin = 0.5*(SQR(FlowVel)-SQR(V[Var_U])-SQR(V[Var_V])-SQR(V[Var_W]));
    V[Var_R] = 1.0 + dKin / (SoundVel*SoundVel);
    V[Var_P] = SoundVel*SoundVel/gam + dKin;
}
//======================================================================================================================
