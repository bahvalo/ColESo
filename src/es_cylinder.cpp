// *********************************************************************************************************************
// *****                                                                                                           *****
// *****                                  Collection of exact solutions (ColESo)                                   *****
// *****                                  Gaussian pulse by cylinder diffraction                                   *****
// *****                                                                                                           *****
// *********************************************************************************************************************

#include "parser.h"
#include "coleso.h"
#include "es_utils.h"
#include "geom_primitive.h"
#ifdef _NOISETTE
#include "lib_base.h"
#endif

//======================================================================================================================
void s_Cylinder::ReadParams(tFileBuffer& FB) {
    tParamManager PM;
    // Cylinder position and radius
    PM.Request(X, "X");
    PM.Request(Y, "Y");
    PM.Request(Radius, "Radius");

    // Gaussian impulse parameters
    PM.Request(Xterm, "Xterm");
    PM.Request(Yterm, "Yterm");
    PM.Request(Aterm, "Aterm");
    PM.Request(Bterm, "Bterm");

    // Integration parameters
    PM.Request(NumK, "NumK"); // default: 200
    PM.Request(NumN, "NumN"); // default: 200
    PM.Request(Kmax, "Kmax"); // default: 6/Bterm
    PM.ReadParamsFromBuffer(FB);
}
//======================================================================================================================


//======================================================================================================================
void s_Cylinder::Init(void) {
    // Discretization by k
    if(IsNaN(Kmax) || Kmax < tiny) Kmax = 6.0 / Bterm;
    GI.Init(NumK);

    // Discretization by radius for intrgrals calculation
    const double R0 = sqrt(SQR(Xterm - X) + SQR(Yterm - Y)); // Distance between cylinder and initial pulse centers
    double rmin = R0 - 7.0 * Bterm;
    double rmax = R0 + 7.0 * Bterm;
    if(rmin < Radius) crash("s_Cylinder::Init error: initial pulse overlaps the cylinder");
   
    // Memory allocation
    Coeff.resize(NumK * NumN, 0.0);

    // For each k, calculation the integral by r
    pprintf("s_Cylinder::Init... ");
    const int n1 = NumN+1;
    double* BesselI = new double[n1*2];
    NativeDouble* BesselOuter = new NativeDouble[NumN];

    const double Alpha = CLN2 / SQR(Bterm);
    double* CoeffOld = new double[NumK * NumN]; // storing the previous result
    for(int ik=0; ik<NumK; ik++) for(int in=0; in<NumN; in++) CoeffOld[ik*NumN+in] = huge;
    double err = huge*NumN*NumK;

    // To integrate by 'r', we will auto-adapt the grid
    // Starting from NumR=32 nodes, we will double NumR until we reach the desired accuracy, but not more than NumR=2048
    for(int NumR = 16; NumR <= 2048 && err > 1e-10; NumR*=2) {
        for(int ik=0; ik<NumK; ik++)
            for(int in=0; in<NumN; in++)
                Coeff[ik*NumN+in] = 0.0;

        double dr = (rmax - rmin) / NumR;
        for(int ir=0; ir<NumR; ir++) {
            printf0("\rs_Cylinder::Init (lasterr=%5.0e): processing %i of %i     \r", err, ir, NumR);
            double r = rmin + ir * dr;

            // Computation of modified Bessel functions
            sf74r_c<double>(2.*Alpha*r*R0, NumN+1, 1e-10, false, BesselI, BesselI+n1);
            double agauss = exp(-Alpha* SQR(r-R0));

            for(int ik=0; ik<NumK; ik++) {
                double k = GI.GN[ik] * Kmax;

                BesselOuterFunction(NativeDouble(k*r), NativeDouble(k*Radius), NumN, BesselOuter, NULL);
                for(int in=0; in<NumN; in++) {
                    double gauss = BesselI[in] * (in==0 ? 1.0 : 2.0) * agauss;
                    double dReA_new = - BesselOuter[in] * gauss * r * k;
                    Coeff[ik*NumN+in] += dReA_new * Aterm * dr;
                }
            }
        }

        // Считаем разницу и сохраняем копию коэффициентов
        err = 0.0;
        for(int i=0; i<NumK*NumN; i++) { err += fabs(Coeff[i] - CoeffOld[i]); CoeffOld[i] = Coeff[i]; }
    }
    delete[] BesselI;
    delete[] BesselOuter;
    delete[] CoeffOld;

    pprintf0(" done (err = %e)\n", err);

    for(int ik=0; ik<NumK; ik++)
        for(int in=0; in<NumN; in++)
            if(IsNaN(Coeff[ik*NumN+in])) crash("s_Cylinder::Init error: NaN detected");
}
//======================================================================================================================

//======================================================================================================================
// Debug output
//======================================================================================================================
void s_Cylinder::PrintCoeffMatrix(FILE* f) {
    for(int ik=0; ik<NumK; ik++) {
        for(int in=0; in<NumN; in++) {
            fprintf(f, "%e ", Coeff[ik*NumN+in]);
        }
        fprintf(f, "\n");
    }
}
//======================================================================================================================

//======================================================================================================================
void s_Cylinder::PointValue(double t, const double* coord, double* V) const {
    double x = coord[0] - X;
    double y = coord[1] - Y;
    double r = sqrt(x*x+y*y);
    double Phi_loc = GetAngle(x, y) - GetAngle(Xterm - X, Yterm - Y);

    if(r < 0.999999*Radius) {
        for(int ivar=Var_R; ivar<=Var_P; ivar++) V[ivar] = 0.0;
        return;
    }
    if(r < Radius) r = Radius;

    double *cosn = GimmeMem<double>(NumN*2, "s_Cylinder");
    double *sinn = cosn + NumN;
    for(int n=0; n<NumN; n++) { cosn[n] = cos(Phi_loc*n); sinn[n] = sin(Phi_loc*n); }

    NativeDouble *Phi = GimmeMem<NativeDouble>(NumN, "s_Cylinder");
    NativeDouble *dPhi = GimmeMem<NativeDouble>(NumN, "s_Cylinder");

    double summ[3] = {0.,0.,0.}; // p, u_r, u_phi
    for(int ik=0; ik<NumK; ik++) {
        double k = GI.GN[ik] * Kmax;
        BesselOuterFunction(NativeDouble(k*r), NativeDouble(k*Radius), NumN, Phi, dPhi);

        double sumn[3] = {0.,0.,0.};
        for(int n=0; n<NumN; n++) {
            double phi = -Phi[n];
            double dphi = -dPhi[n];
            sumn[0] += Coeff[ik*NumN+n] * cosn[n] * phi;
            sumn[1] -= Coeff[ik*NumN+n] * cosn[n] * dphi;
            sumn[2] += Coeff[ik*NumN+n] * sinn[n] * phi * n;
        }
        summ[0] += sumn[0] * cos(k*t) * Kmax * GI.GC[ik];
        summ[1] += sumn[1] * sin(k*t) * Kmax * GI.GC[ik];
        summ[2] += sumn[2] * sin(k*t) / (k*r) * Kmax * GI.GC[ik];
    }
    
    V[Var_R] = summ[0];
    V[Var_U] = summ[1] * cos(Phi_loc) + summ[2] * sin(Phi_loc);
    V[Var_V] = - summ[1] * sin(Phi_loc) + summ[2] * cos(Phi_loc);
    V[Var_W] = 0.0;
    V[Var_P] = summ[0];
    if(IsNaN(summ[0]) || IsNaN(summ[1])|| IsNaN(summ[2])) crash("s_Cylinder::PointValue error: NaN detected");

    FreeMem(cosn);
    FreeMem(Phi);
    FreeMem(dPhi);
}
//======================================================================================================================
