// *********************************************************************************************************************
// *****                                                                                                           *****
// *****                                  Collection of exact solutions (ColESo)                                   *****
// *****                      Blasius profile and compressible modification of Blasius profile                      *****
// *****                                                                                                           *****
// *********************************************************************************************************************

#include "base_parser.h"
#include "coleso.h"

// *********************************************************************************************************************
// *****                                                                                                           *****
// *****                                      Blasius boundary layer                                               *****
// *****                                                                                                           *****
// *********************************************************************************************************************

//======================================================================================================================
void s_Blasius::ReadParams(tFileBuffer& FB) {
    tParamManager PM;
    PM.Request(Rey, "re");
    PM.Request(gam, "gam");
    PM.Request(FlowVel, "FlowVel");
    PM.Request(SoundVel, "SoundVel");
    PM.ReadParamsFromBuffer(FB);
}
//======================================================================================================================

//======================================================================================================================
// Calculating the solution of Blasius problem: 
// 2 f''' + f'' f = 0 with boundary conditions f(0) = f'(0) = 0, f(+infty) = 1.
// In fact, we know that f''(0) = 0.3320573362151946 and solve the Cauchy problem
// F' = G(F), F = (f, f', f''), G(f) = (f', f'', -f f''/2)
// Output: UV0 = f', UV1 = 1/2*(eta*f' - f), where eta is the independent variable
//======================================================================================================================
const double s_Blasius::alpha = 0.3320573362151946;
void s_Blasius::Init() {
    // Parameters of the numerical solution. Chosen in a way to get the solution with double precision
    const double L_sim = 12.5;  // size of region
    const int N_sim = 4096;     // number of intervals
    const double h = L_sim / N_sim; // mesh step
    d_eta = h;

    UV0.resize(N_sim + 1);
    UV1.resize(N_sim + 1);
    double F[3] = {0., 0., alpha}; // initial data
    UV0[0] = UV1[0] = 0.0;

    for (int i = 1; i <= N_sim; ++i){
        double k1[3] = {h*F[1], h*F[2], -0.5*h*F[0]*F[2]};
        double k2[3] = {h*(F[1]+0.5*k1[1]), h*(F[2]+0.5*k1[2]), -0.5*h*((F[0]+0.5*k1[0])*(F[2]+0.5*k1[2]))};
        double k3[3] = {h*(F[1]+0.5*k2[1]), h*(F[2]+0.5*k2[2]), -0.5*h*((F[0]+0.5*k2[0])*(F[2]+0.5*k2[2]))};
        double k4[3] = {h*(F[1]+    k3[1]), h*(F[2]+    k3[2]), -0.5*h*((F[0]+    k3[0])*(F[2]+    k3[2]))};
        for(int j=0; j<3; j++) F[j] += C1_6*(k1[j]+k4[j]) + C1_3*(k2[j]+k3[j]);

        UV0[i] = F[1];
        UV1[i] = 0.5*(i*h*F[1] - F[0]);
    }
}
//======================================================================================================================


//======================================================================================================================
// Interpolation of function with values V[0], V[1], V[2], V[3], V[4] (at points 0, 1, 2, 3, 4) to point x
// using 4-th order Lagrange polynomial
//======================================================================================================================
static double FivePointInterpolation(const double* V, double x) {
    double sum = 0.0;
    for(int i=0; i<5; i++) {
        double m = 1.0;
        for(int j=0; j<5; j++) {
            if(j==i) continue;
            m *= (x - j) / (i - j);
        }
        sum += V[i] * m;
    }
    return sum;
}
//======================================================================================================================

//======================================================================================================================
// Blasius profile: obtaining the solution in the given point using interpolation
//======================================================================================================================
void s_Blasius::PointValue(double /*t*/, const double* coor, double* V) const{
    if(UV0.size()<5 || UV0.size()!=UV1.size()) crash("s_Blasius: Init not done");
    V[Var_R] = 1.0;
    V[Var_U] = V[Var_V] = V[Var_W] = 0.0;
    V[Var_P] = SoundVel*SoundVel / gam;

    double x = coor[0], y = fabs(coor[1]);
    if (y < 1e-50) return;

    const double eta = (x < 1e-10) ? huge : sqrt(Rey) * y / sqrt(x);
    const int i = int(eta / d_eta);
    if(i < 3) { // Use 2 terms of Taylor series
        double eta3 = eta*eta*eta;
        double f  = 0.5 * alpha*eta*eta * (1.0 - alpha*eta3 / 120.0);
        double df = alpha*eta* (1.0 - alpha*eta3 / 48.0);
        V[Var_U] = df;
        V[Var_V] = 0.5*(eta*df - f);
    }
    else if(i > int(UV0.size() - 4)) { // Free stream values
        V[Var_U] = 1.;
        V[Var_V] = UV1[UV0.size() - 1];
    }
    else { // Interpolation
        double dd = eta / d_eta - i; // fractional part, 0 <= dd <= 1
        // Interpolation over ((i-2) - (i-1) - (i) - (i+1) - (i+2)) and ((i-1) - (i) - (i+1) - (i+2) - (i+3))
        V[Var_U] = FivePointInterpolation(&(UV0[i-2]),2.0+dd)*(1.0-dd) + FivePointInterpolation(&(UV0[i-1]),1.0+dd)*dd;
        V[Var_V] = FivePointInterpolation(&(UV1[i-2]),2.0+dd)*(1.0-dd) + FivePointInterpolation(&(UV1[i-1]),1.0+dd)*dd;
    }

    V[Var_U] *= FlowVel;
    V[Var_V] *= FlowVel / sqrt(Rey * x);
    if(coor[1]<0) V[Var_V] = - V[Var_V];
}
//======================================================================================================================



