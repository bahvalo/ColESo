// *********************************************************************************************************************
// *****                                                                                                           *****
// *****                                  Collection of exact solutions (ColESo)                                   *****
// ****                                           Planar steady vortices                                           *****
// *****                                                                                                           *****
// *********************************************************************************************************************

#include "parser.h"
#include "coleso.h"
#include "geom_primitive.h"
#ifdef _NOISETTE
#include "lib_base.h"
#endif

//======================================================================================================================
// Считывание общих параметров для всех вихрей
//======================================================================================================================
void tVortexFunction::ReadParams(tFileBuffer& FB) {
    tParamManager PM;
    // Центр вихря при t=0
    PM.Request(r0[0], "X");
    PM.Request(r0[1], "Y");
    PM.Request(r0[2], "Z");
    // Фоновый поток
    PM.Request(v0[0], "vx");
    PM.Request(v0[1], "vy");
    PM.Request(v0[2], "vz");
    PM.Request(gam, "gam");
    PM.Request(SoundVel, "SoundVel");
    // Некоторым образом определяемое число Маха вихря
    PM.Request(Mach, "Mach");
    // Радиус вихря
    PM.Request(R, "R");
    // Ось вихря (0 - X, 1 - Y, 2 - Z)
    PM.Request(axis, "axis");
    PM.ReadParamsFromBuffer(FB);
}
//======================================================================================================================


//======================================================================================================================
// Проверка общих параметров
//======================================================================================================================
void tVortexFunction::Init() {
    if(IsNaN(Mach) || Mach < 0) crash("tVortexFunction::CheckParams: wrong Mach = %e", Mach);
    if(IsNaN(R) || R <= 0) crash("tVortexFunction::CheckParams: wrong R = %e", R);
    if(axis < 0 || axis > 2) crash("tVortexFunction::CheckParams: wrong axis %i", axis);
    if(IsNaN(r0[0])||IsNaN(r0[1])||IsNaN(r0[2])||IsNaN(v0[0])||IsNaN(v0[1])||IsNaN(v0[2]))
        crash("tVortexFunction::CheckParams: NaN detected");
}
//======================================================================================================================


//======================================================================================================================
// Общая процедура для вычисления решения
//======================================================================================================================
void tVortexFunction::PointValue(double t, const double* coord, double* V) const {
    int axis_x = (axis + 1) % 3;
    int axis_y = (axis + 2) % 3;
    double x = coord[axis_x]-r0[axis_x]-t*v0[axis_x];
    double y = coord[axis_y]-r0[axis_y]-t*v0[axis_y];
    RoundToCentre(x, Per[axis_x]);
    RoundToCentre(y, Per[axis_y]);

    double r = sqrt(x*x+y*y);
    double u_over_r;
    Profile(r, V[Var_R], u_over_r, V[Var_P]);
    V[Var_U+axis_x] = v0[axis_x] - y * u_over_r;
    V[Var_U+axis_y] = v0[axis_y] + x * u_over_r;
    V[Var_U+axis]   = v0[axis];
}
//======================================================================================================================


//======================================================================================================================
//  Вихрь с финитным профилем скорости
//======================================================================================================================

//======================================================================================================================
void s_FiniteVortex::ReadParams(tFileBuffer& FB) {
    tVortexFunction::ReadParams(FB);
    tParamManager PM;
    PM.Request(Per[0], "PerX");
    PM.Request(Per[1], "PerY");
    PM.Request(Per[2], "PerZ");
    PM.Request(mode, "mode");
    PM.Request(deg, "deg");
    PM.ReadParamsFromBuffer(FB);
}
//======================================================================================================================

//======================================================================================================================
void s_FiniteVortex::Profile(double r, double& rho, double& u_over_r, double& p) const {
    rho = 1.0; p = SoundVel*SoundVel / gam; u_over_r = 0.0;
    double a = R;
    if(mode) if(r > 3.0*a) return;
    if(!mode) if(r > 2.0*a) return;

    double C;
    if(mode) {
        C = Mach / (4.0*a*a*a);
        u_over_r = C*SQR(3*a-r);
    }
    else {
        C = Mach / pow(a, 2*deg);
        u_over_r = C*pow(r*(2*a-r), deg-1)*(2*a-r);
    }

    double f = C*C;
    if(mode) f*=pow(3*a-r, 5)*(C1_6*r + 0.1*a);
    else {
        // Beware: BinomialCoeff0(...) is based on factorials and does not work for big 'deg' values
        if(deg>=5 || deg<=0) crash("deg = %i not realized", deg);
        double sum = 0.0;
        for(int i=0; i<=2*deg; i++)
            sum += BinomialCoeff0(2*deg, i)*pow(2*a,i)*pow(-1.0,2*deg-i)*(pow(2*a,4*deg-i)-pow(r,4*deg-i))/(4*deg-i);
        f*=sum;
    }
    double pfunc = 1.0 - (gam-1.0)*f;
    if(pfunc < tiny) crash("s_FiniteVortex::PointValue error: P<0");
    p = pow(pfunc, gam/(gam-1.0)) / gam;
    rho = pow(gam*p, 1.0/gam);
}
//======================================================================================================================



//======================================================================================================================
// Vortex with u_phi(r) = Gamma * r / (2*Pi*(r^2 + R^2)), Gamma = 2*Pi*(SoundVel*Mach)*R
// It is assumed that SoundVel = 1
//======================================================================================================================

//======================================================================================================================
void s_Vortex_BG::Profile(double r, double& rho, double& u_over_r, double& p) const {
    double Gamma = Pi2*Mach*R;
    u_over_r = Gamma / (Pi2*(r*r + R*R));
    double g1g = (gam-1.0)/gam;
    double buf = pow(1.0/gam, g1g) - g1g*pow(gam, 1.0/gam)*Gamma*Gamma/(8.0*PiNumber*PiNumber * (r*r+R*R));
    if(buf<0) crash("s_Vortex_BG::PointValue internal error: buf = %e", buf);
    p = pow(buf, gam/(gam-1.0));
    rho = pow(gam*p, 1.0/gam);
}
//======================================================================================================================



//======================================================================================================================
//  Rankine vortex
//======================================================================================================================

//======================================================================================================================
void s_RankineVortex::Profile(double r, double& rho, double& u_over_r, double& p) const {
    double J;
    if(r < R) {
        J = Mach*Mach*(1 - r*r/(2.*R*R));
        u_over_r = Mach/R;
    } else {
        J = Mach*Mach*R*R / (2*r*r);
        u_over_r = Mach*R/(r*r);
    }
    double buf = 1 - (gam-1)*J / (SoundVel*SoundVel);
    rho = pow(buf,1/(gam-1));
    p = pow(rho,gam) * SoundVel*SoundVel/gam;
}
//======================================================================================================================


//======================================================================================================================
//  Vortex with Gaussian profile for circulation
//======================================================================================================================

const double s_GaussianVortex::alpha0 = 1.2564312086261696770;
const double s_GaussianVortex::Mach2A = 1.3979525473159165448;

//======================================================================================================================
void s_GaussianVortex::Profile(double r, double& rho, double& u_over_r, double& p) const {
    double rr = r*r;
    const double inv_gam = 1.0/gam;
    const double alpha = alpha0 / (R*R);
    double arr = alpha * rr;
    double A = Mach * Mach2A * R;

    u_over_r = (rr > tiny*tiny) ? (1.0 - exp(-arr)) / rr : alpha; // V_theta(r) / r
    u_over_r *= A;

    double J = 0.0;
    if(rr > tiny*tiny) J = alpha*(sf13d_c(2.*arr, 2) - sf13d_c(arr, 2)) - SQR(1. - exp(-arr)) / (2.*rr);
    else J = -alpha*CLN2;
    J *= A*A * (gam - 1.0)/gam * pow(gam / SQR(SoundVel), inv_gam);
    J += pow(gam/SQR(SoundVel), -(gam-1.)/gam);

    p = pow(J, gam/(gam-1.));
    rho = pow(gam/SQR(SoundVel) * p, 1/gam);
}
//======================================================================================================================
