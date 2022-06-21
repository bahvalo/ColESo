// *********************************************************************************************************************
// *****                                                                                                           *****
// *****                                  Collection of exact solutions (ColESo)                                   *****
// *****                                           Miscellaneous cases:                                            *****
// *****                  shock wave by wall reflection, shock wave profile, smooth simple wave,                   *****
// *****                  flow around a sphere with Re << 1, acoustics between coaxial cylinders,                   *****
// *****                        Couette flow between parallel plates and coaxial cylinders                         *****
// *****                                                                                                           *****
// *********************************************************************************************************************

#include "base_parser.h"
#include "coleso.h"
#include "geom_primitive.h"
#include "es_specfunc.h"

//======================================================================================================================
void s_ShockRefl::ReadParams(tFileBuffer& FB) {
    tParamManager PM;
    PM.Request(gam, "gam"); // Specific ratio
    PM.Request(WallPosition, "WallPosition"); // Coordinate of the wall
    PM.Request(MUL[0], "RhoL"); // Values left to the shock
    PM.Request(MUL[1], "UL");   // Values left to the shock
    PM.Request(MUL[2], "PL");   // Values left to the shock
    PM.Request(MUR[0], "RhoR"); // Values right to the shock
    //PM.Request(MUR[1], "UR"); // MUR[1] = 0.0
    PM.Request(MUR[2], "PR");   // Values right to the shock
    PM.ReadParamsFromBuffer(FB);
}
//======================================================================================================================

//======================================================================================================================
void s_ShockRefl::Init() {
    MUR[Var_U] = 0.0;
    if(IsNaN(MUL[0]) || IsNaN(MUL[1]) || IsNaN(MUL[2]) || IsNaN(MUR[0]) || IsNaN(MUR[2]))
        crash("s_ShockRefl error: initial data is not set");
    s_Riemann SR(gam, MUL, MUR);
    if(!SR.IsShock(NULL)) crash("s_ShockRefl error: initial data is not a shock wave");
}
//======================================================================================================================

//======================================================================================================================
void s_ShockRefl::PointValue(double t, const double* coor, double* V) const {
    if(t < tiny) {
        const double* VV = coor[0]>DF_Coor ? MUR : MUL;
        V[0]=VV[0]; V[1]=VV[1]; V[2]=V[3]=0.0; V[4]=VV[2];
        return;
    }
    double rho1, u1, p1, rho2, u2, p2, rho3, u3, p3;
    rho1 = MUR[Var_R]; u1 = MUR[Var_U]; p1 = MUR[2];
    rho2 = MUL[Var_R]; u2 = MUL[Var_U]; p2 = MUL[2];
    double ShockVelInc = (rho2*u2-rho1*u1)/(rho2-rho1);    // Скорость ударной волны, падающей на стенку
    double Trefl = (WallPosition - DF_Coor) / ShockVelInc; // Время отражения

// Параметры за фронтом отражённой волны (Станюкович К.П. Неустановившиеся движения сплошной среды, страницы 263-265)
    u3 = 0.0;
    p3 = p2 + 2.0*gam*p2*(p2-p1)/((gam-1.0)*p2+(gam+1.0)*p1);
    rho3 = ((gam-1.0)/gam+p1/(gam*p2))/rho2;
    rho3 = 1.0/rho3;
    double ShockVelRel = (rho3*u3-rho2*u2)/(rho3-rho2);   // Скорость ударной волны, отражённой от стенки

    double MURR[3];
    const double* VV = coor[0]>DF_Coor+t*ShockVelInc ? MUR : MUL;
    if(t > Trefl) {
        MURR[0] = rho3; MURR[1] = u3; MURR[2] = p3;
        VV = coor[0]>WallPosition+(t-Trefl)*ShockVelRel ? MURR : MUL;
    }
    V[0]=VV[0]; V[1]=VV[1]; V[2]=V[3]=0.0; V[4]=VV[2];
}
//======================================================================================================================



//======================================================================================================================
// Simple wave (modified case from Keldysh Institute preprint 2013-53 by Ladonkina et al.)
//======================================================================================================================
void s_SimpleWave::ReadParams(tFileBuffer& FB) {
    double background_flow = FlowVel * sign_FlowVel;
    tParamManager PM;
    PM.Request(x0, "x0");             // perturbation center
    PM.Request(l, "l");               // perturbation width
    PM.Request(a, "a");               // perturbation strength
    PM.Request(gam, "gam");           // specific ratio
    PM.Request(SoundVel, "SoundVel"); // sound speed in background media
    PM.Request(background_flow, "FlowVel"); // background flow; 0.0 by default; -sqrt(10) in preprint
    PM.ReadParamsFromBuffer(FB);

    FlowVel = fabs(background_flow); // FlowVel can be accessed elsewhere and must be non-negative
    sign_FlowVel = SIGN(background_flow);
}
//======================================================================================================================

//======================================================================================================================
void s_SimpleWave::Init() {
    if(gam<=1.0 || l<=0.0 || SoundVel<=0.0) crash("s_SimpleWave: wrong parameters");
    if(FlowVel<0.0) { FlowVel=-FlowVel; sign_FlowVel=-sign_FlowVel; }
    shock_formation_time = l / (a*SoundVel) * (gam-1.0)/(gam+1.0) * exp(1.0)/sqrt(27.0);
    if(log) pprintf("s_SimpleWave: shock formation time = %25.15f\n", shock_formation_time);
}
//======================================================================================================================

//======================================================================================================================
void s_SimpleWave::PointValue(double t, const double* coor, double* V) const {
    const double background_flow = FlowVel * sign_FlowVel; // background flow = reference frame flow
    const double u0 = 2.*SoundVel/(gam-1.0); // velocity scale
    const double p0 = SoundVel*SoundVel/gam; // background pressure
    const double lambdamin = - u0 * gam;
    const double lambdamax = 0.0;
    double x = coor[0] - (x0 + background_flow*t);

    double rho = 1.0, p = p0, u = 0.0;
    double xl = -l+t*lambdamin, xr = l+t*lambdamax;
    if(x > xl && x < xr) {
        for(int k=0; k<50; k++) {
            double xc = 0.5*(xl+xr);
            double buf = 1.0;
            if(fabs(xc)<l) buf += a*exp(- 2.0*xc*xc/(l*l - xc*xc));
            rho = pow(buf, 2.0/(gam-1.0)); // buf = rho^((gam-1)/2)
            p = rho*buf*buf*p0; // p = p0*rho^gam
            u = u0 * (1.0 - buf);
            double xtilde = xc + t*(u - sqrt(gam*p/rho));
            if(xtilde < x) xl=xc;
            else xr = xc;
        }
        if(fabs(xr - xl) > 1e-10) crash("s_SimpleWave::PointValue: process diverged");
    }

    V[Var_R] = rho;
    V[Var_P] = p;
    V[Var_U] = u + background_flow;
    V[Var_V] = V[Var_W] = 0.0;
}
//======================================================================================================================


//======================================================================================================================
void s_ViscShock::ReadParams(tFileBuffer& FB) {
    tParamManager PM;
    PM.Request(MUL[0], "RhoL"); // Physical values at -inf
    PM.Request(MUL[1], "UL");
    PM.Request(MUL[2], "PL");
    PM.Request(MUR[0], "RhoR"); // Physical values at +inf
    PM.Request(MUR[1], "UR");
    PM.Request(MUR[2], "PR");
    PM.Request(gam, "gam");     // Specific ratio
    PM.Request(mode, "mode");   // mode=0: constant kinematic viscosity; mode=1: constant dynamic viscosity
    PM.Request(Rey, "Rey");     // "Reynolds number": nu=1/Rey if mode==0 and mu=1/Rey if mode==1
    PM.Request(Pr, "Prandtl");  // Prandtl number: 1e50 or 0.75
    PM.RequestOption(SteadyTest, "SteadyTest"); // Check that shock velocity is zero

    PM.Request(x0, "x0");       // Shock position
    #ifdef _NOISETTE
        PM.RequestOption(AutodetectShockPosition, "AutodetectShockPosition");
    #endif
    PM.ReadParamsFromBuffer(FB);
}
//======================================================================================================================

//======================================================================================================================
// Check that the data as x->infty prescribes a shock and determine shock speed
//======================================================================================================================
void s_ViscShock::Init() {
    s_Riemann SR(gam, MUL, MUR);
    if(!SR.IsShock(&ushock))
        crash("s_ViscShock error: data at infinity do not prescribe a shock");
    if(MUL[1]<MUR[1]) crash("s_ViscShock error: entropy inequality fail");
    if(SteadyTest) {
        if(fabs(ushock)>1e-13)
            crash("s_ViscShock error: ushock = %e but steady test is activated", ushock);
        ushock = 0.0;
    }
    pprintf("s_ViscShock: shock speed = %e\n", ushock);

    if(mode!=0 && mode!=1) crash("s_ViscShock error: unknown mode %i", mode);

    if(Pr > 5e49) {} // Zero heat conductivity mode
    else if(fabs(Pr - 0.75) < 1e-13) {} // Pr=0.75 mode
    else crash("s_ViscShock error: not implemented for Pr=%e. Should be 0.75 or 1e50", Pr);
}
//======================================================================================================================

//======================================================================================================================
void s_ViscShock::PointValue(double t, const double* coor, double* V) const {
    V[Var_U] = V[Var_V] = V[Var_W] = 0.0;
    double x = coor[Coor_X] - (x0 + ushock*t);
    double ul = MUL[1] - ushock, ur = MUR[1] - ushock;
    // ul and ur have same signs; ul > ur

    double xi = x;
    if(mode==1) {
        // if ur<0 the we use the transformation x <--> -x, 'l' <--> 'r', u <--> -u
        const double _x = ur>0 ? x : -x;
        const double mu_over_A = ur>0 ? 1.0/(Rey*MUR[0]*ur) : -1.0/(Rey*MUL[0]*ul);
        const double _ur = ur>0 ? ur : -ul;

        double a = 0.375*(gam+1)*fabs(ul-ur);
        if(Pr<1) a /= gam; // Pr = 0.75
        double xiL = MUL[0]*Rey*_x;
        double xiR = MUR[0]*Rey*_x;
        if(xiR < xiL) SWAP(xiL, xiR);

        double xC;
        for(int iter=0; iter<50; iter++) {
            double xiC = 0.5*(xiL+xiR);

            double aux = 0.0;
            if(xiC>0) aux = log(0.5*(1 + exp(-a*xiC)));
            else aux = log(0.5*(1 + exp(a*xiC))) - a*xiC; // avoid floating point overflow if a*xiC>>1
            xC = mu_over_A*(_ur*xiC - fabs(ul-ur)/a * aux);
            if(xC<_x) xiL = xiC;
            else xiR = xiC;
        }
        xi = 0.5*(xiL + xiR);
        if(ur<0) xi = -xi;
        //printf0("xi = %e xC = %e x=%e bnd=%e %e\n", xi, xC, x, MUL[0]*Rey*x, MUR[0]*Rey*x);
        xi /= Rey;
    }

    {
        double a = 0.375*Rey*(gam+1)*fabs(ul-ur);
        if(Pr<1) a /= gam; // Pr = 0.75
        if(xi>0.0) {
            double e = exp(-a*xi);
            V[Var_U] = (ul*e + ur) / (1 + e);
        }
        else {
            double e = exp(a*xi);
            V[Var_U] = (ul + ur*e) / (1 + e);
        }
        V[Var_R] = MUL[0]*ul / V[Var_U];
        if(Pr<1) {
            V[Var_P] = (MUL[2]/MUL[0] + (gam-1.)/(2*gam)*(ul*ul - V[Var_U]*V[Var_U])) * V[Var_R];
        }
        else {
            double C = 0.0;
            if(fabs(a*xi) < 300) C = (MUL[1]-MUR[1])/(exp(0.5*a*xi)+exp(-0.5*a*xi));
            V[Var_P] = MUL[2] + MUL[0]*ul*(ul-V[Var_U]) - V[Var_R]*0.5*(gam+1.0)*C*C;
        }
    }

/* // closed-form expression for Pr=0.75, rhoR/rhoL=2
    if(mode==1) {
        double V_check[5] = {0.,0.,0.,0.,0.};
        if(MUL[0] < MUR[0]) { x=-x; ul=-ur; }
        x *= 0.375*Rey*(gam+1)/gam*MUL[0]*ul;
        if(x>0.0) {
            double se = sqrt(1.0 + 4.*exp(-x));
            V_check[Var_U] = se / (1. + se) * 2;
        }
        else {
            double e = exp(x);
            V_check[Var_U] = 2. + 0.5*e - 0.5*sqrt(e*e + 4*e);
        }
        V_check[Var_U] *= ul;
        V_check[Var_R] = MUL[0]*ul / V_check[Var_U];
        V_check[Var_P] = (MUL[2]/MUL[0] + (gam-1.)/(2*gam)*(ul*ul - V_check[Var_U]*V_check[Var_U])) * V_check[Var_R];
        if(MUL[0] < MUR[0]) {
            V_check[Var_U] *= -1;
            V_check[Var_R] = MUR[0]*ur / V_check[Var_U];
            V_check[Var_P] = (MUR[2]/MUR[0] + (gam-1.)/(2*gam)*(ur*ur - V_check[Var_U]*V_check[Var_U])) * V_check[Var_R];
        }
        printf0("V = %e %e %e; V_check = %e %e %e\n", V[0], V[1], V[4], V_check[0], V_check[1], V_check[4]);
        for(int ivar=0; ivar<5; ivar++) V_check[ivar] -= V[ivar];
        printf0("err = %e %e %e\n", V_check[0], V_check[1], V_check[4]);
    }
*/

    V[Var_U] += ushock;
}
//======================================================================================================================


//======================================================================================================================
void s_ViscSphere::ReadParams(tFileBuffer& FB) {
    tParamManager PM;
    PM.Request(R, "R"); // Sphere radius
    PM.Request(Rey, "Rey"); // Reynolds number
    PM.Request(gam, "gam"); // Specific ratio
    PM.Request(SoundVel, "SoundVel"); // Sound speed in the background flow
    PM.Request(FlowVel, "FlowVel"); // Flow velocity in the background flow
    PM.ReadParamsFromBuffer(FB);
}
//======================================================================================================================

//======================================================================================================================
void s_ViscSphere::PointValue(double, const double* coord, double* V) const {
    V[Var_R] = 1.0;
    V[Var_U] = V[Var_V] = V[Var_W]=0.0;
    V[Var_P] = SQR(SoundVel) / gam;

    double r = sqrt(VDOT(coord, coord));
    if(r<R) r = R;
    double Rr = R/r;

    double u[3] = {1.0, 0.0, 0.0};
    double n[3];
    for(int icoor=0; icoor<3; icoor++) n[icoor] = coord[icoor] / r;
    double un = n[0];

    for(int icoor=0; icoor<3; icoor++) {
        V[icoor+Var_U] = - 0.75*Rr*(u[icoor] + n[icoor]*un) - (Rr*Rr*Rr)*(0.25*u[icoor] - 0.75*n[icoor]*un) + u[icoor];
        V[icoor+Var_U] *= FlowVel;
    }
    V[Var_P] -= 1.5*R/(r*r) * un / Rey * FlowVel*FlowVel;
    if(V[Var_P] <= 0) crash("s_ViscSphere::PointValue error: P=%e\n  (coords = %e %e %e)\n",
        V[Var_P], coord[0], coord[1], coord[2]);
}
//======================================================================================================================


//======================================================================================================================
void s_Coaxial::ReadParams(tFileBuffer& FB) {
    tParamManager PM;
    PM.Request(rmin, "rmin");                   // cylinder radius (may be zero)
    PM.Request(rmax, "rmax");                   // cylinder radius
    PM.Request(AzimuthalMode, "AzimuthalMode"); // azimuthal mode (integer)
    PM.Request(RadialMode, "RadialMode");       // radial mode (integer)
    PM.Request(kz, "kz");                       // axial wave number (double)
    PM.Request(FlowVel, "FlowVel");             // axial flow velocity
    PM.Request(SoundVel, "SoundVel");           // sound speed
    PM.Request(Ampl, "Ampl");                   // amplitude (multiplicator)
    PM.Request(phase, "phase");                 // phase
    PM.Request(CoorAxis, "CoorAxis");           // axis: 0 - X, 1 - Y, 2 - Z
    PM.ReadParamsFromBuffer(FB);
}
//======================================================================================================================

//======================================================================================================================
// Calculating root for kr of equation using Newton method
// Uses initial value for kr as the initial value for the iterative process
// Returns 1 if unconvergence occurs
//======================================================================================================================
int s_Coaxial::FindRoot(int nu, double rmin, double rmax, double& kr) {
    NativeDouble* buf = new NativeDouble[(nu+3)*4];
    NativeDouble* bufN = buf + (nu+3)*2;

    int success = 0;
    const int MAX_ITERS = 100;
    for(int iter = 0; iter<MAX_ITERS; iter++) {
        double x = kr * rmin;
        BesselFunctionsJ(x, nu+2, buf);
        if(BesselFunctionsY(x, nu+2, bufN)!=nu+2) break; // error
        double JL = buf[nu], NL = buf[(nu+3)*2+nu];
        double JpL = - buf[nu+1]  + nu/x*JL;
        double NpL = - bufN[nu+1] + nu/x*NL;

        double y = kr * rmax;
        BesselFunctionsJ(y, nu+2, buf);
        if(BesselFunctionsY(y, nu+2, bufN)!=nu+2) break; // error
        double JR = buf[nu], NR = buf[(nu+3)*2+nu];
        double JpR = - buf[nu+1]  + nu/y*JR;
        double NpR = - bufN[nu+1] + nu/y*NR;

        // function
        double f  = JpL * NpR - JpR * NpL;
        // derivative
        double fk = - f / kr;
        fk -= (rmin*rmin - nu*nu/(kr*kr)) * (JL * NpR - JpR * NL) / rmin;
        fk -= (rmax*rmax - nu*nu/(kr*kr)) * (JpL * NR - JR * NpL) / rmax;
        // offset
        if(fabs(fk) < 1e-50) break; // infinite increment. No convergence
        double df = - f / fk;
        kr += df;
        if(kr < 0.0) break; // got negative value. No convergence

        // We stop iteration process after 2 successive successful iterations
        if(fabs(df) < 1e-10*(1 + fabs(NpR) + fabs(NpL))) success++;
        else success = 0;
        if(success==2) break; // done!
    }
    delete[] buf;

    return (success!=2);
}
//======================================================================================================================


//======================================================================================================================
void s_Coaxial::Init() {
    if(AzimuthalMode < 0 || AzimuthalMode > 50) crash("s_Coaxial: wrong azimuthal mode %i", AzimuthalMode);

    // Будем искать корень, постепенно увеличивая rmin
    // При rmin=0 определим корень из решения задачи об акустике в цилиндре
    const double mc = 0.025; // магическая константа
    double r = 0.0; // внутренний радиус
    kr = BesselPrimeZero(AzimuthalMode, RadialMode, 0) / rmax;
    int fin = (rmin < tiny);
    while(!fin) {
        r += (rmax-r) * mc;
        if(r >= rmin) { r = rmin; fin = 1; }
        if (FindRoot(AzimuthalMode, r, rmax, kr))
            crash("s_Coaxial: iteration process for 'kr' diverged");
    }
#ifdef _NOISETTE
    pprintf("s_Coaxial: kr = %25.15f\n", kr);
#endif
}
//======================================================================================================================

//======================================================================================================================
// Attention! The function calculates V[0..9] (real and imaginary parts)
//======================================================================================================================
void s_Coaxial::PointValue(double t, const double* coor, double* V) const {
    if(IsNaN(kr)) crash("s_Coaxial::PointValue: Init not done");
    const int CoorX = (CoorAxis + 1)%3, CoorY = (CoorAxis + 2)%3; // axes in a plane of the circle
    const double x = coor[CoorX], y = coor[CoorY], z = coor[CoorAxis];

    double r = sqrt(x*x+y*y);
    if(r<rmin) r = rmin;
    double phi = GetAngle(x,y);

    const int NumN = AzimuthalMode+1;
    NativeDouble *Phi = new NativeDouble[NumN*2];
    NativeDouble *dPhi = Phi + NumN;

    BesselOuterFunction(NativeDouble(kr*r), NativeDouble(kr*rmin), NumN, Phi, dPhi);
    double omega = SoundVel * sqrt(kr*kr + kz*kz);
    double _phase = omega*t - kz*(z-FlowVel*t) - AzimuthalMode*phi + phase;
    double c = cos(_phase), s = -sin(_phase);

    V[  Var_P] = c * Phi[AzimuthalMode];
    V[5+Var_P] = s * Phi[AzimuthalMode];
    V[  Var_R] = V[  Var_P] / (SoundVel * SoundVel);
    V[5+Var_R] = V[5+Var_P] / (SoundVel * SoundVel);
    V[  Var_U+CoorAxis] = c * kz * Phi[AzimuthalMode] / omega;
    V[5+Var_U+CoorAxis] = s * kz * Phi[AzimuthalMode] / omega;

    double Re_ur =   s * dPhi[AzimuthalMode] * kr / omega;
    double Im_ur = - c * dPhi[AzimuthalMode] * kr / omega;
    if(r > tiny) {
        double Re_uphi = AzimuthalMode * c * Phi[AzimuthalMode] / (omega*r);
        double Im_uphi = AzimuthalMode * s * Phi[AzimuthalMode] / (omega*r);
        c = cos(phi); s = sin(phi);
        V[  Var_U+CoorX] = Re_ur * c - Re_uphi * s;
        V[  Var_U+CoorY] = Re_ur * s + Re_uphi * c;
        V[5+Var_U+CoorX] = Im_ur * c - Im_uphi * s;
        V[5+Var_U+CoorY] = Im_ur * s + Im_uphi * c;
    }
    else { // regularizing 0/0
        V[Var_U+CoorX] = V[Var_U+CoorY] = V[5+Var_U+CoorX] = V[5+Var_U+CoorY] = 0.0;
        if(AzimuthalMode==1) {
            double mult = 0.5*kr/omega;
            V[  Var_U+CoorX] =  s * mult;
            V[5+Var_U+CoorX] = -c * mult;
            V[  Var_U+CoorY] =  c * mult;
            V[5+Var_U+CoorY] =  s * mult;
        }
    }
    delete[] Phi;
    for(int i=0; i<10; i++) V[i] *= Ampl;
}
//======================================================================================================================


//======================================================================================================================
void s_ConcCyl::ReadParams(tFileBuffer& FB) {
    tParamManager PM;
    PM.Request(CoorAxis, "CoorAxis");        // axis direction
    PM.Request(X[0], "X");                   // a point on the axis
    PM.Request(X[1], "Y");                   // a point on the axis
    PM.Request(X[2], "Z");                   // a point on the axis
    PM.Request(RIn, "RIn");                  // radius of the inner cylinder
    PM.Request(ROut, "ROut");                // radius of the outer cylinder
    PM.Request(OmegaIn, "OmegaIn");          // angular velocity of the inner cylinder
    PM.Request(OmegaOut, "OmegaOut");        // angular velocity of the outer cylinder
    PM.Request(Pr, "Pr");                    // Prandtl number
    PM.Request(gam, "gam");                  // specific ratio
    PM.Request(Twall, "Twall");              // temperature of the outer cylinder
    PM.ReadParamsFromBuffer(FB);
}
//======================================================================================================================

//======================================================================================================================
void s_ConcCyl::Init() {
    if(RIn<=0.0 || ROut<=RIn+1e-10 || Pr<=0.0 || gam<=1.0 || CoorAxis<0 || CoorAxis>2)
        crash("s_ConcCyl: wrong parameters");
}
//======================================================================================================================

//======================================================================================================================
// ATTENTION! Returns the solution for velocity and temperature
//======================================================================================================================
void s_ConcCyl::PointValue(double, const double* coor, double* V) const {
    int CoorX = (CoorAxis + 1)%3, CoorY = (CoorAxis + 2)%3; // оси в плоскости круга
    double x = coor[CoorX]-X[CoorX], y = coor[CoorY]-X[CoorY];

    double invdet = 1.0/(RIn/ROut - ROut/RIn);
    double A = (RIn*OmegaIn/ROut - ROut*OmegaOut/RIn) * invdet;
    double B = RIn*ROut*(OmegaOut-OmegaIn) * invdet;
    double r = sqrt(x*x+y*y);
    if(r<tiny) { V[0]=V[1]=V[2]=V[3]=V[4]=0.0; return; }
    double uphi = A*r + B/r;
    double Temp = Twall + (gam-1.0)/gam*Pr*B*B*(1.0/SQR(ROut) - 1.0/SQR(r) - 2.0/SQR(RIn)*log(r/ROut));
    V[Var_R] = 1.0;
    V[Var_U+CoorX] = - uphi * y/r;
    V[Var_U+CoorY] = uphi * x/r;
    V[Var_U+CoorAxis] = 0.0;
    V[Var_P] = Temp;
}
//======================================================================================================================



//======================================================================================================================
void s_Couette::ReadParams(tFileBuffer& FB) {
    tParamManager PM;
    PM.Request(xL, "xL"); // left wall position
    PM.Request(xR, "xR"); // right wall positions
    PM.Request(vL, "vL"); // velocity (v_y) of the left wall
    PM.Request(vR, "vR"); // velocity (v_y) of the right wall
    PM.Request(condL, "condL"); // type of the condition on the left wall: 0 - adiabatic, 1 - isothermal
    PM.Request(condR, "condR"); // type of the condition on the right wall: 0 - adiabatic, 1 - isothermal
    PM.Request(tL, "tL"); // temperature on the left wall (if condL=1)
    PM.Request(tR, "tR"); // temperature on the right wall (if condR=1)
    PM.Request(ViscType, "ViscType"); // Viscosity coefficient
    // ViscType values: 0: mu=const; 1: mu=mu0/sqrt(T); 2: mu=mu0*T; 3: mu=mu0*sqrt(T); 4: Sutherland law
    PM.Request(pressure, "pressure"); // Pressure (which is constant in space)
    PM.Request(gam, "gam"); // specific ratio
    PM.Request(Pr, "Pr"); // Prandtl number
    PM.Request(TSutherland, "TSutherland"); // Constant in the Sutherland law
    PM.ReadParamsFromBuffer(FB);
}
//======================================================================================================================


//======================================================================================================================
void s_Couette::Init() {
    if(xR-xL<tiny) crash("s_Couette: error: xR <= xL");
    if(ViscType<0 || ViscType>4) crash("s_Couette: wrong ViscType %i", ViscType);
    if(condL!=0 && condL!=1) crash("s_Couette: wrong condL %i", condL);
    if(condR!=0 && condR!=1) crash("s_Couette: wrong condR %i", condR);
    if((condL && tL <= 0.0) || (condR && tR <= 0.0) || (gam<=1.0) || (Pr<=0.0)) crash("s_Couette: wrong parameters");
    if(condL==0 && condR==0) crash("s_Couette: adiabatic BC on both plates set");
    if(ViscType==1) {
        if(!(condL==1 && condR==0)) crash("s_Couette: ViscType==1 is proper only for condL==1 and condR==0");
        double kappa = fabs(vR-vL)/sqrt(tL) * sqrt((gam-1.0)*Pr/(2.0*gam));
        pprintf("s_Couette::Init: kappa = %e, tR = %e\n", kappa, tL*(1.0+kappa*kappa));
    }
    if(ViscType==4) if(IsNaN(TSutherland) || TSutherland<=0.0) crash("s_Couette: wrong TSutherland");

    FlowVel = 0.0;
    //double Temp = condL&&condR ? 0.5*(tL+tR) : (condL ? tL : tR); // характерная температура
    //SoundVel = sqrt(gam*Temp); // do not override, otherwise Sutherland constant will be wrong

    #ifdef _NOISETTE
        if(IsNaN(pressure)) { AutodetectPressure = 1; pressure = Temp; }
    #else
        if(IsNaN(pressure) || pressure <= 0.0) crash("s_ConcCyl: wrong pressure %e", pressure);
    #endif
}
//======================================================================================================================

//======================================================================================================================
void s_Couette::PointValue(double, const double* coor, double* V) const {
    double x = (coor[0]-xL)/(xR-xL);
    double Vy=0, T;
    if(ViscType==0) { // mu=const
        Vy = vL*(1.0-x) + vR*x;
        double Txx = -(gam-1.0)/gam*Pr*SQR(vR-vL);

        double b, c=tL;
        if(condL==1 && condR==1) b = tR-tL-0.5*Txx; // изотермические ГУ с обеих сторон
        else if(condL==1) b = -Txx; // изотермическое ГУ слева, адиабатическое справа
        else if(condR==1) { b = 0.0; c = tR-0.5*Txx; } // адиабатическое ГУ слева, изотермическое справа
        else crash("s_Couette::PointValue error: wrong configuration");

        T = 0.5*Txx*x*x + b*x + c;
    }
    else if(ViscType==1) { // mu=const/sqrt(T)
        if(!(condL==1 && condR==0)) crash("s_Couette::PointValue error: wrong configuration");
        double vdash0 = sqrt(2.*tL*gam/(Pr*(gam-1)));
        double c = atan((vR-vL) / vdash0);
        T = cos(c*(1.0-x))/cos(c);
        T = tL*T*T;
        Vy = vR - vdash0 * sin(c * (1.0-x)) / cos(c);
    }
    else {
        const double varkappa = (gam-1.0)*Pr/gam;
        double c2, c3;
        if(condL==1 && condR==1) { // isothermal b.c. on both walls
            double TL = tL + 0.5*varkappa*vL*vL;
            double TR = tR + 0.5*varkappa*vR*vR;
            // solving system c2*vL+c3=TL, c2*vR+c3=TR
            c2 = (TR - TL) / (vR - vL);
            c3 = TL - c2*vL;
        }
        else if(condL==1) { // isothermal b.c. on the left wall
            c2 = varkappa*vR;
            c3 = tL + 0.5*varkappa*vL*vL - c2*vL;
        }
        else if(condR==1) { // isothermal b.c. on the right wall
            c2 = varkappa*vL;
            c3 = tR + 0.5*varkappa*vR*vR - c2*vR;
        }
        else crash("s_Couette::PointValue error: wrong configuration");

        if(ViscType==2) { // mu=const*T
            double c1, c4; // we put mu0 = 1
            {
                double VL = -C1_6*varkappa*vL*vL*vL + 0.5*c2*vL*vL + c3*vL;
                double VR = -C1_6*varkappa*vR*vR*vR + 0.5*c2*vR*vR + c3*vR;
                // solving system c1*xL+c4=VL, c1*xR+c4=VR. We scale x such that xL=0, xR=1
                c4 = VL;
                c1 = VR-c4;
            }

            // It remains to solve the equation
            double R = c1*x+c4;
            double vl=MIN(vL,vR), vr=MAX(vL,vR);
            for(int i=0; i<50; i++) {
                double vc=0.5*(vl+vr);
                double fc=-C1_6*varkappa*vc*vc*vc + 0.5*c2*vc*vc + c3*vc;
                if(fc<R) { vl=vc; }
                else { vr=vc; }
            }
            Vy = 0.5*(vl+vr);
        }
        else if(ViscType==3) { // mu=const*sqrt(T)
            const double v0 = c2 / varkappa;
            const double W = sqrt(v0*v0 + 2.*c3/varkappa);
            const double ul = (vL-v0)/W, ur = (vR-v0)/W;
            const double Il = asin(ul) + ul*sqrt(1-ul*ul), Ir = asin(ur) + ur*sqrt(1-ur*ur);

            // It remains to solve the equation
            double R = (Ir - Il)*x;
            double vl=MIN(vL,vR), vr=MAX(vL,vR);
            for(int i=0; i<50; i++) {
                double vc=0.5*(vl+vr);
                double uc=(vc-v0)/W;
                double fc=asin(uc) + uc*sqrt(1-uc*uc) - Il;
                if(fc<R) { vl=vc; }
                else { vr=vc; }
            }
            Vy = 0.5*(vl+vr);
        }
        else if(ViscType==4) { // mu=const * T^{3/2} / (T+C)
            const double v0 = c2 / varkappa;
            const double W = sqrt(v0*v0 + 2.*c3/varkappa);
            const double c = W*sqrt(0.5*varkappa);
            const double C = TSutherland;
            const double alpha = sqrt(C/(C+c*c));
            #define II(X) (0.5*c*W*(asin(X) + X*sqrt(1-X*X)) - C*W/c*(asin(X) - alpha*atan(alpha*X/sqrt(1-X*X))))
            const double ul = (vL-v0)/W, ur = (vR-v0)/W;
            const double Il = II(ul);
            const double Ir = II(ur);

            // It remains to solve the equation
            double R = (Ir - Il)*x;
            double vl=MIN(vL,vR), vr=MAX(vL,vR);
            for(int i=0; i<50; i++) {
                double vc=0.5*(vl+vr);
                double uc=(vc-v0)/W;
                double fc=II(uc) - Il;
                if(fc<R) { vl=vc; }
                else { vr=vc; }
            }
            Vy = 0.5*(vl+vr);
            #undef II
        }
        else crash("Internal error");
        T = -0.5*varkappa*Vy*Vy + c2*Vy + c3;
    }

    V[Var_U] = V[Var_W] = 0.0;
    V[Var_V] = Vy;
    V[Var_P] = pressure;
    V[Var_R] = V[Var_P]/T;
}
//======================================================================================================================
