// *********************************************************************************************************************
// *****                                                                                                           *****
// *****                                  Collection of exact solutions (ColESo)                                   *****
// *****                          Exact solution for Riemann problem and Godunov solver                            *****
// *****                                                                                                           *****
// *********************************************************************************************************************

#include "base_parser.h"
#include "coleso.h"

//======================================================================================================================
enum tRiemannProblemType{
    RIEMANN_CUSTOM     =  0,  // general case
    RIEMANN_SOD        =  1,  // Sod shock tube
    RIEMANN_LAX        =  2,  // Lax problem
    RIEMANN_ZHMF       =  3   // Zhmakin-Fursenko problem
};

//======================================================================================================================
int s_Riemann::IsUnsteady() const {
    // check that the fields are a steady solution of Euler equations
    // do not check the entropy inequality
    double tinyloc = 1e-13;
    if(fabs(MUL[0]*MUL[1]-MUR[0]*MUR[1])>tinyloc) return 1; // [rho*u] != 0
    if(fabs(MUL[0]*MUL[1]*MUL[1]+MUL[2]-MUR[0]*MUR[1]*MUR[1]-MUR[2])>tinyloc) return 2; // [rho*u*u+p] != 0
    if(fabs((0.5*MUL[0]*MUL[1]*MUL[1]+gam/(gam-1)*MUL[2])*MUL[1] -
        (0.5*MUR[0]*MUR[1]*MUR[1]+gam/(gam-1)*MUR[2])*MUR[1])>tinyloc) return 3; // [(E+p)*u] != 0
    return 0;
}
//======================================================================================================================

//======================================================================================================================
int s_Riemann::IsShock(double* ushock) const {
    double tinyloc = 1e-13;
    // if the densities are equal check that other variables are also equal. Uniform field can be considered as a shock
    if(fabs(MUL[0] - MUR[0]) < tinyloc) {
        if(fabs(MUL[1] - MUR[1]) > tinyloc || fabs(MUL[2] - MUR[2]) > tinyloc) return 0;
        if(ushock!=NULL) *ushock = 0.5*(MUL[1] + MUR[1]);
        return 1;
    }
    // otherwise determine a shock velocity from the continuity equation
    double us = (MUL[0]*MUL[1] - MUR[0]*MUR[1]) / (MUL[0] - MUR[0]);
    if(ushock!=NULL) *ushock = us;
    // check Rankine--Hugoniot conditions
    if(fabs(MUL[0]*SQR(MUL[1]-us)+MUL[2]-MUR[0]*SQR(MUR[1]-us)-MUR[2])>tinyloc) return 0; // [rho*u*u+p] != 0
    if(fabs((0.5*MUL[0]*SQR(MUL[1]-us)+gam/(gam-1)*MUL[2])*(MUL[1]-us) -
            (0.5*MUR[0]*SQR(MUR[1]-us)+gam/(gam-1)*MUR[2])*(MUR[1]-us))>tinyloc) return 0; // [(E+p)*u] != 0
    // check that this is not a contact discontinuity
    if(fabs(us - MUL[1]) < tinyloc) return 0;
    // check the entropy inequality
    return (MUR[0] > MUL[0]) == (MUL[1]-us>0);
}
//======================================================================================================================

//======================================================================================================================
void s_Riemann::ReadParams(tFileBuffer& FB) {
    int KeyProblem = RIEMANN_CUSTOM; // one can specify a problem by a single key

    tParamManager PM;
    PM.Request(e[0], "ex"); // Shock direction
    PM.Request(e[1], "ey"); // Shock direction
    PM.Request(e[2], "ez"); // Shock direction
    PM.Request(DF_Coor, "DF_Coor"); // Shock position
    PM.RequestParameter(KeyProblem, "KeyProblem", "CUSTOM 0  SOD 1  LAX 2  ZHMF 3", IO_DONTCRASH);
    PM.ReadParamsFromBuffer(FB);

    switch(KeyProblem) {
    case RIEMANN_CUSTOM: break;
    case RIEMANN_SOD: // G. A. Sod, A survey of several finite difference methods for systems of nonlinear hyperbolic conservation law, Journal of Computational Physics, vol. 27, pp. 1-31, 1978.
        MUL[0]=MUL[2]=1.0; MUR[0]=0.125; MUR[2]=0.1; MUL[1]=MUR[1]=0.0;
        if(IsNaN(gam)) gam = 1.4;
        break;
    case RIEMANN_LAX: // P.D Lax, B.Wendroff. System of conservation laws. Comm. Pure Appl. Math., vol. XIII, pp. 217–237, 1960
        MUL[0]=0.445; MUL[1]=0.698; MUL[2]=3.528; MUR[0]=0.5; MUR[1]=0; MUR[2]=0.571;
        if(IsNaN(gam)) gam = 1.4;
        break;
    case RIEMANN_ZHMF: // A. I. Zhmakin, A. A. Fursenko, On a monotonic shock-capturing difference scheme, U.S.S.R. Comput. Maths. Math. Phys. Vol. 20, No. 4, pp. 218-227
        MUL[0]=8.0; MUL[2]=480.0; MUR[0]=MUR[2]=1.0; MUL[1]=MUR[1]=0.0;
        if(IsNaN(gam)) gam = C5_3;
        break;
    default:
        crash("KeyProblem is undefined");
    }

    PM.clear();
    PM.RequestOption(SteadyTest, "SteadyTest");
    PM.Request(gam, "gam"); // Specific ratio
    PM.Request(MUL[0], "RhoL");
    PM.Request(MUL[1], "UL");
    PM.Request(MUL[2], "PL");
    PM.Request(MUR[0], "RhoR");
    PM.Request(MUR[1], "UR");
    PM.Request(MUR[2], "PR");
    PM.ReadParamsFromBuffer(FB);

    // If a specific KeyProblem set, check that the specific ratio corresponds to the problem
    switch(KeyProblem) {
    case RIEMANN_SOD: // G. A. Sod, A survey of several finite difference methods for systems of nonlinear hyperbolic conservation law, Journal of Computational Physics, vol. 27, pp. 1-31, 1978.
    case RIEMANN_LAX: // P.D Lax, Weak solutions of nonlinear hyperbolic equations and their numerical computation, Comm. Pure Appl. Math., vol. 7 , pp. 159–193, 1954
        if(fabs(gam-1.4) > tiny) crash("s_Riemann::ReadParams error: gam=%f but 1.4 needed for Sod and Lax problem", gam);
        break;
    case RIEMANN_ZHMF: // A. I. Zhmakin, A. A. Fursenko, On a monotonic shock-capturing difference scheme, U.S.S.R. Comput. Maths. Math. Phys. Vol. 20, No. 4, pp. 218-227
        if(fabs(gam-C5_3) > tiny) crash("s_Riemann::ReadParams error: gam=%f but 5/3 needed for ZhmF problem", gam);
        break;
    default:;
    }
}
//======================================================================================================================

//======================================================================================================================
void s_Riemann::Init() {
    for(int i=0; i<3; i++) if(IsNaN(MUL[i]) || IsNaN(MUR[i])) crash("s_Riemann: values not set");
    for(int i=0; i<3; i++) if(IsNaN(e[i])) crash("s_Riemann: direction not set");
    { 
        double buf =  sqrt(VDOT(e,e)); 
        if(buf < 1e-50) crash("s_Riemann: direction vector is zero");
        e[0]/=buf; e[1]/=buf; e[2]/=buf;
    }
    for(int i=0; i<3; i++) if(IsNaN(gam)) crash("s_Riemann: specific ratio not set");

    if(SteadyTest) {
        int unst = IsUnsteady();
        if(unst) crash("s_Riemann: Steady flag set but condition %i is not satisfied", unst-1);
    }

    // Specify background field parameters using the background field with lower pressure
    // This background field is for the use in Noisette only. Not used in ColESo itself
    if(MUL[2]>MUR[2]) { FlowVel = MUR[1]; SoundVel = sqrt(gam*MUR[2]/MUR[0]); }
    else { FlowVel = MUL[1]; SoundVel = sqrt(gam*MUL[2]/MUL[0]); }
}
//======================================================================================================================

//======================================================================================================================
void s_Shock::ReadParams(tFileBuffer& FB) {
    double cr = 1.0;

    tParamManager PM;
    PM.Request(gam, "gam"); // Specific ratio
    PM.Request(DF_Coor, "DF_Coor"); // Shock position
    PM.Request(MUR[0], "Rho");
    PM.Request(MUR[1], "U");
    PM.Request(MUR[2], "P");
    PM.Request(cr, "C");
    PM.Request(Mach, "Mach");
    PM.Request(e[0], "ex"); // Shock direction
    PM.Request(e[1], "ey"); // Shock direction
    PM.Request(e[2], "ez"); // Shock direction
    PM.ReadParamsFromBuffer(FB);

    if(PM["C"].GetIsSet() && PM["P"].GetIsSet() && PM["Rho"].GetIsSet()) { // rho, p, c are given
        if(fabs(MUR[2] - cr*cr*MUR[0]/gam) > 1e-10) // Consistency check
            crash("s_Shock: both C, P and Rho given but P!=rho*C^2/gam");
    }
    else { // Two variables including sound speed are given, use them to calculate the third one
        if(PM["C"].GetIsSet() && PM["Rho"].GetIsSet()) MUR[2] = cr*cr*MUR[0]/gam;
        if(PM["C"].GetIsSet() && PM["P"].GetIsSet()) MUR[0] = gam*MUR[2]/(cr*cr);
    }
}
//======================================================================================================================

//======================================================================================================================
void s_Shock::Init() {
    // To this moment, the values in front of the shock and the Mach number of the shock should be set
    if(IsNaN(Mach) || Mach < 1.0) crash("s_Shock error: Mach = %e, should be >= 1", Mach);
    if(MUR[0] < tiny || MUR[2] < tiny) crash("s_Shock error: rhoR = %e, pR = %e", MUR[0], MUR[2]);
    double cr = sqrt(gam*MUR[2]/MUR[0]);

    // "Background flow" parameters as a values in front of the shock (not used in ColESo)
    FlowVel = fabs(MUR[1]); SoundVel = sqrt(gam*MUR[2]/MUR[0]);

    // Parameters behind the shock (MUL)
    double rmc = sqrt(gam*MUR[0]*MUR[2])*Mach;

    MUL[1] = Mach*cr*(1 - 2.0/(gam+1.0) * (gam + 1/(Mach*Mach)));
    MUL[0] = -rmc/MUL[1];
    MUL[2] = rmc*Mach*cr + MUR[2] - MUL[0]*MUL[1]*MUL[1];

    MUL[Var_U] += MUR[Var_U] + Mach*cr;

    pprintf("Field(L): %e %e %e\n", MUL[0], MUL[1], MUL[2]);
    pprintf("Field(R): %e %e %e\n", MUR[0], MUR[1], MUR[2]);
    if(MUL[2]<0 || MUL[0]<0) crash("s_Shock: internal error, density_L = %e, pressure_L = %e", MUL[0], MUL[2]);
}
//======================================================================================================================


//======================================================================================================================
void s_Riemann::PointValue(double t, const double* coor, double* V) const {
    double x = VDOT(coor,e) - DF_Coor;
    if(!IsUnsteady() || t<tiny) {
        double sign = SIGN(x);
        V[Var_R] = 0.5*(MUL[0]+MUR[0]) + 0.5*(MUR[0]-MUL[0])*sign;
        V[Var_W] = 0.5*(MUL[1]+MUR[1]) + 0.5*(MUR[1]-MUL[1])*sign;
        for(int i=0; i<3; i++) V[Var_U+i] = V[Var_W] * e[i];
        V[Var_P] = 0.5*(MUL[2]+MUR[2]) + 0.5*(MUR[2]-MUL[2])*sign;
        return; 
    }

    double lambda = x / t;
    double UAL[Var_N], UAR[Var_N];
    UAL[Var_R]=MUL[0]; UAL[Var_U]=MUL[1]*e[0]; UAL[Var_V]=MUL[1]*e[1]; UAL[Var_W]=MUL[1]*e[2]; UAL[Var_P] = MUL[2];
    UAR[Var_R]=MUR[0]; UAR[Var_U]=MUR[1]*e[0]; UAR[Var_V]=MUR[1]*e[1]; UAR[Var_W]=MUR[1]*e[2]; UAR[Var_P] = MUR[2];
    RMN_PrimVar(UAL, UAR, e, V, lambda, gam, 0);
}
//======================================================================================================================



// *********************************************************************************************************************
// *****                                                                                                           *****
// *****                          Exact solution for Riemann problem and Godunov solver                            *****
// *****                                                                                                           *****
// *********************************************************************************************************************

double f(double x, double p, double r, double gam) {
    double g1=0.5*(gam+1.0);
    double g2=0.5*(gam-1.0);
    double g3=g1/gam;
    double g4=g2/gam;

    double pi=x/p;
    double c=sqrt(gam*p/r);
    if(x >= p) return (x-p)/(r*c*sqrt(g3*pi+g4));
    else return c*(pow(pi,g4)-1.0)/g2;
}

double df(double x, double p, double r, double gam) {
    double g1=0.5*(gam+1.0);
    double g2=0.5*(gam-1.0);
    double g3=g1/gam;
    double g4=g2/gam;
    double g5=gam+1.0;
    double g6=3.0*gam-1.0;

    double pi=x/p;
    double c=sqrt(gam*p/r);
    if(x >= p) return (g5*pi+g6)/(4.*gam*r*c*pow(g3*pi+g4,1.5));
    else return c*pow(pi,g4)/(gam*x);
}

double dnewt(double r1, double p1, double c1, double r2, double p2, double c2, double udif, double gam, int log) {
    double g1=0.5*(gam+1.0);
    double g2=0.5*(gam-1.0);
    double g3=g1/gam;
    double g4=g2/gam;
    double g7=(gam-1.0)/(3.0*gam);

    int iter=0;
    double eps=1e-8;
    double x0=(p1*r2*c2+p2*r1*c1+udif*r1*c1*r2*c2)/(r1*c1+r2*c2);
    if(x0 < 0.0) {
        x0=1.0;
        if(log) printf("Modify iterations\n");
        double z=x0/(p1+p2);
        double a1, a2;
        if(x0 >= p1) a1=sqrt(r1*(g1*x0+g2*p1));
        else a1=g4*r1*c1*(1.0-x0/p1)/(1.0-pow(x0/p1,g4));
        if(x0 >= p2) a2=sqrt(r2*(g1*x0+g2*p2));
        else a2=g4*r2*c2*(1.0-x0/p2)/(1.0-pow(x0/p2,g4));
        double fi=(a2*p1+a1*p2+a1*a2*udif)/(a1+a2);
        double al=g7*(1.0-z)/((1.0-pow(z,g4))*pow(z,g3))-1.0;
        if(al < 0.0) al=0.0;
        double x1=(al*x0+fi)/(1.0+al);
        if(x1<0) x1 = 1e-10; //crash("Godunov solver error. Subroutine dnewt: x1 = %e < 0", x1);
        double ff=f(x1,p1,r1,gam)+f(x1,p2,r2,gam);
        x0=x1;
        iter++;
        if(log) printf("iter=%i, Residual=%e, x=%e\n", iter, fabs(x0-x1),x0);
        if(ff >= udif) crash("dnewt: error, no output value given!\n");
    }
    if(log) printf("Newton iterations\n");
    int flag = 0;
    while(1) {
        double dx = (f(x0,p1,r1,gam)+f(x0,p2,r2,gam)-udif) / (df(x0,p1,r1,gam)+df(x0,p2,r2,gam));
        double x1=x0-dx;
        iter++;
        if(log) printf("iter=%i, Residual=%e, x=%e\n", iter, fabs(dx),x0);
        if(fabs(dx) < eps) { if(flag) return x1; flag = 1; }
        if(iter==1000) crash("Riemann solver: Newton process failed!\n");
        x0=x1;
    }
}

//======================================================================================================================
// Exact solution for Riemann problem
// Input:  UAL, UAR = conservative vars, 
//         normals = (nx, ny, nz) = unit normal
//         lambda = minus reference frame velocity
//         gam = specific ratio
//         Var_NumTurb = number of extra variables moving with fluid velocity
// Output: U0 = physical(!) vars at x/t=lambda
//======================================================================================================================
void RMN_PrimVar(const double* UAL, const double* UAR, const double* normals, double* U0, 
                 double lambda, double gam, int Var_NumTurb) {
    double g1=0.5*(gam+1.0);
    double g2=0.5*(gam-1.0);
    double g3=g1/gam;
    double g4=g2/gam;

    double r1 = UAL[Var_R], p1 = UAL[Var_P], r2 = UAR[Var_R], p2 = UAR[Var_P];
    double u1 = VDOT(UAL+Var_U, normals);
    double u2 = VDOT(UAR+Var_U, normals);

    const int log = 0;

    int lp = 0; // invert x axis
    if(p1 > p2) {
        lp=1;
        u1=-u1;
        u2=-u2;
        SWAP(u1,u2);
        SWAP(r1,r2);
        SWAP(p1,p2);
        lambda=-lambda;
    }

    double udif=u1-u2;
    double c1=sqrt(gam*p1/r1);
    double c2=sqrt(gam*p2/r2);
// -----------------------------------------
    double usw=f(p2,p1,r1,gam)+f(p2,p2,r2,gam);
    double urw=f(p1,p1,r1,gam)+f(p1,p2,r2,gam);
    double uvc=-(c1+c2)/g2;

    double uc; // Скорость контактного разрыва
    double un0; // Скорость при x=0
// -----------------------------------------

    if(udif >= usw) {
//  shock wave-shock wave configuration
        double pc=dnewt(r1,p1,c1,r2,p2,c2,udif,gam,log);
        double a1=sqrt(r1*(g1*pc+g2*p1));
        double a2=sqrt(r2*(g1*pc+g2*p2));
               uc=(a1*u1+a2*u2+p1-p2)/(a1+a2);
        double dl=u1-a1/r1;
        double dr=u2+a2/r2;
        double rl=r1*a1/(a1-r1*(u1-uc));
        double rr=r2*a2/(a2+r2*(u2-uc));
        if(dl>lambda)      { U0[Var_R] = r1; un0 = u1; U0[Var_P] = p1; }
        else if(uc>lambda) { U0[Var_R] = rl; un0 = uc; U0[Var_P] = pc; }
        else if(dr>lambda) { U0[Var_R] = rr; un0 = uc; U0[Var_P] = pc; }
        else               { U0[Var_R] = r2; un0 = u2; U0[Var_P] = p2; }
    }
// -----------------------------------------
    else if(udif >= uvc && udif <= urw) {
//  rarefaction wave-rarefaction wave configuration
        double pc=p1*pow((u1-u2-uvc)/(urw-uvc), 1.0/g4);
        {
            double a1, a2;
            if(fabs(pc-p1) < 1e-7) a1 = r1*c1*(1.-0.5*g3*(1-pc/p1));
            else a1=g4*r1*c1*(1.0-pc/p1)/(1.0-pow(pc/p1,g4));
            if(fabs(pc-p2) < 1e-7) a2 = r2*c2*(1.-0.5*g3*(1-pc/p2));
            else a2=g4*r2*c2*(1.0-pc/p2)/(1.0-pow(pc/p2,g4));
            uc=(a1*u1+a2*u2+p1-p2)/(a1+a2);
        }
        double cc1=c1+g2*(u1-uc);
        double cc2=c2-g2*(u2-uc);
        double rl=gam*pc/SQR(cc1);
        double rr=gam*pc/SQR(cc2);

        if(u1-c1>lambda)       { U0[Var_R] = r1; un0 = u1; U0[Var_P] = p1; }
        else if(uc-cc1>lambda) {
            double par1=g2*u1/c1+1.0;
            un0=(lambda-(u1-c1))/g1+u1;
            double ama = -un0/(lambda-un0);
            double ama1=par1/(g2*ama+1.0);
            U0[Var_R]=r1*pow(ama1,1.0/g2);
            U0[Var_P]=p1*pow(ama1,1.0/g4);
        }
        else if(uc>lambda)     { U0[Var_R] = rl; un0 = uc; U0[Var_P] = pc; }
        else if(uc+cc2>lambda) { U0[Var_R] = rr; un0 = uc; U0[Var_P] = pc; }
        else if(u2+c2>lambda)  {
            double par2=g2*u2/c2-1.0;
            un0=(lambda-(u2+c2))/g1+u2;
            double ama = un0/(lambda-un0);
            double ama1=par2/(g2*ama-1.0);
            U0[Var_R]=r2*pow(ama1,1.0/g2);
            U0[Var_P]=p2*pow(ama1,1.0/g4);
        }
        else              { U0[Var_R] = r2; un0 = u2; U0[Var_P] = p2; }
    }
// -----------------------------------------
    else if(udif >= urw && udif <= usw) {
//  shock wave-rarefaction wave configuration
        double pc=dnewt(r1,p1,c1,r2,p2,c2,udif,gam,log);
        double a1=sqrt(r1*(g1*pc+g2*p1));
        {
            double a2;
            if(fabs(pc-p2) < 1e-7) a2 = r2*c2*(1.-0.5*g3*(1-pc/p2));
            else a2=g4*r2*c2*(1.0-pc/p2)/(1.0-pow(pc/p2,g4));
            uc=(a1*u1+a2*u2+p1-p2)/(a1+a2);
        }
        double dl=u1-a1/r1;
        double rl=r1*a1/(a1-r1*(u1-uc));
        double dr2=u2+c2;
        double cc2=c2-g2*(u2-uc);
        double dr1=uc+cc2;
        double rr=gam*pc/SQR(cc2);

        if(dl>lambda)       { U0[Var_R] = r1; un0 = u1; U0[Var_P] = p1; }
        else if(uc>lambda)  { U0[Var_R] = rl; un0 = uc; U0[Var_P] = pc; }
        else if(dr1>lambda) { U0[Var_R] = rr; un0 = uc; U0[Var_P] = pc; }
        else if(dr2>lambda) {
            double par2=g2*u2/c2-1.0;
            un0=(lambda-dr2)/g1+u2;
            double ama = un0/(lambda-un0);
            double ama1=par2/(g2*ama-1.0);
            U0[Var_R]=r2*pow(ama1,1.0/g2);
            U0[Var_P]=p2*pow(ama1,1.0/g4);
        }
        else           { U0[Var_R] = r2; un0 = u2; U0[Var_P] = p2; }
    }
    else crash("RMN_PrimVar: error: vacuum detected");

    if(lp) { un0 = -un0; uc = -uc; lambda = -lambda; } // restoring signs of the velocities if we changed x <-> -x

    // pressure and density are already set
    // now set the values of velocity and extra variables using the values from the left or from the right
    const double* UU = (uc-lambda > 0 ? UAL : UAR);
    for(int ivar=0; ivar<Var_NumTurb; ivar++) U0[Var_N+ivar] = UU[Var_N+ivar];
    for(int ivar=Var_U; ivar<=Var_W; ivar++) U0[ivar] = UU[ivar];

    // it remains to corrcet the normal component or the velocity
    double dun0 = un0 - VDOT(UU+Var_U, normals);
    U0[Var_U] += dun0*normals[0];
    U0[Var_V] += dun0*normals[1];
    U0[Var_W] += dun0*normals[2];
}
//======================================================================================================================
