// *********************************************************************************************************************
// *****                                                                                                           *****
// *****                                  Collection of exact solutions (ColESo)                                   *****
// *****                                    Planar waves with different forms;                                     *****
// *****                                     acoustic-shock interaction (1D);                                      *****
// *****                                      sine wave in 3D box (viscous)                                        *****
// *****                                                                                                           *****
// *********************************************************************************************************************

#include "coleso.h"
#include "es_utils.h"
#include "base_parser.h"
#include "geom_primitive.h"
#include <complex>

//======================================================================================================================
void s_PlanarSinus::PointValue(double t, const double* coord, double* uex) const {
    double C = sqrt(gam * PRef / RhoRef);

    // Wave front velocity projected on the normal direction
    double Cn = VDOT(URef, Dir) + C;
    double R[3] = { coord[0] - Xterm, coord[1] - Yterm, coord[2] - Zterm };
    double Rn = VDOT(R, Dir);
    double Rt2 = VDOT(R, R) - Rn * Rn;
    double phase = 0.0, expr = 0.0;

    if (Rt2 <= Rmax * Rmax) { // if Rmax (channel radius) is set, check that the point is inside the channel
        phase = t - Rn / Cn + Phase; // phase = delta(T)
        if (phase > 0 || AllowNegativePhase) expr = Aterm * sin(Pi2 * Freq * phase);
    }
    // Aterm is pressure amplitude by definition
    uex[Var_P] = expr;
    uex[Var_R] = expr / (C * C);
    for (int ivar = Var_U; ivar <= Var_W; ivar++)
        uex[ivar] = expr * Dir[ivar - 1] / (RhoRef * C);
}
//======================================================================================================================


//======================================================================================================================
void s_PlanarGauss::PointValue(double t, const double* coord, double* uex) const {
    double C = sqrt(gam * PRef / RhoRef);

    // Wave front velocity projected on the normal direction
    double Cn = VDOT(URef, Dir) + C;
    double R[3] = { coord[0] - Xterm, coord[1] - Yterm, coord[2] - Zterm };
    double Rn = VDOT(R, Dir);
    double Rt2 = VDOT(R, R) - Rn * Rn;
    double phase = 0.0, expr = 0.0;

    if (Rt2 <= Rmax * Rmax) { // if Rmax (channel radius) is set, check that the point is inside the channel
        phase = Rn - Cn * t; // phase is delta(x)
        if (NormalDistance < 0.5 * huge) RoundToCentre(phase, NormalDistance);
        expr = Aterm * exp(-log(2.0) * phase * phase / Bterm / Bterm);
    }

    // Aterm is pressure amplitude by definition
    uex[Var_P] = expr;
    uex[Var_R] = expr / (C * C);
    for (int ivar = Var_U; ivar <= Var_W; ivar++)
        uex[ivar] = expr * Dir[ivar - 1] / (RhoRef * C);
}
//======================================================================================================================



//======================================================================================================================
void s_Planar::ReadParams(tFileBuffer& FB) {
    tParamManager PM;
    PM.Request(Aterm, "Aterm"); // Amplitude of pressure pulsation
    PM.Request(Xterm, "Xterm");
    PM.Request(Yterm, "Yterm");
    PM.Request(Zterm, "Zterm");
    PM.Request(Dir[0], "DirX");
    PM.Request(Dir[1], "DirY");
    PM.Request(Dir[2], "DirZ");
    PM.Request(Rmax, "Rmax"); // If set, solution is set to zero outsize the infinite circular cylinder 
                              // with radius Rmax and axis direction DirX,DirY,DirZ

    // Background field
    PM.Request(RhoRef, "RhoRef");
    PM.Request(URef[0], "URef");
    PM.Request(URef[1], "VRef");
    PM.Request(URef[2], "WRef");
    PM.Request(PRef, "PRef");
    PM.Request(gam, "gam"); // specific ratio

    PM.ReadParamsFromBuffer(FB);
}
//======================================================================================================================

//======================================================================================================================
void s_PlanarSinus::ReadParams(tFileBuffer& FB) {
    tParamManager PM;
    s_Planar::ReadParams(FB);
    PM.Request(Freq, "Freq");
    PM.Request(Phase, "Phase"); // В нулевой фазе синус равен нулю
    PM.RequestOption(AllowNegativePhase, "AllowNegativePhase");
    PM.ReadParamsFromBuffer(FB);
}
//======================================================================================================================

//======================================================================================================================
void s_PlanarGauss::ReadParams(tFileBuffer& FB) {
    tParamManager PM;
    s_Planar::ReadParams(FB);
    PM.Request(Bterm, "Bterm");
    PM.Request(NormalDistance, "NormalDistance");
    PM.ReadParamsFromBuffer(FB);
}
//======================================================================================================================

//======================================================================================================================
void s_Planar::Init() {
    // Normalizing direction
    double rr = 0.0;
    for(int icoor=0; icoor<3; icoor++) rr += SQR(Dir[icoor]);
    if(rr < tiny) crash("s_Planar::Init: wave direction can't be zero");
    for(int icoor=0; icoor<3; icoor++) Dir[icoor] /= sqrt(rr);
}
//======================================================================================================================


// *********************************************************************************************************************
// *****                                                                                                           *****
// *****                     Planar wave with a profile prescribed by a Chebyshev polynomial                        *****
// *****                                                                                                           *****
// *********************************************************************************************************************


//======================================================================================================================
void s_Polynom::PointValue(double T, const double* coord, double* uex) const {
    double x, p1, p2, p3;
    if(Degree==0)
        p3 = 1.0;
    else {
        x = (coord[0]-T*SoundSpeed - 0.5*(Xmin+Xmax)) / (Xmax - Xmin);
        p2 = 1.0;
        p3 = x;
        for(int i=1; i<Degree; i++) {
            p1 = p2;
            p2 = p3;
            p3 = 2.*x*p2 - p1;
        }
    }

    uex[Var_R] = Aterm * p3;
}
//======================================================================================================================

//======================================================================================================================
void s_Polynom::ReadParams(tFileBuffer& FB) {
    tParamManager PM;
    PM.RequestParameter(Aterm, "aterm");
    PM.RequestParameter(Xmin, "Xmin");
    PM.RequestParameter(Xmax, "Xmax");
    PM.RequestParameter(Degree, "Degree");
    PM.ReadParamsFromBuffer(FB);
}
//======================================================================================================================


// *********************************************************************************************************************
// *****                                                                                                           *****
// *****                        Propagation of an acoustic perturbation through a shock                            *****
// *****                                                                                                           *****
// *********************************************************************************************************************

//======================================================================================================================
void s_AcousticShock::ReadParams(tFileBuffer& FB) {
    tParamManager PM;
    // Initial Gaussian acoustic pulse
    PM.Request(Aterm, "Aterm");
    PM.Request(Bterm, "Bterm");
    PM.Request(Xterm, "Xterm");
    // Background field
    PM.Request(gam, "gam"); // Specific ratio
    PM.Request(DF_Coor, "DF_Coor"); // Shock position
    PM.Request(MUL[0], "RhoL");
    PM.Request(MUL[1], "UL");
    PM.Request(MUL[2], "PL");
    PM.Request(MUR[0], "RhoR");
    PM.Request(MUR[1], "UR");
    PM.Request(MUR[2], "PR");
    PM.ReadParamsFromBuffer(FB);
}
//======================================================================================================================

//======================================================================================================================
void s_AcousticShock::Init(void) {
    // Check that the fields do form a shock
    s_Riemann SR;
    SR.gam = gam;
    SR.MUL[0]=MUL[0]; SR.MUL[1]=MUL[1]; SR.MUL[2]=MUL[2]; 
    SR.MUR[0]=MUR[0]; SR.MUR[1]=MUR[1]; SR.MUR[2]=MUR[2];
    if(!SR.IsShock(&VShock)) crash("s_AcousticShock: background field is not a shock");

    if(gam*MUL[2] > MUL[0]*SQR(MUL[1]-VShock)) crash("s_AcousticShock: left field must be supersonic");
    if(gam*MUR[2] < MUR[0]*SQR(MUR[1]-VShock)) crash("s_AcousticShock: left field must be subsonic");
    if(Xterm>=DF_Coor) crash("s_AcousticShock: initial wave position is right to the shock");
}
//======================================================================================================================

//======================================================================================================================
void FillPrim(double *A, double *ref, double _gam) {
    double gam_gam1 = _gam / (_gam - 1.0);
    double &R = ref[0], &U = ref[1], &P = ref[2];
    double UV2 = 0.5*U*U;
    A[0]=U;      A[1]=R;                      A[2]=0.0;
    A[3]=U*U;    A[4]=2.0*R*U;                A[5]=1.0;
    A[6]=U*UV2;  A[7]=gam_gam1*P+R*UV2+R*U*U; A[8]=gam_gam1*U;
}
//======================================================================================================================

//======================================================================================================================
// Sonving a linear 3x3 system using the Cramer rule
// Tolerance -- минимальный допустимый детерминант
//======================================================================================================================
static int SolveSLAUCramer3(const double* M, const double* B, double* X, double tolerance) {
    double det = M[0]*(M[4]*M[8]-M[5]*M[7]) + M[1]*(M[5]*M[6]-M[3]*M[8]) + M[2]*(M[3]*M[7]-M[4]*M[6]);
    if(fabs(det) < tolerance) return 1; // degenerate
    det = 1.0/det;
    X[0] = (B[0]*(M[4]*M[8]-M[5]*M[7]) + M[1]*(M[5]*B[2]-B[1]*M[8]) + M[2]*(B[1]*M[7]-M[4]*B[2]))*det;
    X[1] = (M[0]*(B[1]*M[8]-M[5]*B[2]) + B[0]*(M[5]*M[6]-M[3]*M[8]) + M[2]*(M[3]*B[2]-B[1]*M[6]))*det;
    X[2] = (M[0]*(M[4]*B[2]-B[1]*M[7]) + M[1]*(B[1]*M[6]-M[3]*B[2]) + B[0]*(M[3]*M[7]-M[4]*M[6]))*det;
    return 0;
}
//======================================================================================================================

//======================================================================================================================
void s_AcousticShock::GetConstitutor(double* UR) const {
    double AL[9], AR[9];
    double vec[3], Qref1[3], Qref2[3], UL[3];
    double mUL[3] = {MUL[0], MUL[1]-VShock, MUL[2]};
    double mUR[3] = {MUR[0], MUR[1]-VShock, MUR[2]};
    double CL = sqrt(gam*MUL[2]/MUL[0]); // Sound speed
    double CR = sqrt(gam*MUR[2]/MUR[0]); // Sound speed
    double gam1 = gam - 1.0;

    // One acoustic wave from the left
    UL[0] = 1.0 / (CL*CL);
    UL[1] = 1.0 / (CL*MUL[0]);
    UL[2] = 1.0;
    // Calculation what is from the right
    FillPrim(AL, mUL, gam);
    vec[0] = AL[0]*UL[0] + AL[1]*UL[1] + AL[2]*UL[2]; // Vec = AL*UL
    vec[1] = AL[3]*UL[0] + AL[4]*UL[1] + AL[5]*UL[2]; // Vec = AL*UL
    vec[2] = AL[6]*UL[0] + AL[7]*UL[1] + AL[8]*UL[2]; // Vec = AL*UL

    double num = 0.0;
    num += vec[0] * (0.5*CR*mUR[1] + 0.25*gam1*SQR(mUR[1]));
    num += vec[1] * (-0.5 * (CR + mUR[1]*gam1));
    num += vec[2] * ( 0.5 * gam1);

    double denom = (0.5*gam1*mUR[1] - 0.5*(gam-3.0)*mUL[1] + CR) * mUR[1] * (mUR[0]-mUL[0]);
    double D = -2.0 * num / denom;

    Qref1[0]=mUL[0]; Qref1[1]=mUL[0]*mUL[1]; Qref1[2]=0.5*mUL[0]*mUL[1]*mUL[1]+mUL[2]/gam1;
    Qref2[0]=mUR[0]; Qref2[1]=mUR[0]*mUR[1]; Qref2[2]=0.5*mUR[0]*mUR[1]*mUR[1]+mUR[2]/gam1;
    for(int i=0; i<3; i++) vec[i] += D * (Qref2[i] - Qref1[i]);

    FillPrim(AR, mUR, gam);
    double vvec[3];
    if(SolveSLAUCramer3(AR, vec, vvec, 1e-50)) crash("s_AcousticShock: determinant = 0");

    // Decompose into modes
    UR[0] = 0.5*(-vvec[1]*mUR[0]*CR + vvec[2]);
    UR[1] = 0.5*( vvec[1]*mUR[0]*CR + vvec[2]);
    UR[2] = vvec[0]*CR*CR - vvec[2];

    if(fabs(UR[0])>1e-10) crash("s_AcousticShock: left acoustic wave appeared");
}
//======================================================================================================================

//======================================================================================================================
void s_AcousticShock::PointValue(double t, const double* coord, double* uex) const {
    double CL = sqrt(gam*MUL[2]/MUL[0]); // Sound speed
    double CR = sqrt(gam*MUR[2]/MUR[0]); // Sound speed

    double XShock = DF_Coor + t * VShock; // Shock position
    if(MUR[0]<MUL[0]) crash("s_AcousticShock::PointValue error: RhoL=%e, RhoR=%e", MUL[0], MUR[0]);

    uex[Var_V] = 0.0;
    uex[Var_W] = 0.0;

    if(coord[0]<=XShock || fabs(VShock)<tiny) { // Supersonic domain
        double phase = (coord[0] - Xterm - t * (CL + MUL[1])) / Bterm;
        double expconst = Aterm*exp(-log(2.0)*SQR(phase));
        uex[Var_R] = expconst / (CL*CL);
        uex[Var_U] = expconst / (CL * MUL[0]);
        uex[Var_P] = expconst;
    }
    else { // Subsonic domain: right acoustic wave and entropy wave. Zero value for left acoustic wave
        double UR[3];
        GetConstitutor(UR);
        // For each mode, calculate the time this mode was at the shock position
        // 1. Right acoustic wave
        double tim = t - (coord[0]-XShock) / (CR + MUR[1] - VShock);
        double phase = (DF_Coor + tim * VShock) - Xterm - tim * (CL + MUL[1]);
        double expconst = Aterm*exp(-log(2.0)*SQR(phase/Bterm));
        uex[Var_R] = expconst * UR[1] / (CR*CR);
        uex[Var_U] = expconst * UR[1] / (CR*MUR[0]);
        uex[Var_P] = expconst * UR[1];
        // 2. Entropy wave
        tim = t - (coord[0]-XShock) / (MUR[1] - VShock);
        phase = (DF_Coor + tim * VShock) - Xterm - tim * (CL + MUL[1]);
        expconst = Aterm*exp(-log(2.0)*SQR(phase/Bterm));
        uex[Var_R] += expconst * UR[2] / (CR*CR);
    }
}
//======================================================================================================================


// *********************************************************************************************************************
// *****                                                                                                           *****
// *****                           4 pulses of different forms (scalar case)                                       *****
// *****                                                                                                           *****
// *********************************************************************************************************************

static double F(double x, double alpha, double a) {
    return sqrt(MAX(1.0 - SQR(alpha*(x-a)), 0.0));
}

//======================================================================================================================
void s_4peak::PointValue(double T, const double* coord, double* uex) const {
    double mul = 2.0 / (Xmax - Xmin);
    double x = (coord[0]-T*SoundSpeed - 0.5*(Xmin+Xmax)) * mul;
    RoundToCentre(x, Period*mul);
    
    double f = 0.0;
    double a = 0.5, z = -0.7, delta = 0.005, alpha = 10.0, beta = log(2.0)/(36.0*delta*delta);
    if(x>=-0.8 && x<=-0.6) f = C1_6*(exp(-beta*SQR(x-(z-delta))) + exp(-beta*SQR(x-(z+delta))) + 4.*exp(-beta*SQR(x-z)));
    if(x>=-0.4 && x<=-0.2) f = 1.0;
    if(x>=0.0 && x<=0.2) f = 1.0 - fabs(10.0*(x-0.1));
    if(x>=0.4 && x<=0.6) f = C1_6*(F(x, alpha, a-delta) + F(x, alpha, a+delta) + 4.0*F(x, alpha, a));

    uex[Var_R] = Aterm * f;
}
//======================================================================================================================

//======================================================================================================================
void s_4peak::ReadParams(tFileBuffer& FB) {
    tParamManager PM;
    PM.Request(Aterm, "Aterm");
    PM.Request(Xmin, "Xmin");
    PM.Request(Xmax, "Xmax");
    PM.Request(Period, "Period");
    PM.Request(SoundSpeed, "SoundSpeed");
    PM.ReadParamsFromBuffer(FB);
}
//======================================================================================================================

