// *********************************************************************************************************************
// *****                                                                                                           *****
// *****                                  Collection of exact solutions (ColESo)                                   *****
// *****                    Wave propagation from 2D and 3D Gaussian pulses; entropy/vortex wave                   *****
// *****                                                                                                           *****
// *********************************************************************************************************************

#include "base_parser.h"
#include "geom_primitive.h" 
#include "coleso.h"
#include "es_specfunc.h"
#ifdef EXTRAPRECISION_COLESO
#include "extraprecision.h"
#endif


// *********************************************************************************************************************
// *****                                                                                                           *****
// *****                                Wave propagation from a 2D Gaussian pulse                                  *****
// *****                                                                                                           *****
// *********************************************************************************************************************

//======================================================================================================================
template<typename fpv>
void s_Gaussian2D<fpv>::ReadParams(tFileBuffer& FB) {
    // Reading pulse parameters
    tSpaceForm<fpv>::Read(FB, true /* allow periodics */);

    tParamManager PM;
    // Background flow velocity
    PM.Request(FlowVelX, "FlowVelX");
    PM.Request(FlowVelY, "FlowVelY");
    PM.ReadParamsFromBuffer(FB);
}
//======================================================================================================================

//======================================================================================================================
template<typename fpv>
void s_Gaussian2D<fpv>::Init(void){
    tSpaceForm<fpv>::Init();
    if(tSpaceForm<fpv>::Form!=tSpaceForm<fpv>::FORM_GAUSSIAN) crash("s_Gaussian2D: non-Gaussian form");
    GI.Init(); // init of Gaussian quadrature rules
}
//======================================================================================================================


//======================================================================================================================
// Computation of the solution by the direct integration of what is given by the Fourier transformation
//======================================================================================================================
template<typename fpv> 
tFixBlock<fpv,2> s_Gaussian2D_func(fpv omega, void* args) {
    const fpv& omegar = omega*((fpv*)args)[0];
    const fpv& omegat = omega*((fpv*)args)[1];
    tFixBlock<fpv,2> f;
    f[0] = omega * BesselJ0(omegar) * cos(omegat);
    f[1] = omega * BesselJ1(omegar) * sin(omegat);
    return f;
}

template<typename fpv>
void s_Gaussian2D<fpv>::CalcDirect(fpv t, fpv r, fpv& p, fpv& ur) const {
    fpv mult = sqrt(fpv(2.0)*GetLn2<fpv>()) * tSpaceForm<fpv>::InvBTerm;
    fpv args[2] = {r*mult, t*mult};
    tFixBlock<fpv,2> f = GI.Integrate(0.0, 1e50, 1.0, 0.0, 0, 0, 1.0, MAX(args[0],args[1]), &s_Gaussian2D_func<fpv>, (void*)args);
    p = f[0];
    ur = f[1];
}
//======================================================================================================================


//======================================================================================================================
// Computation of the solution using the Parseval identity for the Bessel-Fourier transformation
//======================================================================================================================
template<typename fpv> 
tFixBlock<fpv,3> s_Gaussian2D_funcBF(fpv xi, void* args) {
    const fpv& r = ((fpv*)args)[0];
    const fpv& t = ((fpv*)args)[1];
    const fpv xx = 1.0 - xi;
    fpv I[2];
    BesselI<fpv>(r*t*xx, 2, get_eps<fpv>()*100.0, false, I); // I0(x)*exp(-x), I1(x)*exp(-x)
    const fpv mult = xx / sqrt(1.0 + xx);
    I[0] *= mult; I[1] *= mult;
    tFixBlock<fpv,3> f;
    f[0] = I[0];       // J_{0,1}
    f[1] = I[1]*xx;    // J_{1,2}
    f[2] = I[0]*xx*xx; // J_{0,3}
    return f;
}

template<typename fpv>
void s_Gaussian2D<fpv>::CalcViaBesselFourier(fpv t, fpv r, fpv& p, fpv& ur) const {
    if(t < r*1e-10) crash("CalcViaBesselFourier error: t is too small");
    fpv mult = sqrt(fpv(2.0)*GetLn2<fpv>()) * tSpaceForm<fpv>::InvBTerm;
    r *= mult; t *= mult;
    fpv args[2] = {r, t};
    fpv q = MAX(MAX(fpv(1.0), r), r*t);
    tFixBlock<fpv,3> f = GI.Integrate(0.0, 1.0, t, 1.0-r/t, 0, 1, 1.0, q, &s_Gaussian2D_funcBF<fpv>, (void*)args);
    p  = f[0] - t*t*f[2] + r*t*f[1];
    ur = - t*t*f[1] + r*t*f[0];
}
//======================================================================================================================


//======================================================================================================================
// Computation of the solution using the Parseval identity for the Fourier transformation
//======================================================================================================================
template<typename fpv> 
tFixBlock<fpv,2> s_Gaussian2D_funcF(fpv xi, void* /*args*/) {
    tFixBlock<fpv,2> f;
    f[0] = 1./sqrt(2.+xi);
    f[1] = -f[0]*(1.+xi);
    return f;
}
template<typename fpv>
void s_Gaussian2D<fpv>::CalcViaFourier(fpv t, fpv r, fpv& p, fpv& ur) const {
    fpv mult = sqrt(fpv(2.0)*GetLn2<fpv>()) * tSpaceForm<fpv>::InvBTerm;
    r *= mult;
    t *= mult;
    if(r < 1e-10) crash("CalcViaFourier error: r is too small");
    fpv q = MAX(fpv(0.5), r);
    tFixBlock<fpv,2> f = GI.Integrate(0.0, 1e50, r, t/r-1.0, 1, 1, sqrt(t/r), q, &s_Gaussian2D_funcF<fpv>, NULL);
    tFixBlock<fpv,2> g = GI.Integrate(0.0, 1e50, r, -t/r-1.0, 1, 1, sqrt(t/r), q, &s_Gaussian2D_funcF<fpv>, NULL);
    mult = r / sqrt(GetPiNumber2<fpv>());
    p  = (f[0]+g[0]) * mult;
    ur = (-f[1]+g[1]) * mult;
}
//======================================================================================================================


//======================================================================================================================
template<typename fpv>
void s_Gaussian2D<fpv>::PointValue(fpv t, const fpv* coord, fpv* uex) const{
    for(int ivar=Var_R; ivar<Var_N; ivar++) uex[ivar] = 0.0;
    if(fabs(tSpaceForm<fpv>::Aterm) < tiny) return;

    for(int iPerX = -tSpaceForm<fpv>::MaxPer[0]; iPerX <= tSpaceForm<fpv>::MaxPer[0]; iPerX++)
    for(int iPerY = -tSpaceForm<fpv>::MaxPer[1]; iPerY <= tSpaceForm<fpv>::MaxPer[1]; iPerY++) {
        if(tSpaceForm<fpv>::Checkerboard) if((iPerX+iPerY)&1) continue;
        if(fabs(tSpaceForm<fpv>::PerX) > 0.5*huge && iPerX) continue;
        if(fabs(tSpaceForm<fpv>::PerY) > 0.5*huge && iPerY) continue;
        // Coordinates, taking into account the background flow and periodical b.c.
        fpv x = coord[0] - tSpaceForm<fpv>::r0[0] - FlowVelX*t + tSpaceForm<fpv>::PerX * iPerX;
        fpv y = coord[1] - tSpaceForm<fpv>::r0[1] - FlowVelY*t + tSpaceForm<fpv>::PerY * iPerY;
        fpv r = sqrt(x*x+y*y);

        fpv p, ur;
        if(t < tiny) { // Calculating the solution using the Taylor expansion at t=0
            fpv alpha = - GetLn2<fpv>() * tSpaceForm<fpv>::InvBTerm * tSpaceForm<fpv>::InvBTerm;
            fpv f0 = exp(alpha*r*r); // значение
            p = f0; ur = 0.0;
            if(t > 0.0) {
                fpv f1_r = 2.0*alpha*f0; // first derivative divided by l
                fpv f2 = 2.0*alpha*(1.0 + 2.0*alpha*r)*f0; // second derivative
                p += 0.5*t*t*(f2 + f1_r);
                ur = - t * f1_r * r;
            }
        }
        else if(r > tSpaceForm<fpv>::Bterm) CalcViaFourier(t, r, p, ur);
        else if(t < tSpaceForm<fpv>::Bterm*(2 + 1.2*sizeof(fpv))) CalcDirect(t, r, p, ur);
        else CalcViaBesselFourier(t, r, p, ur);

        uex[Var_R] += p * tSpaceForm<fpv>::Aterm;
        if(r > tiny) {
            uex[Var_U] += ur*x/r * tSpaceForm<fpv>::Aterm;
            uex[Var_V] += ur*y/r * tSpaceForm<fpv>::Aterm;
        }
    }
    uex[Var_P] = uex[Var_R];
}
//======================================================================================================================


// Template classes instantiation
template struct s_Gaussian2D<NativeDouble>;
#ifdef EXTRAPRECISION_COLESO
template struct s_Gaussian2D<dd_real>;
template struct s_Gaussian2D<qd_real>;
#endif



// *********************************************************************************************************************
// *****                                                                                                           *****
// *****                                Wave propagation from a 3D Gaussian pulse                                  *****
// *****                                                                                                           *****
// *********************************************************************************************************************


//======================================================================================================================
void s_Gaussian3D::ReadParams(tFileBuffer& FB) {
    // Reading the pulse parameters
    tSpaceForm<double>::Read(FB, true /* allow periodics */);

    tParamManager PM;
    // Background flow velocity
    PM.Request(FlowVelX, "FlowVelX");
    PM.Request(FlowVelY, "FlowVelY");
    PM.Request(FlowVelZ, "FlowVelZ");
    PM.ReadParamsFromBuffer(FB);
}
//======================================================================================================================

//======================================================================================================================
void s_Gaussian3D::Init(void){
    tSpaceForm<double>::Init();
    if(Form!=tSpaceForm<double>::FORM_GAUSSIAN)  crash("s_Gaussian3D: non-Gaussian form");
}
//======================================================================================================================


//======================================================================================================================
// Increment from a single pulse
//======================================================================================================================
void s_Gaussian3D::CalcValues(double t, double x, double y, double z, double* uex) const {
    // non-dimentionalizing
    const double qqq = sqrt(2.*CLN2);
    double mult = qqq * InvBTerm; 
    t *= mult; x *= mult; y *= mult; z *= mult;
    double r = sqrt(x*x+y*y+z*z);

    double P = 0.0, URr = 0.0; // pressure, Ur/r
    const double X = t * r;
    if(X < 50.0) {
        double A = 0., B = 0.; // Need to compute A = sh(X)/X, B = (sh(X) - X*ch(X)) / X^3
        const double ep = exp(X), em = exp(-X);
        const double ch = 0.5*(ep + em), sh = 0.5*(ep - em);
        if(X > 0.1) {
            A = sh / X;
            B = (sh - X * ch) / (X*X*X);
        }
        else { // using Taylor expansion
            double x2 = X * X;
            A =   ((((x2/110.+1.)*x2/72.+1.)*x2/42.+1.)*x2/20.+1.)*x2/6.+1.;
            B = -(((((x2/130.+1.)*x2/88.+1.)*x2/54.+1.)*x2/28.+1.)*x2/10.+1.)*C1_3;
        }
        const double m = exp(-0.5*(t*t+r*r));
        P = m*(ch - t*t*A);
        URr = m*(t*A + t*t*t*B);
    }
    else if(r > 1.0) { // r > 1.0, t*r > 50.0
        double AexpMinus = 0.5*exp(-0.5*(t-r)*(t-r));
        double AexpPlus  = 0.5*exp(-0.5*(t+r)*(t+r));
        P = AexpMinus * (1.0 - t/r) + AexpPlus * (1 + t/r);
        URr = AexpMinus * (1.0 - t/r + 1.0/(r*r)) -  AexpPlus * (1.0 + t/r + 1.0/(r*r));
        URr /= r;
    }
    // If t*r > 50.0 and r < 1.0 assume P = URr = 0

    P *= Aterm;
    URr *= Aterm;
    uex[Var_R] += P;
    uex[Var_P] += P;
    uex[Var_U] += URr * x;
    uex[Var_V] += URr * y;
    uex[Var_W] += URr * z;
}
//======================================================================================================================


//======================================================================================================================
void s_Gaussian3D::PointValue(double T, const double* coord, double* uex) const {
    for(int ivar=Var_R; ivar<=Var_P; ivar++) uex[ivar] = 0.0;

    if(NAngularPeriods == 1) {
        for(int iPerX = -MaxPer[0]; iPerX <= MaxPer[0]; iPerX++)
        for(int iPerY = -MaxPer[1]; iPerY <= MaxPer[1]; iPerY++)
        for(int iPerZ = -MaxPer[2]; iPerZ <= MaxPer[2]; iPerZ++) {
            if(Checkerboard) if((iPerX+iPerY+iPerZ)&1) continue;
            if(fabs(PerX) > 0.5*huge && iPerX) continue;
            if(fabs(PerY) > 0.5*huge && iPerY) continue;
            if(fabs(PerZ) > 0.5*huge && iPerZ) continue;

            double r[3]={0.,0.,0.};
            r[0] = coord[0] - r0[0] - FlowVelX*T - PerX * iPerX;
            r[1] = coord[1] - r0[1] - FlowVelY*T - PerY * iPerY;
            r[2] = coord[2] - r0[2] - FlowVelZ*T - PerZ * iPerZ;

            CalcValues(T, r[0], r[1], r[2], uex);
        }
    }
    else {
        double e[3] = {0.,0.,1.};
        for(int iPerPhi = 0; iPerPhi < NAngularPeriods; iPerPhi++)
        for(int iPerZ = -MaxPer[2]; iPerZ <= MaxPer[2]; iPerZ++) {
            double rr0[3] = {r0[0], r0[1], r0[2]};
            double angle = Pi2 / NAngularPeriods * iPerPhi;
            RotateVector(rr0, e, angle, rr0); // Выбираем положение источника по периоду
            double r[3] = {coord[0]-rr0[0], coord[1]-rr0[1], coord[2]-rr0[2]- PerZ * iPerZ};
            CalcValues(T, r[0], r[1], r[2], uex);
        }
    }
    for(int ivar=Var_R; ivar<=Var_P; ivar++) if(IsNaN(uex[ivar])) crash("s_Gaussian3D: NaN detected");
}
//======================================================================================================================



// *********************************************************************************************************************
// *****                                                                                                           *****
// *****                         Entropy and vorticity perturbation with Gaussian profile                          *****
// *****                   (the solution is transported by a background flow without evolution)                    *****
// *****                                                                                                           *****
// *********************************************************************************************************************

//======================================================================================================================
void s_EntropyVortex::ReadParams(tFileBuffer& FB) {
    tParamManager PM;
    PM.Request(PerX, "PerX");
    PM.Request(PerY, "PerY");
    PM.Request(Xterm, "Xterm");
    PM.Request(Yterm, "Yterm");
    PM.Request(Aterm1, "Aterm1");
    PM.Request(Aterm2, "Aterm2");
    PM.Request(Bterm, "Bterm");
    PM.Request(FlowVelX, "FlowVelX");
    PM.Request(FlowVelY, "FlowVelY");
    PM.ReadParamsFromBuffer(FB);
}
//======================================================================================================================


//======================================================================================================================
void s_EntropyVortex::PointValue(double t, const double* coor, double* V) const{
    double x = coor[0] - (FlowVelX * t) - Xterm;
    double y = coor[1] - (FlowVelY * t) - Yterm;
    double alpha = log(2.0) / SQR(Bterm);
    RoundToCentre(x, PerX);
    RoundToCentre(y, PerY);
    double exprr = exp(-alpha * (x*x+y*y));
    V[Var_R] = Aterm1*exprr;
    V[Var_U] = 2 * Aterm2 * (-alpha) * y * exprr;
    V[Var_V] = 2 * Aterm2 * alpha * x * exprr;
    V[Var_W] = V[Var_P] = 0.0;
}
//======================================================================================================================
