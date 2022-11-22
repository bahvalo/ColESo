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
void s_Gaussian2D<fpv>::CalcViaDirect_GaussLaguerre(const fpv& t, const fpv& r, fpv& p, fpv& ur) const {
    p = ur = 0.0;
    #if 0 // use the Gauss - Laguerre quardature formulas
    for(int i=0; i<GLAI.GR; i++) {
        fpv xi = GLAI.GN[i], wi = GLAI.GC[i];
        fpv omega = sqrt(2.*xi);
        p  += wi * BesselJ0(omega*r) * cos(omega*t);
        ur += wi * BesselJ1(omega*r) * sin(omega*t);
    }
    #else // use the Gauss - Legendre quardature formulas
    for(int i=0; i<GLI.GR; i++) {
        fpv omega = GLI.GN[i]*H;
        fpv wi = GLI.GC[i]*exp(-0.5*omega*omega)*omega;
        p  += wi * BesselJ0(omega*r) * cos(omega*t);
        ur += wi * BesselJ1(omega*r) * sin(omega*t);
    }
    p *= H;
    ur *= H;
    #endif
}


template<typename fpv>
void s_Gaussian2D<fpv>::CalcViaBesselFourier_GaussLedengre(const fpv& t, const fpv& r, fpv& p, fpv& ur) const {
    p = ur = 0.0;
    const fpv a = 1.0 - (r+H)/t, b = 1.0;
    for(int i=0; i<GLI.GR; i++) {
        fpv xi = a + (b-a)*GLI.GN[i], wi = GLI.GC[i];
        const fpv txx = t*(1.0-xi); // (1-GLI.GN[i])*(r+H)
        fpv I[2];
        BesselI<fpv>(r*txx, 2, get_eps<fpv>() /*unused*/, false, I); // I0(x)*exp(-x), I1(x)*exp(-x)
        wi *= exp(-0.5*(r-txx)*(r-txx)) * txx / sqrt(xi*(2.-xi));
        p  += wi*((1.-txx*txx)*I[0] + r*txx*I[1]);
        ur += wi*(r*I[0]-txx*I[1]);
    }
    p  *= (b-a)/t;
    ur *= (b-a);
}

template<typename fpv>
void s_Gaussian2D<fpv>::CalcViaFourier_GaussJacobi(const fpv& t, const fpv& r, fpv& p, fpv& ur) const {
    p = ur = 0.0;
    const fpv inv_r = 1.0 / r;
    const fpv b = (t+H)/r - 1.0;
    if(b<0.0) return;

    for(int i=0; i<GJI.GR; i++) {
        fpv xi = b*GJI.GN[i], wi = GJI.GC[i];
        fpv tmp = r - t + r*xi;
        wi *= exp(-0.5*tmp*tmp) / sqrt(2.+xi);
        p  += wi * tmp;
        //ur += wi*(tmp  + 1./(r*(2.+xi)));
        ur += wi*((xi+1.)*tmp + inv_r) / SQR(xi+1.);
    }
    const fpv m = sqrt(b / GetPiNumber2<fpv>());
    p  *= m;
    ur *= m;
}

template<typename fpv>
void s_Gaussian2D<fpv>::CalcViaFourier_uniform(const fpv& t, const fpv& r, fpv& p, fpv& ur) const {
    const NativeDouble Hhat = sqrt(1.16666666*SQR(NativeDouble(H)) + 4.22 + 1.011*log(1.16666666*SQR(NativeDouble(H)) + 4.22));
    const int n = int(Hhat*Hhat / GetPiNumber<NativeDouble>())+1;

    const fpv inv_r = 1.0/r;
    p = ur = 0.0;
    fpv h = Hhat / (n + 0.5);
    if(Hhat+t<r) crash("Internal error");
    for(int i=1; i<=n; i++) {
        fpv eta = i*h;
        fpv U = (t+eta)*inv_r, V = (t-eta)*inv_r;
        fpv u = (U-1.)*(U+1.), v = (V-1.)*(V+1.);
        fpv sqrt_u = sqrt(u),  sqrt_v = sqrt(v);
        fpv w = exp(-0.5*eta*eta) * eta*eta;
        p += w / (u*sqrt_v + v*sqrt_u);
        ur += w / (u*V*sqrt_v + v*U*sqrt_u);
    }
    fpv mult = -4.*t*h*inv_r*inv_r*inv_r / sqrt(GetPiNumber2<fpv>());
    p  *= mult;
    ur *= mult;
}


template<typename fpv>
void s_Gaussian2D<fpv>::CalcAsymptSeries(const fpv& t, const fpv& r, fpv& p, fpv& ur) const {
    fpv eps = get_eps<fpv>();
    const int M_max = int(-log(eps));
    const fpv r2 = r*r;
    const fpv inv_t = fpv(1.0) / t;
    const fpv inv_t2 = inv_t*inv_t;
    eps *= inv_t; // when lg(t)>>1, we want the relative error to be <= epsilon

    {
        const fpv a1 = (0.234375*r2 - 0.75)*r2 + 1.0;
        const fpv a3 = (0.15625*r2 - 0.25)*r2;
        const fpv a5 = 0.015625*r2*r2;
        p = ((-15.*a1*inv_t2 - 3.*a1 + 15.*a3)*inv_t2 -a1 + 3.*a3 - 15.*a5)*inv_t2;
        fpv m = ((-a1*inv_t2 + a3)*inv_t2 - a5)*inv_t2*inv_t2 * 105.0;
        for(int l=4; l<M_max; l++) {
            p += m;
            m *= fpv(2*l+1)*inv_t2;
            if(fabs(m)*M_max < eps) break;
        }
    }
    {
        const fpv b0 = ((0.0390625*r2 - 0.1875)*r2 + 0.5)*r;
        const fpv b2 = ((0.1171875*r2 - 0.375)*r2 + 0.5)*r;
        const fpv b4 = (0.0390625*r2 - 0.0625)*r2*r;
        const fpv b6 = r2*r2*r/fpv(384.0);
        ur = (((15*b0*inv_t2 + 3*b0 - 15*b2)*inv_t2 + b0-3*b2+15*b4)*inv_t2 + b0-b2+3*b4-15*b6)*inv_t;
        fpv m = (((b0*inv_t2 - b2)*inv_t2 + b4)*inv_t2 - b6)*inv_t2*inv_t * 105.0;
        for(int l=4; l<M_max; l++) {
            ur += m;
            m *= fpv(2*l+1)*inv_t2;
            if(fabs(m)*M_max < eps) break;
        }
    }
}

template<typename fpv>
void s_Gaussian2D<fpv>::get_solution(fpv t, fpv r, fpv& p, fpv& ur) const {
    if(GLI.GN==NULL || GJI.GN==NULL || GLAI.GN==NULL) crash("Init not done");

    if(t < get_eps<fpv>()) { // Calculating the solution using the Taylor expansion at t=0
        fpv f0 = exp(-0.5*r*r); // значение
        p = f0; ur = 0.0;
        if(t > 0.0) {
            fpv f1_r = -f0; // first derivative divided by l
            fpv f2 = (r-1.0)*f0; // second derivative
            p += 0.5*t*t*(f2 + f1_r);
            ur = - t * f1_r * r;
        }
        return;
    }
    if(t-r > 1.152*H) {
        if(r > R1) CalcViaFourier_uniform(t,r,p,ur);
        else if(t > 1.3*H) CalcAsymptSeries(t,r,p,ur);
        else CalcViaBesselFourier_GaussLedengre(t,r,p,ur);
    }
    else {
        const fpv alphaH = 1.05*H;
        if(t < r-alphaH) { p=ur=0.0; }
        else if(t+r < alphaH) CalcViaDirect_GaussLaguerre(t,r,p,ur);
        else if(r <= R2) CalcViaBesselFourier_GaussLedengre(t,r,p,ur);
        else CalcViaFourier_GaussJacobi(t,r,p,ur);
    }
}
//======================================================================================================================

//======================================================================================================================
template<typename fpv>
s_Gaussian2D<fpv>::s_Gaussian2D() : tSpaceForm<fpv>() {
    // Default parameters
    FlowVelX=FlowVelY=FlowVelZ=0.0;
    SoundSpeed=1.0;
    // For the correct normalizing of the inital pulse amplitude (if requested)
    tSpaceForm<fpv>::numCoords = 2;

    // Parameters of numerical integration and choosing a method to calculate the solution
    H = sqrt(-2.*log(0.5*NativeDouble(get_eps<fpv>())));
    R1 = pow(15.0*NativeDouble(get_eps<fpv>()), 1./6.);
    R2 = 5.*pow(NativeDouble(get_eps<fpv>()), 0.1);
    gr = num_points_default();
}

template<typename fpv>
void s_Gaussian2D<fpv>::ReadParams(tFileBuffer& FB) {
    // Reading pulse parameters
    tSpaceForm<fpv>::Read(FB, true /* allow periodics */);

    tParamManager PM;
    // Background flow velocity
    PM.Request(FlowVelX, "FlowVelX");
    PM.Request(FlowVelY, "FlowVelY");
    PM.Request(SoundSpeed, "SoundSpeed");
    PM.ReadParamsFromBuffer(FB);
}

template<typename fpv>
void s_Gaussian2D<fpv>::Init(void){
    tSpaceForm<fpv>::Init();
    if(tSpaceForm<fpv>::Form!=tSpaceForm<fpv>::FORM_GAUSSIAN) crash("s_Gaussian2D: non-Gaussian form");

    // Initialization of Gaussian quadrature rules
    if(gr<0) gr = num_points_default();
    GLI.Init(gr, GI_LEGENDRE);
    GJI.Init(gr, GI_JACOBI1);
    GLAI.Init(gr, GI_LAGUERRE);
}

template<typename fpv>
void s_Gaussian2D<fpv>::PointValue(fpv tt, const fpv* coord, fpv* uex) const{
    if(GLI.GN==NULL || GJI.GN==NULL || GLAI.GN==NULL) crash("Init not done");
    for(int ivar=Var_R; ivar<Var_N; ivar++) uex[ivar] = 0.0;
    const fpv mult = sqrt(fpv(2.0)*GetLn2<fpv>()) * tSpaceForm<fpv>::InvBTerm;
    const tSpaceForm<fpv>& SF = *this;

    for(int iPerX = -SF.MaxPer[0]; iPerX <= SF.MaxPer[0]; iPerX++)
    for(int iPerY = -SF.MaxPer[1]; iPerY <= SF.MaxPer[1]; iPerY++) {
        if(SF.Checkerboard) if((iPerX+iPerY)&1) continue;
        // Coordinates, taking into account the background flow and periodical b.c.
        fpv x = coord[0] - SF.r0[0] - FlowVelX*tt + SF.PerX * iPerX;
        fpv y = coord[1] - SF.r0[1] - FlowVelY*tt + SF.PerY * iPerY;
        // Scaling coordinates and time to the case with the initial pulse exp(-0.5*r^2)
        fpv r = sqrt(x*x+y*y) * mult;
        fpv t = tt * mult * SoundSpeed;

        fpv p, ur;
        get_solution(t, r, p, ur);

        uex[Var_R] += p;
        if(r > 1e-100) {
            uex[Var_U] += ur*x*mult/r;
            uex[Var_V] += ur*y*mult/r;
        }
    }
    uex[Var_P] = uex[Var_R];
    for(int ivar=0; ivar<5; ivar++) uex[ivar] *= SF.Aterm;
    uex[1]*=SoundSpeed; uex[2]*=SoundSpeed; uex[3]*=SoundSpeed; uex[4]*=SoundSpeed*SoundSpeed;
}
//======================================================================================================================

// Template classes instantiation
template struct s_Gaussian2D<NativeDouble>;
#ifdef EXTRAPRECISION_COLESO
template struct s_Gaussian2D<dd_real>;
template struct s_Gaussian2D<qd_real>;
#endif


// Verification subroutine
#ifdef EXTRAPRECISION_COLESO
#define real_type1 double //dd_real
#define real_type2 dd_real //qd_real
#define minval 1e-20
#define maxval 1e10
#define threshold 2e-15
void CheckGaussian() {
    s_Gaussian2D<real_type1> S1;
    s_Gaussian2D<real_type2> S2;
    S1.Init();
    S2.Init();

    const double mult = 1.01; // multiplication step
    double max_err = 0.0;
    real_type2 r, t;
    for(r=0.0; r<=maxval; r*=mult) {
        for(t=0.0; t<=maxval; t*=mult) {
            real_type1 p1,ur1; S1.get_solution(double(t),double(r),p1,ur1);
            real_type2 p2,ur2; S2.get_solution(t,r,p2,ur2);
            double err = fabs(double(p1-p2)) + fabs(double(ur1-ur2));
            if(err>max_err) max_err = err;
            if(err>threshold)
                printf("% e % e  % 25.15e % 25.15e\n", double(r), double(t), double(p2), err);
            if(t<0.999*minval) t=minval/mult;
        }
        if(r<0.999*minval) r=minval/mult;
    }
    printf("max_err = %e\n", max_err); // 1.926437e-15
}
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
    PM.Request(SoundSpeed, "SoundSpeed");
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

            CalcValues(T*SoundSpeed, r[0], r[1], r[2], uex);
        }
    }
    else {
        double e[3] = {0.,0.,1.};
        for(int iPerPhi = 0; iPerPhi < NAngularPeriods; iPerPhi++)
        for(int iPerZ = -MaxPer[2]; iPerZ <= MaxPer[2]; iPerZ++) {
            double rr0[3] = {r0[0], r0[1], r0[2]};
            double angle = Pi2 / NAngularPeriods * iPerPhi;
            RotateVector(rr0, e, angle, rr0); // Выбираем положение источника по периоду
            double r[3] = {coord[0]-rr0[0]-FlowVelX*T, coord[1]-rr0[1]-FlowVelY*T, coord[2]-rr0[2]-FlowVelZ*T- PerZ * iPerZ};
            CalcValues(T*SoundSpeed, r[0], r[1], r[2], uex);
        }
    }
    uex[1]*=SoundSpeed; uex[2]*=SoundSpeed; uex[3]*=SoundSpeed; uex[4]*=SoundSpeed*SoundSpeed;
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
