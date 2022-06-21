// *********************************************************************************************************************
// *****                                                                                                           *****
// *****                                            Special functions                                              *****
// *****                                                                                                           *****
// *********************************************************************************************************************

#include "es_specfunc.h"
#include "es_utils.h"
#include "base_const.h"
#ifdef EXTRAPRECISION_COLESO
#include "extraprecision.h"
#endif
#include <math.h>
#include <complex>

//======================================================================================================================
// Bessel functions J_0(x), J_1(x), N_0(x), N_1(x)
//======================================================================================================================
// Double-precision
NativeDouble cephes_j0(NativeDouble);
NativeDouble cephes_j1(NativeDouble);
NativeDouble cephes_y0(NativeDouble);
NativeDouble cephes_y1(NativeDouble);
template<> NativeDouble BesselJ0(NativeDouble x) { return cephes_j0(x); }
template<> NativeDouble BesselJ1(NativeDouble x) { return cephes_j1(x); }
template<> NativeDouble BesselN0(NativeDouble x) { return cephes_y0(x); }
template<> NativeDouble BesselN1(NativeDouble x) { return cephes_y1(x); }

// Double-double-precision
// ** Defined in es_specfunc_dd.cpp **

// Quad-double-precision
#ifdef EXTRAPRECISION_COLESO
template<> qd_real BesselJ0(qd_real x) {
    qd_real re=0.0, im=0.0;
    BesselJComplex<qd_real>(0, x, 0.0, &re, &im);
    return re;
}
template<> qd_real BesselJ1(qd_real x) {
    qd_real re=0.0, im=0.0;
    BesselJComplex<qd_real>(1, x, 0.0, &re, &im);
    return re;
}
template<> qd_real BesselN0(qd_real) { crash("BesselN0<qd_real> not implemented"); }
template<> qd_real BesselN1(qd_real) { crash("BesselN1<qd_real> not implemented"); }
#endif


//======================================================================================================================
// Calculation of Bessel functions J_n(x) for n=0,...,nmax-1
// ATTENTION! Output array J must be of size S=nmax+MIN(nmax,65)
// Buffer array 'buf' must be of size S; if buf==NULL, we assume that J is of size 2*S
//======================================================================================================================
void BesselFunctionsJ(NativeDouble x, int nmax, NativeDouble *J, NativeDouble* buf) {
    if(buf==NULL) buf = J + nmax*2;
    if(x<0 || IsNaN(x)) crash("BesselFunctions: x = %e", x);
    J[0] = BesselJ0(x);
    if(nmax==1) return;
    J[1] = BesselJ1(x);
    if(nmax==2) return;

    // Forward recursion when applicable: J[2], ..., J[n1]
    const int n1 = MIN(nmax-1, int(x));
    const NativeDouble two_over_x = 2./x;
    for(int i=1; i<n1; i++) J[i+1] = i*two_over_x*J[i] - J[i-1];
    if(n1 == nmax-1) return;

    // Estimating
    int n2 = (nmax>=65) ? nmax+64 : nmax*2-1;
    J[n2] = 0.0;
    if(nmax<=32) {
        // Calculating J[n2] using Taylor series at x=0
        const double r1 = x * .5;
        const double r2 = -0.25*x*x;
        double u = 1.0;
        for(int i=n2; i>=1; --i) u *= r1 / NativeDouble(i);
        double p = u;
        int k = (x >= 10.0) ? 2*int(x) : 20; // note that x<nmax-1
        for (int i = 1; i <= k; ++i) {
            u *= r2 / (NativeDouble) (i*(n2+i));
            p += u;
        }
        J[n2] = p;
    }

    // Solving recurrent relations as three-diagonal system with the boundary conditions J[n1] and J[n2]
    buf[n1] = 0.0;
    for(int i=n1+1; i<n2; ++i) {
        buf[i] = 1.0 / (two_over_x * NativeDouble(i) - buf[i-1]);
        J[i] = buf[i] * J[i-1];
    }

    for(int i=n2-1; i>n1; i--) J[i] += buf[i] * J[i+1];
}
//======================================================================================================================

//======================================================================================================================
// Calculation of Bessel functions J_n(x) for n=0,...,nmax-1
// Output array Y must be of size nmax
// Returns n+1 where n is the maximal index such that Y_n(x) is successfully calculated
//======================================================================================================================
int BesselFunctionsY(NativeDouble x, int nmax, NativeDouble *Y) {
    if(x<0 || IsNaN(x)) crash("BesselFunctionsY: x = %e", x);
    // for x<1e-100, we will not calculate Y
    if(x<1e-100) { for(int i=0; i<nmax; i++) Y[i]=-1e100; return 0; }

    Y[0] = BesselN0(x);
    if(nmax==1) return nmax;
    Y[1] = BesselN1(x);
    if(nmax==2) return nmax;

    // Forward recursion
    const NativeDouble two_over_x = 2./x;
    const NativeDouble critical_value = (x<2.0) ? 1e300*x : 2e300;
    for(int i=1; i<nmax-1; i++) {
        Y[i+1] = i*two_over_x*Y[i] - Y[i-1];
        if(fabs(Y[i+1]) > critical_value) {
            for(int j=i+2; j<nmax; j++) Y[i]=-critical_value;
            return i+2;
        }
    }
    return nmax;
}
//======================================================================================================================



//======================================================================================================================
// Radial part of eigenfunctions of the Laplace operator in a circle exterior (r > y)
// with the Neumann conditions at r = rc
// These functions are given by  u_{nu,k}(r) = phi(nu, k*r, k*y), r >= y,
// phi(nu,x,y) = (- N'_nu(y) J_nu(x) + J'_nu(y) N_nu(x)) / sqrt((J'_nu(y))^2 + (N'_nu(y))^2), x >= y.
// This subroutine calculates the values of phi(nu,x,y) and dphi/dx(nu,x,y) for given x,y and each nu=0,...,Nmax-1
//======================================================================================================================
void BesselOuterFunction(NativeDouble x, NativeDouble y, int Nmax, NativeDouble* phi, NativeDouble* dphi) {
    if(x<y) crash("BesselOuterFunction error: x < y");
    if(phi==NULL && dphi==NULL) crash("BesselOuterFunction error: output arrays are not allocated");

    if(x < tiny) { // x=y=0. No inner radius and calculating values at x=0
        if(phi!=NULL)  { for(int i=0; i<Nmax; i++) phi[i] = 0.0;  phi[0] = 1.0; }
        if(dphi!=NULL) { for(int i=0; i<Nmax; i++) dphi[i] = 0.0; if(Nmax>=2) dphi[1] = 1.0; }
        return;
    }

    // Allocating memory and set pointers
    const int bufsize = (Nmax+1)*2;
    NativeDouble* buf = new NativeDouble[bufsize*4 + (Nmax+1)*3];
    NativeDouble* Jx = buf;
    NativeDouble* Nx = buf+bufsize;
    NativeDouble* Jy = buf+bufsize*2;
    NativeDouble* Ny = buf+bufsize*3;

    // Calculating Bessel functions of the first kind (up to n=Nmax inclusively) and second kind (as many as we can)
    BesselFunctionsJ(x, Nmax+1, Jx, Nx); // use Nx as buffer
    int numfuncx = BesselFunctionsY(x, Nmax+1, Nx);

    if(y < tiny) { // no inner radius => calculating Bessel functions of the first kind
        for(int i=0; i<Nmax; i++) {
            if(phi!=NULL)  phi[i] = Jx[i];
            if(dphi!=NULL) dphi[i] = -Jx[i+1] + i*Jx[i]/x;
        }
        delete[] buf;
        return;
    }

    BesselFunctionsJ(y, Nmax+1, Jy, Ny); // use Ny as buffer
    int numfunc = BesselFunctionsY(y, Nmax+1, Ny);
    // Returned values (numfunc*) are the number of BesselY functions actually calculated
    // Internal critical error
    if(numfuncx<=0 || numfunc<=0) crash("BesselOuterFunction error: returned %i, %i values", numfuncx, numfunc);

    // так как y<=x, то количество функций, влезающих в арифметику, для y не больше, чем для x. На всякий случай эта проверка
    if(numfunc > numfuncx) numfunc = numfuncx;

    // Вычисляем искомую функцию в области, где функции Бесселя 2-го рода вычислены корректно
    for(int i=0; i<numfunc-1; i++) {
        NativeDouble jx = Jx[i];
        NativeDouble djx = -Jx[i+1] + i*Jx[i]/x;
        NativeDouble djy = -Jy[i+1] + i*Jy[i]/y;
        NativeDouble nx = Nx[i];
        NativeDouble dnx = -Nx[i+1] + i*Nx[i]/x;
        NativeDouble dny = -Ny[i+1] + i*Ny[i]/y;
        NativeDouble num =  -dny*jx + djy*nx;
        NativeDouble dnum = -dny*djx + djy*dnx;
        NativeDouble denom = fabs(dny);
        // если |N_n'(y)| > 1e20, то ф-ей 1-го рода пренебрегаем. Тем самым избегаем overflow от возведения в квадрат
        if(denom < 1e20) denom = sqrt(djy*djy + dny*dny);
        denom = 1.0/denom;
        if(phi)  phi[i] = num * denom;
        if(dphi) dphi[i] = dnum * denom;
    }
    if(numfunc==Nmax+1) { delete[] buf; return; } // all done

    // Set pointers
    NativeDouble* Ax = buf+bufsize*4;
    NativeDouble* Ay = buf+bufsize*4+(Nmax+1);
    NativeDouble* B  = buf+bufsize*4+(Nmax+1)*2;

    // Calculating A_n(y) = N_n(y) / N_{n-1}(y)
    // For n=numfunc-1 do this directly
    Ay[numfunc-1] = Ny[numfunc-1] / Ny[numfunc-2];
    // For n>numfunc-1 do this recursively
    for(int i=numfunc-1; i<Nmax; i++)
        Ay[i+1] = 2*i/y - 1.0/Ay[i];

    // Calculating A_n(x) = N_n(x) / N_{n-1}(x)
    // For n=numfunc-1,...,numfuncx-1 do this directly
    for(int i=numfunc-1; i<=numfuncx-1; i++) {
        if(fabs(Nx[i-1]) < 1) Ax[i] = 1.23456789e100; // this will not be used
        else Ax[i] = Nx[i] / Nx[i-1];
    }
    // For n>numfuncx-1 do this recursively
    for(int i=numfuncx-2; i<Nmax; i++)
        Ax[i+1] = 2*i/x - 1.0/Ax[i];

    // Calculating B_n(x,y) = N_n(x) / N_n(y)
    // For n=numfunc-1 do this directly. Note: here Nx and Ny are correct and |Ny| >> 1
    B[numfunc-1] = Nx[numfunc-1] / Ny[numfunc-1];
    // For n>numfunc-1 do this recursively
    for(int i=numfunc; i<Nmax; i++) {
        if(fabs(Nx[i]) < 1e100) B[i] = 0.0; // since i>=numfunc, there holds |Ny|>1e250. If |Nx| is small enough, ratio is ~0
        else B[i] = Ax[i]*B[i-1]/Ay[i];
    }

    // Calculating the functions we need for n>=numfunc
    for(int i=numfunc-1; i<Nmax; i++) {
        NativeDouble jx = Jx[i], djx = -Jx[i+1] + i*Jx[i]/x;
        NativeDouble djy = -Jy[i+1] + i*Jy[i]/y;
        if(phi) phi[i] = -jx +  djy*B[i]*(i/y - Ay[i+1]);
        if(dphi) dphi[i] = -djx + djy*(y/x)*B[i]*(i - x*Ax[i+1]) / (i - y*Ay[i+1]);
    }
    delete[] buf;
}
//======================================================================================================================



//======================================================================================================================
// Bessel function of the first kind and index l at a point x
// Also calculates its derivatives up to the order 5 (if the corresponding pointers are nonzero)
//======================================================================================================================
NativeDouble BesselJ(int l, NativeDouble x, NativeDouble* dbdx, NativeDouble* d2bdx2,
                     NativeDouble* d3bdx3, NativeDouble* d4bdx4, NativeDouble* d5bdx5) {
    if(x<tiny) {
        if(dbdx) {
            if(l==1) *dbdx = 0.5;
            else *dbdx = 0.0;
        }
        if(d2bdx2) {
            if(l==0) *d2bdx2 = -0.5;
            else if(l==2) *d2bdx2 = 0.25;
            else *d2bdx2 = 0.0;
        }
        if(d3bdx3) {
            if(l==1) *d3bdx3 = -0.375;
            else if(l==3) *d3bdx3 = 0.125;
            else *d3bdx3 = 0.0;
        }
        if(d4bdx4) {
            if(l==0) *d4bdx4 = 0.375;
            else if(l==2) *d4bdx4 = -0.25;
            else if(l==4) *d4bdx4 = 0.0625;
            else *d4bdx4 = 0.0;
        }
        if(d5bdx5) {
            if(l==1) *d5bdx5 = 0.3125;
            else if(l==3) *d5bdx5 = -0.15625; // -5/32
            else if(l==5) *d5bdx5 = -0.03125; // -1/32
            else *d5bdx5 = 0.0;
        }
        return (l==0) ? 1.0 : 0.0;
    }

    NativeDouble* BesselJ = new NativeDouble[(l+2)*4];
    BesselFunctionsJ(x, l+2, BesselJ);
    NativeDouble bj = BesselJ[l];
    if(dbdx) *dbdx = -BesselJ[l+1] + l*BesselJ[l]/x;
    if(d2bdx2) *d2bdx2 = BesselJ[l+1]/x + BesselJ[l]*(-1.0 + l*(l-1)/(x*x));
    if(d3bdx3) *d3bdx3 = BesselJ[l+1]*(1.0-(2+l*l)/(x*x)) + BesselJ[l]*((1-l)/x + (l*l*l-3*l*l+2*l)/(x*x*x));
    if(d4bdx4) *d4bdx4 = BesselJ[l+1]*(-2+6*(1+l*l)/(x*x))/x + BesselJ[l]*(1 + (-2*l*l-3+2*l)/(x*x) + (-6*l+11*l*l-6*l*l*l+l*l*l*l)/(x*x*x*x));
    if(d5bdx5) *d5bdx5 =-BesselJ[l+1]*(1.0+(-7-2*l*l)/(x*x)+(35*l*l+24+l*l*l*l)/(x*x*x*x))
        - BesselJ[l]*((-l+2)/x + (-12*l*l+7*l+2*l*l*l-12)/(x*x*x) + (50*l*l-l*l*l*l*l-24*l-35*l*l*l+10*l*l*l*l)/(x*x*x*x*x));
    delete[] BesselJ;
    return bj;
}
//======================================================================================================================


//======================================================================================================================
// Printing for test purpose
//======================================================================================================================
void PrintBesselDerivatives(int l) {
    FILE* J0 = fopen("besselderivs.dat", "wt");
    for(int i = 0; i<100000; i++) {
        NativeDouble x = NativeDouble(i) / 100;
        NativeDouble bj, dbj, d2bj, d3bj, d4bj, d5bj;
        bj = BesselJ(l, x, &dbj, &d2bj, &d3bj, &d4bj, &d5bj);
        fprintf(J0, "%f %e %e %e %e %e %e\n", x, bj, dbj, d2bj, d3bj, d4bj, d5bj);
    }
    fclose(J0);
}
//======================================================================================================================


//======================================================================================================================
// k-th zero of J'_n(x), where k=RadialMode, n=AngularMode
// To this moment, only n=0,...,8 are supported
// It was checked that zeroes coincides with the ones given by Maple. The Maple code is the following:
// f := proc (m, n) options operator, arrow; fsolve(BesselJ(m+1, x)*x-m*BesselJ(m, x), x = BesselJZeros(m, n) .. BesselJZeros(m, n+1)) end proc
// for i to 20 do f(1, i) end do
//======================================================================================================================
NativeDouble BesselPrimeZero(int AngularMode, int RadialMode, int log) {
    // Define an initial value for the iterative process
    // First calculate an asymptotic expression for large values of RadialMode
    NativeDouble Nu = PiNumber * (RadialMode + 0.25 + 0.5*AngularMode);

    if(AngularMode > 50) crash("BesselPrimeZero: not implemented");
    if(AngularMode > 8) {
        if(RadialMode>50) crash("BesselPrimeZero: not implemented");
        if(RadialMode==0) {
            double a = pow(NativeDouble(AngularMode), NativeDouble(C1_3));
            // Abramovitz, Stegun, 9.5.16
            Nu = AngularMode + 0.8086165*a + 0.07249/a - 0.05097/AngularMode + 0.0094*a/(AngularMode*AngularMode);
        }
        else {
            // Recursive formula by J.-P. Moreau
            // http://jean-pierre.moreau.pagesperso-orange.fr/Cplus/mjyzo_cpp.txt
            Nu = BesselPrimeZero(AngularMode, RadialMode-1, log);
            Nu += 3.1416+(0.4955+0.0915*AngularMode-0.000435*AngularMode*AngularMode)/RadialMode;
        }
    }
    else {
        // For some specific values, the Newton process starting from this value converges to another root,
        // so we modify the initial value
        if(RadialMode==0) {
            if(AngularMode==0) Nu = 0.0;
            else Nu = AngularMode * 1.05 + 1.0;
        }
        if(RadialMode==1) {
            if(AngularMode==0) Nu = 3.8317;
            else Nu = -0.001*AngularMode*AngularMode + AngularMode * 1.15 + 5.0;
        }
        if(RadialMode==2 && AngularMode==7) Nu = 16.52936588;
        if(RadialMode==3 && AngularMode==7) Nu = 19.94185337;
        if(RadialMode==2 && AngularMode==8) Nu = 17.77401;
        if(RadialMode==3 && AngularMode==8) Nu = 21.22906;
        if(RadialMode==4 && AngularMode==8) Nu = 24.58720;
    }

    // Precise the zero using Newton's method
    if(log) pprintf("Bessel prime zero: %23.16e; ", Nu);

    unsigned int maxit = 20;
    for(unsigned int i=0; i<maxit; i++) {
        NativeDouble dfdx = 0.0, d2fdx2 = 0.0;
        BesselJ(AngularMode, Nu, &dfdx, &d2fdx2, NULL, NULL, NULL); // this is a double-precision subroutine
        double dNu = dfdx / d2fdx2;
        Nu -= NativeDouble(dNu); // Precise the zero
        // Zero is precise enough, so we will make not more than 2 iterations (4 for DD, 8 for QD)
        if(fabs(dNu) < 1e-14) maxit = MIN(maxit, i+1+(sizeof(NativeDouble)>>2));
    }

    if(log) pprintf("precised: %23.16e\n", Nu);

    return Nu;
}
//======================================================================================================================


//======================================================================================================================
// Auxiliary subroutine for calculating Bessel functions of complex argument
// Input:
// z      = argument
// N      = desirable number of functions
// eps    = desirable accuracy (for example, 1e-16)
// Output: 
// n      = the index to start the backward recursion from
// ncalc  = maximal index with reliable value of J_n(x)
//======================================================================================================================
void InitBackwardRecursion(const std::complex<NativeDouble>& z, int N, NativeDouble eps, int& ncalc, int& n) {
    const NativeDouble inv_eps = 1.0/eps;
    const NativeDouble magZ = abs(z); // |z|
    const std::complex<NativeDouble> zinv = 1.0 / z;

    // 1e306 is near the largest value representable in double precision
    const NativeDouble very_huge = 1e306 * eps;

    // Let p_n follow the same recursion relations as Bessel functions, 
    // with the initial conditions p_{n-1}=1 and p_n=2*n/z, where n = [|z|]+1
    std::complex<NativeDouble> p, pold, plast;
    n = int(magZ) + 1;
    plast = 1.;
    p = NativeDouble(n*2) * zinv;

    int overflow_flag = 0;
    // If (N > [|z|]+2), calculate p_n for n=[|z|]+2, ..., N-1. Check for possible overflow
    for (++n; n < N; ++n) {
        pold = plast;
        plast = p;
        p = NativeDouble(n*2)*plast*zinv - pold;
        if (abs(p) >= very_huge) { overflow_flag = 1; break; }
    }
    n--;

    // Main mode
    if(!overflow_flag) {
        // Find 'n' such that |p_n|>=test
        NativeDouble threshold = sqrt(2.*inv_eps) * MAX(abs(p), abs(plast));
        threshold = MAX(threshold,inv_eps);

        ncalc = N; // values for all n=0,...,N-1 will be calculated
        // Calculate p_n until |p_n|>=test
        int m = 0;
        while(1) {
            ++n;
            pold = plast;
            plast = p;
            p = NativeDouble(n*2) * plast*zinv - pold;
            if (norm(p/threshold) < 1.) continue;
            if (m == 1) return;

            // We increase 'test' and check the condition again
            m = 1;
            NativeDouble w = abs(p) / abs(plast);
            NativeDouble q = NativeDouble(n + 1) / magZ;
            if (w + 1./w > q*2.) w = q + sqrt(q*q - 1.);
            threshold /= sqrt(w - 1./w);
            if (norm(p/threshold) >= 1.) return;
        }
    }
    // Now consider the case when  |J_n(x)|=inf  for some n<N
    else {
        ncalc = n + 1; // save
        // Here 'n' is the index where overflow occurred
        // Above we compared |p| with 'very_huge', so we can slightly increase 'n'
        // To do this, divide all data by very_huge and continue the recursion
        p /= very_huge;
        plast /= very_huge;
        std::complex<NativeDouble> psave = p;
        std::complex<NativeDouble> plastsave = plast;

        do {
            ++n;
            pold = plast;
            plast = p;
            p = NativeDouble(2*n) * plast*zinv - pold;
        } while (abs(p) <= inv_eps);
        // Last value of unnormalized 'p' is equal to +inf, so returning to the previous one
        --n;

        // Calculate the threshold
        NativeDouble threshold4; // threshold^4
        {
            NativeDouble w = sqrt(norm(plast) / norm(pold));
            NativeDouble q = NativeDouble(n+1) / magZ;
            if (w + 1./w > q*2.) w = q + sqrt(q*q-1.);
            NativeDouble buf = (1. - 1./(w*w)) * .5 * eps;
            threshold4 = (norm(plast) * buf) * (norm(pold) * buf);
        }
        // Returning to saved values and continue the recursion to check the number of functions
        // we calculated with the requested accuracy
        plast = plastsave;
        p = psave;
        for (; ncalc<=N && ncalc<=n; ++ncalc) {
            pold = plast;
            plast = p;
            p = NativeDouble(n*2) * plast*zinv - pold;
            if (norm(p)*norm(plast) > threshold4) break;
        }
        --ncalc;
    }
}
//======================================================================================================================

// exp(z) is double-precision; we need higher-precision function
template <typename fpv> std::complex<fpv> A_exp(std::complex<fpv> z);

//======================================================================================================================
// This routine calculates Bessel functions J_n of complex argument and integer order using backward recursion method
// ATTENTION: Not suitable for large arguments
// Input: z = argument; N = 1 + highest order to be calculated
//        n_extra_steps = number of additional steps for the backward recursion (for test purpose only)
// Output: J = \{J_0 .. J_{ncalc-1}\}
// Returning value: ncalc (number of functions actually calculated). The user should check that ncalc=N
//======================================================================================================================
template<class fpv>
int BesselFunctionsComplex(std::complex<fpv> z, int N, std::complex<fpv> *J, int n_extra_steps) {
    if(norm(z) > 1e8 || fabs(z.imag()) > 700.0) { // exp(700) ~ 1.01e304
        for(int i=0; i<N; i++) J[i]=0.0;
        return 0; // error
    }

    // Use two-term series for small z
    const NativeDouble SQRT_EPS = (sizeof(fpv)==8) ? 1e-8 : ((sizeof(fpv)==16) ? 1e-16 : 1e-32);
    if (norm(z) < SQRT_EPS) {
        std::complex<fpv> X22 = fpv(0.25)*z*z;
        std::complex<fpv> Xi = std::complex<fpv>(1.0); // Xi = (z/2)^n / n!
        J[0] = std::complex<fpv>(1.0) - X22;
        for (int n = 1; n < N; ++n) {
            Xi *= z / fpv(n*2);
            J[n] = Xi - Xi*X22 / fpv(n+1);
        }
        return N; // all values have been calculated
    }


    // During the recursion, we will calculate the 0.5-multiplied right-hand side of the identity
    // exp(i*z) = J_0(z) + 2 * \sum_{k=1}^{\infty} i^k J_k(z)
    // (see Abramowitz, Stegun, 9.1.47 and 9.1.48)
    // If (ImZ > 0), the signs of the terms in the sum alternates, which may lead to a huge relative error
    // Thus if ImZ>0, remember this and put ImZ=-ImZ
    const int sign_y = (z.imag()>=0.0) ? 1 : -1;
    if(sign_y==1) z = std::complex<fpv>(z.real(), -z.imag());

    std::complex<fpv> invz = std::complex<fpv>(1.0) / z;

    std::complex<fpv> sum = fpv(0.);
    // Macro for the following operation:
    // sum += i^j * V
    #define INCREMENT(j, V) \
    do { \
        if(!(j)) sum += fpv(0.5)*(V); \
        else switch((j)&3) { \
        case 0: sum += V; break; \
        case 1: sum += std::complex<fpv>(-(V).imag(), (V).real()); break; \
        case 2: sum -= V; break; \
        case 3: sum -= std::complex<fpv>(-(V).imag(), (V).real()); break; \
        } \
    } while(0)

    // Now we find 'n' (the index we start the backward recursion from) 
    // and 'ncalc' (maximal index with reliable value of J_n(x))
    // For the case if higher-precision arithmetic model is used, 
    // first argument is converted to double-precision
    int n, ncalc;
    InitBackwardRecursion(std::complex<NativeDouble>(NativeDouble(z.real()), NativeDouble(z.imag())),
        N, SQRT_EPS*SQRT_EPS, ncalc, n);
    n += n_extra_steps; // for test purpose only

    // If initial function index is less than the maximal index of function to calculate (i.e. there holds n < N-1),
    // then we store the initial values and assign zero to all following functions we're unable to calculate
    for (int l = n+1; l < N; ++l) J[l] = 0.;

    // Normalization factor (if needed)
    const fpv very_small = 1e-75;
    const fpv very_huge = fpv(1.0)/(very_small*very_small);

    // If initial function index is greater than the maximal index of function to calculate (i.e. n > N-1),
    // then we start the backward recursion not storing the results
    // Even if n = N-1, do one iteration of this form
    std::complex<fpv> C_n_plus_1(0., 0.); // C_{n+1} = 0
    std::complex<fpv> C_n(1., 0.); // C_n = 1
    for (; n >= N-1 && n >= 1; --n) {
        INCREMENT(n+1, C_n_plus_1); // add the term to the sum
        std::complex<fpv> C_n_plus_2 = C_n_plus_1;
        C_n_plus_1 = C_n;
        // Calculating C_{n-1}, then --n will be called before exiting the loop
        C_n = fpv(n*2)*invz*C_n_plus_1 - C_n_plus_2;

        if(norm(sum) > very_huge) {
            C_n *= very_small;
            C_n_plus_1 *= very_small;
            sum *= very_small;
        }
    }

    // Now we have two first terms, namely, n-th and (n+1)-th, in the recursion. Add them to the sum
    // if N=1, then the array size is just one
    if(N>1) J[n+1] = C_n_plus_1;
    J[n] = C_n;
    INCREMENT(n+1, C_n_plus_1);
    INCREMENT(n, C_n);

    // Proceed with the recursion till n=0
    for(; n>=1; --n) {
        J[n-1] = fpv(n*2)*invz*J[n] - J[n+1];
        INCREMENT(n-1, J[n-1]);
    }
    sum *= fpv(2.0);

    // Now we need to multiply all the values by exp(i*z)/sum
    std::complex<fpv> a = A_exp(std::complex<fpv>(fpv(0.0),fpv(1.0)) * z) / sum;
    for (n=0; n<N; ++n) J[n] *= a;

    if(sign_y==1) for (n=0; n<N; ++n) J[n] = std::complex<fpv>(J[n].real(), -J[n].imag());
    return ncalc;
}
//======================================================================================================================

//======================================================================================================================
// Debug subroutine: calculate the relative difference 
// between the values of J_n(z) obtained with different n_extra_steps
//======================================================================================================================
template<class fpv>
NativeDouble RelativeErr(int n, std::complex<double> z, int n_extra) {
    std::complex<double>* b = new std::complex<double>[n+1];
    int ncalc = BesselFunctionsComplex(z, n+1, b, 0);
    if(ncalc != n+1) crash("BesselJComplex(%i, %e, %e) failed", n, double(z.real()), double(z.imag()));
    std::complex<double> V1 = b[n];
    ncalc = BesselFunctionsComplex(z, n+1, b, n_extra);
    if(ncalc != n+1) crash("BesselJComplex(%i, %e, %e) failed", n, double(z.real()), double(z.imag()));
    std::complex<double> V2 = b[n];
    if(0) {
        printf("J[%5i] = % .16e +i* % .16e\n", n, V1.real(), V1.imag());
        printf("         = % .16e +i* % .16e\n", V2.real(), V2.imag());
    }
    delete[] b;
    return NativeDouble(abs(V1-V2) / (abs(V1) + abs(V2) + 1e-100));
}
//======================================================================================================================


//======================================================================================================================
// Compute the Bessel function of the first kind and index 'l' at z=x+iy
// Also computes its derivative, if the corresponding pointers are nonzero
//======================================================================================================================
template<typename fpv>
void BesselJComplex(int l, fpv x, fpv y,
             fpv* ReBesselJ, fpv* ImBesselJ, fpv* ReDBesselJ, fpv* ImDBesselJ) {

    // Заплатка для больших действительных (или близких к действительной оси) z
    // Оставляем один член разложения; из-за этого на действительной оси ошибка вычисления ~ 1e-7
    if(fabs(y) < 50 && fabs(x) > 9000.0) {
        std::complex<double> z; z = std::complex<double>(double(x), double(y));
        double pi4 = 0.25*GetPiNumber<double>(); // pi/4
        std::complex<double> bj = sqrt(double(1.0) / (double(2.0)*pi4*z)) * cos(z - double(2*l+1)*pi4);
        *ReBesselJ = bj.real();
        *ImBesselJ = bj.imag();
        std::complex<double> dbj = - sqrt(double(1.0) / (double(2.0)*pi4*z)) * sin(z - double(2*l+1)*pi4);
        if(ReDBesselJ) *ReDBesselJ = dbj.real();
        if(ImDBesselJ) *ImDBesselJ = dbj.imag();
        return;
    }

    int nb = l+2;
    std::complex<fpv> z(x, y);
    std::complex<fpv>* b = GimmeMem< std::complex<fpv> >(nb);
    int ncalc = BesselFunctionsComplex(z, nb, b, 0);
    if(ncalc != nb) crash("BesselJComplex(%i, %e, %e) failed", l, double(x), double(y));
    *ReBesselJ = b[l].real();
    *ImBesselJ = b[l].imag();

    if(ReDBesselJ || ImDBesselJ) {
        if(fabs(x)<get_eps<fpv>() && fabs(y)<get_eps<fpv>()) {
            if(ReDBesselJ) *ReDBesselJ = (l==1) ? 0.5 : 0.0;
            if(ImDBesselJ) *ImDBesselJ = 0.0;
        }
        else {
            std::complex<fpv> dbj = - b[l+1] + fpv(l) * b[l] / z;
            if(ReDBesselJ) *ReDBesselJ = dbj.real();
            if(ImDBesselJ) *ImDBesselJ = dbj.imag();
        }
    }
    FreeMem(b);
}
//======================================================================================================================

template void BesselJComplex<NativeDouble>(int l, NativeDouble x, NativeDouble y, NativeDouble* ReBesselJ, NativeDouble* ImBesselJ, NativeDouble* ReDBesselJ, NativeDouble* ImDBesselJ);

#ifdef EXTRAPRECISION_COLESO
template void BesselJComplex<dd_real>(int l, dd_real x, dd_real y, dd_real* ReBesselJ, dd_real* ImBesselJ, dd_real* ReDBesselJ, dd_real* ImDBesselJ);
template void BesselJComplex<qd_real>(int l, qd_real x, qd_real y, qd_real* ReBesselJ, qd_real* ImBesselJ, qd_real* ReDBesselJ, qd_real* ImDBesselJ);
#endif


//======================================================================================================================
// Calculation of modified Bessel functions I_n(x), n = 0, ..., nmax-1
// Backward recursion is used: I_{p-1}(x) = (2p/x)*I_p(x) + I_{p+1}(x), starting with I_m=1, I_{m+1}=0.
// We increase 'm' unless this stops affecting the result (up to 'eps')
// Input: x, nmax, eps, s
// If s==1, then we calculate Ip(x), and if s==0, then we calculate exp(-x) Ip(x)
// Output: t (of size nmax)
//======================================================================================================================
template<typename fpv>
void BesselI(fpv x, int nmax, fpv eps, int s, fpv *t) {
    if(nmax<=0) return;
    const fpv y = fabs(x);
    if (y < 1e-150) {
        t[0] = 1.;
        for (int p = 1; p < nmax; ++p) t[p] = 0.;
        return;
    }
    const fpv two_over_x = 2.0/x;

    const fpv inv_sum_exact = s ? exp(-x) : fpv(1.);
    if (inv_sum_exact < 1e-100) crash("BesselI error: argument is too much");

    int d; // index increment
    int m; // index to start the recursion from
    if(y <= 1.) { m = 10; d = 3; }
    else if(y <= 10.) { m = int(y * 1.5) + 10; d = int(y * .5) + 3; }
    else { m = int(y * .5) + 25; d = int(y * .25 + 100. / y); }
    if (m < nmax+d) m = nmax+d;

    fpv saved_value = 0.0;
    while(1) {
        // During the recursion, we will calculate the right-hand side of the identity
        // exp(z) = I_0(z) + 2 * \sum_{k=1}^{\infty} I_k(z)
        fpv sum = 0.;
        fpv r0 = 0.;
        fpv r1 = 1.;
        fpv r2;
        const fpv very_huge = 1e100;

        for (int n=m-1; n>=0; n--) {
            r2 = two_over_x*(n+1)*r1 + r0;
            if (n < nmax) t[n] = r2;

            // All the data is up to a multiplicative constant,
            // so we normalize data when necessary to avoid FP overflow
            if(abs(r2) > very_huge) {
                const fpv inv_very_huge = 1.0 / very_huge;
                r1 *= inv_very_huge;
                r2 *= inv_very_huge;
                sum *= inv_very_huge;
                for (int k=n; k<nmax; k++) t[k] *= inv_very_huge;
            }

            if (n != 0) sum += r2 * 2.;
            else sum += r2;

            r0 = r1;
            r1 = r2;
        }
        sum *= inv_sum_exact;
        sum = 1.0 / sum;
        for (int n=0; n<nmax; n++) t[n] *= sum;

        // Convergence?
        if(fabs(t[nmax-1] - saved_value) < eps * fabs(t[nmax-1])) return;

        // If not, save the result and increase m
        saved_value = t[nmax-1];
        m += d;
    }
}
//======================================================================================================================

// Instantiation of the modified Bessel functions
template void BesselI<NativeDouble>(NativeDouble x, int nmax, NativeDouble eps, int s, NativeDouble *t);
#ifdef EXTRAPRECISION_COLESO
template void BesselI<dd_real>(dd_real x, int nmax, dd_real eps, int s, dd_real *t);
template void BesselI<qd_real>(qd_real x, int nmax, qd_real eps, int s, qd_real *t);
#endif
