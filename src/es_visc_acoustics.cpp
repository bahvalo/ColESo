// *********************************************************************************************************************
// *****                                                                                                           *****
// *****                                  Collection of exact solutions (ColESo)                                   *****
// ****                Acoustic wave in planar or circular channel with viscosity and heat transfer                *****
// *****                                                                                                           *****
// *********************************************************************************************************************

#include "base_parser.h"
#include "coleso.h"
#include "es_utils.h"
#include "geom_primitive.h"
#include "es_specfunc.h"
#include <complex>
#include <string.h>
#ifdef EXTRAPRECISION_COLESO
#include "extraprecision.h"
#endif

template<typename fpv> 
inline complex<fpv> ImaginaryUnit() { return complex<fpv>(fpv(0.0), fpv(1.0)); }


//======================================================================================================================
template<typename fpv>
void s_WaveInChannel<fpv>::ReadParams(tFileBuffer& FB) {
    tParamManager PM;
    // Equation parameters
    PM.Request(nu, "nu");             // kinematic viscosity
    PM.Request(gamma, "gamma");       // specific ratio
    PM.Request(Prandtl, "Prandtl");   // Prandtl number
    // Geometry parameters
    PM.Request(form, "form");         // 0 - planar channel, 1 - cylindrical channel
    PM.Request(R, "R");               // radius of the cylinder (for form=1)
    PM.Request(CoorAxis, "CoorAxis"); // axial direction (direction normal to the circle): 0=X, 1=Y, 2=Z
    PM.Request(Ymin, "Ymin");         // for planar channel only: channel position
    PM.Request(Ymax, "Ymax");         // for planar channel only: channel position
    // solution parameters
    PM.Request(k, "k");               // axial wave number (real)
    PM.Request(l, "AzimuthalMode");   // azimuthal wave number (integer, >=0)
    PM.Request(kmode, "RadialMode");  // radial wave number (integer, >=0; used if frequency is not set)
    PM.Request(ReOmega, "ReOmega");   // real part of frequency (initial value for the iterative process)
    PM.Request(ImOmega, "ImOmega");   // imaginary part of frequency (initial value for the iterative process)
    PM.Request(ampl, "ampl");         // total multiplier
    PM.Request(phase, "phase");       // phase shift
    // additional
    PM.Request(_dmumax, "dmumax");    // шаг по мю при нахождении частоты (инит.)
    PM.Request(loglevel, "loglevel"); // уровень вывода отчёта об инициализации
    PM.RequestOption(ControlRootsJump, "ControlRootsJump"); // enhanced control of switching from one mode to another
    PM.ReadParamsFromBuffer(FB);
}
//======================================================================================================================


//======================================================================================================================
// Structure of auxiliary data calculated at the initialization step
//======================================================================================================================
template<typename fpv>
struct s_WaveInChannel_PrivateData {
    complex<fpv> Omega; // main parameter in need: eigenfrequency
    complex<fpv> varkappa_vort; // inverse boundary layer width for the vortex velocity field
    complex<fpv> varkappa_plus, varkappa_minus; // inverse boundary layer width for the two potential velocity fields
    complex<fpv> v1a, v1b; // some coefficients
    complex<fpv> Ynu_a2R, Ynu_sqlpR, Ynu_sqlmR; // some coefficients (form=1)
    complex<fpv> v2, alpha, beta; // some coefficients (form=1)
    fpv BZ; // zero of the Bessel function we have started from (for debug purpose, form=1)

    s_WaveInChannel_PrivateData() { Omega=varkappa_vort=varkappa_plus=varkappa_minus=v1a=v1b=Ynu_sqlpR=Ynu_sqlmR=alpha=beta=BZ=0.0; }
};
//======================================================================================================================


//======================================================================================================================
template<typename fpv>
void s_WaveInChannel<fpv>::Free(void) {
    if(data!=NULL) FreeMemSingle(data);
    data = NULL;
}
//======================================================================================================================


//======================================================================================================================
template<typename fpv>
s_WaveInChannel<fpv>& s_WaveInChannel<fpv>::operator=(const s_WaveInChannel<fpv> &gi) {
    if(this == &gi) return *this;
    if(data!=NULL) FreeMem(data);
    void *dstptr = (void*)this; // в s_WaveInChannel есть не подовые свойства (возможно сам fvp)
    void *srcptr = (void*)&gi;  // копировать по memcpy опасно! (заглушка ворнинга)
    memcpy(dstptr, srcptr, sizeof(s_WaveInChannel<fpv>));
    if(gi.data!=NULL) {
        data = GimmeMemSingle< s_WaveInChannel_PrivateData<fpv> >("s_WaveInChannel");
        *data = *(gi.data);
    }
    return *this;
}
//======================================================================================================================



//======================================================================================================================
// ---- Functions of a complex argument for arbitrary precision --------------------------------------------------------
//======================================================================================================================

//----------------------------------------------------------------------------------------------------------------------
template <typename fpv>
complex<fpv> A_exp(complex<fpv> z) {
//----------------------------------------------------------------------------------------------------------------------
    fpv e = exp(z.real());
    return complex<fpv>(e * cos(z.imag()), e * sin(z.imag()));
}
//----------------------------------------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------------------------------------
template <>
complex<double> A_exp(complex<double> z) {
//----------------------------------------------------------------------------------------------------------------------
    return exp(z);
}
//----------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------------------
template <class fpv>
complex<fpv> A_sqrt(complex<fpv> z) {
//----------------------------------------------------------------------------------------------------------------------
    complex<fpv> _I = ImaginaryUnit<fpv>();
    fpv phi = GetAngle(z.real(), z.imag());
    if(phi > GetPiNumber<fpv>()) phi -= GetPiNumber2<fpv>(); // выбираем нужную ветвь корня
    phi *= 0.5;
    fpv r = sqrt(sqrt(z.real()*z.real() + z.imag()*z.imag()));
    if(r*cos(phi)<0) crash("phi = %e", double(phi));
    return r*cos(phi) + _I*r*sin(phi);
}
//----------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------------------
template <>
complex<double> A_sqrt(complex<double> z) {
//----------------------------------------------------------------------------------------------------------------------
    return sqrt(z);
}
//----------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------------------
template<class fpv>
fpv cabs(complex<fpv> V) {
//----------------------------------------------------------------------------------------------------------------------
    return sqrt(SQR(V.real()) + SQR(V.imag()));
}
//======================================================================================================================


//======================================================================================================================
// ---- Auxiliary arithmetic subroutines -------------------------------------------------------------------------------
//======================================================================================================================

//======================================================================================================================
// cos(a*z) / cos(z) for a<1
//======================================================================================================================
template<typename fpv>
complex<fpv> cbb(fpv a, complex<fpv> z) {
    if(fabs(z.imag()) < 10.0) return cos(a*z) / cos(z);
    if(z.imag() < 0.0) z = -z;
    if(a < 0.0) a = -a;
    complex<fpv> iz = complex<fpv>(-z.imag(), z.real());
    return A_exp((1.0-a)*iz) * (fpv(1.0) + A_exp(fpv(2.0)*a*iz)) / (fpv(1.0) + A_exp(fpv(2.0)*iz));
}
//======================================================================================================================

//======================================================================================================================
// sin(a*z) / cos(z) for a<1
//======================================================================================================================
template<typename fpv>
complex<fpv> sbb(fpv a, complex<fpv> z) {
    if(fabs(z.imag()) < 10.0) return sin(a*z) / cos(z);
    complex<fpv> mult = complex<fpv>(0.0, 1.0);
    if(z.imag() < 0.0) { z = -z; mult = -mult; }
    if(a < 0.0) { a = -a; mult = -mult; }
    complex<fpv> iz = complex<fpv>(-z.imag(), z.real());
    return mult * A_exp((fpv(1.0)-a)*iz) * (fpv(1.0) - A_exp(fpv(2.0)*a*iz)) / (fpv(1.0) + A_exp(fpv(2.0)*iz));
}
//======================================================================================================================

//======================================================================================================================
// BesselJ(l,a*z) / BesselJ(l,z) for 0 <= a <= 1
//======================================================================================================================
template<class fpv>
complex<fpv> zbb(int l, fpv a, complex<fpv> z) {
    complex<fpv> _I = ImaginaryUnit<fpv>();
    if(fabs(z.imag()) * (1.0-a) > 240.0) {
        return fpv(0.0);
    }
    if(fabs(z.imag()) * a > 240.0) {
        int inv = 0;
        if(z.imag()>0) {
            z = z.real() - _I*z.imag();
            inv = 1;
        }

        complex<fpv> num, denom; num = denom = 0.0;
        complex<fpv> z1k, z2k; z1k = z2k = 1.0;
        const complex<fpv> mult1 = _I / (a*z), mult2 = _I / z;
        fpv ak = 1.0, bk = 1.0;
        for(int m=0; m<2*(int)sizeof(fpv); m++) { // ряд, возникающий в асимптотике бесселевских функций
            num += ak * z1k;
            denom += ak * z2k;
            z1k *= mult1;
            z2k *= mult2;
            bk = ak * (m+0.5);
            ak *= (4.0*l*l - (2*m+1)*(2*m+1)) / (8.0 * (m+1));
            bk += ak;
        }
        complex<fpv> res = num / denom * sqrt(1. / a) * A_exp(_I * (a-1.)*z);
        if(inv) res = res.real() - _I * res.imag();
        return res;
    }

    fpv bj1Re, bj1Im, bj2Re, bj2Im;
    BesselJComplex(l, z.real()*a, z.imag()*a, &bj1Re, &bj1Im);
    BesselJComplex(l, z.real(),   z.imag(),   &bj2Re, &bj2Im);
    if(bj2Re*bj2Re + bj2Im*bj2Im < 1e-300) crash("zbb: floating point error");
    return (bj1Re + _I*bj1Im) / (bj2Re + _I*bj2Im);
}
//======================================================================================================================

//======================================================================================================================
// Y(l,z) = z * BesselJ'(l,z) / BesselJ(l,z)
//======================================================================================================================
template<class fpv>
void zdbb(int l, complex<fpv> z, complex<fpv>* f, complex<fpv>* df) {
    if(f==NULL) crash("zdbb: NULL pointer given for f");
    complex<fpv> _I = ImaginaryUnit<fpv>();
    if(cabs(z) < 1e-50) {
        *f = fpv(l);
    }
    else if(fabs(z.imag()) > 50.0) {
        complex<fpv> sigma = _I / z;
        int inv = 0;
        if(z.imag() > 0) {
            sigma = -sigma.real() + _I*sigma.imag();
            inv = 1;
        }
        complex<fpv> num, denom; num = denom = 0.0;
        complex<fpv> sigmak; sigmak = 1.0;
        fpv ak = 1.0, bk = 1.0;
        for(int m=0; m<20; m++) { // ряд, возникающий в асимптотике бесселевских функций
            num += bk * sigmak;
            denom += ak * sigmak;
            sigmak *= sigma;
            bk = ak * (m+0.5);
            ak *= fpv(4*l*l - (2*m+1)*(2*m+1)) / fpv(8*(m+1));
            bk += ak;
        }
        *f = _I * z * num / denom;
        if(inv) *f = f->real() - _I*f->imag();
    }
    else if(fabs(z.real()) > 9000.0) {
        fpv pi4 = atan(fpv(1.0)); // pi/4
        *f = - z * tan(z - fpv(2*l+1)*pi4); // осторожно: возможно деление на 0
    }
    else {
        fpv bjRRe, bjRIm, bjxRRe, bjxRIm; // Bessel function value and derivative at x=R
        BesselJComplex(l, z.real(), z.imag(), &bjRRe, &bjRIm, &bjxRRe, &bjxRIm);
        if(bjRRe*bjRRe + bjRIm*bjRIm < 1e-300) crash("zdbb: floating point error");
        *f = z * (bjxRRe + _I*bjxRIm) / (bjRRe + _I*bjRIm);
    }
    if(df) {
        if(cabs(z) < 1e-50) *df = 0.0;
        else *df = (fpv(l*l) - z*z - *f**f) / z;
    }
}
//======================================================================================================================


//======================================================================================================================
// Calculating point value of the solution
// ATTENTION: s_WaveInChannel is a tComplexPulsFunction, it returns 10 components
//======================================================================================================================
template<typename fpv>
void s_WaveInChannel<fpv>::PointValue(fpv t, const fpv* coor, fpv* V) const {
    if(data==NULL) crash("s_WaveInChannel: Init not done");
    complex<fpv> _I = ImaginaryUnit<fpv>();
    complex<fpv> rho, p, ux, uy, uz;
    complex<fpv> mult = ampl;

    complex<fpv> omega = data->Omega;
    complex<fpv> CC4_3 = fpv(4.0)/fpv(3.0);
    complex<fpv> buf = fpv(1.0) / (fpv(1.0) + CC4_3*_I*nu*omega*gamma);

    int CoorX=0, CoorY=1, CoorZ=2;
    if(form==0) {
        const fpv H = 0.5*(Ymax - Ymin);
        const fpv y = coor[1] - 0.5*(Ymax + Ymin);

        complex<fpv> kappa_plus = cbb(y/H, H*data->varkappa_plus); // cos(y * varkappa_plus) / cos(H * varkappa_plus);
        complex<fpv> kappa_minus = cbb(y/H, H*data->varkappa_minus); // cos(y * varkappa_minus) / cos(H * varkappa_minus);
        complex<fpv> kappa_vort = cbb(y/H, H*data->varkappa_vort); // cos(y * varkappa_vort) / cos(H * varkappa_vort);
        complex<fpv> kappa_plus_dash = - data->varkappa_plus * sbb(y/H, H*data->varkappa_plus); // sin(y * varkappa_plus) / cos(H * varkappa_plus);
        complex<fpv> kappa_minus_dash = - data->varkappa_minus * sbb(y/H, H*data->varkappa_minus); // sin(y * varkappa_minus) / cos(H * varkappa_minus);
        complex<fpv> kappa_vort_dash = - data->varkappa_vort * sbb(y/H, H*data->varkappa_vort); // sin(y * varkappa_vort) / cos(H * varkappa_vort);
        complex<fpv> Phi = data->v1a * kappa_plus - data->v1b * kappa_minus;
        complex<fpv> PhiDash = data->v1a * kappa_plus_dash - data->v1b * kappa_minus_dash;
        complex<fpv> Eps = - (kappa_plus - kappa_minus);

        p = _I*omega*buf * (- Phi + CC4_3*nu*Eps);
        rho = gamma*p - Eps;
        ux = -_I*k*Phi;
        uy = PhiDash;
        uz = fpv(0.0);

        if(nu/Prandtl > 1e-50) {
            complex<fpv> As = _I*k*(data->v1a - data->v1b) / (data->varkappa_vort);
            ux += data->varkappa_vort*As*kappa_vort;
            uy += _I*k*As*kappa_vort_dash;
        }
        mult = ampl * A_exp(_I*(data->Omega*t - k*coor[0] + phase));
    }
    if(form == 1) {
        CoorZ = CoorAxis; CoorX = (CoorAxis + 1)%3; CoorY = (CoorAxis + 2)%3; // оси в плоскости круга
        fpv x = coor[CoorX], y = coor[CoorY], z = coor[CoorZ];
        fpv r = sqrt(x*x+y*y);
        for(int i=0; i<10; i++) V[i] = 0.0;

        if(r > R*1.0000000001) crash("s_WaveInChannel::PointValue error: r=%.16e > R=%.16e\n", double(r), double(R));
        if(r>R) r = R;

        const fpv tinylocr = (sizeof(fpv)==32) ? 1e-70 : 1e-25;

        complex<fpv>* Null = (complex<fpv>*) NULL;
        complex<fpv> Ynu_a2r, Ynu_sqlpr, Ynu_sqlmr;
        zdbb(l, data->varkappa_vort*r, &Ynu_a2r, Null);
        zdbb(l, data->varkappa_plus*r, &Ynu_sqlpr, Null);
        zdbb(l, data->varkappa_minus*r, &Ynu_sqlmr, Null);

        complex<fpv> kappa    = zbb(l, r/R, data->varkappa_vort*R);
        complex<fpv> Znu_sqlp = zbb(l, r/R, data->varkappa_plus*R);
        complex<fpv> Znu_sqlm = zbb(l, r/R, data->varkappa_minus*R);
        // Если теплопроводность отсутствует, считаем, что ГУ по температуре удовлетворять не нужно)
        if(Prandtl > fpv(9e49)*nu) Znu_sqlm = fpv(0.0);

        complex<fpv> Phi = (gamma-1.0) * ( - data->v1a * Znu_sqlp + data->v1b * Znu_sqlm);
        complex<fpv> Eps = (gamma-1.0) * (               Znu_sqlp -             Znu_sqlm);
        p = _I*omega*buf * (- Phi + CC4_3*nu*Eps);
        rho = gamma*p - Eps;
        uz = _I*k*Phi + data->alpha*kappa*(-_I*data->varkappa_vort);

        complex<fpv> ur = fpv(0.0), uphi = fpv(0.0);
        if(r>tinylocr) {
            ur = (gamma-1.0) / r * ( - data->v1a*Ynu_sqlpr*Znu_sqlp + data->v1b*Ynu_sqlmr*Znu_sqlm);
            uphi = _I*fpv(l)/r * Phi;

            if(cabs(data->alpha) > 1e-100) {
                complex<fpv> varkappa_vort = data->varkappa_vort;
                ur += data->alpha*kappa*(-k/(varkappa_vort*r)*(fpv(l) - Ynu_a2r));
                uphi += data->alpha*kappa*(_I*k/(varkappa_vort*r)*(fpv(l) - Ynu_a2r));
            }

            if(cabs(data->beta) > 1e-100) {
                ur += data->beta*kappa*(-_I*fpv(l)/r);
                uphi += data->beta*kappa*(Ynu_a2r / r);
            }
        }
        else if(l==1) {
            complex<fpv> InvJ1_plus, InvJ1_minus; InvJ1_plus = InvJ1_minus = 0.0; // 1 / J_nu(varkappa_+ * R), // 1 / J_nu(varkappa_- * R)
            if(fabs(data->varkappa_plus.imag()) < 100.0) {
                fpv bj1Re, bj1Im;
                BesselJComplex(l, data->varkappa_plus.real()*R, data->varkappa_plus.imag()*R, &bj1Re, &bj1Im);
                InvJ1_plus = fpv(1.0) / (bj1Re + _I * bj1Im);
            }            
            if(fabs(data->varkappa_minus.imag()) < 100.0) {
                fpv bj1Re, bj1Im;
                BesselJComplex(l, data->varkappa_minus.real()*R, data->varkappa_minus.imag()*R, &bj1Re, &bj1Im);
                InvJ1_minus = fpv(1.0) / (bj1Re + _I * bj1Im);
            }            
            ur   = fpv(0.5)*(gamma-1.0) * ( - data->v1a * data->varkappa_plus * InvJ1_plus + data->v1b * data->varkappa_minus * InvJ1_minus);

            if(cabs(data->beta) > 1e-100) {
                complex<fpv> InvJ1_a2; InvJ1_a2 = 0.0;
                if(fabs(data->varkappa_vort.imag()) < 100.0) {
                    fpv bj1Re, bj1Im;
                    BesselJComplex(l, data->varkappa_vort.real()*R, data->varkappa_vort.imag()*R, &bj1Re, &bj1Im);
                    InvJ1_a2 = fpv(1.0) / (bj1Re + _I * bj1Im);
                }            
                ur -= _I*data->beta*fpv(0.5)*data->varkappa_vort * InvJ1_a2;
            }
            uphi = _I * ur;
        }

        if(r<tinylocr) { // предполагаем, что выполняется ur + iuphi = 0; иное означало бы особенность в нуле
            ux = ur;
            uy = uphi;
            // при l=0 домножать на exp(i*l*phi) не нужно, а при l>1 rho=p=uz = 0
        }
        else {
            complex<fpv> phi = GetAngle(x, y);
            complex<fpv> mult0 = A_exp(_I*fpv(l  )*phi);
            complex<fpv> multp = A_exp(_I*fpv(l+1)*phi);
            complex<fpv> multm = A_exp(_I*fpv(l-1)*phi);

            ux = fpv(0.5)*multp*( ur + _I*uphi) + fpv(0.5)*multm*(ur - _I*uphi);
            uy = fpv(0.5)*multp*(-_I*ur + uphi) + fpv(0.5)*multm*(_I*ur + uphi);
            p   *= mult0;
            rho *= mult0;
            uz  *= mult0;
        }

        mult = ampl * A_exp(_I*(data->Omega*t + k*z + phase));
    }

    uz *= mult;
    ux *= mult;
    uy *= mult;
    rho *= mult;
    p *= mult;

    V[Var_R]         = rho.real();
    V[Var_P]         = p.real();
    V[Var_U+CoorX]   = ux.real(); // in RZ: radial
    V[Var_U+CoorY]   = uy.real(); // in RZ: azimuthal
    V[Var_U+CoorZ]   = uz.real();
    V[Var_R+5]       = rho.imag();
    V[Var_P+5]       = p.imag();
    V[Var_U+CoorX+5] = ux.imag();
    V[Var_U+CoorY+5] = uy.imag();
    V[Var_U+CoorZ+5] = uz.imag();
}
//======================================================================================================================


//======================================================================================================================
// Вычисление функции и её производной при решении методом Ньютона уравнения для определения Omega
// Вместе с этим пересчитывается структура data
//======================================================================================================================
template<typename fpv>
complex<fpv> CalcFun(const s_WaveInChannel<fpv>& S, 
                     const s_WaveInChannel_PrivateData<fpv>* data_prev,  
                     complex<fpv> omega, complex<fpv> mu,
                     complex<fpv> *dfun, complex<fpv> *dfundmu) {
    const fpv CC1 = 1.0;
    const fpv CC2 = 2.0;
    const fpv CC4_3 = fpv(4.0)/3.0;
    complex<fpv> _I = ImaginaryUnit<fpv>();

    complex<fpv> buf = CC1 / (CC1 + CC4_3*_I*mu*omega*S.gamma);
    complex<fpv> a = omega*omega*S.gamma * buf;
    complex<fpv> b = -_I*omega * buf;
    complex<fpv> c = S.Prandtl*omega*omega*(S.gamma-1.0) / mu * buf;
    complex<fpv> d = -_I*omega * S.Prandtl/mu * (CC1 + CC4_3*_I*mu*omega) * buf;
    complex<fpv> discr = SQR(fpv(0.5)*(a-d)) + b*c;
    complex<fpv> sqdiscr = A_sqrt(discr);
    if((sqdiscr / d).real() > 0.0) sqdiscr = -sqdiscr; // выбираем корень
    complex<fpv> lambda_plus  = fpv(0.5)*(a+d) + sqdiscr;
    complex<fpv> lambda_minus = fpv(0.5)*(a+d) - sqdiscr;

    if(cabs(lambda_plus)*100.0 < cabs(sqdiscr)) { // случай, близкий к нетепловодному. Нужно уточнить значение корня lambda_plus
        complex<fpv> m = a*a/(d*d)-fpv(2.0)*a/d + fpv(4.0)*b*c/(d*d); // discr / (0.25*d*d) - 1;
        if(cabs(m) > 0.9) crash("s_WaveInChannel::CalcFun internal error in lambda_plus series");
        // Теперь lambda_plus = 0.5*a + 0.5*d*(1 - sqrt(1 + m))
        lambda_plus = 0.0;
        complex<fpv> ak; ak = 1.0;
        for(int kk=1; cabs(ak) > 1e-100; kk++) {
            ak = ak * m * fpv((1.5-kk)/kk);
            lambda_plus += ak;
        }
        lambda_plus = fpv(0.5)*(a - d*lambda_plus);
    }

    complex<fpv> varkappa_plus  = A_sqrt(lambda_plus  - S.k*S.k);
    complex<fpv> varkappa_minus = A_sqrt(lambda_minus - S.k*S.k);
    complex<fpv> varkappa_vort  = A_sqrt(-_I*omega/mu - S.k*S.k);

    // На всякий случай сохраняем те же знаки, что были на предыдущей итерации
    if(data_prev!=NULL) {
        if(varkappa_vort.real()*data_prev->varkappa_vort.real() + varkappa_vort.imag()*data_prev->varkappa_vort.imag() < 0.0) varkappa_vort = -varkappa_vort;
        if(varkappa_plus.real() *data_prev->varkappa_plus.real()  +  varkappa_plus.imag()*data_prev->varkappa_plus.imag()  < 0.0) varkappa_plus  = -varkappa_plus;
        if(varkappa_minus.real()*data_prev->varkappa_minus.real() + varkappa_minus.imag()*data_prev->varkappa_minus.imag() < 0.0) varkappa_minus = -varkappa_minus;
    }
    complex<fpv> v1a = (lambda_minus - a)/c;
    complex<fpv> v1b = b/(lambda_plus - d);

    s_WaveInChannel_PrivateData<fpv>& data = *(S.data);
    data.Omega = omega;
    data.varkappa_plus = varkappa_plus;
    data.varkappa_minus = varkappa_minus;
    data.varkappa_vort = varkappa_vort;
    data.v1a = v1a;
    data.v1b = v1b;

    complex<fpv> fun; // returning value

    // Auxiliary data for cylindrical channel
    complex<fpv> Ynu_a2R, Ynu_sqlpR, Ynu_sqlmR;
    complex<fpv> dYnu_a2R0, dYnu_sqlpR0, dYnu_sqlmR0;
    complex<fpv> v1, v2, buffun;

    if(S.form==1) { // Cylindrical channel
        zdbb(S.l, varkappa_plus*S.R, &Ynu_sqlpR, &dYnu_sqlpR0);
        zdbb(S.l, varkappa_minus*S.R, &Ynu_sqlmR, &dYnu_sqlmR0);
        zdbb(S.l, varkappa_vort*S.R, &Ynu_a2R, &dYnu_a2R0);

        v1 = v1a - v1b;
        v2 = v1a * Ynu_sqlpR - v1b * Ynu_sqlmR;
        buffun = S.k*S.k / (varkappa_vort*varkappa_vort) * (Ynu_a2R*Ynu_a2R - fpv(S.l*S.l)) - fpv(S.l*S.l);
        fun = v1 * buffun + Ynu_a2R * v2;

        // Запоминаем все насчитанные данные (не насчитаны только alpha и beta)
        data.Ynu_a2R = Ynu_a2R;
        data.Ynu_sqlmR = Ynu_sqlmR;
        data.Ynu_sqlpR = Ynu_sqlpR;
    }
    else if(S.form==0) { // Planar channel
        fpv H = 0.5*(S.Ymax - S.Ymin);
        complex<fpv> s1 = v1a * varkappa_plus*tan(H*varkappa_plus);
        complex<fpv> s2 = -v1b * varkappa_minus * tan(H*varkappa_minus);
        complex<fpv> s3 = S.k*S.k*(v1a-v1b) * tan(H*varkappa_vort)/varkappa_vort;
        fun = s1 + s2 + s3;
    }
    else { // free space
        fun = lambda_plus - S.k * S.k;
    }

    // считаем производные по omega и mu
    for(int iii=0; iii<2; iii++) {
        complex<fpv>* DFUN = (iii==0) ? dfun : dfundmu;
        if(DFUN==NULL) continue;
        
        complex<fpv> dbuf, da, db, dc, dd, da2;
        if(iii==0) { // производные по omega
            dbuf = - CC4_3*_I*mu*S.gamma * buf*buf; // buf = CC1 / (CC1 + CC4_3*_I*mu*omega*S.gamma);
            da = CC2*omega*S.gamma * buf + omega*omega*S.gamma * dbuf; // a = omega*omega*S.gamma * buf;
            db = -_I * buf - _I*omega * dbuf; // b = -_I*omega * buf;
            dc = S.Prandtl*(S.gamma-1.0)/mu * (CC2*omega* buf + omega*omega*dbuf); // c = S.Prandtl*omega*omega*(S.gamma-1.0) / mu * buf;
            // d = -_I*omega * S.Prandtl/mu * (CC1 + CC4_3*_I*mu*omega) * buf;
            dd = -_I       * S.Prandtl/mu * (CC1 + CC4_3*_I*mu*omega) * buf
                 -_I*omega * S.Prandtl/mu * (      CC4_3*_I*mu      ) * buf 
                 -_I*omega * S.Prandtl/mu * (CC1 + CC4_3*_I*mu*omega) * dbuf;
            da2 = (-_I/mu) / (CC2 * varkappa_vort); // varkappa_vort = A_sqrt(-_I*omega/mu - k*k);
        }
        else { // производные по mu
            dbuf = -CC4_3*_I*omega*S.gamma * buf*buf;
            da = omega*omega*S.gamma * dbuf;
            db = -_I*omega * dbuf;
            dc = S.Prandtl*omega*omega*(S.gamma-1.0) * (dbuf / mu - buf / (mu*mu));
            dd = -_I*omega * S.Prandtl * ((CC1 + CC4_3*_I*mu*omega) * (dbuf / mu - buf / (mu*mu)) + CC4_3*_I*omega * (buf /mu));
            da2 = (_I*omega/(mu*mu)) / (CC2 * varkappa_vort); // varkappa_vort = A_sqrt(-_I*omega/mu - k*k);
        }

        // дальнейшие операции одинаковые для обоих дифференцирований
        complex<fpv> ddiscr = fpv(0.5)*(a-d)*(da-dd) + b*dc + db*c; // discr = SQR(fpv(0.5)*(a-d)) + b*c;
        complex<fpv> dsqdiscr = ddiscr / (CC2*sqdiscr); // sqdiscr = A_sqrt(discr);
        complex<fpv> dlambda_plus  = fpv(0.5)*(da+dd) + dsqdiscr;
        complex<fpv> dlambda_minus = fpv(0.5)*(da+dd) - dsqdiscr;
        if(cabs(lambda_plus)*100.0 < cabs(sqdiscr)) { // случай, близкий к нетепловодному. Уточняем производную от lambda_plus
            complex<fpv> m = a*a/(d*d)-fpv(2.0)*a/d + fpv(4.0)*b*c/(d*d); // discr / (0.25*d*d) - 1;
            complex<fpv> dm = CC2*a/(d*d) * (da - a*dd/d) - CC2*(da/d - a*dd/(d*d)) + CC2*CC2*((db*c+dc*b)/(d*d) - CC2*b*c*dd/(d*d*d));
            if(cabs(m) > 0.9) crash("s_WaveInChannel::CalcFun internal error in lambda_plus series");
            // Теперь lambda_plus = 0.5*a + 0.5*d*(1 - sqrt(1 + m))
            dlambda_plus = 0.0;
            complex<fpv> _lambda_plus; _lambda_plus = 0.0;
            complex<fpv> ak; ak = 1.0;
            complex<fpv> dak; dak = 0.0;
            for(int kk=1; cabs(ak) > 1e-100; kk++) {
                dak = (dak * m + ak * dm) * fpv((1.5-kk)/kk);
                ak = ak * m * fpv((1.5-kk)/kk);
                _lambda_plus += ak;
                dlambda_plus += dak;
            }
            dlambda_plus = fpv(0.5)*(da - dd*_lambda_plus - d*dlambda_plus);
        }
        complex<fpv> dvarkappa_plus = dlambda_plus  / (CC2*varkappa_plus); // varkappa_plus = A_sqrt(lambda_plus - k*k);
        complex<fpv> dvarkappa_minus = dlambda_minus  / (CC2*varkappa_minus); // varkappa_minus = A_sqrt(lambda_minus - k*k);

        complex<fpv> dv1a = ((dlambda_minus - da)*c - dc*(lambda_minus-a))/(c*c); // v1a = (lambda_minus - a)/c;
        complex<fpv> dv1b = (db*(lambda_plus - d) - (dlambda_plus - dd)*b)/(SQR(lambda_plus - d)); // v1b = b/(lambda_plus - d);

        if(S.form==1) { // Cylindrical channel
            complex<fpv> dYnu_a2R   = dYnu_a2R0 * S.R * da2;
            complex<fpv> dYnu_sqlpR = dYnu_sqlpR0 * S.R * dvarkappa_plus;
            complex<fpv> dYnu_sqlmR = dYnu_sqlmR0 * S.R * dvarkappa_minus;

            complex<fpv> dv1 = dv1a - dv1b; // v1 = v1a - v1b;
            complex<fpv> dv2 = dv1a * Ynu_sqlpR + v1a*dYnu_sqlpR - dv1b * Ynu_sqlmR - v1b*dYnu_sqlmR; // v2 = v1a * Ynu_sqlpR - v1b * Ynu_sqlmR;

            // complex<fpv> buffun = k*k / (varkappa_vort*varkappa_vort) * (Ynu_a2R*Ynu_a2R - fpv(l*l)) - fpv(l*l);
            complex<fpv> dbuffun = S.k*S.k * (-CC2*da2 / (varkappa_vort*varkappa_vort*varkappa_vort) * (Ynu_a2R*Ynu_a2R - fpv(S.l*S.l)) + CC2/(varkappa_vort*varkappa_vort) * Ynu_a2R * dYnu_a2R);
            // complex<fpv> fun = v1 * buffun + Ynu_a2R * v2;
            *DFUN = dv1 * buffun + v1 * dbuffun + dYnu_a2R * v2 + Ynu_a2R * dv2;
        }
        else if(S.form==0) { // Planar channel
            fpv H = 0.5*(S.Ymax - S.Ymin);
            complex<fpv> TanPlus  = tan(H*varkappa_plus);
            complex<fpv> TanMinus = tan(H*varkappa_minus);
            complex<fpv> TanV     = tan(H*varkappa_vort);

            complex<fpv> ds1 = dv1a*varkappa_plus*TanPlus
                                + v1a*dvarkappa_plus*TanPlus
                                + v1a*varkappa_plus*(CC1 + TanPlus*TanPlus) * H*dvarkappa_plus;
            complex<fpv> ds2 = -dv1b * varkappa_minus * TanMinus
                - v1b*dvarkappa_minus*(TanMinus + varkappa_minus*(CC1+TanMinus*TanMinus));
            complex<fpv> ds3 = (dv1a-dv1b) * TanV / varkappa_vort
                + (v1a-v1b)*da2* (TanV / (varkappa_vort*varkappa_vort) + (CC1 + TanV*TanV) / varkappa_vort);
            ds3 *= (S.k*S.k);

            *DFUN = ds1 + ds2 + ds3;
        }
        else { // free space
            *DFUN = dlambda_plus;
        }
    }
    return fun;
}
//======================================================================================================================


//======================================================================================================================
// Initialization of the solution for wave in planar or cylindrical channel, 
// generally with viscosity and heat conductivity
//======================================================================================================================
template<typename fpv>
void s_WaveInChannel<fpv>::Init() {
    // Checking parameters
    if(nu<0.) crash("s_WaveInChannel::Init: wrong nu=%e", double(nu));
    if(gamma<=1.) crash("s_WaveInChannel::Init: wrong gamma=%e", double(gamma));
    if(Prandtl<=0.) crash("s_WaveInChannel::Init: wrong Prandtl=%e", double(Prandtl));
    if(form!=0 && form!=1) crash("s_WaveInChannel::Init: wrong form=%i", form);
    if(form==1) {
        if(R<=0.) crash("s_WaveInChannel::Init: wrong R=%e", double(R));
        if(CoorAxis<0 || CoorAxis>2) crash("s_WaveInChannel::Init: wrong CoorAxis=%i", CoorAxis);
        if(l<0) crash("s_WaveInChannel::Init: wrong AzimuthalMode=%i", l);
    }
    if(!IsNaN(_dmumax)) if(_dmumax<=0.) crash("s_WaveInChannel::Init: wrong dmumax=%e", double(_dmumax));
    if((form==1 || form==0) && kmode<0) crash("s_WaveInChannel::Init: wrong RadialMode=%i", kmode);

    // Allocating data, solving the equation for omega and printing the root found
    int errcode = SolveEquation();
    if(errcode) crash("s_WaveInChannel::Init: error solving main equation, errcode = %i", errcode);

    if(loglevel>=1) {
        pprintf("s_WaveInChannel: omega = % .14e + i * % .14e\n", double(data->Omega.real()), double(data->Omega.imag()));
        pprintf("s_WaveInChannel: kap_v = % .14e + i * % .14e\n", double(data->varkappa_vort.real()), double(data->varkappa_vort.imag()));
        pprintf("s_WaveInChannel: kap_+ = % .14e + i * % .14e\n", double(data->varkappa_plus.real()), double(data->varkappa_plus.imag()));
        pprintf("s_WaveInChannel: kap_- = % .14e + i * % .14e\n", double(data->varkappa_minus.real()), double(data->varkappa_minus.imag()));
    }

    // Checking that the solution satisfies the boundary condition
    if(form==1) {
        int inviscid = nu<1e-50;
        int CoorR = (CoorAxis+1)%3;
        fpv coor[3] = {0.,0.,0.};
        coor[CoorR] = R;
        fpv V[10];
        PointValue(0.0, coor, V);
        for(int ireim=0; ireim<2; ireim++) {
            fpv* pV = V + ireim*5;
            for(int i=(inviscid?CoorR+1:1); i<=(inviscid?CoorR+1:3); i++) if(fabs(pV[i])>1e-10)
                crash("Error: nonzero velocities at boundary: \n  uz   = (% e + i*% e)\n  ur   = (% e + i*% e)\n  uphi = (% e + i*% e)",
                double(V[1]), double(V[6]), double(V[2]), double(V[7]), double(V[3]), double(V[8]));
        }
    }
}
//======================================================================================================================


//======================================================================================================================
// Solving algebraic equation using Newton method
// Input: mu, omega, tolerance, MAX_ITERS; prev (obsolete)
// Output: omega, dfun, dfundmu
// Returns nonzero if error
//======================================================================================================================
template<typename fpv>
int NewtonProcess(s_WaveInChannel<fpv>& S, fpv mu, complex<fpv>& omega, fpv tolerance, int MAX_ITERS, 
                                        const s_WaveInChannel_PrivateData<fpv>* prev, complex<fpv>& dfun, complex<fpv>& dfundmu) {
    fpv residual = 1e300;
    fpv coeff = 1.0;
    for(int iter=0; iter<MAX_ITERS; iter++) {
        complex<fpv> fun;
        fun = CalcFun<fpv>(S, prev, omega, mu, &dfun, &dfundmu);

        if(S.loglevel>=3){
            pprintf("trying om = %e %e...\n", double(omega.real()), double(omega.imag()));
            pprintf("       varkappa_vort = %e %e...\n", double(S.data->varkappa_vort.real()), double(S.data->varkappa_vort.imag()));
            pprintf("       xp = %e %e...\n", double(S.data->varkappa_plus.real()), double(S.data->varkappa_plus.imag()));
            pprintf("       xm = %e %e...\n", double(S.data->varkappa_minus.real()), double(S.data->varkappa_minus.imag()));
            pprintf("fun, dfun, domega = %e %e %e\n", double(cabs(fun)), double(cabs(dfun)), double(cabs(fun/dfun)));
        }

        if(cabs(fun / dfun) < tolerance) return 0; // Success
        if(cabs(fun) > residual) return 1; // Ньютоновский процесс потерял монотонность

        residual = cabs(fun);
        complex<fpv> omega_new = omega - coeff * fun / dfun;

        if(omega.real() < 0.0) { // если пытаемся пересечь ось, то останавливаемся на ней
            omega += (omega_new - omega) * omega.real() / (omega.real() - omega_new.real());
        }
        else omega = omega_new;
    }
    return 2; // No success for MAX_ITERS iterations
}
//======================================================================================================================

//======================================================================================================================
// Solving algebraic equation using successive increasing of viscosity coefficient from zero up to the required value
// Returns nonzero if error
//======================================================================================================================
template<typename fpv>
int s_WaveInChannel<fpv>::SolveEquation(void) {
    if(data==NULL) data = GimmeMemSingle< s_WaveInChannel_PrivateData<fpv> >( "s_WaveInChannel");

    // Начальное приближение итерационного процесса -- omega для невязкой задачи
    complex<fpv> _I(fpv(0.0), fpv(1.0));
    if(!IsNaN(ReOmega) && !IsNaN(ImOmega)) {
        data->Omega = complex<fpv>(ReOmega, ImOmega);
    }
    else if(form==1) {
        if(l<0) crash("s_WaveInChannel::SolveEquation: wrong azimuthal mode %i", l);
        if(kmode<0) crash("s_WaveInChannel::SolveEquation: wrong radial mode %i", kmode);
        data->BZ = BesselPrimeZero(l, kmode, loglevel>1); // double precision
        if(IsNaN(data->BZ) || data->BZ<l) crash("s_WaveInChannel::Init: internal error, Zero = %e < l = %i", double(data->BZ), l);
        data->Omega = sqrt(SQR(data->BZ / R) + k*k);
    }
    else if(form==0) {
        const double H = 0.5*(Ymax - Ymin);
        data->Omega = sqrt(k*k + SQR(kmode*PiNumber/H));
    }
    else {
        data->Omega = fabs(k);
    }
    if(loglevel>=1) pprintf("omega(initial) = %e %e\n", double(data->Omega.real()), double(data->Omega.imag()));

    if(form==1 || form==0) {
        if(form==1) data->varkappa_plus = data->BZ / R;
        else data->varkappa_plus = kmode*PiNumber/(0.5*(Ymax - Ymin));
        data->varkappa_vort = data->varkappa_minus = 1e60; // скорость затухания вязкого и теплопроводного погранслоёв бесконечны
        data->v1a = -_I / (data->Omega * (fpv(gamma)-1.0));
    }

    // Исходя из линейного приближения оцениваем максимально возможный коэффициент вязкости, при котором Re(omega)>0
    // Нарушение этого условия, как правило, не позволяет найти корень уравнения
    // Но эта оценка очень грубая. Реальное значение может отличаться на 1-2 порядка
    //pprintf("Estimated maximal visc = %e\n", 2.0*R/BZ*SQR(BZ*BZ-l*l)*pow(1.0+k*k*R*R/BZ/BZ,1.5)/(BZ*BZ+k*k*R*R));

    // В цикле увеличиваем коэффициент вязкости
    complex<fpv> omega = data->Omega;
    fpv mu_prev = 0.0;
    s_WaveInChannel_PrivateData<fpv> prev = *data; // данные с предыдущей итерации
    const fpv dmumax = IsNaN(_dmumax) ? MAX(0.01 * nu, 1e-3 / (1.0 + 1.0/Prandtl)) : _dmumax;
    const fpv DMU_MIN = IsNaN(_dmumin) ? (sizeof(fpv)==8 ? 1e-5*nu : fpv(1e-10)) : _dmumin;
    fpv dmu = dmumax; // шаг по коэффициенту вязкости
    fpv mu = 0.0; // текущий коэффициент вязкости
    // если nu=0, то формально решение уже найдено. Но если точность выше двойной, запускаем 1 итерацию для уточнения нуля функции Бесселя
    int finishflag = (nu < 1e-50) && (IsNaN(ReOmega) && IsNaN(ImOmega)) && sizeof(fpv)==8;
    const fpv tolerance = (sizeof(fpv)==32) ? 1e-45 : 1e-12;
    while(!finishflag) {
        // Вычисляем новое значение коэффициента вязкости
        mu = mu_prev + dmu;
        if(mu >= nu) { mu = nu; finishflag = 1; }

        if(loglevel) pprintf("Try mu = %.7e ... ", double(mu));

        // Запускаем Ньютоновский процесс
        complex<fpv> dfundmu, dfun;
        int errcode = NewtonProcess<fpv>(*this, mu, omega, tolerance, 1000, &prev, dfundmu, dfun);

        if(errcode==0 && ControlRootsJump) {
            // процесс сошёлся! Проверяем, что сошёлся куда надо
            complex<fpv> domega = - dmu * dfundmu / dfun;
            complex<fpv> err = (omega - prev.Omega) / (domega) - fpv(1.0);
            if(cabs(err) > 1e-4 / sqrt(nu - mu_prev + 1e-10) && mu>1e-4) {
                errcode = 3;
                if(loglevel>=1) {
                    complex<fpv> omega_taylor = prev.Omega + domega;
                    pprintf("Wrong convergence:\n");
                    pprintf("omega_found  = %20.15e %20.15e\n", double(omega.real()), double(omega.imag()));
                    pprintf("omega_Taylor = %20.15e %20.15e\n", double(omega_taylor.real()), double(omega_taylor.imag()));
                }
            }
        }

        if(loglevel>=1) {
            pprintf("  om = %.13e %.13e\n", double(omega.real()), double(omega.imag()));
        }
        if(loglevel>=2 || (loglevel>=1 && errcode)) {
            pprintf("  vv = %.13e %.13e\n", double(data->varkappa_vort.real()), double(data->varkappa_vort.imag()));
            pprintf("  v+ = %.13e %.13e\n", double(data->varkappa_plus.real()), double(data->varkappa_plus.imag()));
            pprintf("  v- = %.13e %.13e\n", double(data->varkappa_minus.real()), double(data->varkappa_minus.imag()));
        }

        // Ньютоновский процесс потерял монотонность либо стагнирует
        if(errcode) {
            if(dmu < DMU_MIN) return 1;

            dmu *= 0.5; // уменьшаем шаг по nu
            if(loglevel) pprintf("Fallback: now dmu = %e\n", double(dmu));
            *data = prev;
            omega = prev.Omega;
            finishflag = 0;
            continue;
        }

        // success
        mu_prev = mu;
        prev = *data;

        if(dmu<dmumax) {
            dmu *= 1.1;
            if(dmu>dmumax) dmu=dmumax;
        }
    }

    // Copying data
    ReOmega = data->Omega.real();
    ImOmega = data->Omega.imag();

    // Досчитываем переменные, которые были не нужны ранее
    if(form==1) {
        data->v2 = data->v1a * (data->Ynu_sqlpR+fpv(l)) - data->v1b * (data->Ynu_sqlmR+fpv(l));
        if(!(nu < 1e-50)) {
            data->alpha = - fpv(k) / data->varkappa_vort * (fpv(gamma) - 1.0) * (data->v1a - data->v1b);
            data->beta = _I * (fpv(gamma) - 1.0)/(data->Ynu_a2R + fpv(l)) * data->v2;
        }
    }
    return 0;
}
//======================================================================================================================


// Сообщаем компилятору о желании иметь функции нужных нам типов
template struct s_WaveInChannel<double>;
#ifdef EXTRAPRECISION_COLESO
template struct s_WaveInChannel<dd_real>;
template struct s_WaveInChannel<qd_real>;
#endif





// *********************************************************************************************************************
// *****                                                                                                           *****
// *****                                             s_VortexInCylinder                                            *****
// *****           simple partial case of the general solution (zero azimuthal mode, no heat conductivity)         *****
// *****                                                                                                           *****
// *********************************************************************************************************************

//======================================================================================================================
void s_VortexInCylinder::ReadParams(tFileBuffer& FB) {
    tParamManager PM;
    PM.Request(nu, "nu"); // Kinematic viscosity
    PM.Request(R, "R"); // Cylinder radius
    PM.Request(CoorAxis, "CoorAxis"); // axial direction
    PM.Request(NumModes, "NumModes"); // number of terms in series
    PM.Request(profile, "profile"); // initial data for u_phi(x), where x = r/R
    PM.ReadParamsFromBuffer(FB);
}
//======================================================================================================================

//======================================================================================================================
void s_VortexInCylinder::Init() {
    if(CoorAxis<0 || CoorAxis>2) crash("s_VortexInCylinder::Init error: wrong CoorAxis");
    if(NumModes<=0) crash("s_VortexInCylinder::Init error: NumModes should be > 0");
    if(nu <= 0.0 || R <= 0.0) crash("s_VortexInCylinder::Init error: wrong parameters");
    
    tGaussIntegrator<double> GI;
    GI.Init(128);

    // Parsing the string with profile function
    double x; // argument of function
    tFormula F;
    F.AddArg("x", 0.0, &x);
    F.Parse(profile.c_str(), true /*log*/);

    data.resize(NumModes*2);
    for(int imode=0; imode<NumModes; imode++) {
        // Определяем нуль функции Бесселя
        double bz = 3.83170597020751 + PiNumber * imode; // начальное приближение
        for(int i=0; i<4; i++) { // уточняем методом Ньютона
            double bj0 = BesselJ0(bz), bj1 = BesselJ1(bz);
            bz *= 1.0 + bj1 / (bj1 - bj0 * bz);
        }

        // Разлагаем начальные данные
        double sum = 0.0;
        for(int j=0; j<GI.GR; j++) {
            x = GI.GN[j];
            double f = F.Evaluate(); // default: f = x*sin(PiNumber*x)
            sum += GI.GC[j] * x*BesselJ1(x*bz)*f;
        }
        sum /= 0.5 * SQR(BesselJ0(bz));
        if(IsNaN(sum)) crash("s_VortexInCylinder::Init: NaN detected");

        data[imode] = bz / R;
        data[imode+NumModes] = sum;
    }
}
//======================================================================================================================

//======================================================================================================================
void s_VortexInCylinder::PointValue(double t, const double* coor, double* V) const {
    int CoorX = (CoorAxis + 1)%3, CoorY = (CoorAxis + 2)%3; // оси в плоскости круга

    V[Var_R] = V[Var_U] = V[Var_V] = V[Var_W] = V[Var_P] = 0.0;
    double r = sqrt(SQR(coor[CoorX]) + SQR(coor[CoorY]));
    if(r<tiny) return;

    double uphi = 0.0;
    for(int imode=0; imode<NumModes; imode++) {
        double omega = data[imode];
        uphi += data[imode+NumModes] * exp(-nu*omega*omega*t) * BesselJ1(omega*r);
    }
    
    V[Var_U + CoorX] = -uphi*coor[CoorY]/r;
    V[Var_U + CoorY] =  uphi*coor[CoorX]/r;
}
//======================================================================================================================



// *********************************************************************************************************************
// *****                                                                                                           *****
// ****                               s_SinusVisc (single mode in free space)                                      *****
// *****                                                                                                           *****
// *********************************************************************************************************************

//======================================================================================================================
void s_SinusVisc::ReadParams(tFileBuffer& FB) {
    tParamManager PM;
    PM.Request(mu, "mu"); // viscosity
    PM.Request(gamma, "gamma"); // specific ratio
    PM.Request(Prandtl, "Prandtl"); //Prandtl number
    PM.Request(Ampl, "Ampl"); // amplitude of acoustic mode
    PM.Request(AmplVX, "AmplVX"); // amplitudes of vortex mode (no-heat-conductive case only)
    PM.Request(AmplVY, "AmplVY"); // amplitudes of vortex mode (no-heat-conductive case only)
    PM.Request(AmplVZ, "AmplVZ"); // amplitudes of vortex mode (no-heat-conductive case only)
    PM.Request(kx, "kx"); // components of the wave vector
    PM.Request(ky, "ky"); // components of the wave vector
    PM.Request(kz, "kz"); // components of the wave vector
    PM.ReadParamsFromBuffer(FB);
}
//======================================================================================================================

//======================================================================================================================
void s_SinusVisc::Init(void){
    if(SQR(kx) + SQR(ky) + SQR(kz) < 1e-30) crash("s_SinusVisc error: kx=ky=kz=0");

    if(Prandtl < 5e49 && mu > 1e-50) { // heat-conductive case
        if(fabs(AmplVX) + fabs(AmplVY) + fabs(AmplVZ) > tiny)
            crash("s_SinusVisc: vortex modes are allowed only for no-heat-conductive case (P>1e50)");

        s_WaveInChannel<double> S;
        S.nu = mu; S.gamma = gamma; S.Prandtl = Prandtl;
        S.form = -1; // free space
        S.k = sqrt(SQR(kx) + SQR(ky) + SQR(kz));
        S.ReOmega = ReOmega; S.ImOmega = ImOmega; // complex frequency. If set, kmode is ignored
        S.loglevel = 0; // уровень вывода отчёта об инициализации

        int errcode = S.SolveEquation();
        if(errcode) crash("s_SinusVisc::Init: error solving equation for omega, errcode = %i", errcode);
        ReOmega = S.ReOmega; ImOmega = S.ImOmega;
        ReCoeff = S.data->v1a.real(); ImCoeff = S.data->v1a.imag();
    }
}
//======================================================================================================================

//======================================================================================================================
void s_SinusVisc::PointValue(double t, const double* coor, double* V) const {
    if(Prandtl < 5e49 && mu > 1e-50) { // heat-conductive case
        if(IsNaN(ReOmega) || IsNaN(ImOmega) || IsNaN(ReCoeff) || IsNaN(ImCoeff)) crash("s_SinusVisc: init not done");
        double kr = kx*coor[0] + ky*coor[1] + kz*coor[2];

        complex<double> _I = ImaginaryUnit<double>();
        complex<double> omega = ReOmega + _I * ImOmega;
        complex<double> Eps = 1.0;
        complex<double> W = - Eps * (ReCoeff + _I * ImCoeff);
        complex<double> p = _I*omega*(C4_3*mu*Eps - W) / (double(1.) + double(C4_3)*_I*omega*mu*gamma);
        complex<double> rho = gamma * p - Eps;
        complex<double> mult = Ampl * exp(_I*(omega*t + kr /*+ phase*/));
        V[Var_R] = (mult * rho).real();
        V[Var_U] = (_I * mult * W).real() * kx;
        V[Var_V] = (_I * mult * W).real() * ky;
        V[Var_W] = (_I * mult * W).real() * kz;
        V[Var_P] = (mult * p).real();
    }
    else {
        double absk = sqrt(SQR(kx) + SQR(ky) + SQR(kz));
        double reOmega = absk, imOmega = 0.0, ReKappa = -absk, ImKappa = 0.0, _Ampl = Ampl / absk;
        if(mu > 0.0) {
            double buf = C2_3*mu*absk;
            if(buf <= 1.0) { // вязкость малая или нулевая
                reOmega = absk * sqrt(1.0 - buf*buf);
                imOmega = absk * buf;
                ReKappa = - reOmega;
                ImKappa = imOmega;
            }
            else {
                reOmega = 0.0;
                imOmega = absk * (buf - sqrt(buf*buf - 1.0)); // берём медленно затухающую волну
                ReKappa = 0.0;
                ImKappa = absk * (buf + sqrt(buf*buf - 1.0));
            }
            _Ampl *= exp(- imOmega * t);
        }
        double kr = kx*coor[0] + ky*coor[1] + kz*coor[2];
        double c = cos(reOmega*t+kr), s = sin(reOmega*t+kr);
        V[Var_R] = _Ampl*(ReKappa*c - ImKappa*s);
        V[Var_U] = _Ampl*kx*c;
        V[Var_V] = _Ampl*ky*c;
        V[Var_W] = _Ampl*kz*c;
        V[Var_P] = V[Var_R];

        // Добавляем вихревые моды
        _Ampl = exp(-mu*absk*absk*t) * cos(kr);
        V[Var_U] += _Ampl*(AmplVY*kz - AmplVZ*ky);
        V[Var_V] += _Ampl*(AmplVZ*kx - AmplVX*kz);
        V[Var_W] += _Ampl*(AmplVX*ky - AmplVY*kx);
    }
}
//======================================================================================================================
