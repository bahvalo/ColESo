#include "base_lib.h"
#include "base_const.h"
#include <math.h>
#include "es_utils.h"
#ifdef EXTRAPRECISION_COLESO
#include "extraprecision.h"
#endif

//======================================================================================================================
//
//                 Integration of int f(x) dx or int f(x) / sqrt(x) dx using compound Gauss formulae
//
//======================================================================================================================


//======================================================================================================================
// Calculation using compound Gauss formula of the integral
// int exp(-alpha^2 (x-x0)^2/2) (x-x0)^k h(x) x^gamma dx, x = xmin..xmax,
// where k=0,1,... (not big), alpha>0, gamma=0 (mode=0) or -1/2 (mode=1),
// h -- analytic function on [xmin..xmax] with Linf-norm M and convergence raduis >=r at each point, q = 1/r.
// If mode = 1, xmin should be zero
//======================================================================================================================
template<typename fpv> template<unsigned int N>
tFixBlock<fpv,N> tCompoundGaussIntegrator<fpv>::Integrate(fpv xmin, fpv xmax, fpv alpha, fpv x0, int k, int mode,
                                             fpv M, fpv q, tFixBlock<fpv,N> (*func)(fpv, void*), void* args) const {
    static const fpv four_over_e = 1.471517764685769286382; // extra precision for estimates is not in need
    if(mode && fabs(xmin) > 1e-50) crash("tCompoundGaussIntegrator: mode=1 but xmin!=0");
    if(alpha<0.0) alpha = -alpha;
    if(alpha<1e-50 || q<0.0) crash("tCompoundGaussIntegrator: alpha=0 or q<0");
    if(GLI.GN==NULL) crash("tCompoundGaussIntegrator: Init not done");
    const fpv inv_alpha = 1.0 / alpha;
    q = fpv(MAX(q, alpha) * reducer);

    // Gaussian truncating
    fpv beta = log(11.0*M*(mode ? sqrt(alpha) : alpha) / get_eps<fpv>());
    fpv H = 1.41;
    if(beta>1.0) {
        beta = sqrt(beta);
        H *= (beta + log(beta) / beta);
    }
    fpv xminus = x0 - H*inv_alpha;
    fpv xplus  = x0 + H*inv_alpha;
    xplus = MIN(xplus, xmax);
    xminus = MAX(xminus, xmin);

    // Gauss -- Jacobi quadrature -- if in need
    tFixBlock<fpv,N> sum;
    sum = 0.0;
    if(mode) {
        fpv xminusmin = four_over_e/q;
        if(xminus<xminusmin) {
            if(GJI.GN==NULL) crash("tCompoundGaussIntegrator: Init not done");
            xminus = xminusmin;
            if(xminus > xmax) xminus = xmax;
            for(int j=0; j<GJI.GR; j++) {
                fpv x = GJI.GN[j]*xminus;
                fpv tmp = alpha * (x-x0);
                sum += func(x, args) * exp(-0.5*tmp*tmp) * pow(x-x0, double(k)) * GJI.GC[j];
            }
            sum *= sqrt(xminus);
        }
    }
    if(xplus <= xminus) return sum;

    // Gauss -- Legendre compound quadrature
    int NumK = int(q*(xplus - xminus)*0.68) + 1;
    fpv dx = (xplus - xminus) / NumK;
    for(int i=0; i<NumK; i++) {
        for(int j=0; j<GLI.GR; j++) {
            fpv x = xminus + (i + GLI.GN[j]) * dx;
            fpv tmp = alpha * (x-x0);
            fpv f = exp(-0.5*tmp*tmp) * pow(x-x0, double(k));
            if(mode) f /= sqrt(x);
            sum += func(x, args) * (f * GLI.GC[j] * dx);
        }
    }
    return sum;
}

//======================================================================================================================
// Instantiating functions of the template class
//======================================================================================================================
#define INSTANTIATE(T) \
    template struct tGaussIntegrator<T>; \
    template struct tCompoundGaussIntegrator<T>; \
    template tFixBlock<T,1> tCompoundGaussIntegrator<T>::Integrate<1>(T xmin, T xmax, T alpha, T x0, int k, int mode, \
                                             T M, T q, tFixBlock<T,1> (*func)(T, void*), void* args) const;   \
    template tFixBlock<T,2> tCompoundGaussIntegrator<T>::Integrate<2>(T xmin, T xmax, T alpha, T x0, int k, int mode, \
                                             T M, T q, tFixBlock<T,2> (*func)(T, void*), void* args) const;   \
    template tFixBlock<T,3> tCompoundGaussIntegrator<T>::Integrate<3>(T xmin, T xmax, T alpha, T x0, int k, int mode, \
                                             T M, T q, tFixBlock<T,3> (*func)(T, void*), void* args) const;

INSTANTIATE(NativeDouble)
#ifdef EXTRAPRECISION_COLESO
INSTANTIATE(dd_real)
INSTANTIATE(qd_real)
#endif
