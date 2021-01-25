// *********************************************************************************************************************
// *****                                                                                                           *****
// *****                                            Special functions                                              *****
// *****                                       Source: Cephes Math Library                                         *****
// *****                                                                                                           *****
// ***** Cephes Math Library Release 2.8: June, 2000. Copyright 1984, 1987, 1988, 1992, 2000 by Stephen L. Moshier *****
// *********************************************************************************************************************

#include <math.h>
#ifdef double
    #undef double
#endif

// Common subroutines and constants
#define PI     3.14159265358979323846       /* pi */
#define PIO4   7.85398163397448309616E-1    /* pi/4 */
#define SQ2OPI 7.9788456080286535587989E-1  /* sqrt( 2/pi ) */
#define TWOOPI 6.36619772367581343075535E-1 /* 2/pi */
#define THPIO4 2.35619449019234492885       /* 3*pi/4 */
#define MACHEP 1.11022302462515654042E-16   /* 2**-53 */
#define MAXLOG 7.09782712893383996732E2     /* log(MAXNUM) */
#define SQRTH  7.07106781186547524401E-1    /* sqrt(2)/2 */

static const unsigned char NaN_char[8] = {0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF, 0xFF};
static const double* NaN = (const double*) NaN_char;

//======================================================================================================================
// Evaluates polynomial of degree N:  y  =  C0  + C1 x + C2 x^2  +...+ CN x^N
// Coefficients are stored in reverse order:  coef[0] = CN  , ..., coef[N] = C0
// The function p1evl() assumes that coef[N] = 1.0 and is omitted from the array.
// Its calling arguments are otherwise the same as polevl().
//======================================================================================================================
template<typename fpv> inline fpv polevl(fpv x, const fpv* coef, int N){
    const fpv *p = coef;
    fpv ans = *p++;
    int i = N;

    do ans = ans * x  +  *p++;
    while( --i );
    return ans;
}

//======================================================================================================================
// The same for two fraction of two polynomials (~30% faster than separate computations of numerator and denominator)
//======================================================================================================================
template<typename fpv> inline fpv div2polevl(fpv x, const fpv* coefN, const fpv* coefD, int N){
    const fpv *pN = coefN;
    const fpv *pD = coefD;
    fpv ansN = *pN++;
    fpv ansD = *pD++;
    int i = N;

    do {
        ansN = ansN * x  +  *pN++;
        ansD = ansD * x  +  *pD++;
    } while( --i );
    return ansN / ansD;
}

//======================================================================================================================
// Evaluate polynomial when coefficient of x^N is 1.0.
// Otherwise same as polevl.
//======================================================================================================================
template<typename fpv> inline fpv p1evl(fpv x, const fpv* coef, int N){
    const fpv *p = coef;
    fpv ans = x + *p++;
    int i = N-1;

    do ans = ans * x  + *p++;
    while( --i );
    return ans;
}


static double PP0[7] = {
  7.96936729297347051624E-4, 8.28352392107440799803E-2, 1.23953371646414299388E0, 5.44725003058768775090E0,
  8.74716500199817011941E0,  5.30324038235394892183E0,  9.99999999999999997821E-1};
static double PQ0[7] = {
  9.24408810558863637013E-4, 8.56288474354474431428E-2, 1.25352743901058953537E0, 5.47097740330417105182E0,
  8.76190883237069594232E0,  5.30605288235394617618E0,  1.00000000000000000218E0};
static double QP0[8] = {
 -1.13663838898469149931E-2,-1.28252718670509318512E0, -1.95539544257735972385E1,-9.32060152123768231369E1,
 -1.77681167980488050595E2, -1.47077505154951170175E2, -5.14105326766599330220E1,-6.05014350600728481186E0};
static double QQ0[7] = {
   6.43178256118178023184E1, 8.56430025976980587198E2,  3.88240183605401609683E3, 7.24046774195652478189E3,
   5.93072701187316984827E3, 2.06209331660327847417E3,  2.42005740240291393179E2};
static double QQ0_[8] = {
   1.0,                       6.43178256118178023184E1, 8.56430025976980587198E2,  3.88240183605401609683E3,
   7.24046774195652478189E3,  5.93072701187316984827E3, 2.06209331660327847417E3,  2.42005740240291393179E2};
static double YP0[8] = {
   1.55924367855235737965E4,-1.46639295903971606143E7,  5.43526477051876500413E9,-9.82136065717911466409E11,
  8.75906394395366999549E13,-3.46628303384729719441E15, 4.42733268572569800351E16,-1.84950800436986690637E16};
static double YQ0[7] = {
   1.04128353664259848412E3, 6.26107330137134956842E5,  2.68919633393814121987E8, 8.64002487103935000337E10,
  2.02979612750105546709E13, 3.17157752842975028269E15, 2.50596256172653059228E17};
static double DR1 = 5.78318596294678452118E0; //  5.783185962946784521175995758455807035071
static double DR2 = 3.04712623436620863991E1; // 30.47126234366208639907816317502275584842

static double RP0[4] = {
    -4.79443220978201773821E9, 1.95617491946556577543E12, -2.49248344360967716204E14, 9.70862251047306323952E15};
static double RQ0[8] = {
     4.99563147152651017219E2, 1.73785401676374683123E5, 4.84409658339962045305E7, 1.11855537045356834862E10,
    2.11277520115489217587E12, 3.10518229857422583814E14, 3.18121955943204943306E16, 1.71086294081043136091E18};


//======================================================================================================================
// Cephes Math Library Release 2.8:  June, 2000. Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
// Returns Bessel function of order zero of the argument.
// The domain is divided into the intervals [0, 5] and (5, infinity).
// In the first interval the following rational approximation is used:
// (x - r_1^2) (x - r_2^2) P_3 (x) / Q_8 (x)
// where r_1 and r_2 are zeros of the function.
// In the second interval, the Hankel asymptotic expansion
// is employed with two rational functions of degree 6/6 and 7/7.
// ACCURACY: Absolute error:
// arithmetic   domain     # trials      peak         rms
//    IEEE      0, 30       60000       4.2e-16     1.1e-16
//======================================================================================================================
double cephes_j0(double x) {
    if( x < 0 )    x = -x;
    if( x <= 5.0 ){
        double z = x * x;
        if( x < 1.0e-5 )
            return( 1.0 - 0.25*z );

        double p = (z - DR1) * (z - DR2);
        p *= polevl( z, RP0, 3)/p1evl( z, RQ0, 8 );
        return p;
    }

    double w = 5.0/x;
    double q = w*w;
    double p = div2polevl( q, PP0, PQ0, 6 );
           q = div2polevl( q, QP0, QQ0_, 7 );
    double xn = x - PIO4;
           p = p * cos(xn) - w * q * sin(xn);
    return( p * SQ2OPI / sqrt(x) );
}

//======================================================================================================================
// Cephes Math Library Release 2.8:  June, 2000. Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
// Returns Bessel function of the second kind, of order zero, of the argument.
// The domain is divided into the intervals [0, 5] and (5, infinity).
// In the first interval a rational approximation R(x) is employed to compute
// y0(x)  = R(x)  +   2 * log(x) * j0(x) / PI.
// Thus a call to j0() is required.
// In the second interval, the Hankel asymptotic expansion
// is employed with two rational functions of degree 6/6 and 7/7.
// ACCURACY: Absolute error, when y0(x) < 1; else relative error:
// arithmetic   domain     # trials      peak         rms
//    IEEE      0, 30       30000       1.3e-15     1.6e-16
//======================================================================================================================
double cephes_y0(double x) {
    double w, z, p, q, xn;
    if( x <= 5.0 ) {
        if( x <= 0.0 ) return *NaN;
        z = x * x;
        // calculation of w(x) =  y0(x)  -  2 * log(x) * j0(x) / PI,
        // w(0) = 2 * ( log(0.5) + EUL ) / PI = 0.073804295108687225.
        w = polevl( z, YP0, 7) / p1evl( z, YQ0, 7 );
        w += TWOOPI * log(x) * cephes_j0(x);
        return( w );
    }
    w = 5.0/x;
    z = 25.0 / (x * x);
    p = polevl( z, PP0, 6)/polevl( z, PQ0, 6 );
    q = polevl( z, QP0, 7)/p1evl( z, QQ0, 7 );
    xn = x - PIO4;
    p = p * sin(xn) + w * q * cos(xn);
    return( p * SQ2OPI / sqrt(x) );
}


static double RP1[4] =
    {-8.99971225705559398224E8, 4.52228297998194034323E11,-7.27494245221818276015E13, 3.68295732863852883286E15};
static double RQ1[8] = {
    6.20836478118054335476E2, 2.56987256757748830383E5, 8.35146791431949253037E7, 2.21511595479792499675E10,
    4.74914122079991414898E12, 7.84369607876235854894E14, 8.95222336184627338078E16, 5.32278620332680085395E18};
static double PP1[7] = {
    7.62125616208173112003E-4, 7.31397056940917570436E-2, 1.12719608129684925192E0, 5.11207951146807644818E0,
    8.42404590141772420927E0, 5.21451598682361504063E0, 1.00000000000000000254E0};
static double PQ1[7] = {
    5.71323128072548699714E-4, 6.88455908754495404082E-2, 1.10514232634061696926E0, 5.07386386128601488557E0,
 8.39985554327604159757E0, 5.20982848682361821619E0, 9.99999999999999997461E-1};
static double QP1[8] = {
    5.10862594750176621635E-2, 4.98213872951233449420E0, 7.58238284132545283818E1, 3.66779609360150777800E2,
    7.10856304998926107277E2, 5.97489612400613639965E2, 2.11688757100572135698E2, 2.52070205858023719784E1};
static double QQ1[7] = {
    7.42373277035675149943E1, 1.05644886038262816351E3, 4.98641058337653607651E3, 9.56231892404756170795E3,
    7.99704160447350683650E3, 2.82619278517639096600E3, 3.36093607810698293419E2};
static double YP1[6] = {
    1.26320474790178026440E9,-6.47355876379160291031E11, 1.14509511541823727583E14,
    -8.12770255501325109621E15, 2.02439475713594898196E17,-7.78877196265950026825E17};
static double YQ1[8] = {
    5.94301592346128195359E2, 2.35564092943068577943E5, 7.34811944459721705660E7, 1.87601316108706159478E10,
    3.88231277496238566008E12, 6.20557727146953693363E14, 6.87141087355300489866E16, 3.97270608116560655612E18};
static double Z1 = 1.46819706421238932572E1;
static double Z2 = 4.92184563216946036703E1;

//======================================================================================================================
// Cephes Math Library Release 2.8:  June, 2000. Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
// Returns Bessel function of order one of the argument
//======================================================================================================================
double cephes_j1(double x) {
    double w, z, p, q, xn;
    w = x;
    if( x < 0 )    w = -x;

    if( w <= 5.0 ){
        z = x * x;
        w = polevl( z, RP1, 3 ) / p1evl( z, RQ1, 8 );
        w = w * x * (z - Z1) * (z - Z2);
        return( w );
    }

    w = 5.0/x;
    z = w * w;
    p = polevl( z, PP1, 6)/polevl( z, PQ1, 6 );
    q = polevl( z, QP1, 7)/p1evl( z, QQ1, 7 );
    xn = x - THPIO4;
    p = p * cos(xn) - w * q * sin(xn);
    return( p * SQ2OPI / sqrt(x) );
}

//======================================================================================================================
// Cephes Math Library Release 2.8:  June, 2000. Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
// Returns Bessel function of the second kind of order one of the argument
//======================================================================================================================
double cephes_y1(double x) {
    double w, z, p, q, xn;
    if( x <= 5.0 ){
        if( x <= 0.0 ) return *NaN;
        z = x * x;
        w = x * (polevl( z, YP1, 5 ) / p1evl( z, YQ1, 8 ));
        w += TWOOPI * ( cephes_j1(x) * log(x)  -  1.0/x );
        return( w );
    }

    w = 5.0/x;
    z = w * w;
    p = polevl( z, PP1, 6)/polevl( z, PQ1, 6 );
    q = polevl( z, QP1, 7)/p1evl( z, QQ1, 7 );
    xn = x - THPIO4;
    p = p * sin(xn) + w * q * cos(xn);
    return( p * SQ2OPI / sqrt(x) );
}


//======================================================================================================================
// Cephes Math Library Release 2.8:  June, 2000. Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
// Fresnel integrals C(x) = int_0^x cos(pi/2 t**2) dt, S(x) = int_0^x sin(pi/2 t**2) dt
// The integrals are evaluated by a power series for x < 1. 
// For x >= 1 auxiliary functions f(x) and g(x) are employed such that
// C(x) = 0.5 + f(x) sin( pi/2 x**2 ) - g(x) cos( pi/2 x**2 )
// S(x) = 0.5 - f(x) cos( pi/2 x**2 ) - g(x) sin( pi/2 x**2 )
// Acuracy: Arithmetic  function   domain     # trials      peak         rms
//            IEEE       S(x)      0, 10       10000       2.0e-15     3.2e-16
//            IEEE       C(x)      0, 10       10000       1.8e-15     3.3e-16
//======================================================================================================================

/* S(x) for small x */
static double sn[6] = {
-2.99181919401019853726E3,
 7.08840045257738576863E5,
-6.29741486205862506537E7,
 2.54890880573376359104E9,
-4.42979518059697779103E10,
 3.18016297876567817986E11,
};
static double sd[6] = {
/* 1.00000000000000000000E0,*/
 2.81376268889994315696E2,
 4.55847810806532581675E4,
 5.17343888770096400730E6,
 4.19320245898111231129E8,
 2.24411795645340920940E10,
 6.07366389490084639049E11,
};

/* C(x) for small x */
static double cn[6] = {
-4.98843114573573548651E-8,
 9.50428062829859605134E-6,
-6.45191435683965050962E-4,
 1.88843319396703850064E-2,
-2.05525900955013891793E-1,
 9.99999999999999998822E-1,
};
static double cd[7] = {
 3.99982968972495980367E-12,
 9.15439215774657478799E-10,
 1.25001862479598821474E-7,
 1.22262789024179030997E-5,
 8.68029542941784300606E-4,
 4.12142090722199792936E-2,
 1.00000000000000000118E0,
};

/* Auxiliary function f(x) */
static double fn[10] = {
  4.21543555043677546506E-1,
  1.43407919780758885261E-1,
  1.15220955073585758835E-2,
  3.45017939782574027900E-4,
  4.63613749287867322088E-6,
  3.05568983790257605827E-8,
  1.02304514164907233465E-10,
  1.72010743268161828879E-13,
  1.34283276233062758925E-16,
  3.76329711269987889006E-20,
};
static double fd[10] = {
/*  1.00000000000000000000E0,*/
  7.51586398353378947175E-1,
  1.16888925859191382142E-1,
  6.44051526508858611005E-3,
  1.55934409164153020873E-4,
  1.84627567348930545870E-6,
  1.12699224763999035261E-8,
  3.60140029589371370404E-11,
  5.88754533621578410010E-14,
  4.52001434074129701496E-17,
  1.25443237090011264384E-20,
};

/* Auxiliary function g(x) */
static double gn[11] = {
  5.04442073643383265887E-1,
  1.97102833525523411709E-1,
  1.87648584092575249293E-2,
  6.84079380915393090172E-4,
  1.15138826111884280931E-5,
  9.82852443688422223854E-8,
  4.45344415861750144738E-10,
  1.08268041139020870318E-12,
  1.37555460633261799868E-15,
  8.36354435630677421531E-19,
  1.86958710162783235106E-22,
};
static double gd[11] = {
/*  1.00000000000000000000E0,*/
  1.47495759925128324529E0,
  3.37748989120019970451E-1,
  2.53603741420338795122E-2,
  8.14679107184306179049E-4,
  1.27545075667729118702E-5,
  1.04314589657571990585E-7,
  4.60680728146520428211E-10,
  1.10273215066240270757E-12,
  1.38796531259578871258E-15,
  8.39158816283118707363E-19,
  1.86958710162783236342E-22,
};

void fresnl(double xxa, double& ss, double& cc) {
    double x = fabs(xxa);
    double x2 = x * x;
    if( x2 < 2.5625 ) {
        double t = x2 * x2;
        ss = x * x2 * polevl( t, sn, 5)/p1evl( t, sd, 6 );
        cc = x * polevl( t, cn, 5)/polevl(t, cd, 6 );
    }
    else if( x > 36974.0 ) {
        cc = 0.5;
        ss = 0.5;
    }
    else {
        // Asymptotic power series auxiliary functions for large argument
        double t = PI * x2;
        double u = 1.0/(t * t);

        t = 1.0/t;
        double f = 1.0 - u * polevl( u, fn, 9)/p1evl(u, fd, 10);
        double g = t * polevl( u, gn, 10)/p1evl(u, gd, 11);

        t = 0.5*PI * x2;
        double c = cos(t);
        double s = sin(t);
        t = PI * x;
        cc = 0.5  +  (f * s  -  g * c)/t;
        ss = 0.5  -  (f * c  +  g * s)/t;
    }

    if( xxa < 0.0 )    {
        cc = -cc;
        ss = -ss;
    }
}
//======================================================================================================================




/*                            gamma.c
 *
 *    Gamma function
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, gamma();
 * extern int sgngam;
 *
 * y = gamma( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns gamma function of the argument.  The result is
 * correctly signed, and the sign (+1 or -1) is also
 * returned in a global (extern) variable named sgngam.
 * This variable is also filled in by the logarithmic gamma
 * function lgam().
 *
 * Arguments |x| <= 34 are reduced by recurrence and the function
 * approximated by a rational function of degree 6/7 in the
 * interval (2,3).  Large arguments are handled by Stirling's
 * formula. Large negative arguments are made positive using
 * a reflection formula.
 *
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    IEEE    -170,-33      20000       2.3e-15     3.3e-16
 *    IEEE     -33,  33     20000       9.4e-16     2.2e-16
 *    IEEE      33, 171.6   20000       2.3e-15     3.2e-16
 *
 * Error for arguments outside the test range will be larger
 * owing to error amplification by the exponential function.
 *
 */
/*                            lgam()
 *
 *    Natural logarithm of gamma function
 *
 *
 *
 * SYNOPSIS:
 *
 * double x, y, lgam();
 * extern int sgngam;
 *
 * y = lgam( x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Returns the base e (2.718...) logarithm of the absolute
 * value of the gamma function of the argument.
 * The sign (+1 or -1) of the gamma function is returned in a
 * global (extern) variable named sgngam.
 *
 * For arguments greater than 13, the logarithm of the gamma
 * function is approximated by the logarithmic version of
 * Stirling's formula using a polynomial approximation of
 * degree 4. Arguments between -33 and +33 are reduced by
 * recurrence to the interval [2,3] of a rational approximation.
 * The cosecant reflection formula is employed for arguments
 * less than -33.
 *
 * Arguments greater than MAXLGM return MAXNUM and an error
 * message.  MAXLGM = 2.035093e36 for DEC
 * arithmetic or 2.556348e305 for IEEE arithmetic.
 *
 *
 *
 * ACCURACY:
 *
 *
 * arithmetic      domain        # trials     peak         rms
 *    DEC     0, 3                  7000     5.2e-17     1.3e-17
 *    DEC     2.718, 2.035e36       5000     3.9e-17     9.9e-18
 *    IEEE    0, 3                 28000     5.4e-16     1.1e-16
 *    IEEE    2.718, 2.556e305     40000     3.5e-16     8.3e-17
 * The error criterion was relative when the function magnitude
 * was greater than one but absolute when it was less than one.
 *
 * The following test used the relative error criterion, though
 * at certain points the relative error could be much higher than
 * indicated.
 *    IEEE    -200, -4             10000     4.8e-16     1.3e-16
 *
 */

/*                            gamma.c    */
/*    gamma function    */

/*
Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1989, 1992, 2000 by Stephen L. Moshier
*/

static double P[] = {
  1.60119522476751861407E-4,
  1.19135147006586384913E-3,
  1.04213797561761569935E-2,
  4.76367800457137231464E-2,
  2.07448227648435975150E-1,
  4.94214826801497100753E-1,
  9.99999999999999996796E-1
};
static double Q[] = {
-2.31581873324120129819E-5,
 5.39605580493303397842E-4,
-4.45641913851797240494E-3,
 1.18139785222060435552E-2,
 3.58236398605498653373E-2,
-2.34591795718243348568E-1,
 7.14304917030273074085E-2,
 1.00000000000000000320E0
};
#define MAXGAM 171.624376956302725
static double LOGPI = 1.14472988584940017414;

/* Stirling's formula for the gamma function */
static double STIR[5] = {
 7.87311395793093628397E-4,
-2.29549961613378126380E-4,
-2.68132617805781232825E-3,
 3.47222221605458667310E-3,
 8.33333333333482257126E-2,
};
#define MAXSTIR 143.01608
static double SQTPI = 2.50662827463100050242E0;

/* Gamma function computed by Stirling's formula.
 * The polynomial STIR is valid for 33 <= x <= 172.
 */
static double stirf(double x) {
double y, w, v;

w = 1.0/x;
w = 1.0 + w * polevl( w, STIR, 4 );
y = exp(x);
if( x > MAXSTIR )
    { /* Avoid overflow in pow() */
    v = pow( x, 0.5 * x - 0.25 );
    y = v * (v / y);
    }
else
    {
    y = pow( x, x - 0.5 ) / y;
    }
y = SQTPI * y * w;
return( y );
}


double gamma(double x) {
double p, q, z;
int i;

int sgngam = 1;
q = fabs(x);

if( q > 33.0 )
    {
    if( x < 0.0 )
        {
        p = floor(q);
        if( p == q ) return *NaN;
        i = int(p);
        if( (i & 1) == 0 )
            sgngam = -1;
        z = q - p;
        if( z > 0.5 )
            {
            p += 1.0;
            z = q - p;
            }
        z = q * sin( PI * z );
        if( z == 0.0 ) return *NaN;
        z = fabs(z);
        z = PI/(z * stirf(q) );
        }
    else
        {
        z = stirf(x);
        }
    return( sgngam * z );
    }

z = 1.0;
while( x >= 3.0 )
    {
    x -= 1.0;
    z *= x;
    }

while( x < 0.0 )
    {
    if( x > -1.E-9 )
        goto small;
    z /= x;
    x += 1.0;
    }

while( x < 2.0 )
    {
    if( x < 1.e-9 )
        goto small;
    z /= x;
    x += 1.0;
    }

if( x == 2.0 )
    return(z);

x -= 2.0;
p = polevl( x, P, 6 );
q = polevl( x, Q, 7 );
return( z * p / q );

small:
if( x == 0.0 ) return *NaN;
else
    return( z/((1.0 + 0.5772156649015329 * x) * x) );
}



/* A[]: Stirling's formula expansion of log gamma
 * B[], C[]: log gamma function between 2 and 3
 */
static double A[] = {
 8.11614167470508450300E-4,
-5.95061904284301438324E-4,
 7.93650340457716943945E-4,
-2.77777777730099687205E-3,
 8.33333333333331927722E-2
};
static double B[] = {
-1.37825152569120859100E3,
-3.88016315134637840924E4,
-3.31612992738871184744E5,
-1.16237097492762307383E6,
-1.72173700820839662146E6,
-8.53555664245765465627E5
};
static double C[] = {
/* 1.00000000000000000000E0, */
-3.51815701436523470549E2,
-1.70642106651881159223E4,
-2.20528590553854454839E5,
-1.13933444367982507207E6,
-2.53252307177582951285E6,
-2.01889141433532773231E6
};
/* log( sqrt( 2*pi ) ) */
static double LS2PI  =  0.91893853320467274178;
#define MAXLGM 2.556348e305

/* Logarithm of gamma function */
double lgam(double x) {
double p, q, u, w, z;
//int i;

//int sgngam = 1;
if( x < -34.0 )
    {
    q = -x;
    w = lgam(q); /* note this modifies sgngam! */
    p = floor(q);
    if( p == q ) return *NaN;
    //i = int(p);
    //if( (i & 1) == 0 )
    //    sgngam = -1;
    //else
    //    sgngam = 1;
    z = q - p;
    if( z > 0.5 )
        {
        p += 1.0;
        z = p - q;
        }
    z = q * sin( PI * z );
    if( z == 0.0 ) return *NaN;
    z = LOGPI - log( z ) - w;
    return( z );
    }

if( x < 13.0 )
    {
    z = 1.0;
    p = 0.0;
    u = x;
    while( u >= 3.0 )
        {
        p -= 1.0;
        u = x + p;
        z *= u;
        }
    while( u < 2.0 )
        {
        if( u == 0.0 ) return *NaN;
        z /= u;
        p += 1.0;
        u = x + p;
        }
    if( z < 0.0 )
        {
        //sgngam = -1;
        z = -z;
        }
    //else
    //   sgngam = 1;
    if( u == 2.0 )
        return( log(z) );
    p -= 2.0;
    x = x + p;
    p = x * polevl( x, B, 5 ) / p1evl( x, C, 6);
    return( log(z) + p );
    }

if( x > MAXLGM ) return *NaN;
q = ( x - 0.5 ) * log(x) - x + LS2PI;
if( x > 1.0e8 )
    return( q );

p = 1.0/(x*x);
if( x >= 1000.0 )
    q += ((   7.9365079365079365079365e-4 * p
        - 2.7777777777777777777778e-3) *p
        + 0.0833333333333333333333) / x;
else
    q += polevl( p, A, 4 ) / x;
return( q );
}


/*                            expn.c
 *
 *        Exponential integral En
 *
 *
 *
 * SYNOPSIS:
 *
 * int n;
 * double x, y, expn();
 *
 * y = expn( n, x );
 *
 *
 *
 * DESCRIPTION:
 *
 * Evaluates the exponential integral
 *
 *                 inf.
 *                   -
 *                  | |   -xt
 *                  |    e
 *      E (x)  =    |    ----  dt.
 *       n          |      n
 *                | |     t
 *                 -
 *                  1
 *
 *
 * Both n and x must be nonnegative.
 *
 * The routine employs either a power series, a continued
 * fraction, or an asymptotic formula depending on the
 * relative values of n and x.
 *
 * ACCURACY:
 *
 *                      Relative error:
 * arithmetic   domain     # trials      peak         rms
 *    DEC       0, 30        5000       2.0e-16     4.6e-17
 *    IEEE      0, 30       10000       1.7e-15     3.6e-16
 *
 */

/*                            expn.c    */

/* Cephes Math Library Release 2.8:  June, 2000
   Copyright 1985, 2000 by Stephen L. Moshier */

#define EUL 0.57721566490153286060
#define BIG  1.44115188075855872E+17

double expn( int n, double x ) {
double ans, r, t, yk, xk;
double pk, pkm1, pkm2, qk, qkm1, qkm2;
double psi, z;
int i, k;
static double big = BIG;

if( n < 0 ) return *NaN;
if( x < 0 ) return *NaN;

if( x > MAXLOG )
    return( 0.0 );

if( x == 0.0 )
    {
    if( n < 2 ) return *NaN;
    else
        return( 1.0/(n-1.0) );
    }

if( n == 0 )
    return( exp(-x)/x );

/*                            expn.c    */
/*        Expansion for large n        */

if( n > 5000 )
    {
    xk = x + n;
    yk = 1.0 / (xk * xk);
    t = n;
    ans = yk * t * (6.0 * x * x  -  8.0 * t * x  +  t * t);
    ans = yk * (ans + t * (t  -  2.0 * x));
    ans = yk * (ans + t);
    ans = (ans + 1.0) * exp( -x ) / xk;
    goto done;
    }

if( x > 1.0 )
    goto cfrac;

/*                            expn.c    */

/*        Power series expansion        */

psi = -EUL - log(x);
for( i=1; i<n; i++ )
    psi = psi + 1.0/i;

z = -x;
xk = 0.0;
yk = 1.0;
pk = 1.0 - n;
if( n == 1 )
    ans = 0.0;
else
    ans = 1.0/pk;
do
    {
    xk += 1.0;
    yk *= z/xk;
    pk += 1.0;
    if( pk != 0.0 )
        {
        ans += yk/pk;
        }
    if( ans != 0.0 )
        t = fabs(yk/ans);
    else
        t = 1.0;
    }
while( t > MACHEP );
k = int(xk);
t = n;
r = n - 1;
ans = (pow(z, r) * psi / gamma(t)) - ans;
goto done;

/*                            expn.c    */
/*        continued fraction        */
cfrac:
k = 1;
pkm2 = 1.0;
qkm2 = x;
pkm1 = 1.0;
qkm1 = x + n;
ans = pkm1/qkm1;

do
    {
    k += 1;
    if( k & 1 )
        {
        yk = 1.0;
        xk = n + (k-1)/2;
        }
    else
        {
        yk = x;
        xk = k/2;
        }
    pk = pkm1 * yk  +  pkm2 * xk;
    qk = qkm1 * yk  +  qkm2 * xk;
    if( qk != 0 )
        {
        r = pk/qk;
        t = fabs( (ans - r)/r );
        ans = r;
        }
    else
        t = 1.0;
    pkm2 = pkm1;
    pkm1 = pk;
    qkm2 = qkm1;
    qkm1 = qk;
if( fabs(pk) > big )
        {
        pkm2 /= big;
        pkm1 /= big;
        qkm2 /= big;
        qkm1 /= big;
        }
    }
while( t > MACHEP );

ans *= exp( -x );

done:
return( ans );
}


//======================================================================================================================
// Cephes Math Library Release 2.8:  June, 2000. Copyright 1984, 1987, 1988, 1992, 2000 by Stephen L. Moshier
// Returns erf(x) = 2/sqrt(pi) \int_0^x exp(-t*t) dt.
// For 0 <= |x| < 1, erf(x) = x * P4(x**2)/Q5(x**2); otherwise erf(x) = 1 - erfc(x).
// Acuracy (relative error): Arithmetic  domain     # trials      peak         rms
//                              IEEE     0,1         30000       3.7e-16     1.0e-16
//
// Returns erfc(x) = 1.0-erf(x)
// For small x, erfc(x) = 1 - erf(x); otherwise rational approximations are computed.
// Acuracy (relative error): Arithmetic  domain     # trials      peak         rms
//                              IEEE     0,26.6417   30000       5.7e-14     1.5e-14
//======================================================================================================================

static double PP[] = {
 2.46196981473530512524E-10,
 5.64189564831068821977E-1,
 7.46321056442269912687E0,
 4.86371970985681366614E1,
 1.96520832956077098242E2,
 5.26445194995477358631E2,
 9.34528527171957607540E2,
 1.02755188689515710272E3,
 5.57535335369399327526E2
};
static double QQ[] = {
/* 1.00000000000000000000E0,*/
 1.32281951154744992508E1,
 8.67072140885989742329E1,
 3.54937778887819891062E2,
 9.75708501743205489753E2,
 1.82390916687909736289E3,
 2.24633760818710981792E3,
 1.65666309194161350182E3,
 5.57535340817727675546E2
};
static double R[] = {
 5.64189583547755073984E-1,
 1.27536670759978104416E0,
 5.01905042251180477414E0,
 6.16021097993053585195E0,
 7.40974269950448939160E0,
 2.97886665372100240670E0
};
static double S[] = {
/* 1.00000000000000000000E0,*/
 2.26052863220117276590E0,
 9.39603524938001434673E0,
 1.20489539808096656605E1,
 1.70814450747565897222E1,
 9.60896809063285878198E0,
 3.36907645100081516050E0
};
static double T[] = {
 9.60497373987051638749E0,
 9.00260197203842689217E1,
 2.23200534594684319226E3,
 7.00332514112805075473E3,
 5.55923013010394962768E4
};
static double U[] = {
/* 1.00000000000000000000E0,*/
 3.35617141647503099647E1,
 5.21357949780152679795E2,
 4.59432382970980127987E3,
 2.26290000613890934246E4,
 4.92673942608635921086E4
};

double cephes_erf(double x);
double cephes_erfc(double a) {
    double p,q,x,y,z;
    if( a < 0.0 ) x = -a;
    else x = a;

    if( x < 1.0 ) return( 1.0 - cephes_erf(a) );
    z = -a * a;
    if( z < -MAXLOG ) return (a<0.0 ? 2.0 : 0.0);
    z = exp(z);

    if( x < 8.0 ) {
        p = polevl( x, PP, 8 );
        q = p1evl( x, QQ, 8 );
    }
    else {
        p = polevl( x, R, 5 );
        q = p1evl( x, S, 6 );
    }
    y = (z * p)/q;

    if( a < 0 ) y = 2.0 - y;
    if( y == 0.0 ) return (a<0.0 ? 2.0 : 0.0);
    return(y);
}

double cephes_erf(double x) {
    if( fabs(x) > 1.0 ) return( 1.0 - cephes_erfc(x) );
    double z = x * x;
    double y = x * polevl( z, T, 4 ) / p1evl( z, U, 5 );
    return( y );
}
