// *********************************************************************************************************************
// *****                                                                                                           *****
// *****                                            Special functions                                              *****
// *****       Sources: Cephes Math Library, PORT3 subroutine library, Numerical analysis library of SRCC MSU      *****
// *****                                                                                                           *****
// *********************************************************************************************************************

#include "pointfuncs.h"
#include "es_utils.h"
#ifdef _NOISETTE
#include "lib_base.h"
#endif
#ifdef EXTRAPRECISION_COLESO
#include "extraprecision.h"
#endif
#include <math.h>

//#pragma warning(disable:6297)  //Arithmetic overflow: 32-bit value is shifted, then cast to 64-bit value... 


//======================================================================================================================
// Cephes Math Library -- common subroutines and constants
//======================================================================================================================
#define PIO4 .78539816339744830962
#define SQ2OPI .79788456080286535588
#define TWOOPI 0.63661977236758134307553505349006
#define THPIO4 2.35619449019234492885 // 3*pi/4

//======================================================================================================================
// Evaluates polynomial of degree N:  y  =  C0  + C1 x + C2 x^2  +...+ CN x^N
// Coefficients are stored in reverse order:  coef[0] = CN  , ..., coef[N] = C0
// The function p1evl() assumes that coef[N] = 1.0 and is omitted from the array.
// Its calling arguments are otherwise the same as polevl().
// Cephes Math Library Release 2.1:  December, 1988. Copyright 1984, 1987, 1988 by Stephen L. Moshier
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


static NativeDouble PP0[7] = {
  7.96936729297347051624E-4, 8.28352392107440799803E-2, 1.23953371646414299388E0, 5.44725003058768775090E0, 
  8.74716500199817011941E0,  5.30324038235394892183E0,  9.99999999999999997821E-1};
static NativeDouble PQ0[7] = {
  9.24408810558863637013E-4, 8.56288474354474431428E-2, 1.25352743901058953537E0, 5.47097740330417105182E0,
  8.76190883237069594232E0,  5.30605288235394617618E0,  1.00000000000000000218E0};
static NativeDouble QP0[8] = {
 -1.13663838898469149931E-2,-1.28252718670509318512E0, -1.95539544257735972385E1,-9.32060152123768231369E1,
 -1.77681167980488050595E2, -1.47077505154951170175E2, -5.14105326766599330220E1,-6.05014350600728481186E0};
static NativeDouble QQ0[7] = {
   6.43178256118178023184E1, 8.56430025976980587198E2,  3.88240183605401609683E3, 7.24046774195652478189E3,
   5.93072701187316984827E3, 2.06209331660327847417E3,  2.42005740240291393179E2};
static NativeDouble QQ0_[8] = {
   1.0,                       6.43178256118178023184E1, 8.56430025976980587198E2,  3.88240183605401609683E3, 
   7.24046774195652478189E3,  5.93072701187316984827E3, 2.06209331660327847417E3,  2.42005740240291393179E2};
static NativeDouble YP0[8] = {
   1.55924367855235737965E4,-1.46639295903971606143E7,  5.43526477051876500413E9,-9.82136065717911466409E11,
  8.75906394395366999549E13,-3.46628303384729719441E15, 4.42733268572569800351E16,-1.84950800436986690637E16};
static NativeDouble YQ0[7] = {
   1.04128353664259848412E3, 6.26107330137134956842E5,  2.68919633393814121987E8, 8.64002487103935000337E10,
  2.02979612750105546709E13, 3.17157752842975028269E15, 2.50596256172653059228E17};
static NativeDouble DR1 = 5.78318596294678452118E0; //  5.783185962946784521175995758455807035071
static NativeDouble DR2 = 3.04712623436620863991E1; // 30.47126234366208639907816317502275584842

static NativeDouble RP0[4] = {
    -4.79443220978201773821E9, 1.95617491946556577543E12, -2.49248344360967716204E14, 9.70862251047306323952E15};
static NativeDouble RQ0[8] = {
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
template<> NativeDouble BesselJ0(NativeDouble x) {
    if( x < 0 )    x = -x;
    if( x <= 5.0 ){
        NativeDouble z = x * x;
        if( x < 1.0e-5 )
            return( 1.0 - 0.25*z );

        NativeDouble p = (z - DR1) * (z - DR2);
        p *= polevl( z, RP0, 3)/p1evl( z, RQ0, 8 );
        return p;
    }

    NativeDouble w = 5.0/x;
    NativeDouble q = w*w;
    NativeDouble p = div2polevl( q, PP0, PQ0, 6 );
           q = div2polevl( q, QP0, QQ0_, 7 );
    NativeDouble xn = x - PIO4;
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
template<> NativeDouble BesselN0(NativeDouble x) {
    NativeDouble w, z, p, q, xn;
    if( x <= 5.0 ) {
        if( x <= 0.0 ) { MakeNaN(x); return x; }
        z = x * x;
        // calculation of w(x) =  y0(x)  -  2 * log(x) * j0(x) / PI,
        // w(0) = 2 * ( log(0.5) + EUL ) / PI = 0.073804295108687225.
        w = polevl( z, YP0, 7) / p1evl( z, YQ0, 7 );
        w += TWOOPI * log(x) * BesselJ0(x);
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


static NativeDouble RP1[4] = 
    {-8.99971225705559398224E8, 4.52228297998194034323E11,-7.27494245221818276015E13, 3.68295732863852883286E15};
static NativeDouble RQ1[8] = {
    6.20836478118054335476E2, 2.56987256757748830383E5, 8.35146791431949253037E7, 2.21511595479792499675E10,
    4.74914122079991414898E12, 7.84369607876235854894E14, 8.95222336184627338078E16, 5.32278620332680085395E18};
static NativeDouble PP1[7] = {
    7.62125616208173112003E-4, 7.31397056940917570436E-2, 1.12719608129684925192E0, 5.11207951146807644818E0,
    8.42404590141772420927E0, 5.21451598682361504063E0, 1.00000000000000000254E0};
static NativeDouble PQ1[7] = {
    5.71323128072548699714E-4, 6.88455908754495404082E-2, 1.10514232634061696926E0, 5.07386386128601488557E0,
 8.39985554327604159757E0, 5.20982848682361821619E0, 9.99999999999999997461E-1};
static NativeDouble QP1[8] = {
    5.10862594750176621635E-2, 4.98213872951233449420E0, 7.58238284132545283818E1, 3.66779609360150777800E2,
    7.10856304998926107277E2, 5.97489612400613639965E2, 2.11688757100572135698E2, 2.52070205858023719784E1};
static NativeDouble QQ1[7] = {
    7.42373277035675149943E1, 1.05644886038262816351E3, 4.98641058337653607651E3, 9.56231892404756170795E3,
    7.99704160447350683650E3, 2.82619278517639096600E3, 3.36093607810698293419E2};
static NativeDouble YP1[6] = {
    1.26320474790178026440E9,-6.47355876379160291031E11, 1.14509511541823727583E14,
    -8.12770255501325109621E15, 2.02439475713594898196E17,-7.78877196265950026825E17};
static NativeDouble YQ1[8] = {
    5.94301592346128195359E2, 2.35564092943068577943E5, 7.34811944459721705660E7, 1.87601316108706159478E10,
    3.88231277496238566008E12, 6.20557727146953693363E14, 6.87141087355300489866E16, 3.97270608116560655612E18};
static NativeDouble Z1 = 1.46819706421238932572E1;
static NativeDouble Z2 = 4.92184563216946036703E1;

//======================================================================================================================
// Cephes Math Library Release 2.8:  June, 2000. Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
// Returns Bessel function of order one of the argument
//======================================================================================================================
template<> NativeDouble BesselJ1(NativeDouble x) {
    NativeDouble w, z, p, q, xn;
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
template<> NativeDouble BesselN1(NativeDouble x){
    NativeDouble w, z, p, q, xn;
    if( x <= 5.0 ){
        if( x <= 0.0 ) { MakeNaN(x); return x; }
        z = x * x;
        w = x * (polevl( z, YP1, 5 ) / p1evl( z, YQ1, 8 ));
        w += TWOOPI * ( BesselJ1(x) * log(x)  -  1.0/x );
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
// Взято из Библиотеки Численного Анализа Научно-Исследовательского Вычислительного Центра МГУ
// Описание метода: В.Н.Родин, Стандартные программы вычисления полных эллиптических интегралов, интеграла вероятностей, 
//      гамма - функции и функции Бесселя, Сб. "Численный анализ на ФОРТРАНе", вып. 8, Изд - во МГУ, 1974 г.
// Вход:
//      x - значение аргумента
//      N - максимальный порядок функции плюс один
//      nmax = MAX(N, 2*[x]+4) - размер массивов fun1, fun2, work
//      fun1 = значения функции BesselJ(n,x) для n = 0..Size()-1
//      fun2 = значения функции BesselY(n,x) для n = 0..Size()-1
//      work = рабочий массив
// Вычисляемое значение:
//      При x>1 вычисляются значения Jn(x) и Nn(x)
//      При x<=1 вычисляются:
//          in (x) = n! (x / 2)^(-n) Jn (x) ,
//          k0 (x) = N0 (x) ,
//          kn (x) = 2 (x / 2)^n Nn (x) / (n - 1)!
//======================================================================================================================
int sf33r(NativeDouble x, int n, int nmax, NativeDouble *fun1, NativeDouble *fun2, NativeDouble *work) {
    static NativeDouble sys039 = 0.57721566490153286060; // Euler's gamma constant
    static NativeDouble sys015 = 2.0/PiNumber;
    static NativeDouble sys005 = 0.5*PiNumber;

    int i1, i2;

    int i;
    int j, k;
    NativeDouble p, _r, u, v;
    int n1, n2;
    NativeDouble p0, p1, q0, q1, q2, r1, r2, r3, r4, r5, u1, u2, u3, u4, u5,
            v1, v2, v3, v4, v5, v6, v7;
    int na, nc1, nc2;

    if (n < 2) return -1;
    i1 = n; i2 = ((int) (x) << 1) + 4;
    if (nmax < i1 || nmax < i2 || nmax < 8) return -2;
    if (x <= 0.0) return -3;

// Computing MAX
    nc1 = (int) MAX(x + 1.0, 3.0);
    nc2 = (nc1 << 1) + 1;
    if (n - nc1 <= 0)
        n1 = n2 = n;
    else {
        n1 = nc1;
        if (n - nc2 <= 0)
            n2 = nc2 + 1;
        else
            n2 = n;
    }

    if(x<=1.0) {
        r1 = x * .5;
        r2 = -(r1 * r1);
        u1 = 1.0;
        u2 = 1.0;
        u3 = 1.0;
        u4 = 1.0;
        p0 = 1.0;
        p1 = 1.0;
        q0 = -((NativeDouble)log(r1) + sys039);
        v1 = q0;
        q1 = v1 + .5;
        v2 = q1;
        for (i = 1; i <= 8; ++i) {
            u1 = u1 * r2 / (NativeDouble) ((n - 1 + i) * i);
            p0 += u1;
            u2 = u2 * r2 / (NativeDouble) ((n - 2 + i) * i);
            p1 += u2;
            r3 = 1.0 / (NativeDouble) i;
            r4 = 1.0 / (NativeDouble) (i + 1);
            u3 = u3 * r2 * (r3 * r3);
            v1 += r3;
            q0 += u3 * v1;
            u4 = u4 * r2 * r3 * r4;
            v2 += (r3 + r4) * .5;
            q1 += u4 * v2;
        }
        fun1[n-1] = p0;
        fun1[n-2] = p1;
        r3 = -sys015;
        fun2[0] = r3 * q0;
        fun2[1] = r3 * (1.0 - r2 * 2.0 * q1);
        if (n == 2) return 0;

        fun2[2] = fun2[1] + r2 * fun2[0] * 2;
        i1 = n;
        for (i = 3; i < i1; ++i) {
            fun2[i] = fun2[i - 1] + fun2[i - 2] * r2 / (NativeDouble) ((i - 1) * (i - 2));
        }
        na = n - 2;
        i1 = na;
        for (j = 1; j <= i1; ++j) {
            i = na - j;
            fun1[i] = fun1[i + 1] + fun1[i + 2] * r2 /
                    (NativeDouble) ((i+1)*(i+2));
        }
        return 0;
    }

    if (x >= 10.0) {
        r1 = .125 / x;
        r2 = -(r1 * r1);
        p0 = 1.0;
        q0 = 1.0;
        p1 = 1.0;
        q1 = 1.0;
        u1 = 1.0;
        u2 = 1.0;
        u3 = 1.0;
        u4 = 1.0;
        for (i = 1; i <= 9; ++i) {
            v1 = (NativeDouble) (i * 2);
            v2 = v1 * 2.0;
            v3 = SQR((v2 - 3.0) * (v2 - 1.0));
            v4 = v2 * v2 - 1.0;
            v5 = v4 * (v2 - 3.0);
            v6 = v5 * (v2 + 3.0);
            v5 *= v2 - 5.0;
            v4 = v4 * v4;
            v7 = r2 / ((v1 - 1.0) * v1);
            v2 = r2 / ((v1 + 1.0) * v1);
            u1 *= v3 * v7;
            p0 += u1;
            u2 *= v4 * v2;
            q0 += u2;
            u3 *= v5 * v7;
            p1 += u3;
            u4 *= v6 * v2;
            q1 += u4;
        }
        q0 = r1 * q0;
        q1 = r1 * 3.0 * q1;
        u1 = sys005;
        u2 = sqrt((NativeDouble)(u1 * x));
        u3 = x + u1 * .5;
        u4 = sin(u3) / u2;
        u5 = cos(u3) / u2;
        fun1[0] = p0 * u4 - q0 * u5;
        fun1[1] = q1 * u4 - p1 * u5;
        fun2[0] = -(p0 * u5 + q0 * u4);
        fun2[1] = -(p1 * u4 + q1 * u5);
    }
    else {
        r1 = x * .5;
        r2 = -(r1 * r1);
        u1 = 1.0;
        u2 = 1.0;
        p0 = 1.0;
        p1 = 1.0;
        q0 = -((NativeDouble)log(r1) + sys039);
        v1 = q0;
        q1 = v1 + .5;
        v2 = q1;
        for (i = 1; i <= 20; ++i) {
            r3 = 1.0 / (NativeDouble) i;
            r4 = 1.0 / (NativeDouble) (i + 1);
            u1 = u1 * r2 * (r3 * r3);
            p0 += u1;
            v1 += r3;
            q0 += u1 * v1;
            u2 = u2 * r2 * r3 * r4;
            p1 += u2;
            v2 += (r3 + r4) * .5;
            q1 += u2 * v2;
        }
        fun1[0] = p0;
        fun1[1] = r1 * p1;
        r2 = -sys015;
        fun2[0] = r2 * q0;
        fun2[1] = (.5 / r1 + r1 * q1) * r2;
    }
    r1 = 2.0 / x;
    for (i = 2; i < n1; ++i) {
        r2 = r1 * (NativeDouble) (i - 1);
        fun1[i] = r2 * fun1[i - 1] - fun1[i - 2];
        fun2[i] = r2 * fun2[i - 1] - fun2[i - 2];
    }
    if (n1 == n) return 0;
    r1 = x * .5;
    r2 = r1 * r1;
    int IMAX = n2-1;
    q2 = 1.0; u = 1.0; v = 1.0; r3 = r1; // to avoid compiler warnings
    int ok = 0;
    while(!ok) { // делаем несколько раз, поскольку нужно знать правильный IMAX в знаменателе
        q2 = 1.0;
        u = 1.0;
        v = 1.0;
        r3 = r1;
        int imaxnew = IMAX;
        ok = 1;
        for (i = 1; i < IMAX; ++i) {
            v += 1.0 / NativeDouble(i + 1);
            u *= r2 / NativeDouble(i * (IMAX - i));
            q2 += u;
            r3 *= r1 / NativeDouble(i);
            // Внимание - отсечка! Чтобы r3 не выходило за разрядную сетку типа NativeDouble, выходим из цикла
            NativeDouble val = fabs(q2/r3);
            // отмечаем, что в след. попытку нужно будет ограничиться этой величиной, но пока ещё не начинаем заново
            if(val>1e200 && imaxnew==IMAX) imaxnew = i;
            // если уж совсем плохо, то начинаем заново
            if(val>1e250) { IMAX = imaxnew; ok = 0; break; }
        }
    }
    u = 1.0;
    p = 1.0;
    r2 = -r2;
    q1 = v * 0.5 - (NativeDouble)log(r1) - sys039;
    if (x - 10.0 >= 0.0) k = nc2 - 1;
    else k = 20;
    for (i = 1; i <= k; ++i) {
        r4 = 1.0 / (NativeDouble) i;
        r5 = 1.0 / (NativeDouble) (IMAX + i);
        u = u * r2 * r4 * r5;
        p += u;
        v += (r4 + r5) * .5;
        q1 += u * v;
    }
    q2 *= 0.5 / r3;
    r3 /= NativeDouble(IMAX);
    p *= r3;
    q1 *= r3;
    r2 = -sys015;

    _r = 2.0 / x;
    work[n1-1] = 0.0;
    for (i = n1; i < n2-1; ++i) {
        work[i] = 1.0 / (_r * NativeDouble(i) - work[i-1]);
        fun1[i] = work[i] * fun1[i-1];
        fun2[i] = work[i] * fun2[i-1];
    }

    fun1[n2-1] = (IMAX==n2-1 ? p : 0.0);
    for(i = n2-2; i >= n1; i--)
        fun1[i] += work[i] * fun1[i+1];

    for(i = n2-1; i > IMAX; i--)
        fun2[i] = -1e100;

    fun2[IMAX] = r2 * (q1 + q2);
    for(i = IMAX-1; i >= n1; i--)
        fun2[i] += work[i] * fun2[i+1];

    if(IMAX!=n2-1 && IMAX<n-1) return IMAX;
    else return 0;
} // sf33r_c
//======================================================================================================================


//======================================================================================================================
// Вычисление функции Бесселя 1-го и 2-го рода целых индексов
//======================================================================================================================
void BesselFunctions(NativeDouble x, int n, int nmax, NativeDouble *fun1, NativeDouble *fun2, NativeDouble *work) {
    if(n==1) {
        fun1[0] = BesselJ0(x);
        fun2[0] = BesselN0(x);
        return;
    }
    int i = sf33r(x, n, nmax, fun1, fun2, work);
    if(i<0) crash("BesselFunctions error: returned %i value\n", i);
    for(i=0; i<n; i++) if(fun2[i]<-1e100) fun2[i]=-1e100;
    if (x > 1.0) return;

    NativeDouble factorial_im1 = 1.;
    NativeDouble factorial_i = 1.;
    for(i=0; i<n; i++) {
        if(i > 30) {
            fun1[i] = 0.0;
            fun2[i] = -huge;
        }
        else {
            if(i) factorial_i *= i;
            fun1[i] *= pow(0.5*x, i) / factorial_i;
            if(i) fun2[i] *= pow(0.5*x, -i) * 0.5 * factorial_im1;
            factorial_im1 = factorial_i;
        }
    }
}
//======================================================================================================================


//======================================================================================================================
// Вычисление радиальной части собственных функций оператора Лапласа во внешности круга радиуса rc
// с граничным условием 2-го рода на границе круга
// Функции определяются формулой u_{nu,k}(r) = phi(nu, k*r, k*r_c), r >= rc,
// phi(nu,x,y) = (- N'_nu(y) J_nu(x) + J'_nu(y) N_nu(x)) / sqrt((J'_nu(y))^2 + (N'_nu(y))^2), x >= y.
// Процедура вычисляет значения функций phi(nu,x,y) и dphi/dx(nu,x,y) при заданных x,y и всех nu=0,...,Nmax-1
//======================================================================================================================
void BesselOuterFunction(NativeDouble x, NativeDouble y, int Nmax, NativeDouble* phi, NativeDouble* dphi) {
    if(x<y) crash("BesselOuterFunction error: x < y");
    if(phi==NULL && dphi==NULL) crash("BesselOuterFunction error: output arrays are not allocated");

    int nmax = MAX(Nmax, 2);
    const int bufsize = MAX(MAX(nmax+1, 2*MAX(int(x),int(y))+4),8);

    if(x < tiny) { // внутреннего радиуса нет и вычисляется значение в нуле
        for(int i=0; i<Nmax; i++) phi[i] = dphi[i] = 0.0;
        phi[0] = 1.0;
        if(Nmax>=1) dphi[1] = 1.0;
        return;
    }

    // Выделяем массив один куском и выставляем указатели
    NativeDouble* buf = new NativeDouble[bufsize*7];
    NativeDouble* Jx = buf;
    NativeDouble* Nx = buf+bufsize;
    NativeDouble* Jy = buf+bufsize*3;
    NativeDouble* Ny = buf+bufsize*4;

    // Вычисляем стандартными программами все функции Бесселя 1-го рода, а функции 2-го рода -- те, которые удаётся посчитать
    // Если flag*=1, то вместо функций Бесселя 1-го и 2-го рода вычисляются масштабированные величины (см. описание sf33r)
    int flagx = (x <= 1.0);
    int flagy = (y <= 1.0);
    int numfuncx = sf33r(x, nmax+1, bufsize, buf, buf+bufsize, buf+bufsize*2); // значения при аргументе x

    if(y < tiny) { // внутреннего радиуса нет - вычисляем функции Бесселя 1-го рода
        if(numfuncx==0) numfuncx = nmax+1;
        if(flagx) {
            NativeDouble multJ = 1.0; // множитель для функции 1-го рода
            for(int i=1; i<numfuncx; i++) {
                multJ *= (0.5*x) / NativeDouble(i);
                Jx[i] *= multJ;
            }
            for(int i=numfuncx; i<Nmax; i++) Jx[i] = 0.0;
        }
        for(int i=0; i<Nmax; i++) {
            if(phi)  phi[i] = Jx[i];
            if(dphi) dphi[i] = -Jx[i+1] + i*Jx[i]/x;
        }
        delete[] buf;
        return;
    }

    int numfuncy = sf33r(y, nmax+1, bufsize, buf+bufsize*3, buf+bufsize*4, buf+bufsize*5); // значения при аргументе y
    // Возвращённые значения numfunc* -- число фактически вычисленных значений
    // Внутренняя критическая ошибка
    if(numfuncx<0 || numfuncy<0) crash("BesselOuterFunction error: returned %i, %i values", numfuncx, numfuncy);
    if(numfuncx==0) numfuncx = nmax+1;
    if(numfuncy==0) numfuncy = nmax+1;

    // Проводим денормализацию функций Бесселя. Были значения:
    //    in (x) = n! (x / 2)^(-n) Jn (x) ,
    //    k0 (x) = N0 (x) ,
    //    kn (x) = 2 (x / 2)^n Nn (x) / (n - 1)!
    if(flagx) {
        NativeDouble multJ = 1.0; // множитель для функции 1-го рода
        NativeDouble multN = 0.5; // множитель для функции 2-го рода
        NativeDouble C2_x = 1.0/(0.5*x);
        int numfuncy_new = numfuncx;
        for(int i=1; i<numfuncx; i++) {
            multJ *= (0.5*x) / NativeDouble(i);
            multN *= C2_x;
            if(i!=1) multN *= (i-1);

            Jx[i] *= multJ;
            Nx[i] *= multN;
            // Избегаем floating point overflow
            if(Nx[i] <= -1e250) { multN = 0.0; numfuncy_new = i+1; }
        }
        numfuncx = numfuncy_new;
    }
    if(flagy) {
        NativeDouble multJ = 1.0; // множитель для функции 1-го рода
        NativeDouble multN = 0.5; // множитель для функции 2-го рода
        NativeDouble C2_y = 1.0/(0.5*y);
        int numfuncy_new = numfuncy;
        for(int i=1; i<numfuncy; i++) {
            multJ *= (0.5*y) / NativeDouble(i);
            multN *= C2_y;
            if(i!=1) multN *= (i-1);

            Jy[i] *= multJ;
            Ny[i] *= multN;
            // Избегаем floating point overflow
            if(Ny[i] <= -1e250) { multN = 0.0; numfuncy_new = i+1; }
        }
        numfuncy = numfuncy_new;
    }

    // For large arguments and small indices, forward recursion is preferable then the values given by sf33r.
    // Recalculating Bessel functions of the second kind
    int flg1 = 10*x > Nmax, flg2 = 10*y > Nmax;
    if(numfuncx==nmax+1 && numfuncy==nmax+1 && (flg1 || flg2)) {
        Nx[0] = BesselN0(x);
        Nx[1] = BesselN1(x);
        Ny[0] = BesselN0(y);
        Ny[1] = BesselN1(y);
        for(int i=1; i<Nmax; i++) {
            if(flg1) Nx[i+1] = 2*i*Nx[i]/x - Nx[i-1];
            if(flg2) Ny[i+1] = 2*i*Ny[i]/y - Ny[i-1];
        }

        /*
        Jx[0] = BesselJ0(x);
        Jx[1] = BesselJ1(x);
        Jy[0] = BesselJ0(y);
        Jy[1] = BesselJ1(y);
        Jx[2] = 2./x*Jx[1] - Jx[0];
        Jy[2] = 2./y*Jy[1] - Jy[0];
        Jx[3] = 4./x*Jx[2] - Jx[1];
        Jy[3] = 4./y*Jy[2] - Jy[1];
        */
    }

    // так как y<=x, то количество функций, влезающих в арифметику, для y не больше, чем для x. На всякий случай эта проверка
    if(numfuncy > numfuncx) numfuncy = numfuncx;
    // и ещё отсечка
    int numfunc = MIN(MIN(numfuncx, numfuncy), Nmax+1);

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
    if(numfunc==Nmax+1) { delete[] buf; return; } // всё уже вычислено

    // Выставляем указатели для вспомогательных переменных
    NativeDouble* Ax = buf+bufsize*2;
    NativeDouble* Ay = buf+bufsize*5;
    NativeDouble* B  = buf+bufsize*6;
    
    // Вычисляем значение Ay по функциям Бесселя 2-го рода максимально доступного индекса
    Ay[numfuncy-2] = Ny[numfuncy-2] / Ny[numfuncy-3];
    // Рекуррентно вычисляем бОльшие значения
    for(int i=numfuncy-2; i<Nmax; i++)
        Ay[i+1] = 2*i/y - 1.0/Ay[i];

    // Вычисляем значение Ax по функциям Бесселя 2-го рода максимально доступного индекса
    // Ещё нужно что-то сделать, если при numfuncy <= i < numfuncx
    for(int i=numfuncy-2; i<=numfuncx-2; i++) {
        if(fabs(Nx[i-1]) < 1e-10) Ax[i] = 1.23456789e100; // это значение не понадобится
        else Ax[i] = Nx[i] / Nx[i-1];
    }

    // Рекуррентно вычисляем бОльшие значения
    for(int i=numfuncx-2; i<Nmax; i++)
        Ax[i+1] = 2*i/x - 1.0/Ax[i];

    // Вычисляем значение B по функциям Бесселя 2-го рода максимально доступного индекса
    B[numfuncy-2] = Nx[numfuncy-2] / Ny[numfuncy-2];  // здесь Nx и Ny уже определены, причём знаменатель уже далёк от нуля
    for(int i=numfuncy-1; i<Nmax; i++) {
        if(fabs(Nx[i]) < 1e100) B[i] = 0.0; // знаменатель уже 1e250. Если числитель маленький, то частное равно 0
        else B[i] = Ax[i]*B[i-1]/Ay[i]; // рекуррентно вычисляем бОльшие значения
    }
        
    // Вычисляем искомую функцию для недостающих значений аргумента
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
// Вычисление выражений H_n(y) / H'_n(x), 
// где H_n - функция Ганкеля индекса n, H'_n - её производная, y >= x > 0 - некоторые числа.
// Значения вычисляются для всех индексов n от 0 до N-1
// Вход: x, y, N. Выход: Re[0..N-1], Im[0..N-1]
//======================================================================================================================
void HankelFractions(NativeDouble x, NativeDouble y, int N, NativeDouble *Re, NativeDouble* Im) {
    if(y<x) crash("HankelFractions error: y < x");
    if(Im!=NULL) crash("HankelFractions error: Imaginary part calculation not realized");

    int n = MAX(N, 2);
    const int bufsize = MAX(n+1, 2*MAX(int(x),int(y))+4);
    NativeDouble* buf = GimmeMem<NativeDouble>(bufsize*6);
    int flagx = (x <= 1.0);
    int flagy = (y <= 1.0);
    int err1 = sf33r(x, n+1, bufsize, buf, buf+bufsize, buf+bufsize*2);
    int err2 = sf33r(y, n, bufsize, buf+bufsize*3, buf+bufsize*4, buf+bufsize*5);
    if(err1<0 || err2<0) crash("BesselFunctions error: returned %i, %i values", err1, err2);

    // Если не смогли вычислить все функции, обрезаем хвост
    if(err1) n = MIN(n, err1-1);
    if(err2) n = MIN(n, err2-1);

    // dH_0/dx = -H_1; dH_n/dx = -H_{n+1} + H_n*n/x
    NativeDouble* Jx = buf;
    NativeDouble* Nx = buf+bufsize;
    NativeDouble* Jy = buf+bufsize*3;
    NativeDouble* Ny = buf+bufsize*4;

    // Значение функций Бесселя 1-го рода обычно много меньше, чем у функций 2-го рода
    // Поэтому меняем их нормировку на принятую у функций 2-го рода
    // Нужно умножить на (x/2)^k / k! для денормализации
    // и затем ещё на 2(x/2)^k (k-1)! для задания новой нормировки (искл. случай k=0)
    // Итого при k=0 ничего не делаем; при k=1,... домножаем на 2/k (x/2)^{2k}
    if(flagx) {
        NativeDouble m = 0.25*x*x;
        NativeDouble mm = 1.0;
        for(int k=1; k<=n; k++) {
            mm *= m;
            Jx[k] *= mm * 2. / k;
        }
    }
    if(flagy) {
        NativeDouble m = 0.25*y*y;
        NativeDouble mm = 1.0;
        for(int k=1; k<n; k++) {
            mm *= m;
            Jy[k] *= mm * 2. / k;
        }
    }

    // Заменяем функцию Ганкеля в знаменателе на её производную (что требуется по условию)
    const NativeDouble invx = 1.0/x;
    if(flagx) { // а здесь ещё и учитываем нормировку
        Jx[0] = -Jx[1]*invx;
        Nx[0] = -Nx[1]*invx;
        for(int k=1; k<n; k++) {
            Jx[k] = (-2.*Jx[k+1] + Jx[k])*NativeDouble(k)*invx;
            Nx[k] = (-2.*Nx[k+1] + Nx[k])*NativeDouble(k)*invx;
        }
    }
    else {
        Jx[0] = -Jx[1];
        Nx[0] = -Nx[1];
        for(int k=1; k<n; k++) {
            Jx[k] = -Jx[k+1] + NativeDouble(k)*Jx[k]*invx;
            Nx[k] = -Nx[k+1] + NativeDouble(k)*Nx[k]*invx;
        }
    }

    // Делим числитель (представленный комплексным числом) на знаменатель (тоже комплексный)
    for(int k=0; k<MIN(N,n); k++) {
        Re[k] = (Jy[k]*Jx[k] + Ny[k]*Nx[k]) / (Jx[k]*Jx[k] + Nx[k]*Nx[k]);
    }
    for(int k=MIN(N,n); k<N; k++) {
        Re[k] = 0.0; // не смогли вычислить
    }
    N = MIN(N,n);

    // Если 1<x<=y, то нормировки нет
    // Если x<=y<=1, то у числителя и знаменателя множители отличаются на (x/y)^n
    if(flagx && flagy) {
        NativeDouble x_y = x / y;
        NativeDouble mult = 1.0;
        for(int k=1; k<N; k++) {
            mult *= x_y;
            Re[k] *= mult;
        }
    }
    // В случае x<=1<y нужно учесть нормировку, т. е. домножить на 2 (x / 2)^n / (n - 1)!
    if(flagx && !flagy) {
        NativeDouble mult = x;
        for(int k=1; k<N; k++) {
            Re[k] *= mult;
            mult *= 0.5*x / NativeDouble(k);
        }
    }
    FreeMem(buf);
}
//======================================================================================================================

//======================================================================================================================
// Вычисление функции Бесселя 1-го рода индекса l в точек x
// Также вычисляет её производные до 5-го порядка, если соответствующие указатели ненулевые
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

    int Nmax = int(x) * 2 + 4;
    if(Nmax < l+3) Nmax = l+3;
    if(Nmax < 8) Nmax = 8;
    NativeDouble* BesselJ = GimmeMem<NativeDouble>(Nmax*3, "BesselJ");
    NativeDouble* BesselY = BesselJ + Nmax;
    NativeDouble* Work    = BesselY + Nmax;
    BesselFunctions(x, l+2, Nmax, BesselJ, BesselY, Work);
    NativeDouble bj = BesselJ[l];
    if(dbdx) *dbdx = -BesselJ[l+1] + l*BesselJ[l]/x;
    if(d2bdx2) *d2bdx2 = BesselJ[l+1]/x + BesselJ[l]*(-1.0 + l*(l-1)/(x*x));
    if(d3bdx3) *d3bdx3 = BesselJ[l+1]*(1.0-(2+l*l)/(x*x)) + BesselJ[l]*((1-l)/x + (l*l*l-3*l*l+2*l)/(x*x*x));
    if(d4bdx4) *d4bdx4 = BesselJ[l+1]*(-2+6*(1+l*l)/(x*x))/x + BesselJ[l]*(1 + (-2*l*l-3+2*l)/(x*x) + (-6*l+11*l*l-6*l*l*l+l*l*l*l)/(x*x*x*x));
    if(d5bdx5) *d5bdx5 =-BesselJ[l+1]*(1.0+(-7-2*l*l)/(x*x)+(35*l*l+24+l*l*l*l)/(x*x*x*x)) 
        - BesselJ[l]*((-l+2)/x + (-12*l*l+7*l+2*l*l*l-12)/(x*x*x) + (50*l*l-l*l*l*l*l-24*l-35*l*l*l+10*l*l*l*l)/(x*x*x*x*x));
    FreeMem(BesselJ);
    return bj;
}
//======================================================================================================================


//======================================================================================================================
// Printing for test purpose
//======================================================================================================================
void PrintBesselFunctions() {
    const int Nmax = 10000;
    NativeDouble *fun1 = GimmeMem<NativeDouble>(Nmax, "BesselFunctions");
    NativeDouble *fun2 = GimmeMem<NativeDouble>(Nmax, "BesselFunctions");
    NativeDouble *work = GimmeMem<NativeDouble>(Nmax, "BesselFunctions");

    FILE* J0 = fopen("Bessels.dat", "wt");
    for(int i = 1; i<100000; i++) {
        NativeDouble x = NativeDouble(i) / 100;
        BesselFunctions(x, 50, Nmax, fun1, fun2, work);
        fprintf(J0, "%f %e %e %e %e %e %e %e %e %e %e\n", x,
            fun1[33], fun1[34], fun1[35], fun1[36], fun1[37],
            fun2[33], fun2[34], fun2[35], fun2[36], fun2[37]);
    }
    fclose(J0);
    FreeMem(fun1);
    FreeMem(fun2);
    FreeMem(work);
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
// Вычисление нуля (номер RadialMode) производной функции Бесселя 1-го рода индекса AngularMode
// Проверялось, что нули совпадают с результатами, полученныим в Maple.
// Код для Maple, выводящий таблицу нулей производной функции Бесселя:
// f := proc (m, n) options operator, arrow; fsolve(BesselJ(m+1, x)*x-m*BesselJ(m, x), x = BesselJZeros(m, n) .. BesselJZeros(m, n+1)) end proc
// for i to 20 do f(1, i) end do
//======================================================================================================================
NativeDouble BesselPrimeZero(int AngularMode, int RadialMode, int log) {
    if(AngularMode > 8) crash("BesselPrimeZero: AngularMode > 8 not realized");

    // Ищем какое-то осмысленное начальное приближение
    NativeDouble Nu = PiNumber * (RadialMode + 0.25 + 0.5*AngularMode); // Асимптотика при больших значениях корня

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

    // Уточняем корень методом Ньютона
    if(log) pprintf("Bessel prime zero: %23.16e; ", Nu);

    unsigned int maxit = 20;
    for(unsigned int i=0; i<maxit; i++) {    
        NativeDouble dfdx = 0.0, d2fdx2 = 0.0;
        BesselJ(AngularMode, Nu, &dfdx, &d2fdx2, NULL, NULL, NULL); // у этой функции двойная точность...
        double dNu = dfdx / d2fdx2;
        Nu -= NativeDouble(dNu); // Уточняем нуль производной функции Бесселя
        if(fabs(dNu) < 1e-14) maxit = MIN(maxit, i+1+(sizeof(NativeDouble)>>2)); // делаем ещё не более 2 итераций (для dd - 4, для qd -- 8)
    }
   
    if(log) pprintf("precised: %23.16e\n", Nu);

    return Nu;
}
//======================================================================================================================


//======================================================================================================================
NativeDouble BesselK0(NativeDouble x, int mode, int *ierr){ // sf10d_c
//======================================================================================================================
    const NativeDouble sys059 = 1.7e308;
    const NativeDouble pp[9] = { 103.0259990294628,2781.556951284448,
            19273.80098600572,52759.23966621326,67868.97902962336,
            44234.46358185698,14773.24772717615,2371.138845157722,
            142.7782534685513 };
    const NativeDouble qq[9] = { 172.3547893760824,3220.845810178108,
            18937.36743618128,47461.83671148535,58017.2702086999,
            36670.51460217709,12017.60312809182,1906.135146578901,
            113.9205640609552 };
    const NativeDouble f[4] = { -1.641445283729907,-296.0165789295884, -17733.78468495299,-403203.4076114548 };
    const NativeDouble g[3] = { -250.6497244587799,29865.71316305403, -1612813.630445819 };
    const NativeDouble sys073 = 2.2227587494851e-162;
    const NativeDouble sys055 = 2.225073858507202e-308;
    const NativeDouble const__ = -.11593151565841;
    const NativeDouble xmax = 177.85287909124;
    const NativeDouble p[5] = { -3.272279992574785,-513.5620533725103, -25470.74686782375,-356272.9668888909,-149269.5386816498 };
    const NativeDouble q[3] = { -228.6835794674987,25325.71801733452, -1287566.524373463 };

    if(ierr) *ierr = 0;
    if (x <= 0.0) { if(ierr) *ierr = 66; return sys059; }
    NativeDouble a = 1.0;
    switch (mode) {
    case 1: 
        if (x > xmax) { if(ierr) *ierr = 67; return 0.0; }
        if (x > 1.0) a = exp(-x);
        break;
    case 2: 
        if (x <= 1.0) a = exp(x);
        break;
    default:
        if(ierr) *ierr = 65;
        return sys059;
    }
    if (x <= 1.0) {
        NativeDouble temp = log(x);
        if (x < sys073) return const__ - temp;
        NativeDouble xx = x * x;
        NativeDouble sump = (((p[0] * xx + p[1]) * xx + p[2]) * xx + p[3]) * xx + p[4];
        NativeDouble sumq = ((xx + q[0]) * xx + q[1]) * xx + q[2];
        NativeDouble sumf = ((f[0] * xx + f[1]) * xx + f[2]) * xx + f[3];
        NativeDouble sumg = ((xx + g[0]) * xx + g[1]) * xx + g[2];
        NativeDouble x1 = (sump / sumq - xx * sumf * temp / sumg - temp) * a;
        if (x1 >= 1.0) return x1 += sys055;
        else return x1;
    }
    else {
        NativeDouble xx = 1.0 / x;
        NativeDouble sump = pp[0];
        for (int i = 2; i <= 9; ++i)
            sump = sump * xx + pp[i - 1];
        NativeDouble sumq = xx;
        for (int i = 1; i <= 8; ++i)
            sumq = (sumq + qq[i - 1]) * xx;
        sumq += qq[8];
        return sump / sumq / sqrt(x) * a;
    }
}
//======================================================================================================================


//======================================================================================================================
NativeDouble BesselK1(NativeDouble x, int mode, int *ierr) { // sf11d_c
//======================================================================================================================
    const NativeDouble sys059 = 1.7e308;
    const NativeDouble xmax = 177.85567464849;
    const NativeDouble sys073 = 2.2227587494851e-162;
    const NativeDouble sys061 = 2.23e-308;
    const NativeDouble p[5] = { .1203176761421961,24.99784339185733,
            1797.13456510212,44333.31008786754,179847.3001635515 };
    const NativeDouble q[4] = { .25,-70.35978938634682,9316.074668016924,
            -553734.3719560826 };
    const NativeDouble f[5] = { .00210761874039124,.4076347371250282,
            30.92138338821924,951.7918567206316,8501.049475393063 };
    const NativeDouble g[2] = { -221.6786554069839,17002.09895078612 };
    const NativeDouble pp[10] = { .06982646013142393,6.943429490437559,
            101.624513525677,519.0226276746315,1207.119908800295,
            1430.398225341416,900.5826147591565,298.6452400790679,
            48.33007694228657,2.958171781643913 };
    const NativeDouble qq[8] = { 30.17103765395105,228.1150283978967,
            673.0124200199513,924.6057876003073,638.5720888284594,
            224.4322524683622,37.67671736737388,2.360279592776385 };

    if(ierr) *ierr = 0;
    if (x < sys061) { if(ierr) *ierr = 66; return sys059; }

    NativeDouble a = 1.0;
    if(mode==1) {
        if (x > xmax) { if(ierr) *ierr = 67; return 0.0; }
        if (x > 1.0) a = exp(-(x));
    }
    else if(mode==2) {
        if (x < 1.0) a = exp(x);
    }
    else { if(ierr)  *ierr = 65; return sys059; }

    if(x <= 1.0) {
        if (x < sys073) return 1.0 / x;
        NativeDouble xx = x * x;
        NativeDouble sump = ((((p[0] * xx + p[1]) * xx + p[2]) * xx + p[3]) * xx + p[4]) * xx + q[3];
        NativeDouble sumq = ((q[0] * xx + q[1]) * xx + q[2]) * xx + q[3];
        NativeDouble sumf = (((f[0] * xx + f[1]) * xx + f[2]) * xx + f[3]) * xx + f[4];
        NativeDouble sumg = (xx + g[0]) * xx + g[1];
        return (xx * log(x) * sumf / sumg + sump / sumq) / x * a;
    }    
    else {
        NativeDouble xx = 1.0 / x;
        NativeDouble sump = pp[0];
        for (int i__ = 2; i__ <= 10; ++i__) {
            sump = sump * xx + pp[i__ - 1];
        }
        NativeDouble sumq = xx;
        for (int i__ = 1; i__ <= 7; ++i__) {
            sumq = (sumq + qq[i__ - 1]) * xx;
        }
        sumq += qq[7];
        return a * sump / sumq / sqrt(x);
    }
}
//======================================================================================================================

//======================================================================================================================
NativeDouble BesselK2(NativeDouble x, int mode, int *ierr) { // sf11d_c
//======================================================================================================================
    int err0, err1;
    NativeDouble K0 = BesselK0(x, mode, &err0);
    NativeDouble K1 = BesselK1(x, mode, &err1);
    if(ierr) *ierr = err0 + err1;
    return 2.0/x*K1 + K0;
}
//======================================================================================================================


//======================================================================================================================
// This routine calculates Bessel functions I and J of complex argument and integer order 
// Based on a FORTRAN subroutine from PORT3 library
// X      Real part of the complex argument
//        If I*S are to be calculated, abs(x) must not exceed EXPARG (see below). 
// Y      Imaginary part of the argument. If J*S are to be calculated, abs(y) must not excced EXPARG. 
// NB     1 + highest order to be calculated. Must be positive.
// IZE    0 if J*S are to be calculated, 1 if I*S are to be calculated. 
// BR[NB] If the routine terminates normally (NCALC==NB), it returns the real part of J(OR I)-SUB-ZERO 
//        through J(OR I)-SUB-NB-MINUS-ONE of Z in this vector. 
// BI[NB] Imaginary analog of BR. 
// NCALC  Before using the results, the user should check that NCALC=NB, 
//        i.e. all orders have been calculatde to the desired accuracy.  See ERROR RETURNS below.

//        EXPLANATION OF MACHINE-DEPENDENT CONSTANTS
// NSIG   DECIMAL SIGNIFICANCE DESIRED.  SHOULD BE SET TO IFIX(ALOG10(2)*NBIT+1), WHERE NBIT IS THE NUMBER OF
//        BITS IN THE MANTISSA OF A DOUBLE PRECISION VARIABLE. SETTING NSIG HIGHER WILL INCREASE CPU TIME WITHOUT 
//        INCREASING ACCURACY, WHILE SETTING NSIG LOWER WILL DECREASE ACCURACY.  IF ONLY SINGLE-PRECISION 
//        ACCURACY IS DESIRED, REPLACE NBIT BY THE NUMBER OF BITS IN THE MANTISSA OF A SINGLE-PRECISION VARIABLE. 
//        THE RELATIVE TRUNCATION ERROR IS LIMITED TO T=.5*10 **-NSIG FOR ORDER GREATER THAN ABS(Z), AND FOR ORDER 
//        LESS THAN ABS(Z) (GENERAL TEST), THE RELATIVE ERROR IS LIMITED TO T FOR FUNCTION VALUES OF MAGNITUDE AT 
//        LEAST 1, AND THE ABSOLUTE ERROR IS LIMITED TO T FOR SMALLER VALUES. 
// NTEN   LARGEST INTEGER K SUCH THAT 10**K IS MACHINE- REPRESENTABLE IN DOUBLE PRECISION. 
// LARGEZ UPPER LIMIT ON THE MAGNITUDE OF Z.  BEAR IN MIND THAT IF ABS(Z)=N, THEN AT LEAST N ITERATIONS OF THE 
//        BACKWARD RECURSION WILL BE EXECUTED. 
// EXPARG LARGEST DOUBLE PRECISION ARGUMENT THAT THE LIBRARY DEXP ROUTINE CAN HANDLE.
//
// IF N.GT.NCALC AND ABS(G(NCALC-1)/G(N-1)).EQ.10**-K, THEN THE LAST K SIGNIFICANT FIGURES OF G(N-1) (=BR(N)+I*BI(N)) ARE
// ERRONEOUS. IF THE USER WISHES TO CALCULATE G(N-1) TO HIGHER ACCURACY, HE SHOULD USE AN ASYMPTOTIC FORMULA FOR LARGE ORDER.
//======================================================================================================================
template<class real>
int db1slc_(real x, real y, int nb, int ize, real *br, real *bi, int *ncalc) {
    int k, l, m, n;
    real pi, pr;
    int npl, imag, nend, magz;
    real sign;
    int nbmz;
    real sumi, test, sumr, poldi, poldr;
    int mlast;
    real tover, zinvi, zinvr, tempai, tempbi, tempci;
    int imrecr;
    real psavei, tempar, tempbr, tempcr, plasti;
    int mrecur;
    real psaver, plastr;
    int nstart;

    const int nsig = 1 + 2 * sizeof(real);
    const int nten = 306;

    const real C10   = 10.0;
    const real C1_10 = 0.1;

    tempar = sqrt(x*x + y*y);
    magz = (int) tempar;
    if(tempar > real(10000) || fabs(y) > 0.5 * nten / log10(2.0)) {
        for(int i=0; i<nb; i++) br[i]=bi[i]=0.0;
        return 1; //crash("db1slc_: values %e %e are outside the limits", x, y);
    }

    --bi;
    --br;

    sign = 1 - (ize << 1);
    *ncalc = nb;
// USE 2-TERM ASCENDING SERIES FOR SMALL Z
    {
        real d__1 = tempar*tempar;
        if (d__1 * d__1 < pow(C1_10, nsig)) { // TWO-TERM ASCENDING SERIES FOR SMALL Z
            tempar = 1.;
            tempai = 0.;
            tempcr = (x*x - y*y) * .25;
            tempci = x * .5 * y;
            br[1] = 1. - sign * tempcr;
            bi[1] = -sign * tempci;
            if (nb == 1) return 0;
            for (n = 2; n <= nb; ++n) {
                tempbr = (tempar * x - tempai * y) / (real) ((real) ((n*2) - 2));
                tempai = (tempar * y + tempai * x) / (real) ((real) ((n*2) - 2));
                tempar = tempbr;
                tempbr = n;
                br[n] = tempar * (1. - sign * tempcr / tempbr) + tempai * tempci / tempbr;
                bi[n] = tempai * (1. - sign * tempcr / tempbr) - tempar * tempci / tempbr;
            }
            return 0;
        }
    }

// INITIALIZE THE CALCULATION OF THE P*S
    nbmz = nb - magz;
    n = magz + 1;
    if (fabs(x) < fabs(y)) {
        zinvi = -1. / (y + x*x / y);
        zinvr = -(x) * zinvi / y;
    }
    else {
        zinvr = 1. / (x + y*y / x);
        zinvi = -(y) * zinvr / x;
    }
    plastr = 1.;
    plasti = 0.;
    pr = sign * (real) ((real) (n*2)) * zinvr;
    pi = sign * (real) ((real) (n*2)) * zinvi;
    test = pow(C10, nsig) * 2.;
    m = 0;
    if (nbmz >= 3) {
    // CALCULATE P*S UNTIL N=NB-1.  CHECK FOR POSSIBLE OVERFLOW
        tover = pow(C10, nten - nsig);
        nstart = magz + 2;
        nend = nb - 1;

        for (n = nstart; n <= nend; ++n) {
            poldr = plastr;
            poldi = plasti;
            plastr = pr;
            plasti = pi;
            pr = sign * ((real) ((real) (n*2)) * (plastr * zinvr - plasti * zinvi) - poldr);
            pi = sign * ((real) ((real) (n*2)) * (plasti * zinvr + plastr * zinvi) - poldi);
            real d__1 = pr / tover;
            real d__2 = pi / tover;
            if (d__1 * d__1 + d__2 * d__2 - 1. <= 0.) continue;
            goto L7;
        }
        n = nend;
    // CALCULATE SPECIAL SIGNIFICANCE TEST FOR NBMZ.GT.2
        real d__1 = fabs(pr), d__2 = fabs(pi);
        tempbi = MAX(d__1,d__2);
        d__1 = pr / tempbi; // Computing 2nd power
        d__2 = pi / tempbi; // Computing 2nd power
        real d__3 = plastr / tempbi; // Computing 2nd power
        real d__4 = plasti / tempbi; // Computing 2nd power
        tempbi *= sqrt(pow(C10, nsig) * 2. * sqrt((d__1 * d__1 + d__2 * d__2) * (d__3 * d__3 + d__4 * d__4)));
        test = MAX(test,tempbi);
    }

// CALCULATE P*S UNTIL SIGNIFICANCE TEST IS PASSED
    while(1) {
        ++n;
        poldr = plastr;
        poldi = plasti;
        plastr = pr;
        plasti = pi;
        pr = sign * ((real) ((real) (n*2)) * (plastr * zinvr - plasti * zinvi) - poldr);
        pi = sign * ((real) ((real) (n*2)) * (plasti * zinvr + plastr * zinvi) - poldi);
        real d__1 = pr / test;
        real d__2 = pi / test;
        if (d__1 * d__1 + d__2 * d__2 < 1.) continue;
        if (m == 1) break;

// CALCULATE STRICT VARIANT OF SIGNIFICANCE TEST, AND CALCULATE P*S UNTIL THIS TEST IS PASSED
        m = 1;
        d__1 = fabs(pr), d__2 = fabs(pi);
        tempbi = MAX(d__1,d__2);
        d__1 = pr / tempbi;
        d__2 = pi / tempbi;
        real d__3 = plastr / tempbi;
        real d__4 = plasti / tempbi;
        tempbr = sqrt((d__1 * d__1 + d__2 * d__2) / (d__3 * d__3 + d__4 * d__4));
        tempbi = (real) ((real) (n + 1)) / tempar;
        if (tempbr + 1. / tempbr > tempbi * 2.) {
            d__1 = tempbi;
            tempbr = tempbi + sqrt(d__1 * d__1 - 1.); // Computing 2nd power
        }
        test /= sqrt(tempbr - 1. / tempbr);
        d__1 = pr / test; // Computing 2nd power
        d__2 = pi / test; // Computing 2nd power
        if (d__1 * d__1 + d__2 * d__2 - 1. >= 0.) break;
    }
    goto L12;

L7:
    nstart = n + 1;
// TO AVOID OVERFLOW, NORMALIZE P*S BY DIVIDING BY TOVER.
// CALCULATE P*S UNTIL UNNORMALIZED P WOULD OVERFLOW.
    pr /= tover;
    pi /= tover;
    plastr /= tover;
    plasti /= tover;
    psaver = pr;
    psavei = pi;
    tempcr = plastr;
    tempci = plasti;
    test = pow(C10, nsig << 1);
    
    do {
        ++n;
        poldr = plastr;
        poldi = plasti;
        plastr = pr;
        plasti = pi;
        pr = sign * ((real) ((real) (n*2)) * (plastr * zinvr - plasti * zinvi) - poldr);
        pi = sign * ((real) ((real) (n*2)) * (plasti * zinvr + plastr * zinvi) - poldi);
    } while (pr*pr + pi*pi <= test);

// CALCULATE BACKWARD TEST, AND FIND NCALC, THE HIGHEST N SUCH THAT THE TEST IS PASSED.
    tempbr = sqrt((plastr * plastr + plasti * plasti) / (poldr * poldr + poldi * poldi));
    tempbi = (real) ((real) n) / tempar;
    if (tempbr + 1. / tempbr > tempbi * 2.) {
        tempbr = tempbi + sqrt(tempbi * tempbi - 1.);
    }
    test = (1. - 1. / (tempbr * tempbr)) * .5 / pow(C10, nsig);
    test = (plastr * plastr + plasti * plasti) * test * ((poldr * poldr + poldi * poldi) * test);
    pr = plastr * tover;
    pi = plasti * tover;
    --n;
    nend = MIN(nb,n);
    for (*ncalc = nstart; *ncalc <= nend; ++(*ncalc)) {
        poldr = tempcr;
        poldi = tempci;
        tempcr = psaver;
        tempci = psavei;
        psaver = sign * ((real) ((real) (n*2)) * (tempcr * zinvr - tempci * zinvi) - poldr);
        psavei = sign * ((real) ((real) (n*2)) * (tempci * zinvr + tempcr * zinvi) - poldi);
        if ((psaver * psaver + psavei * psavei) * (tempcr * tempcr + tempci * tempci) - test > 0.) break;
    }
    --(*ncalc);

// THE COEFFICIENT OF B(N) IN THE NORMALIZATION SUM IS
// M*SQRT(-1)**IMAG, WHERE M=-2,0, OR 2, AND IMAG IS 0 OR 1.
// CALCULATE RECURSION RULES FOR M AND IMAG, AND INITIALIZE THEM.
L12:
    ++n;
    tempbr = ize ? x : y;
    const int ipos = tempbr > 0 ? 1 : (tempbr < 0 ? -1 : 0);
    mrecur = ((ize + 2 + ipos) / 2 << 2) - 3 - ((ize + ipos) << 1);
    k = ipos + 2 + (ize << 1) * (ipos * ipos) - ize;
    l = n - (n / 4 << 2);
    mlast = (k * l / 4 << 3) + 2 - (k * l / 2 << 2);
    if (ipos == 0 && (l == 1 || l == 3)) mlast = 0;
    l = l + 3 - ((l + 3) / 4 << 2);
    m = (k * l / 4 << 3) + 2 - (k * l / 2 << 2);
    if (ipos == 0 && (l == 1 || l == 3)) m = 0;
    imrecr = (1 - ize) * (ipos * ipos);
    imag = imrecr * (l - (l / 2 << 1));

// INITIALIZE THE BACKWARD RECURSION AND THE NORMALIZATION SUM
    tempbr = 0.;
    tempbi = 0.;
    if (fabs(pi) > fabs(pr)) {
        tempai = -1. / (pi + pr * (pr / pi));
        tempar = -(pr * tempai) / pi;
    }
    else {
        tempar = 1. / (pr + pi * (pi / pr));
        tempai = -(pi * tempar) / pr;
    }
    if (imag != 0) {
        sumr = -m * tempai;
        sumi = m * tempar;
    }
    else {
        sumr = m * tempar;
        sumi = m * tempai;
    }
    nend = n - nb;

    if (nend < 0) {
// N.LT.NB, SO STORE BR(N), BI(N), AND SET HIGHER ORDERS ZERO
        br[n] = tempar;
        bi[n] = tempai;
        nend = -nend;
        for (l = 1; l <= nend; ++l) {
            npl = n + l;
            br[npl] = 0.;
            bi[npl] = 0.;
        }
    }
    else {
    // RECUR BACKWARD VIA DIFFERENCE EQUATION CALCULATING (BUT NOT STORING) BR(N) AND BI(N) UNTIL N=NB.
        for (l = 1; l <= nend; ++l) {
            --n;
            tempcr = tempbr;
            tempci = tempbi;
            tempbr = tempar;
            tempbi = tempai;
            pr = (real) ((real) (n*2)) * zinvr;
            pi = (real) ((real) (n*2)) * zinvi;
            tempar = pr * tempbr - pi * tempbi - sign * tempcr;
            tempai = pr * tempbi + pi * tempbr - sign * tempci;
            imag = (1 - imag) * imrecr;
            k = mlast;
            mlast = m;
            m = k * mrecur;
            if (imag != 0) {
                sumr -= m * tempai;
                sumi += m * tempar;
            }
            else {
                sumr += m * tempar;
                sumi += m * tempai;
            }
        }
    // STORE BR(NB), BI(NB)
        br[n] = tempar;
        bi[n] = tempai;
        if (n <= 1) {
        // NB=1.  SINCE 2*TEMPAR AND 2*TEMPAI WERE ADDED TO SUMR AND SUMI RESPECTIVELY, WE MUST SUBTRACT TEMPAR AND TEMPAI
            sumr -= tempar;
            sumi -= tempai;
            goto L35;
        }

    // CALCULATE AND STORE BR(NB-1),BI(NB-1)
        --n;
        pr = (real) ((real) (n*2)) * zinvr;
        pi = (real) ((real) (n*2)) * zinvi;
        br[n] = pr * tempar - pi * tempai - sign * tempbr;
        bi[n] = pr * tempai + pi * tempar - sign * tempbi;
        if (n == 1) goto L34;

        imag = (1 - imag) * imrecr;
        k = mlast;
        mlast = m;
        m = k * mrecur;
        if (imag != 0) {
            sumr -= m * bi[n];
            sumi += m * br[n];
        }
        else {
            sumr += m * br[n];
            sumi += m * bi[n];
        }
    }

    nend = n - 2;
// CALCULATE VIA DIFFERENCE EQUATION AND STORE BR(N),BI(N), UNTIL N=2
    for (l = 1; l <= nend; ++l) {
        --n;
        pr = (real) ((real) (n*2)) * zinvr;
        pi = (real) ((real) (n*2)) * zinvi;
        br[n] = pr * br[n + 1] - pi * bi[n + 1] - sign * br[n + 2];
        bi[n] = pr * bi[n + 1] + pi * br[n + 1] - sign * bi[n + 2];
        imag = (1 - imag) * imrecr;
        k = mlast;
        mlast = m;
        m = k * mrecur;
        if (imag != 0) {
            sumr -= m * bi[n];
            sumi += m * br[n];
        }
        else {
            sumr += m * br[n];
            sumi += m * bi[n];
        }
    }
// CALCULATE AND STORE BR(1), BI(1)
    br[1] = (br[2] * zinvr - bi[2] * zinvi) * 2. - sign * br[3];
    bi[1] = (br[2] * zinvi + bi[2] * zinvr) * 2. - sign * bi[3];
L34:
    sumr += br[1];
    sumi += bi[1];
// CALCULATE NORMALIZATION FACTOR, TEMPAR +I*TEMPAI
L35:
    if (ize == 1) {
        tempcr = ipos * x;
        tempci = ipos * y;
    }
    else {
        tempcr = ipos * y;
        tempci = -ipos * x;
    }
    tempcr = exp(tempcr);
    tempbr = cos(tempci);
    tempbi = sin(tempci);
    if (fabs(sumr) < fabs(sumi)) {
        tempci = sumr / sumi;
        tempcr = tempcr / sumi / (tempci * tempci + 1.);
        tempar = tempcr * (tempbr * tempci + tempbi);
        tempai = tempcr * (tempbi * tempci - tempbr);
    }
    else {
        tempci = sumi / sumr;
        tempcr = tempcr / sumr / (tempci * tempci + 1.);
        tempar = tempcr * (tempbr + tempbi * tempci);
        tempai = tempcr * (tempbi - tempbr * tempci);
    }
// NORMALIZE
    for (n = 1; n <= nb; ++n) {
        tempbr = br[n] * tempar - bi[n] * tempai;
        bi[n] = br[n] * tempai + bi[n] * tempar;
        br[n] = tempbr;
    }
    return 0;
}
//======================================================================================================================


#include <complex>
using namespace std; 

//======================================================================================================================
// Вычисление функции Бесселя 1-го рода целого индекса l в точке x+iy
// Также вычисляет её 1-ю производную, если соответствующие указатели ненулевые
//======================================================================================================================
template<typename fpv>
void BesselJComplex(int l, fpv x, fpv y, 
             fpv* ReBesselJ, fpv* ImBesselJ, fpv* ReDBesselJ, fpv* ImDBesselJ) {
    // Заплатка для больших действительных (или близких к действительной оси) z
    // Оставляем один член разложения; из-за этого на действительной оси ошибка вычисления ~ 1e-7
    if(fabs(y) < 50 && fabs(x) > 9000.0) {
        double pi4 = 0.25*GetPiNumber<double>(); // pi/4
        complex<double> _I(double(0.0), double(1.0));
        complex<double> z = double(x) + _I*double(y);
        complex<double> bj = sqrt(double(1.0) / (double(2.0)*pi4*z)) * cos(z - double(2*l+1)*pi4);
        *ReBesselJ = bj.real();
        *ImBesselJ = bj.imag();
        complex<double> dbj = - sqrt(double(1.0) / (double(2.0)*pi4*z)) * sin(z - double(2*l+1)*pi4);
        if(ReDBesselJ) *ReDBesselJ = dbj.real();
        if(ImDBesselJ) *ImDBesselJ = dbj.imag();
        return;
    }

    int nb = l+2;
    fpv* br = GimmeMem<fpv>(nb*2);
    fpv* bi = br + nb;
    int ncalc = 0;
    if(db1slc_(x, y, nb, 0, br, bi, &ncalc))
        crash("BesselJ(%i, %e, %e) failed", l, double(x), double(y));
    *ReBesselJ = br[l];
    *ImBesselJ = bi[l];

    if(ReDBesselJ || ImDBesselJ) {
        if(fabs(x)<get_eps<fpv>() && fabs(y)<get_eps<fpv>()) {
            if(ReDBesselJ) *ReDBesselJ = (l==1) ? 0.5 : 0.0;
            if(ImDBesselJ) *ImDBesselJ = 0.0;
        }
        else {
            complex<fpv> _I(fpv(0.),fpv(1.));
            complex<fpv> z = x + _I * y;
            complex<fpv> bjl  = br[l  ] + _I * bi[l  ];
            complex<fpv> bjl1 = br[l+1] + _I * bi[l+1];
            complex<fpv> dbj = - bjl1 + fpv(l) * bjl / z;
            if(ReDBesselJ) *ReDBesselJ = dbj.real();
            if(ImDBesselJ) *ImDBesselJ = dbj.imag();
        }
    }
    FreeMem(br);
}
//======================================================================================================================

template void BesselJComplex<NativeDouble>(int l, NativeDouble x, NativeDouble y, NativeDouble* ReBesselJ, NativeDouble* ImBesselJ, NativeDouble* ReDBesselJ, NativeDouble* ImDBesselJ);

#ifdef EXTRAPRECISION_COLESO
template void BesselJComplex<dd_real>(int l, dd_real x, dd_real y, dd_real* ReBesselJ, dd_real* ImBesselJ, dd_real* ReDBesselJ, dd_real* ImDBesselJ);
template void BesselJComplex<qd_real>(int l, qd_real x, qd_real y, qd_real* ReBesselJ, qd_real* ImBesselJ, qd_real* ReDBesselJ, qd_real* ImDBesselJ);
#endif


//======================================================================================================================
// sf13d.f -- вычисление экспоненциальных интегралов:
// l=1 -- вычисляется интеграл \int\limits_{-\infty}^x e^t/t dt для x>0 и -\int\limits_{-x}^{\infty} e^{-t}/t dt для x<0
// l=2 -- вычисляется интеграл \int\limits_x^{\infty} e^{-t}/t dt для x>0
// l=3 -- вычисляется интеграл e^{-x} \int\limits_{-\infty}^x e^t/t dt для x>0 и -e^{-x} \int\limits_{-x}^{\infty} e^{-t}/t dt для x<0
//======================================================================================================================
NativeDouble sf13d_c(NativeDouble x, int l) {
    static const NativeDouble dexp40 = 2.3538526683702e17;
    static const NativeDouble xmin = -175.047544754057132;
    static const NativeDouble xmax = 179.859659356407516;
    static const NativeDouble x0 = .3725074107813666351995962;
    static const NativeDouble x01 = .3725074107805994572117924;
    static const NativeDouble x02 = 7.671772501993941448417297e-13;
    static const NativeDouble a[6] = { -.5772156649015328655494272,
            .7541643136630166166511912,.1298492329273732287520104,
            .0240681355683977412501795,.001320843092096085028691598,
            6.577393997532644860800309e-5 };
    static const NativeDouble b[6] = { 1.,.4258991938115898184813445,
            .0797794718410228254068528,.008302084760987716791080792,
            4.864271383930163900576048e-4,1.306551958228488786542564e-5 };
    static const NativeDouble c__[8] = { 8.677459548384437804144743e-8,
            .9999955193013903009813247,11.84831055549458445064203,
            45.59306442533898362512446,69.9279451291003023,
            42.5202034768840779,8.83671808803843927826449,
            .4013776649406647217821486 };
    static const NativeDouble d__[8] = { 1.,12.84819353791566509670475,
            56.44335695618032033848976,106.6451837699138813775333,
            89.73110971252897982708421,31.49718491704407341558179,
            3.795590037621222379016217,.09088045691888692323434639 };
    static const NativeDouble e[8] = { -.9999999999999734101585602,
            -34.4061995006684895,-427.532671201988539,-2396.0194324749054,
            -6168.85210055476351,-6576.09698748021179,
            -2106.077371426332888404432,-14.89908499729481694551225 };
    static const NativeDouble f[8] = { 1.,36.40619950064597887262609,
            494.3450702099036675463139,3190.27237489543301762751,
            10337.07530858409791107989,16324.14535577835067670094,
            11149.77528710966180369721,2378.138991021602237196926 };
    static const NativeDouble p0[6] = { 152.5388359511120022204979,
            340.2682862739600295753913,259.7386446160080026857031,
            77.8709658676071185823275,7.460510544921463704781672,
            .05574718225325587389606951 };
    static const NativeDouble q0[6] = { 152.5388359511120022204979,
            416.5377042495160253565701,417.161218090395379931579,
            185.7403824840773474136315,34.90362129565328430658155,2. };
    static const NativeDouble p1[9] = { 5.531977362081978988328501,
            206.3133336244559714600651,14272.34409068235072481913,
            36857.38950923287211480784,4493416.458218791522085663,
            -1684497.821007959544658662,352960804.7950283885002138,
            -125194997.4431756809353824,2997849734.461850762367254 };
    static const NativeDouble q1[9] = { 25.62890624999999999999998,
            -1512.618411191137056448496,42688.55000903745440155035,
            -747877.2860127960884710789,8857915.400539994239807139,
            -72490357.19651194289326679,401413891.4734782576560978,
            -1398366755.614922940731048,1279632488.038081765174863 };
    static const NativeDouble p2[9] = { -2.469409834483613064293194,
            -36.77831134783114919173387,23.27302338390391156508484,
            7.894722092944572056083529,-19.41329675144307032041979,
            5.886582407532811034300834,4.142039348134066489848237,
            5.731167034972717733154468,.9989576665165518232214481 };
    static const NativeDouble q2[8] = { 2.639830073180245983976988,
            965.4052174292802988020413,-8.387670841896406992432844,
            307.2972365373511820507703,52.31655687345585903358371,
            341.3652125243755222072648,-199.1496002312352011642817,
            1.146252532490162012734913 };
    static const NativeDouble p3[10] = { -1.647721172463463057411559,
            -18.60092121726437852657908,-10.00641913989284992325678,
            -21.05740799548039987598713,-.9134835699998743036021053,
            -33.23612519739351256475861,24.95487730402059156631367,
            26.52575818452800149316314,-1.845086232391278890929697,
            .9999933106160570628340167 };
    static const NativeDouble q3[9] = { 97.92403599217290022238568,
            64.03800405352415836546239,59.32515296193498244292641,
            253.2935526652900151134416,44.29413178337928158612163,
            1192.83242396860100598133,199.1004470817736304866228,
            -10.93556195391090990476357,1.001533852045342953829276 };
    static const NativeDouble p4[10] = { 175.3388012654659995348538,
            -223.1276707776323995346957,-18.19496649298688950580072,
            -27.97985286243054048327395,-7.631477016202536045597071,
            -15.28566236369296005292994,-7.068109778950294019850275,
            -5.000066404131310093816864,-3.000000003209813081994639,
            1.000000000000011102230246 };
    static const NativeDouble q4[9] = { 39784.59771674147032172193,
            3.972771091004144938807485,137.7903902357479992701883,
            117.1792205020865011988464,70.48318471804246954093281,
            -12.01877635471547001166926,-7.992435957763396991992976,
            -2.999998940403250102804122,1.999999999990480947786863 };

    /* Local variables */
    NativeDouble frac, sump, sumq, ret_val;
    int i;
    NativeDouble t, w, y;
    NativeDouble px[9], qx[9], xx0, xmx0;

    int ierr = 0;
    if (l < 1 || l > 3) crash("Errcode: 66"); // ret_val = sys059;
    if (x == 0) crash("Errcode: 67"); // ret_val = (l == 2) ? sys059 : -sys059;

    y = x;
    if(l!=2 && x>0.) {
        if (x >= 24.) {
            if (x >= xmax && l < 3) {
                crash("Errcode: 68"); // ret_val = sys059;
            }
            y = 1. / x;
            frac = 0.;
            for (i = 1; i <= 9; ++i) {
                frac = q4[i - 1] / (p4[i - 1] + x + frac);
            }
            frac = p4[9] + frac;
            ret_val = y + y * y * frac;
            if (l == 3)
                return ret_val;
            if (x > 40.)
                return ret_val * exp(x - 40.) * dexp40;
            else
                return ret_val * exp(x);
        }
        else if (x >= 12.) {
            frac = 0.;
            for (i = 1; i <= 9; ++i) {
                frac = q3[i - 1] / (p3[i - 1] + x + frac);
            }
            ret_val = (p3[9] + frac) / x;
            if (l != 3) ret_val *= exp(x);
            return ret_val;
        }
        else if (x >= 6.) {
            frac = 0.;
            for (i = 1; i <= 8; ++i) {
                frac = q2[i - 1] / (p2[i - 1] + x + frac);
            }
            ret_val = (p2[8] + frac) / x;
            if (l != 3) ret_val *= exp(x);
            return ret_val;
        }
        else {
            t = x + x;
            t = t / 3.0 - 2.0;
            px[0] = 0.;
            qx[0] = 0.;
            px[1] = p1[0];
            qx[1] = q1[0];
            for (i = 2; i <= 8; ++i) {
                px[i] = t * px[i - 1] - px[i - 2] + p1[i - 1];
                qx[i] = t * qx[i - 1] - qx[i - 2] + q1[i - 1];
            }
            sump = 0.5 * t * px[8] - px[7] + p1[8];
            sumq = 0.5 * t * qx[8] - qx[7] + q1[8];
            frac = sump / sumq;
            xmx0 = x - x01 - x02;
            if (abs(xmx0) < .037) {
                y = xmx0 / x0;
                sump = ((((p0[5] * y + p0[4]) * y + p0[3]) * y + p0[2]) * y + p0[2]) * y 
                        + p0[0];
                sumq = ((((q0[5] * y + q0[4]) * y + q0[3]) * y + q0[2]) * y + q0[1]) * y 
                        + q0[0];
                ret_val = (sump / (sumq * x0) + frac) * xmx0;
                if (l == 3) ret_val = exp(-(x)) * ret_val;
                return ret_val;
            }
            else {
                xx0 = x / x0;
                ret_val = log(xx0) + xmx0 * frac;
                if (l == 3) ret_val = exp(-(x)) * ret_val;
                return ret_val;
            }
        }
    }

    if(x < 0.) { y=-x; if(l==2) ierr = 1; }
    w = 1.0 / y;
    if (y > 4.0) {
        if (-abs(x) < xmin && l < 3) return 0.0; //crash("Errcode: 69");
        ret_val = -w * (w * (((((((e[7] * w + e[6]) * w + e[5]) * w + e[4]) * w + 
                e[3]) * w + e[2]) * w + e[1]) * w + e[0]) / (((((((f[7] * w + f[6])
                * w + f[5]) * w + f[4]) * w + f[3]) * w + f[2]) * w + f[1]) * w 
                + f[0]) + 1.);
        if (l == 3) return ret_val;
        ret_val *= exp(-y);
    }
    else if (y > 1.0) {
        ret_val = -(((((((c__[7] * w + c__[6]) * w + c__[5]) * w + c__[4]) * w + 
                c__[3]) * w + c__[2]) * w + c__[1]) * w + c__[0]) / (((((((d__[7] 
                * w + d__[6]) * w + d__[5]) * w + d__[4]) * w + d__[3]) * w +
                d__[2]) * w + d__[1]) * w + d__[0]);
        if (l == 3) return ret_val;
        ret_val *= exp(-y);
    }
    else {
        ret_val = log(y) - (((((a[5] * y + a[4]) * y + a[3]) * y + a[2]) * y +
                a[1]) * y + a[0]) / (((((b[5] * y + b[4]) * y + b[3]) * y + b[2])
                * y + b[1]) * y + b[0]);
        if (l == 3) ret_val *= exp(y);
    }

    if (l == 2) ret_val = -ret_val;
    if (ierr == 1) crash("Errcode: %i", ierr);
    return ret_val;
} /* sf13d_c */
//======================================================================================================================


//======================================================================================================================
// Функция sf34r_c вычисляет значения функции ошибок (интеграла вероятностей)
// errf(x) = 2/sqrt(pi) \int_0^x exp(-t*t) dt.
// J.F.Hart, E.W.Cheney, C.L.Lawson, Computer Approximations, Wiley, New York, 1968
//======================================================================================================================
NativeDouble errf(NativeDouble x) {
    NativeDouble t = x*x;
    if (t >= .2275) {
        t = abs(x);
        if (t >= 6.5) return SIGN(t);
        NativeDouble s = exp(-t * t) * ((((((t * .56418976853 + 6.80289887975) * t + 
                38.7114283284) * t + 131.126592722) * t + 278.597756334) * t + 
                355.969003302) * t + 224.182771538) / (((((((t + 12.0578408552) * 
                t + 69.1138352376) * t + 238.450272778) * t + 527.553795978) * t 
                + 741.521361689) * t + 608.932172346) * t + 224.182771538);
        return (1. - s) * SIGN(x);
    }
    else {
        return x * 1.1283791671 * (((((((((t * -1.45038522232e-7 + 
                1.45891690009e-6) * t - 1.32275132275e-5) * t + 1.06837606838e-4) 
                * t - 7.57575757576e-4) * t + .00462962962963) * t - 
                .0238095238095) * t + .1) * t - .3333333333333) * t + 1);
    }
}
//======================================================================================================================


//======================================================================================================================
// Функция sf74r_c -- вычисление модифицированных функций Бесселя первого рода In (x) для последовательности целых индексов
// sf74r_c вычисляет значения модифицированных функций Бесселя первого рода  Ip (x) для действительного аргумента  x 
// и последовательности целых индексов  p от  n до 0. Расчеты осуществляются по рекуррентной формуле:
//     Ip-1(x) = (2p/x)*Ip(x) + Ip+1(x) . 
// Ошибка вычисления  Ip (x) не возрастает в pекуppентном процессе, применяемом при убывании  p.
// Вычисления начинаются с задания значений  Ik + 1 = 0,   Ik = 1 и числа  k.
// Величина  k выбирается в зависимости от значений  x и  n. 
// Затем последовательно вычисляются значения функций для p = k - 1, k - 2, ..., 0.
// Для достижения заданной точности  eps рекурсия повторяется с большим  k до тех пор, пока разность между pезультатами 
// двух последних рекурсий не станет меньше  eps.
// 1. Абрамовиц М., Стиган И. Справочник по специальным функциям - М. Наука, 1979.
// 2. Goldstein M., Thaler R.M. Recurrence techniques for the calculation of Bessel functions. - M TAC, 1959, V. 13, N 66.

// Параметры 
// x - заданное значение аргумента  x (тип: вещественный);
// n1 - максимальный порядок последовательности функций   Ip (x), увеличенный на единицу (тип: целый);
// eps - относительная точность вычисления последней функции последовательности  In (x) (тип: вещественный);
// s - если TRUE, то вычисляются функции  Ip(x); если FALSE, то вычисляются функции exp(-x) Ip(x);
// t - вещественный вектоp длины n1 вычисленных значений функций Бесселя порядков от 0 до n
// t1 - вещественный вектоp длины n1, используемый как рабочий.
//======================================================================================================================
template<typename fpv>
void sf74r_c(fpv x, int n1, fpv eps, bool s, fpv *t, fpv *t1) {
    const int n = n1 - 1;
    if (x == 0.) {
        t[0] = 1.;
        for (int p = 1; p < n1; ++p) t[p] = 0.;
        return;
    }

    const fpv c__ = s ? exp(-x) : fpv(1.);
    if (c__ < 1.175495e-38) crash("sf74r_c error: argument is too much");
    const fpv y = fabs(x);
    
    int d, m;
    if (y <= 1.) {
        m = 10;
        d = 3;
    } 
    else if (y <= 10.) {
        m = int(y * 1.5) + 10;
        d = int(y * .5) + 3;
    }
    else {
        m = int(y * .5) + 25;
        d = int(y * .25 + 100. / y);
    }

    if (m <= n) m = n + d;

    bool f = false;
    while(1) {
        fpv sum = 0.;
        fpv r0 = 0.;
        fpv r1 = 1.;
        fpv r2;
        const fpv z__ = 3.4e38 / 4. / m * y;

        for (int ip1 = 1; ip1 <= m; ++ip1) {
            const int ip = m - ip1;
            r2 = (ip+1) * 2. / x * r1 + r0;
            if (ip <= n) t[ip] = r2;
            else t[n] = r2;

            if(abs(r2) > z__) {
                r1 /= z__;
                r2 /= z__;
                sum /= z__;
                for (int k1 = ip; k1 <= n; ++k1) t[k1] /= z__;
            }

            if (ip != 0) sum += r2 * 2.;
            else sum += r2;

            r0 = r1;
            r1 = r2;
        }
        sum *= c__;
        for (int ip1 = 0; ip1 < n1; ++ip1) t[ip1] /= sum;

        if (f) {
            if(fabs(t[n] - t1[n]) < eps * fabs(t[n1])) return;
        }

        f = true;
        for (int ip = 0; ip < n1; ++ip) t1[ip] = t[ip];
        m += d;
    }
}
//======================================================================================================================

template void sf74r_c<NativeDouble>(NativeDouble x, int n1, NativeDouble eps, bool s, NativeDouble *t, NativeDouble *t1);
#ifdef EXTRAPRECISION_COLESO
template void sf74r_c<dd_real>(dd_real x, int n1, dd_real eps, bool s, dd_real *t, dd_real *t1);
template void sf74r_c<qd_real>(qd_real x, int n1, qd_real eps, bool s, qd_real *t, qd_real *t1);
#endif


//======================================================================================================================
// Fresnel integrals C(x) = int_0^x cos(pi/2 t**2) dt, S(x) = int_0^x sin(pi/2 t**2) dt
// The integrals are evaluated by a power series for x < 1. For x >= 1 auxiliary functions f(x) and g(x) are employed such that
// C(x) = 0.5 + f(x) sin( pi/2 x**2 ) - g(x) cos( pi/2 x**2 )
// S(x) = 0.5 - f(x) cos( pi/2 x**2 ) - g(x) sin( pi/2 x**2 )
// Acuracy: Arithmetic  function   domain     # trials      peak         rms
//            IEEE       S(x)      0, 10       10000       2.0e-15     3.2e-16
//            IEEE       C(x)      0, 10       10000       1.8e-15     3.3e-16
// Cephes Math Library Release 2.8:  June, 2000. Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
//======================================================================================================================

/* S(x) for small x */
static NativeDouble sn[6] = {
-2.99181919401019853726E3,
 7.08840045257738576863E5,
-6.29741486205862506537E7,
 2.54890880573376359104E9,
-4.42979518059697779103E10,
 3.18016297876567817986E11,
};
static NativeDouble sd[6] = {
/* 1.00000000000000000000E0,*/
 2.81376268889994315696E2,
 4.55847810806532581675E4,
 5.17343888770096400730E6,
 4.19320245898111231129E8,
 2.24411795645340920940E10,
 6.07366389490084639049E11,
};

/* C(x) for small x */
static NativeDouble cn[6] = {
-4.98843114573573548651E-8,
 9.50428062829859605134E-6,
-6.45191435683965050962E-4,
 1.88843319396703850064E-2,
-2.05525900955013891793E-1,
 9.99999999999999998822E-1,
};
static NativeDouble cd[7] = {
 3.99982968972495980367E-12,
 9.15439215774657478799E-10,
 1.25001862479598821474E-7,
 1.22262789024179030997E-5,
 8.68029542941784300606E-4,
 4.12142090722199792936E-2,
 1.00000000000000000118E0,
};

/* Auxiliary function f(x) */
static NativeDouble fn[10] = {
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
static NativeDouble fd[10] = {
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
static NativeDouble gn[11] = {
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
static NativeDouble gd[11] = {
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

void fresnl(NativeDouble xxa, NativeDouble& ss, NativeDouble& cc) {
    NativeDouble x = fabs(xxa);
    NativeDouble x2 = x * x;
    if( x2 < 2.5625 ) {
        NativeDouble t = x2 * x2;
        ss = x * x2 * polevl( t, sn, 5)/p1evl( t, sd, 6 );
        cc = x * polevl( t, cn, 5)/polevl(t, cd, 6 );
    }
    else if( x > 36974.0 ) {
        cc = 0.5;
        ss = 0.5;
    }
    else {
        // Asymptotic power series auxiliary functions for large argument
        NativeDouble t = PiNumber * x2;
        NativeDouble u = 1.0/(t * t);
        
        t = 1.0/t;
        NativeDouble f = 1.0 - u * polevl( u, fn, 9)/p1evl(u, fd, 10);
        NativeDouble g = t * polevl( u, gn, 10)/p1evl(u, gd, 11);

        t = 0.5*PiNumber * x2;
        NativeDouble c = cos(t);
        NativeDouble s = sin(t);
        t = PiNumber * x;
        cc = 0.5  +  (f * s  -  g * c)/t;
        ss = 0.5  -  (f * c  +  g * s)/t;
    }

    if( xxa < 0.0 )    {
        cc = -cc;
        ss = -ss;
    }
}
//======================================================================================================================
