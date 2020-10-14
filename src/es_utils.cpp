// *********************************************************************************************************************
// *****                                                                                                           *****
// *****                                  Collection of exact solutions (ColESo)                                   *****
// *****                       Auxiliary classes: impulse form, numerical quadratures etc.                         *****
// *****                                                                                                           *****
// *********************************************************************************************************************

#include "parser.h"
#include "es_utils.h"
#include "pointfuncs.h"
#include "geom_primitive.h" 
#ifdef _NOISETTE
#include "lib_base.h"
#endif
#ifdef EXTRAPRECISION_COLESO
#include "extraprecision.h"
#endif
#include "es_jacobi_rule.hpp" // CGQF subroutine computes knots and weights of a Gauss quadrature formula


//======================================================================================================================
//
//                                                Структура для формы импульса
//
//======================================================================================================================

//----------------------------------------------------------------------------------------------------------------------
// Считывание параметров
//----------------------------------------------------------------------------------------------------------------------
template<typename fpv>
void tSpaceForm<fpv>::Read(tFileBuffer& FB, int CanBePeriodic) {
    tParamManager PM;
    PM.RequestParameter(Form, "Form", SpaceFormTypes, IO_DONTCRASH);
    PM.ReadParamsFromBuffer(FB);

    if(Form==tSpaceForm::FORM_CONST) return;

    PM.clear();
    // Амплитуда и полуширина импульса (по умлочанию равны 1)
    PM.Request(Aterm, "Aterm");
    PM.Request(Bterm, "Bterm");
    PM.RequestOption(NormalizeForm, "NormalizeForm");
    // для нормировки нужно знать размерность (numCoords), но её не считываем из файла

    // Положение центра импульса (по умолчанию -- 0,0,0)
    PM.Request(r0[Coor_X], "CoorX");
    if(numCoords>=2) PM.Request(r0[Coor_Y], "CoorY");
    if(numCoords>=3) PM.Request(r0[Coor_Z], "CoorZ");
    // Временные синонимы (оставить что-то одно - TODO)
    PM.Request(r0[Coor_X], "Xterm");
    if(numCoords>=2) PM.Request(r0[Coor_Y], "Yterm");
    if(numCoords>=3) PM.Request(r0[Coor_Z], "Zterm");

    if(CanBePeriodic) {
        // Пространственный период (по умолчанию -- huge) и число доп. копий (по умолчанию -- 0)
        PM.Request(PerX, "PerX");
        if(numCoords>=2) PM.Request(PerY, "PerY");
        if(numCoords>=3) PM.Request(PerZ, "PerZ");
        PM.Request(MaxPer[0], "NumPeriodsX");
        if(numCoords>=2) PM.Request(MaxPer[1], "NumPeriodsY");
        if(numCoords>=3) PM.Request(MaxPer[2], "NumPeriodsZ");
        if(numCoords>=2) PM.Request(NAngularPeriods, "NAngularPeriods");
        // Оставляем только чётные суммы iperx+ipery+iperz
        if(numCoords>=2) PM.RequestOption(Checkerboard, "Checkerboard");
    }
    PM.ReadParamsFromBuffer(FB);
}

//----------------------------------------------------------------------------------------------------------------------
// Вычисление множителя, переводящего Linf-норму функции в L1-норму
//----------------------------------------------------------------------------------------------------------------------
template<typename fpv>
fpv tSpaceForm<fpv>::AmpliduteLinf2L1(fpv bterm) {
    const fpv Pi = GetPiNumber<fpv>();
    const fpv NormalizingConst2D = (0.5 * (Pi*Pi - 4.0)) / Pi;
    const fpv NormalizingConst3D = (fpv(2.)/fpv(3.) * (Pi*Pi - 6.0)) / Pi;

    switch ( Form ) {
    case tSpaceForm::FORM_GAUSSIAN:
        return pow((SQR(Bterm)*Pi)/GetLn2<fpv>(), 0.5*numCoords);
    case tSpaceForm::FORM_COS2:
        switch(numCoords) {
        case 1: return bterm;
        case 2: return bterm * bterm * NormalizingConst2D;
        case 3: return bterm * bterm * bterm * NormalizingConst3D; break;
        default: crash("tSpaceForm::AmpliduteLinf2L1: wrong numCoords %i", numCoords);
        }
    default: crash("tSpaceForm::AmpliduteLinf2L1: form %i is undefined!", Form);
    }
}

//----------------------------------------------------------------------------------------------------------------------
// Вычисление InvBTerm=1/Bterm и нормировка амплитуды
//----------------------------------------------------------------------------------------------------------------------
template<typename fpv>
void tSpaceForm<fpv>::Init() {
    if(Form==tSpaceForm<fpv>::FORM_CONST) return;

    if(Bterm < get_eps<fpv>()) crash("tSourceStruct::Init: bterm = %e too small", double(Bterm));
    InvBTerm = 1.0/Bterm;

    if(NormalizeForm) {
        Aterm /= AmpliduteLinf2L1(Bterm);
        NormalizeForm = 0;
    }
}


//----------------------------------------------------------------------------------------------------------------------
// Вычисление значения для одного импульса
//----------------------------------------------------------------------------------------------------------------------
template<typename fpv>
fpv tSpaceForm<fpv>::SpaceForm_rr(fpv rr, fpv* df) const {
    if(Form == tSpaceForm::FORM_CONST) { 
        if(df!=NULL) *df = 0.0;
        return Aterm;
    }

    switch ( Form ) {
    case tSpaceForm::FORM_GAUSSIAN: 
    {
        fpv alpha = GetLn2<fpv>()*InvBTerm*InvBTerm;
        fpv expr = Aterm * exp(-alpha*rr);
        if(df!=NULL) *df = expr * (-alpha);
        return expr;
    }
    case tSpaceForm::FORM_COS2:
    {
        if(df!=NULL) *df = dSpaceForm_rr(rr);
        if(rr > Bterm*Bterm) return 0.0;
        fpv Expr = cos(0.5*GetPiNumber<fpv>()*sqrt(rr)*InvBTerm);
        return Aterm*Expr*Expr;
    }
    default: crash("tSpaceForm::SpaceForm: form %i is undefined!", Form); break;
    }
}


//----------------------------------------------------------------------------------------------------------------------
// Вычисление производной от функции формы импульса по квадрату расстояния до центра
//----------------------------------------------------------------------------------------------------------------------
template<typename fpv>
fpv tSpaceForm<fpv>::dSpaceForm_rr(fpv rr) const {
    if(Form == tSpaceForm::FORM_CONST) return 0.0;

    switch ( Form ) {
    case tSpaceForm::FORM_GAUSSIAN:
    {
        fpv alpha = GetLn2<fpv>()*InvBTerm*InvBTerm;
        return - Aterm * exp(-alpha*rr) * alpha;
    }
    case tSpaceForm::FORM_COS2:
    {
        if(rr > Bterm*Bterm) return 0.0;
        fpv buf = 0.5*GetPiNumber<fpv>()*InvBTerm;
        if(rr<get_eps<fpv>()) return -Aterm * buf*buf;
        fpv r = sqrt(rr);
        return -Aterm * sin(2.0*buf*r) * 0.5*buf/r;
    }
    default: crash("tSpaceForm::dSpaceForm_rr: form %i is undefined!", Form); break;
    }
}

//----------------------------------------------------------------------------------------------------------------------
// Вычисление производной от функции формы импульса по расстоянию до центра
//----------------------------------------------------------------------------------------------------------------------
template<typename fpv>
fpv tSpaceForm<fpv>::dSpaceForm(fpv r) const {
    switch ( Form ) {
    case tSpaceForm::FORM_CONST: return 0.0;
    case tSpaceForm::FORM_GAUSSIAN: {
        // return 2.*r*dSpaceForm_rr(r*r);
        fpv alpha = GetLn2<fpv>()*InvBTerm*InvBTerm;
        fpv buf = 2.*alpha*r;
        return Aterm * exp(-alpha*r*r) * (-buf);
        }
    case tSpaceForm::FORM_COS2:
    {
        if(r > Bterm) return 0.0;
        fpv buf = GetPiNumber<fpv>()*InvBTerm;
        return -0.5*Aterm * sin(buf*r) * buf;
    }
    default: crash("tSpaceForm::dSpaceForm: form %i is undefined!", Form); break;
    }
}

//----------------------------------------------------------------------------------------------------------------------
// Вычисление второй производной от функции формы импульса по расстоянию до центра
//----------------------------------------------------------------------------------------------------------------------
template<typename fpv>
fpv tSpaceForm<fpv>::ddSpaceForm(fpv r) const {
    switch ( Form ) {
    case tSpaceForm::FORM_CONST: return 0.0;
    case tSpaceForm::FORM_GAUSSIAN: 
    {
        fpv alpha = GetLn2<fpv>()*InvBTerm*InvBTerm;
        fpv buf = 2.*alpha*r;
        return Aterm * exp(-alpha*r*r) * (buf*buf - 2.*alpha);
    }
    case tSpaceForm::FORM_COS2:
    {
        if(r > Bterm) return 0.0;
        fpv buf = GetPiNumber<fpv>()*InvBTerm;
        return -0.5*Aterm * cos(buf*r) * buf*buf;
    }
    default: crash("tSpaceForm::dSpaceForm: form %i is undefined!", Form); break;
    }
}

//----------------------------------------------------------------------------------------------------------------------
// Вычисление значения для решётки импульсов
//----------------------------------------------------------------------------------------------------------------------
template<typename fpv>
fpv tSpaceForm<fpv>::SpaceForm(const fpv* Coords) const {
    fpv r[3] = {0.0, 0.0, 0.0};

    fpv sum = 0.0;
    const fpv dphi = GetPiNumber2<fpv>() / NAngularPeriods;
    for(int iPer0 = -MaxPer[0]; iPer0 <= MaxPer[0]; iPer0++) {
    for(int iPer1 = -MaxPer[1]; iPer1 <= MaxPer[1]; iPer1++) {
    for(int iPer2 = -MaxPer[2]; iPer2 <= MaxPer[2]; iPer2++) {
        if(Checkerboard) if((iPer0+iPer1+iPer2)&1) continue;
        for(int iperPhi = 0; iperPhi<NAngularPeriods; iperPhi++) {
            fpv c[3] = {Coords[0],Coords[1],Coords[2]};
            if(iperPhi) RotateVector2D(c, dphi*iperPhi);
            r[0] = c[0] - r0[0];
            r[1] = c[1] - r0[1];
            r[2] = c[2] - r0[2];
            if(MaxPer[0] || MaxPer[1] || MaxPer[2]) {
                r[0] += PerX * iPer0;
                r[1] += PerY * iPer1;
                r[2] += PerZ * iPer2;
            }
            else {
                RoundToCentre(r[0], PerX);
                RoundToCentre(r[1], PerY);
                RoundToCentre(r[2], PerZ);
            }
            fpv expr = SpaceForm_rr(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
            if(IsNaN(expr)) crash("tSpaceForm::SpaceForm error: NaN detected");

            sum += expr;
        }
    }}}
    return sum;
}
//======================================================================================================================



//======================================================================================================================
//
//                                          Определение узлов и весов квадратур Гаусса
//
//======================================================================================================================

//======================================================================================================================
// Вычисление значения полинома Лежандра, определённого на отрезке [0,1], в точке x
//======================================================================================================================
template<typename real>
real _LegendrePolynomial(real x, int order, real* dp) {
    if(order<=0) return 1.0;
    real y = 2.0*x-1.0; // переводим из [0,1] в [-1,1]
    real p1=0.0, p2=1.0;
    real p3=y;
    for(int n=1; n<=order-1; n++) {
        p1 = p2;
        p2 = p3;
        p3 = ((2*n+1)*y*p2 - n*p1)/(n+1);
    }
    if(dp!=NULL) *dp = 2.0*order*(p2-y*p3)/(1.0-y*y);
    return p3;
}

double LegendrePolynomial(double x, int order, double* dp) { return _LegendrePolynomial<double>(x, order, dp); }
//======================================================================================================================

//======================================================================================================================
// Вычисление положения i-го узла k-точечной гауссовой квадратуры на (0,1)
//======================================================================================================================
template<typename real>
real GetGaussPoint(int i, int k) {
    if(i<0 || i>=k) crash("GetGaussPoint error: i = %i, k = %i", i, k);
    if(k==1) return 0.5;
    // Найдём корни многочлена Лежандра ньютоновским процессом
    // Источник: ru.wikipedia.org/wiki/Многочлены_Лежандра
    real x = cos(PiNumber*real(NativeDouble(4*i+3))/(4*k+2));
    x = 0.5 + 0.5 * x; // переходим к [0,1]
    for(int j=0; j<100; j++) { // будем делать 100 итераций Ньютоновского процесса
        real dp;
        real p = _LegendrePolynomial<real>(x, k, &dp);
        x -= p / dp;
    }
    return x;
}
//======================================================================================================================


//======================================================================================================================
// Заполнение массива коэффициентов для квадратур Гаусса int f(x) dx, x=0..1, и int f(x)/sqrt{x} dx, x=0..1
//======================================================================================================================
template<typename fpv>
void GaussPointsInit(int NumPoints, fpv* Nodes, fpv* Weights, tGItype gi_type) {
    if(NumPoints<=0) return;
    if(gi_type == GI_LEGENDRE) {
        if(NumPoints==1) { Nodes[0] = 0.5; Weights[0] = 1.0; return; }

        // Заполним положения квадратурных точек
        for(int i=0; i < NumPoints; i++) {
            Nodes[i] = GetGaussPoint<fpv>(i, NumPoints);
        }
        // Веса квадратуры определяются по формуле: 
        // [СПРАВОЧНИК ПО СПЕЦИАЛЬНЫМ ФУНКЦИЯМ Под редакцией М. АБРАМОВИЧА и И. СТИГАН. Формула 25.4.29]
        for(int i=0; i < NumPoints; i++) {
            fpv dp;
            _LegendrePolynomial<fpv>(Nodes[i], NumPoints, &dp);
            fpv w = 1.0/(Nodes[i] * (1.0-Nodes[i]) * dp * dp);
            Weights[i] = w;
        }
        return;
    }
    #ifdef ES_JACOBI_RULE_HPP
    if(gi_type == GI_JACOBI1) {
        //  Compute the Gauss quadrature formula int_{-1}^1 f(x) / sqrt(1+x) dx.
        cdgqf<fpv>(NumPoints, 4 /* Gauss-Jacobi */, 0.0 /*alpha*/, -0.5 /*beta*/, Nodes, Weights );

        //  Scale the quadrature formula to [0,1]
        fpv p = sqrt(fpv(0.5));
        for (int k = 0; k < NumPoints; k++ ) { Nodes[k] = 0.5 + 0.5 * Nodes[k]; Weights[k] *= p; }
        return;
    }
    #endif
    crash("GaussPointsInit: wrong quadrature type %i", gi_type);
}
//======================================================================================================================

//----------------------------------------------------------------------------------------------------------------------
// Нанотехнологии
//----------------------------------------------------------------------------------------------------------------------

#ifdef EXTRAPRECISION_HEADER
int IsNaN(dd_real x) {
    return IsNaN(x[0]) || IsNaN(x[1]); 
}
int IsNaN(qd_real x) {
    return IsNaN(x[0]) || IsNaN(x[1]) || IsNaN(x[2]) || IsNaN(x[3]); 
}
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
template<typename fpv> template<int N>
tFixArray<fpv,N> tCompoundGaussIntegrator<fpv>::Integrate(fpv xmin, fpv xmax, fpv alpha, fpv x0, int k, int mode, 
                                             fpv M, fpv q, tFixArray<fpv,N> (*func)(fpv, void*), void* args) const {
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
    tFixArray<fpv,N> sum;
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
    template void GaussPointsInit<T>(int NumPoints, T* Nodes, T* Weights, tGItype); \
    template struct tGaussIntegrator<T>; \
    template struct tCompoundGaussIntegrator<T>; \
    template tFixArray<T,1> tCompoundGaussIntegrator<T>::Integrate<1>(T xmin, T xmax, T alpha, T x0, int k, int mode, \
                                             T M, T q, tFixArray<T,1> (*func)(T, void*), void* args) const;   \
    template tFixArray<T,2> tCompoundGaussIntegrator<T>::Integrate<2>(T xmin, T xmax, T alpha, T x0, int k, int mode, \
                                             T M, T q, tFixArray<T,2> (*func)(T, void*), void* args) const;   \
    template tFixArray<T,3> tCompoundGaussIntegrator<T>::Integrate<3>(T xmin, T xmax, T alpha, T x0, int k, int mode, \
                                             T M, T q, tFixArray<T,3> (*func)(T, void*), void* args) const;   \
    template struct tSpaceForm<T>;

INSTANTIATE(NativeDouble)
#ifdef EXTRAPRECISION_COLESO
INSTANTIATE(dd_real)
INSTANTIATE(qd_real)
#endif
