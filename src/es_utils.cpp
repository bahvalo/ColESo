// *********************************************************************************************************************
// *****                                                                                                           *****
// *****                                  Collection of exact solutions (ColESo)                                   *****
// *****                       Auxiliary classes: impulse form, numerical quadratures etc.                         *****
// *****                                                                                                           *****
// *********************************************************************************************************************

#include "base_parser.h"
#include "es_utils.h"
#include "pointfuncs.h"
#include "geom_primitive.h"
#include "spaceform.h"
#include "gaussrule.h"
#ifdef EXTRAPRECISION_COLESO
#include "extraprecision.h"
#endif
#include "es_jacobi_rule.hpp" // CGQF subroutine computes knots and weights of a Gauss quadrature formula


//======================================================================================================================
//
//                                                Pulse form structure
//
//======================================================================================================================

//----------------------------------------------------------------------------------------------------------------------
// Reading parameters
//----------------------------------------------------------------------------------------------------------------------
template<typename fpv>
void tSpaceForm<fpv>::Read(tFileBuffer& FB, int CanBePeriodic) {
    tParamManager PM;
    PM.RequestParameter(Form, "Form", SpaceFormTypes, IO_DONTCRASH);
    PM.ReadParamsFromBuffer(FB);

    if(Form==tSpaceForm::FORM_CONST) return;

    PM.clear();
    // Амплитуда и полуширина импульса (по умолчанию равны 1)
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

    if(Bterm < get_eps<fpv>()) crash("Bterm = %e too small", double(Bterm));
    InvBTerm = 1.0/Bterm;

    if(fabs(PerX) > 0.5*huge) MaxPer[0] = 0;
    if(fabs(PerY) > 0.5*huge) MaxPer[1] = 0;
    if(fabs(PerZ) > 0.5*huge) MaxPer[2] = 0;

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
//                                      Weights and points of Gaussian quadrature rules
//
//======================================================================================================================

template<typename fpv>
void GaussPointsInit(int NumPoints, fpv* Nodes, fpv* Weights, tGItype gi_type) {
    if(NumPoints<=0) return;
    #ifdef ES_JACOBI_RULE_HPP
        if(gi_type == GI_LEGENDRE) { // int_0^1 f(x) dx
            cdgqf<fpv>(NumPoints, 1 /* Gauss-Legendre */, 0.0 /*unused*/, -0.5 /*unused*/, Nodes, Weights );

            //  Scale the quadrature formula to [0,1]
            for (int k = 0; k < NumPoints; k++ ) { Nodes[k] = 0.5 + 0.5 * Nodes[k]; Weights[k] *= 0.5; }
            return;
        }
        if(gi_type == GI_JACOBI1) { //  int_{-1}^1 f(x) / sqrt(1+x) dx.
            cdgqf<fpv>(NumPoints, 4 /* Gauss-Jacobi */, 0.0 /*alpha*/, -0.5 /*beta*/, Nodes, Weights );

            //  Scale the quadrature formula to [0,1]
            fpv p = sqrt(fpv(0.5));
            for (int k = 0; k < NumPoints; k++ ) { Nodes[k] = 0.5 + 0.5 * Nodes[k]; Weights[k] *= p; }
            return;
        }
        if(gi_type == GI_LAGUERRE) { // int_0^{\infty} f(x) exp(-x) dx.
            cdgqf<fpv>(NumPoints, 5 /* Generalized Laguerre */, 0.0 /*alpha*/, 0.0 /*beta*/, Nodes, Weights );
            return;
        }
        crash("Wrong quadrature type %i", gi_type);
    #else
        // backup version - just for a case
        if(gi_type != GI_LEGENDRE) crash("ES_JACOBI_RULE is not available");
        GaussLegendrePointsInit<fpv>(NumPoints, Nodes, Weights);
    #endif
}
//======================================================================================================================


//----------------------------------------------------------------------------------------------------------------------
// Nanotechlology
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
// Instantiating template classes and functions
//======================================================================================================================
#define INSTANTIATE(T) \
    template void GaussPointsInit<T>(int NumPoints, T* Nodes, T* Weights, tGItype); \
    template struct tGaussIntegrator<T>; \
    template struct tSpaceForm<T>;

INSTANTIATE(NativeDouble)
#ifdef EXTRAPRECISION_COLESO
INSTANTIATE(dd_real)
INSTANTIATE(qd_real)
#endif
