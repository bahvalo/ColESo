//======================================================================================================================
// Gaussian quadrature rules for a segment
//======================================================================================================================
#pragma once
#ifndef GAUSSRULE_HEADER
#define GAUSSRULE_HEADER

#include <math.h>
#include "base_lib.h"

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

inline double LegendrePolynomial(double x, int order, double* dp=NULL) { 
    return _LegendrePolynomial<double>(x, order, dp); 
}

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

template<typename fpv>
void GaussLegendrePointsInit(int NumPoints, fpv* Nodes, fpv* Weights) {
    if(NumPoints<=0) return;
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
}

#endif