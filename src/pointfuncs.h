//======================================================================================================================
// Абстрактный класс функции, возвращающей значения по координатам
//======================================================================================================================
#pragma once
#ifndef _POINTFUNCS_H
#define _POINTFUNCS_H

#include "base_lib.h" // MakeNaN

enum tFuncType {
    FUNC_SCALAR          = 0, // scalar function
    FUNC_PULSATIONS      = 1, // physical pulsations (rho', u', v', w', p')
    FUNC_PHYSICAL        = 2, // physical variables (rho, u, v, w, p)
    FUNC_CONSERVATIVE    = 3, // conservative variables (rho, rho*u, rho*v, rho*w, E) or its pulsations
    FUNC_TEMPVEL         = 4, // velocities and temperature (1, u, v, w, p/rho)
    FUNC_PULSATIONS_COMPLEX = 5, // complex physical pulsations in cyl. coord. system ( Re(rho', uz', ur', uphi', p'), Im(...) )
    FUNC_PULSCONS_COMPLEX = 6, // complex conservative pulsations in cyl. coord. system
    FUNC_VECTOR          = 7  // vector field (three components)
};

#define FuncTypeNames "Scalar 0   Puls 1   Phys 2   Cons 3   TempVel 4   Complex 5   ComplexCons 6   Vector 7"

// Основная структура для функции, принимающей на вход t и x[3] и возвращающей несколько (от 1 до 10) значений
struct tPointFunction {
    // Описание функции
    virtual const char* description() const { return NULL; }
    // Тип возвращаемого значения: скаляр, пульсации, полные переменные и т. д.
    virtual tFuncType Type() const = 0;
    // Количество переменных, вычисляемой функцией
    virtual int NumVars() const {
        switch(Type()) {
        case FUNC_SCALAR: return 1;
        case FUNC_VECTOR: return 3;
        case FUNC_TEMPVEL: return 5;
        case FUNC_PULSATIONS: case FUNC_PHYSICAL: case FUNC_CONSERVATIVE: return 5;
        case FUNC_PULSATIONS_COMPLEX: case FUNC_PULSCONS_COMPLEX: return 10;
        default: return 0;
        }
    }
    // Максимальное количество переменных (для массивов в стеке)
    static const int NumVarsMax = 10;

    // Передача в функцию параметров уравнения и характерных скоростей звука и потока
    virtual void SetEquationConstants(double /*_gam*/, double /*_Rey*/, double /*_Pr*/, 
        double /*_SoundVel*/, double /*_FlowVel*/) {}

    // Вычисление значения в точке
    virtual void PointValue(double t, const double* coor, double* V) const = 0;

    // Освобождение выделенной памяти
    virtual ~tPointFunction(void) {}

    // Считывание параметров из файла
    virtual const char* filename() const { return NULL; }
    virtual void ReadParams(class tFileBuffer&) {}
    virtual void ReadParamsFromFile(const char* fname = NULL); // default body: es_aux.cpp
    // Инициализация (после считывания параметров)
    virtual void Init() {}

private: // Только для Noisette (вне ColESo)
    friend struct s_function;
    virtual void Init_AfterMeshRead() {} // заключительная часть инициализации - после чтения сетки
};


// Структура для функции, принимающей на вход t и x[3] и возвращающей несколько (от 1 до 10) значений,
// для которой есть реализации различной точности
template<typename fpv>
struct tPointFunction_EP : tPointFunction {
    virtual void PointValue(fpv t, const fpv* coor, fpv* V) const = 0; // основная функция
    // поскольку у базового класса есть функция от double, нужно сделать ей реализацию по умолчанию
    // вынести в *.cpp это тело невозможно, так как эта функция отсутствует при fpv = double (см. ниже)
    void PointValue(double t, const double* coor, double* V) const {
        fpv _t = t;
        fpv _coor[3] = {coor[0], coor[1], coor[2]};
        fpv _V[NumVarsMax];
        PointValue(_t, _coor, _V);
        const int numVars = NumVars();
        for(int ivar=0; ivar<numVars; ivar++) V[ivar] = (double)_V[ivar]; 
    }
};
// tPointFunction_EP<double> совпадает с tPointFunction
template<> struct tPointFunction_EP<double> : tPointFunction {};


template<tFuncType type>
struct tTypePointFunc : tPointFunction {
    tFuncType Type() const OVERRIDE { return type; }
};

struct tScalarFunction : tTypePointFunc<FUNC_SCALAR> {};
struct tPulsFunction : tTypePointFunc<FUNC_PULSATIONS> {};
struct tComplexPulsFunction : tTypePointFunc<FUNC_PULSATIONS_COMPLEX> {};

struct tPhysFunction : tTypePointFunc<FUNC_PHYSICAL>{
// Абстрактный класс фонового поля
// По сравнению с абстрактным классом функции имеет параметры уравнений Эйлера,
// а также характерные скорости звука (SoundVel) и потока (FlowVel).
// FlowMach = 1.0/SoundVel; RefMach = FlowVel/SoundVel
    double gam, Rey, Pr; // Параметры уравнений Навье-Стокса
    double SoundVel, FlowVel; // Характерные скорости потока и звука
    tPhysFunction() { MakeNaN(SoundVel); MakeNaN(FlowVel); MakeNaN(gam); MakeNaN(Rey); MakeNaN(Pr); }
    void SetEquationConstants(double _gam, double _Rey, double _Pr, double _SoundVel, double _FlowVel) OVERRIDE {
        gam = _gam; Rey = _Rey; Pr = _Pr; SoundVel = _SoundVel; FlowVel = _FlowVel;
    }
};

struct tArbitraryPointFunction : tPointFunction {
// По сравнению с абстрактным классом функции имеет переменную типа и все параметры от tPhysFunction
    tFuncType FuncType; // тип данных (пульсации/полные)
    double gam, Rey, Pr; // Параметры уравнений Навье-Стокса
    double SoundVel, FlowVel; // Характерные скорости потока и звука

    tArbitraryPointFunction() {
        MakeNaN(SoundVel); MakeNaN(FlowVel); MakeNaN(gam); MakeNaN(Rey); MakeNaN(Pr); 
        FuncType = FUNC_PHYSICAL; 
    }
    virtual ~tArbitraryPointFunction(){}
    virtual tFuncType Type() const OVERRIDE { return FuncType; }
    void SetEquationConstants(double _gam, double _Rey, double _Pr, double _SoundVel, double _FlowVel) OVERRIDE {
        gam = _gam; Rey = _Rey; Pr = _Pr; SoundVel = _SoundVel; FlowVel = _FlowVel;
    }
};

//======================================================================================================================
struct s_Zero : tPulsFunction {
//======================================================================================================================
    const char* description() const OVERRIDE { return "zero"; }
    void PointValue(double, const double*, double* V) const OVERRIDE {
        for(int ivar=0; ivar<NumVars(); ivar++) V[ivar] = 0.0; 
    }
};
//======================================================================================================================


//======================================================================================================================
// Вычисление производной по времени от некоторой PointFunction по схеме 2-го порядка
//======================================================================================================================
template<typename fpv>
void ddt_PointValue(const tPointFunction_EP<fpv>& F, fpv t, const fpv* coor, fpv* V) {
    const fpv dt = (sizeof(fpv)==32) ? 1e-24 : ((sizeof(fpv)==16) ? 1e-12 : 1e-6);
    const int NumVarsMax = tPointFunction_EP<fpv>::NumVarsMax;
    const int NumVars = F.NumVars();
    if(t-0.5*dt >= 0.0 || t < 0.0) {
        fpv UApuls[NumVarsMax], UBpuls[NumVarsMax];
        F.PointValue(t+0.5*dt, coor, UApuls);
        F.PointValue(t-0.5*dt, coor, UBpuls);
        for(int ivar=0; ivar<NumVars; ivar++) 
            V[ivar] = (UApuls[ivar] - UBpuls[ivar]) / dt; 
    }
    else {
        fpv U0[NumVarsMax], U1[NumVarsMax], U2[NumVarsMax];
        F.PointValue(t+dt, coor, U2);
        F.PointValue(t+0.5*dt, coor, U1);
        F.PointValue(t, coor, U0);
        for(int ivar=0; ivar<NumVars; ivar++) 
            V[ivar] = (- U2[ivar] + 4.*U1[ivar] - 3.*U0[ivar]) / dt; 
    }
}
inline void ddt_PointValue(const tPointFunction& F, double t, const double* coor, double* V) {
    ddt_PointValue((const tPointFunction_EP<double>&)F, t, coor, V);
}
//======================================================================================================================

#endif
