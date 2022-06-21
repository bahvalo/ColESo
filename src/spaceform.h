//======================================================================================================================
// EXACT SOLUTION MODULE
//======================================================================================================================
#pragma once
#ifndef SPACEFORM_H
#define SPACEFORM_H

#include "base_config.h"

#define SpaceFormTypes " GAUSSIAN 0  COS2 1  CONST 2 "
// Пространственные формы источника: гауссиан, квадрат косинуса
template<typename fpv>
struct tSpaceForm {
    enum tSpaceFormType {
        FORM_GAUSSIAN = 0,  // гауссиан
        FORM_COS2     = 1,  // квадрат косинуса на периоде
        FORM_CONST    = 2   // константа
    };
    fpv Aterm, Bterm; // Амплитуда и полуширина импульса
    fpv r0[3];        // координаты центра импульса
    fpv PerX, PerY, PerZ; // Период по каждому направлению
    int MaxPer[3];       // Количество дополнительных источников по каждому направлению
    int Checkerboard;    // Оставляем только чётные суммы iperx+ipery+iperz
    int NAngularPeriods; // количество образов при вращении вокруг оси OZ (1 - нет периодики)
    tSpaceFormType Form; // 0 - гауссиан, 1 - квадрат косинуса на периоде, 2 - константа
    int NormalizeForm;   // нормализация на единичный интеграл
    int numCoords;       // количество координат (нужно знать для нормализации; Read не считывает лишние координаты)

protected:
    fpv InvBTerm;
    fpv AmpliduteLinf2L1(fpv bterm); // вычисление множителя, переводящего Linf-норму функции в L1-норму

public:
    tSpaceForm() {
        Aterm=Bterm=1.0; r0[0]=r0[1]=r0[2]=0.0; PerX=PerY=PerZ=1e50; MaxPer[0]=MaxPer[1]=MaxPer[2]=0;
        Checkerboard=0; NAngularPeriods=1; Form=FORM_GAUSSIAN; NormalizeForm=0; numCoords=3;
        InvBTerm = 1.0;
    }

    void Read(class tFileBuffer &FB, int CanBePeriodic = 1); // Считывание параметров из уже открытого файла
    void Init(); // вычисление нормировочной константы
    fpv SpaceForm_rr(fpv rr, fpv* df=NULL) const; // вычисление амплитуды одного импульса и производной, если df!=NULL. На вход квадрат расстояния до центра
    fpv dSpaceForm_rr(fpv rr) const; // вычисление производной от амплитуды одного импульса по квадрату расстояния до центра
    inline fpv SpaceForm(fpv r) const { return SpaceForm_rr(r*r); } // то же, на вход расстояние до центра
    fpv dSpaceForm(fpv r) const; // d(SpaceForm(r))/dr
    fpv ddSpaceForm(fpv r) const; // d^2(SpaceForm(r))/dr^2
    fpv SpaceForm(const fpv* Coords) const; // вычисление значения для решётки импульсов в заданной точке
    // Определение нормированной амплитуды
    inline fpv GetAmplitudeL1(void) { if(NormalizeForm) return Aterm; else return Aterm*AmpliduteLinf2L1(Bterm); }
};

#endif
