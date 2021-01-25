//======================================================================================================================
// EXACT SOLUTION MODULE
//======================================================================================================================
#pragma once
#ifndef ES_UTILS_HEADER
#define ES_UTILS_HEADER

#include "personal.h"


// Вычисление биномиальных коэффициентов
inline int Factorial(int n) {int i, f = 1; for(i=2; i<=n; i++) f*=i; return f;}
inline int BinomialCoeff0(int n, int k){return Factorial(n) / Factorial(k) / Factorial(n-k);}


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

enum tGItype {
    GI_LEGENDRE = 0,  // Gauss-Legendre quadratures for int_0^1 f(x) dx
    GI_JACOBI1  = 1   // Gauss-Jacobi quadratures for int_0^1 f(x)/sqrt(x) dx
};

template<typename fpv>
void GaussPointsInit(int NumPoints, fpv* Nodes, fpv* Weights, tGItype gi_type);

template<typename fpv>
struct tGaussIntegrator {
    int GR;
    fpv *GN, *GC;

    tGaussIntegrator() { GN=GC=NULL; GR=4*sizeof(fpv); }
    tGaussIntegrator(int _GR) { GN=GC=NULL; GR=0; Init(_GR); }
    tGaussIntegrator(const tGaussIntegrator<fpv>& gi) { GN=GC=NULL; *this = gi; }
    inline void Clear() { GN=GC=NULL; GR=0; } // not a destructor - just nullify pointers and data
    inline void Free() { if(GN!=NULL) delete[] GN; GN=GC=NULL; } // do not touch GR
    ~tGaussIntegrator() { Free(); }

    inline void Init(int _GR = 0, tGItype GItype = GI_LEGENDRE) { // выделение памяти и определение узлов и коэффициентов квадратуры
        Free();
        if(_GR==0) _GR = GR;
        if(_GR <= 0 || _GR > 0x7FFF) return;
        GR = _GR;
        GN = new fpv[GR*2];
        GC = GN + GR;
        GaussPointsInit<fpv>(GR, GN, GC, GItype);
    }
    inline tGaussIntegrator<fpv>& operator=(const tGaussIntegrator<fpv> &gi){
        if(this == &gi) return *this;
        Free();
        GR = gi.GR;
        if(GR>0) {
            GN = new fpv[GR*2]; 
            GC = GN + GR;
            for(int i=0; i<GR*2; ++i) GN[i]=gi.GN[i];
        }
        return *this;
    }
};


template<typename fpv, int N> struct tFixArray {
    fpv V[N];
    inline tFixArray<fpv,N>& operator=(fpv x){ for(int i=0; i<N; i++) V[i]=x;  return *this; }
    inline tFixArray<fpv,N>& operator*=(fpv x){ for(int i=0; i<N; i++) V[i]*=x;  return *this; }
    inline tFixArray<fpv,N>& operator+=(const tFixArray<fpv,N>& object){ for(int i=0; i<N; i++) V[i]+=object.V[i]; return *this; }
    inline tFixArray<fpv,N>& operator*=(const tFixArray<fpv,N>& object){ for(int i=0; i<N; i++) V[i]*=object.V[i]; return *this; }
    inline       fpv& operator[](int i){ return V[i]; }
    inline const fpv& operator[](int i)const{ return V[i]; }
};
#pragma warning ( push )
#pragma warning (disable:4701) // warning C4701: potentially uninitialized local variable 'R' used
template<typename fpv, int N> inline tFixArray<fpv,N> operator+(const tFixArray<fpv,N> &a, const tFixArray<fpv,N> &b)
    { tFixArray<fpv,N> R; for(int i=0; i<N; i++) R.V[i]=a.V[i]+b.V[i]; return R; }
template<typename fpv, int N> inline tFixArray<fpv,N> operator-(const tFixArray<fpv,N> &a, const tFixArray<fpv,N> &b)
    { tFixArray<fpv,N> R; for(int i=0; i<N; i++) R.V[i]=a.V[i]-b.V[i]; return R; }
template<typename fpv, int N> inline tFixArray<fpv,N> operator*(const tFixArray<fpv,N> &a, fpv b)
    { tFixArray<fpv,N> R; for(int i=0; i<N; i++) R.V[i]=a.V[i]*b; return R; }
template<typename fpv, int N> inline tFixArray<fpv,N> operator*(fpv b, const tFixArray<fpv,N> &a)
    { tFixArray<fpv,N> R; for(int i=0; i<N; i++) R.V[i]=a.V[i]*b; return R; }
#pragma warning ( pop )


// Calculation using compound Gauss formula of the integral
// int exp(-alpha^2 (x-x0)^2/2) (x-x0)^k h(x) x^gamma dx, x = xmin..xmax,
// where k=0,1,... (not big), alpha>0, gamma=0 (mode=0) or -1/2 (mode=1),  
// h -- analytic function on [0..xmax] with convergence raduis >=r at each point, q = 1/r.
// If mode==1, xmin should be zero
template<typename fpv>
struct tCompoundGaussIntegrator {
    tGaussIntegrator<fpv> GLI, GJI; // Gauss-Legendre and Gauss-Jacobi integrators
    double reducer; // Workaround: multiply number of quadrature nodes by this value (set < 1 to reduce computation time)
    tCompoundGaussIntegrator() { reducer = 1.0; }
    inline void Clear() { GLI.Clear(); GJI.Clear(); } // not a destructor - just nullify pointers and data
    inline void Free() { GLI.Free(); GJI.Free(); } // the same as destructor
    inline void Init(int GR = 0) {
        if(GR<=0) GR = sizeof(fpv)*4;
        GLI.Init(GR, GI_LEGENDRE); 
        GJI.Init(GR, GI_JACOBI1);
    }
    // Automatic integration subroutine
    template<int N> tFixArray<fpv ,N> Integrate(fpv xmin, fpv xmax, fpv alpha, fpv x0, int k, int mode, 
                                            fpv M, fpv q, tFixArray<fpv ,N> (*fun_ptr)(fpv, void*), void*) const;
};

#endif
