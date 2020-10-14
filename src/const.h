//======================================================================================================================
// Основные определения
//======================================================================================================================
//#pragma once // Выключено - мешает компиляции от ПБ для EXTRAPRECISION через универсальный интерфейс
#include "personal.h" // чтобы знать флаг EXTRAPRECISION

#ifndef CONSTANTS_HEADER
#define CONSTANTS_HEADER

//----------------------------------------------------------------------------------------------------------------------
// ПОЗИЦИИ КООРДИНАТ И ПЕРЕМЕННЫХ
//----------------------------------------------------------------------------------------------------------------------
#define Coor_X 0                  //позиция X соординаты
#define Coor_Y 1                  //позиция Y координаты
#define Coor_R 1                  //для ZR-геометрии — позиция R координаты
#define Coor_Z 2                  //позиция Z координаты

// Stiefel (cell-centered) moments 
#define Coor_XX    3
#define Coor_XY    4
#define Coor_XZ    5
#define Coor_YY    6
#define Coor_YZ    7
#define Coor_ZZ    8

#define Var_R  0                  //позиция плотности
#define Var_U  1                  //позиция Х-скорости или X-импульса
#define Var_V  2                  //позиция Y-скорости или Y-импульса
#define Var_W  3                  //позиция Z-скорости или Z-импульса
#define Var_P  4                  //позиция давления
#define Var_E  4                  //позиция полной энергии
#define Var_N  5                  //число базовых переменных без учета турбуля
#define Var_NN 25                 //квадрат предыдушего числа

#define Var_Nu 5                  //позиция Ню       
#define Var_K  5                  //позиция К
#define Var_D  6                  //позиция Эпсилон

#define Var_Reserv 7              //стратегический резерв
#define Var_NumMax 8              //пусть уж кратно 64 байтам


//----------------------------------------------------------------------------------------------------------------------
// ГЛОБАЛЬНЫЕ УТИЛЬ-МАКРОСЫ
//----------------------------------------------------------------------------------------------------------------------
template<typename T>
inline void SWAP(T &a, T &b){ T swapbuf = a; a = b; b = swapbuf;}

#define MIN(a,b) ( ((a)<(b)) ? (a) : (b))
#define MAX(a,b) ( ((a)<(b)) ? (b) : (a))
#define SQR(X) ((X)*(X))
#define SIGN(a) ( ((a)<0.0) ? -1.0 : (((a)>0.0) ? 1.0 : 0.0) )
#define VDOT(a,b) ((a)[0]*(b)[0] + (a)[1]*(b)[1] + (a)[2]*(b)[2] )
#define VELDOT(a,b) (a[Var_U]*b[Var_U] + a[Var_V]*b[Var_V] + a[Var_W]*b[Var_W] ) // 5
#define DIST(v1, v2) (sqrt(SQR(v1[0]-v2[0])+SQR(v1[1]-v2[1])+SQR(v1[2]-v2[2]))) 
//----------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------------------
// constants
//----------------------------------------------------------------------------------------------------------------------

// epsilons
#define tiny    1e-16  // #define DBL_EPSILON   2.2204460492503131e-016 /* smallest such that 1.0+DBL_EPSILON != 1.0 */
#define tinyflt 1e-8   // #define FLT_EPSILON   1.192092896e-07F        /* smallest such that 1.0+FLT_EPSILON != 1.0 */
#define huge    1e+16
#define hugeflt 1e+8  

// Undefs
#define UNDEF -1E+35
#define UNDEF_INT -2123456789
#define IS_UNDEF(X) ( (X < -1E+34) || IsNaN(X) )
#define IS_UNDEF_INT(X) (X==UNDEF_INT)

#ifndef EXTRAPRECISION
#define PiNumber 3.141592653589793238462643383279
#define Pi2      6.283185307179586476925286766558

#define C1_12    0.083333333333333333333333333333
#define C1_6     0.166666666666666666666666666666
#define C1_3     0.333333333333333333333333333333
#define C1_30    0.033333333333333333333333333333
#define C2_3     0.666666666666666666666666666666
#define C4_3     1.333333333333333333333333333333
#define C5_3     1.666666666666666666666666666666
#define C5_6     0.833333333333333333333333333333
#define C11_60   0.183333333333333333333333333333
#define C13_60   0.216666666666666666666666666666
#define C47_60   0.783333333333333333333333333333
#define CLN2     0.693147180559945309417232121458
#endif
//----------------------------------------------------------------------------------------------------------------------

// Nanotechnology
inline int IsNaN(NativeDouble x) {
    unsigned char* X = (unsigned char*) &x;
    if(((X[6] | 0x0F) & X[7] & 0x7F) == 0x7F) return 1;
    return !(((X[6]^0xF0)&0xF7) || ((X[7]^0x7F)&0x7F));
}

template<typename T> 
inline void MakeNaN(T& x) {
    unsigned char* X = (unsigned char*) &x;
    for(unsigned int i=0; i<sizeof(T); i++) X[i]=0xFF;
}

// Минимальное число, представимое в заданной системе и отличное от нуля (с небольшим запасом)
template <typename fpv> inline fpv get_min_value() { return 1e-300; }
template <> inline float get_min_value<float>() { return 1e-37f; }

template <typename fpv> inline fpv get_eps();
template<> inline NativeDouble get_eps() { return 2.2204460492503131e-016; } /* smallest such that 1.0+DBL_EPSILON != 1.0 */

// Числа Pi, 2*Pi и ln(2)
template<typename fpv> inline fpv GetPiNumber();
template<typename fpv> inline fpv GetPiNumber2();
template<typename fpv> inline fpv GetLn2();
template<> inline NativeDouble GetPiNumber() { return 3.141592653589793238462643383279; }
template<> inline NativeDouble GetPiNumber2() { return 6.283185307179586476925286766558; } // 2*Pi
template<> inline NativeDouble GetLn2() { return 0.693147180559945309417232121458; } // ln(2)

//to silence unused warnings - only for virt functions with given interface!
inline void UnuseIt(double x){x = 0.0;}
inline void UnuseIt(char* x){x = 0;}
inline void UnuseIt(int x){x = 0;}

// calculating a^n where 'a' is a double and 'n' is a non-negative const int
template<int n> inline double pow_to_const_int(double a) { return pow_to_const_int<n-1>(a) * a; }
template<> inline double pow_to_const_int<0>(double) { return 1.0; }

#endif
