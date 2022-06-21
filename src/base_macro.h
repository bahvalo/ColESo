//======================================================================================================================
//Basic macro definitions - lightweight header
//======================================================================================================================
#pragma once
#ifndef BASE_MACRO_HEADER
#define BASE_MACRO_HEADER

#include "base_config.h"

// Макрос для определения, принадлежит ли ячейка локальной подобласти - true (или является гало - false)
#define IS_OWN_CELL(i) (i>=0 && i<Num_Cell_Own) // временно тут
#define IS_OWN_SEG(cells) (IS_OWN_CELL(cells[0])|| IS_OWN_CELL(cells[1]))
#define IS_OWN_BCNODE(i) (i>=0 && i<Num_BndNode_Own) // временно тут

// Undefs
#define UNDEF -1E+35
#define UNDEF_INT -2123456789
#define IS_UNDEF(X) ( (X < -1E+34) || IsNaN(X) )
#define IS_UNDEF_INT(X) (X==UNDEF_INT)

//separator lines to print to log
#define LINE_SEPARATOR  "===============================================================================\n"
#define LINE_SEPARATOR_ "-------------------------------------------------------------------------------\n"

#define PRINT_ULL(X) ((unsigned long long)X) // use %llu to print
#define _PRINT_MACRO(X) #X
#define PRINT_MACRO(x) _PRINT_MACRO(x) // turns macro definition into char string 

#define MAX_UINT64             (0xEFFFFFFFFFFFFFFFull)
#define MAX_ALLOWED_BLOCK_SIZE (0xEFFFFFFFl) // 4026531839(3.8Gb)
#define MAX_SIGNED_INT         (0x7FFFFFFF ) // 2147483647

#define MIN(a,b) ( ((a)<(b)) ? (a) : (b))
#define MAX(a,b) ( ((a)<(b)) ? (b) : (a))
#define SQR(X) ((X)*(X))
#define SIGN(a) ( ((a)<0.0) ? -1.0 : (((a)>0.0) ? 1.0 : 0.0) )
#define VDOT(a,b) ((a)[0]*(b)[0] + (a)[1]*(b)[1] + (a)[2]*(b)[2] )
#define VELDOT(a,b) (a[Var_U]*b[Var_U] + a[Var_V]*b[Var_V] + a[Var_W]*b[Var_W] ) // 5
#define DIST(v1, v2) (sqrt(SQR(v1[0]-v2[0])+SQR(v1[1]-v2[1])+SQR(v1[2]-v2[2])))

#endif
