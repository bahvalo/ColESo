//======================================================================================================================
// Distributed file IO for parallel execution
//======================================================================================================================
#pragma once
#ifndef BASE_IO_HEADER
#define BASE_IO_HEADER

#include "base_config.h"

#include <vector>
#include <string>
using namespace std;

//---------------------------------------------------------------------------------------------------------------------
// Interface for reading and writing to files
//---------------------------------------------------------------------------------------------------------------------

// Access modes 
#define IO_READ   0 // access on read
#define IO_WRITE  1 // access on write
#define IO_APPEND 2 // access on append

// File types
#define IO_BINARY 1 // access binary
#define IO_TEXT   0 // access text

// Error behaviour (can be predefined from outside)
#ifndef IO_DONTCRASH
    #define IO_DONTCRASH 0 // dont crash if cant open file 
#endif
#ifndef IO_CRASH
    #define IO_CRASH 1 // crash if cant open file
#endif
#define IO_ANTICRASH -1 // crash if can - для патчёвой затычки, чтобы не те параметры не сували в не те файлы
#define IO_ANTIDONTCRASH -2 // warning if can - для патчёвой затычки, чтобы не те параметры не сували в не те файлы

// Сodes of data types (do not touch!)
#define IO_CHAR   8872
#define IO_INT    1267
#define IO_FLOAT  8867
#define IO_DOUBLE 1272
#define IO_TYPE_NAMES "CHAR 8872  IO_INT 1267  FLOAT 8867  DOUBLE 1272"
#define IO_REAL_TYPE_NAMES "FLOAT 8867  DOUBLE 1272"

#endif
