//======================================================================================================================
// Parser for user input 
//======================================================================================================================
#pragma once
#ifndef BASE_PARSER_HEADER
#define BASE_PARSER_HEADER

#include "base_param.h"
#include <map>

//=====================================================================================================================
// Parser input processing
//=====================================================================================================================

#define GT 1    // ">"
#define GE 2    // ">="
#define LT 3    // "<"
#define LE 4    // "<="
#define UNLIM 5 // unlimited

// Parser value types
enum tParType{
    PARTYPE_UNDEFTYPE  = 0,
    PARTYPE_INT        = 1,  // integer
    PARTYPE_DOUBLE     = 2,  // double precision floating point
    PARTYPE_DOUBLE3    = 3,  // 3 double precision floating point values, space separated
    PARTYPE_WORD       = 4,  // sequence of characters, limited by space character
    PARTYPE_STRING     = 5,  // sequence of characters, limited by new line character
    PARTYPE_BOOL       = 6,  // bool flag of a form +FlagName or -FlagName
    PARTYPE_ENUM       = 7,  // text label of a enum code (uses corresponding dictionary)
    PARTYPE_INTLIST    = 8,  // list of ints till end of line or comment sign
    PARTYPE_WORDLIST   = 9,  // list of words till end of line or comment sign
    PARTYPE_DOUBLELIST = 10, // list of doubles till end of line or comment sign
    PARTYPE_STRINGLIST = 11  // list of sequence of characters till end of line or comment sign
};

//=====================================================================================================================
// Command line classes
//=====================================================================================================================

class tCmdLineArgUnnamed { // Class for command line argument (unnamed)
private:
    string Value; // строка для обработки
    bool Used;    // запрашивался ли аргумент
    inline void SetUsed() { Used = true; }
public:
    tCmdLineArgUnnamed(const string& tmpValue) : Value(tmpValue), Used(false) {}
    tCmdLineArgUnnamed(const tCmdLineArgUnnamed& tmp) : Value(tmp.Value), Used(tmp.Used) {}
    friend class tCmdLine;
};

class tCmdLineArgStandard : public tCmdLineArgUnnamed { // Class for command line argument (-par=value)
private:
    string Name;  // имя
public:
    tCmdLineArgStandard(const string& tmpName, const string& tmpValue) : tCmdLineArgUnnamed(tmpValue), Name(tmpName) {}
    tCmdLineArgStandard(const tCmdLineArgStandard& tmp) : tCmdLineArgUnnamed(tmp), Name(tmp.Name) {}
    friend class tCmdLine;
};

class tCmdLineArgFile : public tCmdLineArgStandard { // Class for command line argument (file.txt:par=value)
private:
    string Fname; // имя файла
public:
    tCmdLineArgFile(const string& tmpFname, const string& tmpName, const string& tmpValue) :
                    tCmdLineArgStandard(tmpName, tmpValue), Fname(tmpFname) {}
    tCmdLineArgFile(const tCmdLineArgFile& tmp) : tCmdLineArgStandard(tmp), Fname(tmp.Fname) {}
    friend class tCmdLine;
};

class tCmdLine {// Class which describes command line
private:
    vector<tCmdLineArgStandard> ListStandardArguments; // список обычных аргументов
    vector<tCmdLineArgUnnamed> ListUnnamedArguments;  // список старорежимных аргументов фиксированного порядка
    vector<tCmdLineArgFile> ListFileArguments; // список оверрайда файлового ввода
    // Добавление аргумента определенного типа (убрано под приват, чтоб не возникало ошибок по формату)
    void AddArg(char* param);              // для отдельного аргумента
    void AddStandardArgument(char* param); // обычные аргументы
    void AddFileArgument(char* param);     // оверрайд файлового ввода
    void AddUnnamedArgument(char* param);  // старорежимные безымянные фиксированные аргументы

public:
    tErrAccumulator ErrAccumulator;

    // Добавление аргументов
    void AddArguments(int argc, char** argv); // для всей ком строки (не забывать убирать экзешник при передаче!)
                                              // т.е. (argc-1,argv+1)
    void AddArgument(char* param);            // для отдельного аргумента

    // Получение аргументов
    string GetNextUnnamedArgument(); // следующий неиспользованный безымянник для личного пользования
    string GetStandardArgument(const string &parName); // обычный аргумент
    // параметр (для указанного файла или с пропущенным именем файла)
    string GetFileParam(const string &fileName, const string &paramName);
    // макроподстановки (для указанного файла или с пропущенным именем файла)
    vector<string> GetFileMacros(const string &fileName, bool globalMacro = false);

    // Печать аргументов
    void PrintAllArguments(void) const; // все аргументы
    void PrintUnusedArguments(void) const; // все неиспользованные аргументы
    void PrintUnusedUnnamedArguments(void) const; // неиспользованные неименованные аргументы
    void PrintUnusedStandardArguments(void) const; // неиспользованные стандартные аргументы
    void PrintUnusedFileArguments(const string &fileName) const; // неиспользованные аргументы для указанного файла

    // Удаление аргументов
    void ClearAllArguments(void);      // из всех списков
    void ClearStandardArguments(void); // обычные аргументы ком. строки
    void ClearFileArguments(void);     // оверрайд файлового ввода
    void ClearUnnamedArguments(void);  // старорежимные безымянные фиксированные аргументы
};

extern tCmdLine CmdLine;

//======================================================================================================================
// Classes for initialization of parameters
//======================================================================================================================

//---------------------------------------------------------------------------------------------------------------------
// Class which describes the limitation of the parameter
//---------------------------------------------------------------------------------------------------------------------
class tLimitation{
private:
    double lowLimit,upLimit;
    int symLowLim,symUpLim;
    tParType parType;
public:
    tLimitation():lowLimit(0.0),upLimit(0.0),symLowLim(UNLIM),symUpLim(UNLIM),parType(PARTYPE_UNDEFTYPE){}
    tLimitation(const tLimitation &aLimitation){
        symLowLim = aLimitation.symLowLim;
        lowLimit  = aLimitation.lowLimit;
        symUpLim  = aLimitation.symUpLim;
        upLimit   = aLimitation.upLimit;
        parType = PARTYPE_UNDEFTYPE;
    }
    void operator=(const tLimitation &aLimitation);
    void SetLimits(const int _symLowLim,const double _lowLim,const int _symUpLim,const double _upLim);
    void SetLimits(const int _symLowLim,const int _lowLim,const int _symUpLim,const int _upLim);
    bool CheckValue(const int value) const;
    bool CheckValue(const double value) const;
};

class tParameter{ // Class which describes the parameter
private:
    tLimitation limit;
    int *ptrInt;        // Pointer to the int parameter
    double *ptrDouble;  // Pointer to the double parameter
    string *ptrString;  // Pointer to the string parameter
    vector<string> *ptrStringArray; // Pointer to the array of strings parameter
    vector<int> *ptrIntArray; // Pointer to the «intlist» parameter
    vector<double> *ptrDoubleArray; // Pointer to the «doublelist» parameter
    int crashIt; //flag to crash or not to crash if error happens
    string parName; // name of parameter
    tParType parType; // type of parameter
    bool isSet; // flag that parameter value has been set
    bool canRedef; // флаг запрещающий\разрешающий переопределение параметра - по умолчанию off
    string Dictionary; // stores text labels for enums

    inline void _init(){
        ptrInt=NULL; ptrDouble=NULL; ptrString=NULL; ptrStringArray = NULL; ptrIntArray=NULL;
        isSet = false; canRedef = false;  crashIt = 0;
        parType=PARTYPE_UNDEFTYPE;
        ptrDoubleArray = NULL;
    }
    void initParamBase(const string& _parName, const tParType _parType, const int _crashIt);
public:
    tParameter(){_init();}
    tParameter(const tParameter &aPar);
    tParameter(int &par,const string& _parName,const tParType _parType,const int _crashIt,const int _symLowLim,
               const int _lowLim,const int _symUpLim,const int _upLim);
    tParameter(double &par,const string& _parName,const tParType _parType,const int _crashIt,const int _symLowLim,
               const double _lowLim,const int _symUpLim,const double _upLim);
    tParameter(double* par,const string& _parName,const tParType _parType,const int _crashIt);
    tParameter(string& par,const string& _parName,const tParType _parType,const int _crashIt);
    tParameter(vector<string>& par,const string& _parName,const tParType _parType,const int _crashIt);
    tParameter(int &par,const string& _parName,const tParType _parType,const int _crashIt, const char *dict);
    tParameter(vector<int> &par,const string& _parName,const tParType _parType,const int _crashIt);
    tParameter(vector<double> &par,const string& _parName,const tParType _parType,const int _crashIt);
    inline bool CheckValue(const int value) const {return limit.CheckValue(value);}
    inline bool CheckValue(const double value) const {return limit.CheckValue(value);}
    inline bool GetIsSet(void) const { return isSet; }
    inline bool GetCrashIt(void)     const { return crashIt > 0 ? true : false ; }
    inline bool GetAntiCrashIt(void) const { return crashIt == -1 ? true : false ; }
    inline bool GetAntiDontCrashIt(void) const { return crashIt == -2 ? true : false ; }
    inline const char* GetParName(void) const { return parName.c_str(); }
    inline tParType GetParType(void) const { return parType; }
    inline void SetRedef(const bool redef) { canRedef=redef; }
    bool GetValue(string value, bool log, tErrAccumulator* ErrAccumulator = NULL);
    void PrintParameter() const;
    bool CheckPtr(const string& par) const;
    static bool CheckParamName(const string& _parName);
};

// Reads list of pairs string and int from file buffer FB to list
void parse_pairs_list(tFileBuffer& FB, vector< pair<string, int> >& list, tErrAccumulator* ErrAccumulator = NULL);

// Reads list of strings from file buffer FB to list
void parse_words_list(tFileBuffer& FB, vector<string>& list);

//---------------------------------------------------------------------------------------------------------------------
// Main parser class which manages initialization of list of parameters
//---------------------------------------------------------------------------------------------------------------------
class tParamManager: public std::vector<tParameter> {
public:
    typedef std::vector<tParameter> tParamsVectorBase; //list of parameters
    tErrAccumulator ErrAccumulator;

    tParamManager(){};

    // Получение пути по умолчанию, ранее запрошенного в переменную path
    // Вынесено в функцию для отдельного вызова до обработки остального ввода
    // (до инициализации каталога по умолчанию некуда писать лог,
    // а при повторном разборе оригинальной CMD будут некорректно обрабатываться безымянные параметры)
    void GetUnnamedPathFromCommandLine(std::string &path);
    void AddParameter(const tParameter &newParameter); //добавляет новый параметр
    void RequestParameter(int &par,const string& _parName,const tParType _parType=PARTYPE_INT,
                          const int _crashIt=IO_CRASH, 
                          const int _symLowLim=UNLIM, const int _lowLim=0,
                          const int _symUpLim =UNLIM, const int _upLim=0);
    template<typename T>
    inline void RequestParameter(T &par,const string& _parName, const char *dictionary, const int _crashIt=IO_CRASH){
        AddParameter(tParameter((int&) par,_parName,PARTYPE_ENUM,_crashIt, dictionary ));
    }
    inline void RequestParameter(double* par, const string& _parName, const tParType _parType=PARTYPE_DOUBLE3,
                                 const int _crashIt=IO_CRASH){
        AddParameter(tParameter(par,_parName,_parType,_crashIt));
    }
    inline void RequestParameter(vector<int>& par, const string& _parName, const tParType _parType=PARTYPE_INTLIST,
                                 const int _crashIt=IO_CRASH){
        AddParameter(tParameter(par,_parName,_parType,_crashIt));
    }
    inline void RequestParameter(vector<double>& par, const string& _parName,
                                 const tParType _parType=PARTYPE_DOUBLELIST, const int _crashIt=IO_CRASH){
        AddParameter(tParameter(par,_parName,_parType,_crashIt));
    }
    inline void RequestParameter(double &par, const string& _parName, const tParType _parType=PARTYPE_DOUBLE,
                                 const int _crashIt=IO_CRASH, 
                                 const int _symLowLim=UNLIM, const double _lowLim=0.0, 
                                 const int _symUpLim =UNLIM, const double _upLim=0.0){
        AddParameter(tParameter(par,_parName,_parType,_crashIt,_symLowLim,_lowLim,_symUpLim,_upLim));
    }
    inline void RequestParameter(string& par, const string& _parName,
                                 const tParType _parType=PARTYPE_WORD, const int _crashIt=IO_CRASH){
        AddParameter(tParameter(par,_parName,_parType,_crashIt));
    }
    inline void RequestParameter(vector<string>& par,const string& _parName,
                                 const tParType _parType=PARTYPE_WORDLIST, const int _crashIt=IO_CRASH){
        AddParameter(tParameter(par,_parName,_parType,_crashIt));
    }

    // Synonyms for optional parameters
    template<typename T> void Request(T &par, const string& _parName);

    template<typename T>
    inline void Request(T &par, const string& _parName, const char *dictionary){
        AddParameter(tParameter((int&) par,_parName,PARTYPE_ENUM,IO_DONTCRASH, dictionary ));
    }

    inline void RequestOption(int &par, const string& _parName){
        RequestParameter(par, _parName, PARTYPE_BOOL, IO_DONTCRASH);
    }

    void ReadParamsFromFile(const string& fileName, bool log = true, bool CMDLine = true);
    void ReadParamsFromBuffer(tFileBuffer &FB, bool log = true, bool CMDLine = true);
    void ReadParamsFromCommandLine(bool log = true);
    int WarnAboutAllAbsentParameters();
    tParameter& operator[] (const char* parname);

    // Extracts parameters pairs <parName,parValue> to map<string, string> (casts parName to lowercase)
    static void ExtractParametersPairs(tFileBuffer &fb, map<string, string>& par);

    void SetRedef(const bool redef); // запретить\разрешить переопределение параметров

private:
    // печать ошибок\предупреждений\информации об обязательных\необязательных\запрещенных параметрах
    void StatusInformation(bool log);
};

inline void changeParValue(map<string, string>& parPair, string parname, string newparvalue) {
    tolowercase(parname);
    map<string, string>::iterator i = parPair.find(parname);
    if(i != parPair.end()) i->second = newparvalue;
}

inline bool CheckIfIntegerParamIsPresent(tFileBuffer &FB, const char *pname){
    int nn;
    tParamManager PM;
    PM.RequestParameter(nn, pname, PARTYPE_INT, IO_DONTCRASH);
    PM.ReadParamsFromBuffer(FB, false/*no log*/,false/*no cmd line*/);
    return PM[pname].GetIsSet(); 
}


//=====================================================================================================================
// Formula parser
//=====================================================================================================================

// internal data structure
struct s_action {  // computing of string expressions
    int cmd; // action code
    int arg1, arg2, arg3; // action parameters
    double number; // stores value if needed
    s_action():cmd(0),arg1(0),arg2(0),arg3(0),number(0.0){}
};

typedef double(*tFuncDouble1Arg)(double);

struct tFormulaBuiltinArgs{ // обратная совместимость со всякими дефолтными аргументами
    int index; // номер узла "index"
    int label; // номер граничной метки "label"
    const double *coo; // координаты узла - "x","y","z"
    double dist; // расстояние до стенки - "dist"
    inline tFormulaBuiltinArgs(){index=-1; label=-1; coo=NULL; dist=0.0;}
};

// crooked freaky formula parser
class tFormula : public vector<s_action> {
private:
    int ntr; // number of threads it can work with
    vector< pair<string, tFuncDouble1Arg> > funcs; // list of built-in functions of 1 argument
    vector< pair<string, double> > consts; // list of constants
    // list of registered arguments
    vector<string> argnames; // names
    vector<int>    argused;  // flags that argument appeared in formula
    // vectors of values and pointers. Outer vector - for threads
    vector< vector<double> >  argvals;  // values, if given
    vector< vector<const double*> > argptrs;  // pointers to values, if given
    
    int FindArg(const char *n)const;
    int FindConst(const char *n)const;
    int FindFunc(const char *n)const;

    // parsing of transformed expression
    int Parse2(vector< pair<int,int> > &P, int Imin, int Imax, tErrAccumulator* ErrAccumulator);
    // auxiliary subroutine
    void Substitute(vector< pair<int,int> > &P, int Open, int Close, int& Imax, unsigned int newcode);
    int BinaryOperations(vector< pair<int,int> > &P, int Imin, int& Imax,
        int op_min, int op_max, int op_min2 = -1, int op_max2 = -1, bool unMinusFlag = false,
        tErrAccumulator* ErrAccumulator = NULL); // auxiliary subroutine

public:
    inline vector<const double*>& _argptrs(int trn = 0) { return argptrs[trn]; }
    inline const vector<const double*>& _argptrs(int trn = 0) const { return argptrs[trn]; }
    inline vector<double>& _argvals(int trn = 0) { return argvals[trn]; }
    inline const vector<double>& _argvals(int trn = 0) const { return argvals[trn]; }

    // указываем в конструкторе число нитей, на которых формула будет вычисляться одновременно
    // если вычисления проводится с явной передачей вектора значений или указателей,
    // то указывать число нитей здесь не нужно
    tFormula(int _ntr = 1);
    tFormula(const tFormula &f);
    void clear(); // reset to defaults
    void SetNumThreads(int _ntr);

    // add/set - returns arg number on success, -1 on fail
    int AddArg(const char *n, double v = 0.0, double* p = NULL, int DoCheck = 1); // add argument with given value
    int AddConst(const char *n, double v); // add constant with given value
    int SetArg(const char *n, double v, int trn = -1); // set/add argument with given value
    int SetArgPtr(const char *n, double *p, int trn = -1); // set/add argument with given pointer so value can be changed by user
    void AddDefaultVars();

    inline bool CheckArgUsed(int i) const { return argused[i]>0; } // check if argument appeared in formula
    bool CheckArgUsed(const char *n) const; // check if argument appeared in formula
    inline bool CheckCoorsUsed()const{ return CheckArgUsed("x") ||CheckArgUsed("y") || CheckArgUsed("z");}

    // parsing string expression
    int Parse(const char* str, bool parseLog = false, tErrAccumulator* ErrAccumulator = NULL);

    // evaluation of expression - general version
    double Evaluate(const vector<const double*>& ptrs, const vector<double>& vals) const; // main definition
    double Evaluate(int trn = 0) const; // версия, автоматически подставляющая внутренние массивы

    //specialized backward-compatibility crap. Variables index, label, x, y, z should be registered, otherwise - crash!
    double Evaluate(tFormulaBuiltinArgs args, int trn = -1); // evaluation of expression for a given vector
    bool Test(tFormulaBuiltinArgs args, int trn = -1); // int(Evaluate)
};

#endif
