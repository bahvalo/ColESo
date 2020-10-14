#pragma once
#ifndef PARSER_H
#define PARSER_H

#include "personal.h"
#include <map>
#include <vector>
#include <string>
using namespace std;

#ifndef IO_DONTCRASH
#define IO_DONTCRASH 0  // dont crash if cant open file
#endif
#ifndef IO_CRASH
#define IO_CRASH     1  // crash if cant open file
#endif

//======================================================================================================================
// Parser error counter: simple error accumulator
//======================================================================================================================
enum tErrQuiet {
    ERR_CRASH = -1, // immediate crash if a error is occured
    ERR_STORE = 0,  // store error messages and call crash in the destructor
    ERR_QUIET = 1   // store error messages and ignore them in the destructor
};
struct tErrAccumulator {
    string ERRMSG; // error accumulator
    string ERRPREF; // prefics for errors
    int ERRCOUNT;// number of errors
    int ERRCHECK;// number of errors at the last check
    tErrQuiet ERRQUIET; // quiet mode: ignore all new errors messages

    tErrAccumulator() { ERRCOUNT=0; ERRCHECK=0; ERRQUIET=ERR_STORE; }
    inline void Clear(){ ERRMSG = ""; ERRCOUNT = 0; }

    // adds a error message
    inline void Add(string sbuf){ 
        ERRMSG += ERRPREF + sbuf; 
        ERRCOUNT++;
        if(ERRQUIET==ERR_CRASH) Check_crash();
    }
    // adds all messages from another error accumulator
    inline void Add(const tErrAccumulator& ErrAccumulator) { 
        if(!ErrAccumulator.ERRCOUNT) return;
        ERRMSG += ErrAccumulator.ERRMSG; 
        ERRCOUNT += ErrAccumulator.ERRCOUNT;
        if(ERRQUIET==ERR_CRASH) Check_crash();
    }
    inline void MoveTo(tErrAccumulator* ErrAccumulator) { 
        if(ErrAccumulator!=NULL && ERRCOUNT>0) {
            ErrAccumulator->ERRMSG += ERRMSG;
            ErrAccumulator->ERRCOUNT += ERRCOUNT; 
            if(ErrAccumulator->ERRQUIET==ERR_CRASH) ErrAccumulator->Check_crash();
        }
        Clear();
    }
    // check, if more errors added since last check
    inline bool Check(){ bool res = (ERRCOUNT>ERRCHECK); ERRCHECK=ERRCOUNT; return res; }
    // call crash(...) with collected messages if there are errors
    void Check_crash() const;

    // Stupid getters & setters
    inline int Count() const{ return ERRCOUNT; } // get number of errors
    inline void IO_error_quiet_ON() { ERRQUIET = ERR_QUIET; } // set quiet mode on
    inline void IO_error_quiet_OFF(){ ERRQUIET = ERR_STORE; } // set quiet mode off
    inline void IO_error_set_pref(string sbuf) { ERRPREF = sbuf + " "; } // set prefix
    inline void IO_error_remove_pref() { ERRPREF = ""; } // remove prefix

    // Automatically call crash() when destroying an instance of tErrAccumulator
    ~tErrAccumulator() { if(ERRQUIET!=ERR_QUIET) Check_crash(); }
};
// accumulation of error messages in input files etc.
void IO_error_add_str(tErrAccumulator* ErrAccumulator, string sbuf); //adds error
void IO_error_add_str(tErrAccumulator& ErrAccumulator, string sbuf); //adds error
// wrapper. Object or pointer ErrAccumulator should be defined in the calling function
#define IO_error_add(...) IO_error_add_str(ErrAccumulator, stringprintf(__VA_ARGS__)+"\n")

#define GT 1    // ">"
#define GE 2    // ">="
#define LT 3    // "<"
#define LE 4    // "<="
#define UNLIM 5 // unlimited

enum tParType{
    PARTYPE_UNDEFTYPE  = 0,
    PARTYPE_INT        = 1, // integer
    PARTYPE_DOUBLE     = 2, // double precision floating point
    PARTYPE_DOUBLE3    = 3, // 3 double precision floating point values, space separated 
    PARTYPE_WORD       = 4, // sequence of characters, limited by space character
    PARTYPE_STRING     = 5, // sequence of characters, limited by new line character
    PARTYPE_BOOL       = 6, // bool flag of a form +FlagName or -FlagName
    PARTYPE_ENUM       = 7, // text label of a enum code (uses corresponding dictionary)
    PARTYPE_INTLIST    = 8, // list of ints till end of line or comment sign
    PARTYPE_WORDLIST   = 9, // list of words till end of line of comment sign
    PARTYPE_DOUBLELIST = 10 // list of doubles till end of line or comment sign
};
// Parser interfaces
void GetNextLine(string &line, FILE* pf, bool fullLine = false);
bool GetNextWord(string &word, string &pline, int lowerCase=1);
bool GetNextElem(string &word, string &pline, int lowerCase = 1);

string& trim_string(string& s, const char *ws = 0);
string& trim_right_string(string& s, const char *ws = 0);
void tolowercase(string& str);

// Returns true if word1 and word2 are equal (literally, in lowercase)
inline bool CompareWords(const string& word1, const string& word2){
    if(word1.size() != word2.size()) return false;
    if(word1.empty() && word2.empty()) return true;
    for (size_t c = 0; c < word1.size(); c++) {
        if (tolower(word1.at(c)) != tolower(word2.at(c))) return false;
    }
    return true;
}

int IsInteger(const char* word);
int IsNumber(const char* word);
int GetIntFromWord(const char* word, tErrAccumulator* ErrAccumulator = NULL);
int ParseIntFromWord(const char* word, tErrAccumulator* ErrAccumulator = NULL);
double GetDoubleFromWord(const char* word, tErrAccumulator* ErrAccumulator = NULL);
double ParseDoubleFromWord(const char* word, tErrAccumulator* ErrAccumulator = NULL);

// Gets integer parameter <par> from <line>
// Returns: true if success (and clears line), false otherwise
bool GetIntFromLine(int &par, string &line, tErrAccumulator* ErrAccumulator = NULL);
bool GetIntIntFromLine(int &par1,int &par2,string &line, tErrAccumulator* ErrAccumulator = NULL);
bool GetIntDoubleFromLine(int &par1,double &par2,string &line, tErrAccumulator* ErrAccumulator = NULL);
int GetDoubleArrayFromLine(double *par,int maxSize,string &line, tErrAccumulator* ErrAccumulator = NULL);
int GetDoubleArrayFromLine(vector<double>& par, int maxSize,string &line, tErrAccumulator* ErrAccumulator = NULL);

string GetEnumName(int val, const char *Dictionary); //finds enum option name for a given int value
int CheckEnumValue(const int value, string Dictionary); //check numerical code for a given Dictionary
int GetEnumValue(const string& word, string Dictionary); //returns numerical code for a given word, UNDEF_INT if error
int ProcessEnum(const string& word, const string& Dictionary); //same as GetEnumValue, admits numerical and literal codes

// Formatted output to std::string conversion
string stringprintf(const char* fmt, ...);

inline int round_noisette(double value) { return int(value + 0.5*(value>0.0 ? 1.0 : -1.0)); }

//======================================================================================================================
// type MapString2String and functions
//======================================================================================================================

typedef map<string, string> MapString2String;

inline void changeParValue(MapString2String& parPair, string parname, string newparvalue) {
    tolowercase(parname);
    MapString2String::iterator i = parPair.find(parname);
    if(i != parPair.end()) i->second = newparvalue;
}

//======================================================================================================================
// Command line classes
//======================================================================================================================

class tCmdLineArgUnnamed { // Class for command line argument (unnamed)
    string Value; // строка для обработки
    bool Used;    // запрашивался ли аргумент

public:
    tCmdLineArgUnnamed(const string& tmpValue) : Value(tmpValue), Used(false) {}
    tCmdLineArgUnnamed(const tCmdLineArgUnnamed& tmp) : Value(tmp.Value), Used(tmp.Used) {}

    friend class tCmdLine;

private:
    inline void SetUsed() { Used = true; }
};

class tCmdLineArgStandard : public tCmdLineArgUnnamed { // Class for command line argument (-par=value)
    string Name;  // имя

public:
    tCmdLineArgStandard(const string& tmpName, const string& tmpValue) : tCmdLineArgUnnamed(tmpValue), Name(tmpName) {}
    tCmdLineArgStandard(const tCmdLineArgStandard& tmp) : tCmdLineArgUnnamed(tmp), Name(tmp.Name) {}

    friend class tCmdLine;
};

class tCmdLineArgFile : public tCmdLineArgStandard { // Class for command line argument (file.txt:par=value)
    string Fname; // имя файла

public:
    tCmdLineArgFile(const string& tmpFname, const string& tmpName, const string& tmpValue) :
        tCmdLineArgStandard(tmpName, tmpValue), Fname(tmpFname) {}
    tCmdLineArgFile(const tCmdLineArgFile& tmp) : tCmdLineArgStandard(tmp), Fname(tmp.Fname) {}

    friend class tCmdLine;
};

class tCmdLine {// Class which describes command line
    vector<tCmdLineArgStandard> ListStandardArguments; // список обычных аргументов
    vector<tCmdLineArgUnnamed> ListUnnamedArguments;  // список старорежимных аргументов фиксированного порядка 
    vector<tCmdLineArgFile> ListFileArguments;     // список оверрайда файлового ввода

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
    void PrintAllArguments(void) const;                          // все аргументы
    void PrintUnusedArguments(void) const;                       // все неиспользованные аргументы
    void PrintUnusedUnnamedArguments(void) const;                // неиспользованные неименнованные аргументы
    void PrintUnusedStandardArguments(void) const;               // неиспользованные стандартные аргументы
    void PrintUnusedFileArguments(const string &fileName) const; // неиспользованные аргументы для указанного файла
    
    // Удаление аргументов
    void ClearAllArguments(void);      // из всех списков 
    void ClearStandardArguments(void); // обычные аргументы ком. строки
    void ClearFileArguments(void);     // оверрайд файлового ввода
    void ClearUnnamedArguments(void);  // старорежимные безымянные фиксированные аргументы 

private: // Добавление аргумента определенного типа (убрано под приват, чтоб не возникало ошибок по формату) 
    void AddArg(char* param);               // для отдельного аргумента
    void AddStandardArgument(char* param); // обычные аргументы
    void AddFileArgument(char* param);     // оверрайд файлового ввода
    void AddUnnamedArgument(char* param);  // старорежимные безымянные фиксированные аргументы 
};

extern tCmdLine CmdLine;

//======================================================================================================================
// Macro substitutions classes
//======================================================================================================================
#define MACROTXT "macro.txt"

class tMacroSubstitution { //macro substitution for parser
    string Name; //word to be replaced
    string Value; //value to replace by

private:
    tMacroSubstitution(void) {};
    tMacroSubstitution(const string &tName, const string &tValue) { Name = tName; Value = tValue; };
    //try to create new macro substitution(return success)
    bool InitMacroSubstitution(const string& line, bool dollars, tErrAccumulator* ErrAccumulator);

    friend class tMacroSubstitutionVector;
};

class tMacroSubstitutionVector : public vector<tMacroSubstitution> { //macro substitutions vector for parser
public:
    tErrAccumulator ErrAccumulator;
    bool AddNew(const string& line, bool log, bool fromMacroFile); //try to create one new (return success)
    bool AddOrReplace(const string& line, bool log); //try to create one new or replace its value(return success)
    bool Evaluate(bool log); //try evaluate value of macro substitutions starting with EVAL(return success)
    int Substitution(string &line, bool log); //find value for expression in line(return error code, 0 - OK)

    static void InitGlobalMacro(string fname, bool CMDLine, bool log); //init global list of macro substitutions
    static tMacroSubstitutionVector& GlobalMacro(void); // return global list of macro substitutions
};

//======================================================================================================================
// Classes for initialization of parameters
//======================================================================================================================
class tLimitation{ // Class which describes the limitation of the parameter
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
    string Dictionary; // stores text labels for enums

    inline void _init(){
        ptrInt=NULL; ptrDouble=NULL; ptrString=NULL; ptrStringArray = NULL; ptrIntArray=NULL;
        isSet=false; crashIt=0;
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
    bool GetValue(string value, bool log, tErrAccumulator* ErrAccumulator = NULL);
    void PrintParameter() const;
    inline bool GetIsSet(void) const { return isSet; }
    inline bool GetCrashIt(void)     const { return crashIt > 0 ? true : false ; }
    inline bool GetAntiCrashIt(void) const { return crashIt == -1 ? true : false ; }
    inline bool GetAntiDontCrashIt(void) const { return crashIt == -2 ? true : false ; }
    inline const char* GetParName(void) const { return parName.c_str(); }
    inline tParType GetParType(void) const { return parType; }
    bool CheckPtr(const string& par) const;

    static bool CheckParamName(const string& _parName);
};

class tFileBuffer {// Class for file buffer 
private:
    bool* Used;                              // flags for lines
    int FirstLine;                           // number of the first line in current zone
    int EnfOfZoneLine;                       // number of the first line after current zone 
    int NextZoneLine;                        // number of the line with the begining of next zone after keyword
    int CurrentLine;                         // last returned line number
    unsigned int Sum;                        // check sum
    string FNameFull;                        // file name
    vector<string> Lines;                    // array of lines
    tMacroSubstitutionVector ListMacroLocal; // macro list for current buffer

public:
    tFileBuffer() : Used(NULL) { Clear(); }
    tFileBuffer(const string& fName, int CrashIt=IO_CRASH, bool log=true) : Used(NULL)
        { Clear(); LoadFile(fName, CrashIt, log); }
    ~tFileBuffer() { Clear(); }

    // общие функции
    inline bool Empty() const { return Lines.empty(); } // проверяет заполнен ли буфер
    string GetFnameShort(void) const; // возвращает только имя файла без пути
    void GetNextLine(string &line); // возвращает следующую строку (пустые и комментарии не считаются)
    inline unsigned int CheckSum(void) { return Sum; } // возвращает сумму всех символов (проверка изменений update)
    inline bool EndOfZone(void) { if(CurrentLine == EnfOfZoneLine || !EnfOfZoneLine) return true; else return false; }
    inline tMacroSubstitutionVector& GetMacroList(void) { return ListMacroLocal; } // возвращает список макроподстановок

    // заполнение буфера данными и очистка
    void AllocText(const string& text); // создает буффер с текстом
    void AllocText(const char* text, size_t length); // создает буффер с текстом
    void AllocText(FILE* file); // создает буффер с текстом
    void Clear(); // удаление внутренних данных
    void CopyActualZoneTo(tFileBuffer& FB) const; // Copies current zone (part of file) to file buffer FB
    int LoadFile(const string& fName, int CrashIt = IO_CRASH, bool log = true); // Returns 1 if file cannot be opened 

    // навигация по зонам буфера
    inline int MakeZone(const string& KeyWord){ return MakeZone(KeyWord.c_str()); } // устанавливает следующую зону
    int MakeZone(const char *KeyWord); // устанавливает начало следующей зоны (Returns 1 if EOF)
    void GotoZone(void); // устанавливает указатель для чтения на начало текущей зоны
    void GotoNextZone(void); // делает текущей следующую зону и устанавливает указатель для чтения на ее начало

    // маркировка использованных строк
    void AllocUsed(void); // создает массив флагов об использовании строк (если его еще нет)
    void MarkPrevLineAsUsed(void); // ставит флаг об использовании предыдущей строке
    int PrintUnknownKeywords(void); // выводит незакомментированные строки в текущей зоне

    tErrAccumulator ErrAccumulator;

private:
    // Переводит файл с bool в формате (+\-)name в формат name 0(1)
    void PatchOldToNewConvert(FILE* pFILE, string &file, const string& fName, bool log);
    void InitMacro(bool log);
};

// Reads list of strings from file buffer FB to list
void parse_words_list(tFileBuffer& FB, vector<string>& list);

// Main parser class which manages initialization of list of parameters
class tParamManager:public std::vector<tParameter> {
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
public: 
    void RequestParameter(int &par,const string& _parName,const tParType _parType=PARTYPE_INT,const int _crashIt=IO_CRASH,
        const int _symLowLim=UNLIM,const int _lowLim=0,const int _symUpLim=UNLIM,const int _upLim=0);
    template<typename T>
    void RequestParameter(T &par,const string& _parName, const char *dictionary, const int _crashIt=IO_CRASH){
        AddParameter(tParameter((int&) par,_parName,PARTYPE_ENUM,_crashIt, dictionary ));
    }
    void RequestParameter(double* par,const string& _parName,const tParType _parType=PARTYPE_DOUBLE3,const int _crashIt=IO_CRASH){
        AddParameter(tParameter(par,_parName,_parType,_crashIt));
    }
    void RequestParameter(vector<int>& par,const string& _parName,const tParType _parType=PARTYPE_INTLIST,const int _crashIt=IO_CRASH){
        AddParameter(tParameter(par,_parName,_parType,_crashIt));
    }
    void RequestParameter(vector<double>& par,const string& _parName,const tParType _parType=PARTYPE_DOUBLELIST,const int _crashIt=IO_CRASH){
        AddParameter(tParameter(par,_parName,_parType,_crashIt));
    }
    void RequestParameter(double &par,const string& _parName,const tParType _parType=PARTYPE_DOUBLE,const int _crashIt=IO_CRASH,
        const int _symLowLim=UNLIM,const double _lowLim=0.0,const int _symUpLim=UNLIM,const double _upLim=0.0){
        AddParameter(tParameter(par,_parName,_parType,_crashIt,_symLowLim,_lowLim,_symUpLim,_upLim));
    }
    void RequestParameter(string& par,const string& _parName,const tParType _parType=PARTYPE_WORD,const int _crashIt=IO_CRASH){
        AddParameter(tParameter(par,_parName,_parType,_crashIt));
    }
    void RequestParameter(vector<string>& par,const string& _parName,const tParType _parType=PARTYPE_WORDLIST,const int _crashIt=IO_CRASH){
        AddParameter(tParameter(par,_parName,_parType,_crashIt));
    }

    // Synonyms for optional parameters
    template<typename T> void Request(T &par, const string& _parName);
    void RequestOption(int &par, const string& _parName){
        RequestParameter(par, _parName, PARTYPE_BOOL, IO_DONTCRASH);
    }

    void ReadParamsFromFile(const string& fileName, bool log = true, bool CMDLine = true);
    void ReadParamsFromBuffer(tFileBuffer &FB, bool log = true, bool CMDLine = true);
    void ReadParamsFromCommandLine(bool log = true);
    int WarnAboutAllAbsentParameters();
    tParameter& operator[] (const char* parname);

    // Extracts parameters pairs <parName,parValue> to MapString2String
    //   (casts parName to lowercase)
    static void ExtractParametersPairs(tFileBuffer &fb, MapString2String& par);

private:
    // печать ошибок\предупреждений\информации об обязательных\необязательных\запрещенных параметрах
    void StatusInformation(bool log);
};

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

private:
    int FindArg(const char *n)const;
    int FindConst(const char *n)const;
    int FindFunc(const char *n)const;

public:
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
    double Evaluate(int index, int label, const double* vec, int trn = -1); // evaluation of expression for a given vector
    int Test(int index, int label, const double* vec, int trn = -1); // int(Evaluate)

private:
    // parsing of transformed expression
    int Parse2(vector< pair<int,int> > &P, int Imin, int Imax, tErrAccumulator* ErrAccumulator); 
    // auxiliary subroutine
    void Substitute(vector< pair<int,int> > &P, int Open, int Close, int& Imax, unsigned int newcode);
    int BinaryOperations(vector< pair<int,int> > &P, int Imin, int& Imax, 
        int op_min, int op_max, int op_min2 = -1, int op_max2 = -1, bool unMinusFlag = false,
        tErrAccumulator* ErrAccumulator = NULL); // auxiliary subroutine
};
//=====================================================================================================================


#endif
