//======================================================================================================================
// Basic processing of user input params 
//======================================================================================================================
#pragma once
#ifndef BASE_PARAMS_HEADER
#define BASE_PARAMS_HEADER

#include "base_lib.h"
#include "base_io.h"

//----------------------------------------------------------------------------------------------------------------------
// Input error counter: simple error accumulator
//----------------------------------------------------------------------------------------------------------------------

enum tErrQuiet {
    ERR_CRASH = -1, // immediate crash if a error is occurred
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

    // Getters & setters
    inline int Count() const{ return ERRCOUNT; } // get number of errors
    inline void IO_error_quiet_ON() { ERRQUIET = ERR_QUIET; } // set quiet mode on
    inline void IO_error_quiet_OFF(){ ERRQUIET = ERR_STORE; } // set quiet mode off
    inline void IO_error_set_pref(string sbuf) { ERRPREF = sbuf + " "; } // set prefix
    inline void IO_error_remove_pref() { ERRPREF = ""; } // remove prefix

    // Automatically call crash() when destroying an instance of tErrAccumulator
    ~tErrAccumulator() { if(ERRQUIET!=ERR_QUIET) Check_crash(); }
};

// accumulation of error messages in input files etc.
void IO_error_add_str(tErrAccumulator* ErrAccumulator, string sbuf);
void IO_error_add_str(tErrAccumulator& ErrAccumulator, string sbuf);

// wrapper. Object or pointer ErrAccumulator should be defined in the calling function
#define IO_error_add(...) IO_error_add_str(ErrAccumulator, stringprintf(__VA_ARGS__)+"\n")


//----------------------------------------------------------------------------------------------------------------------
// Input params processing
//----------------------------------------------------------------------------------------------------------------------

// Params processing interfaces
bool IsCommentLine(const char* pline);
void GetNextLine(string &line, FILE* pf, bool fullLine = false);
bool GetNextWord(string &word, string &pline, int lowerCase=1);
bool GetNextElem(string &word, string &pline, int lowerCase = 1);

// String processing
string& trim_string(string& s, const char *ws = 0);
string& trim_right_string(string& s, const char *ws = 0);
void tolowercase(string& str);
void split(const string& source, char delim, vector<string>& elems); // Split <source> to vector of strings <elems>
void ParseString(const string& S, vector<string>& VS, char comma); // the same operation (?)

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
string EnumToDefines(const char* Dictionary, const string prefix=""); //converts dictionary (enums) to defines (for OpenCL)

// Formatted output to std::string conversion
string stringprintf(const char* fmt, ...);

inline int round_noisette(double value){ return int(value + 0.5*(value>0.0 ? 1.0 : -1.0)); }


//---------------------------------------------------------------------------------------------------------------------
// Macro substitutions for parameters
//---------------------------------------------------------------------------------------------------------------------
#define MACROTXT "macro.txt"

class tMacroSubstitution { //macro substitution for parser
    string Name; //word to be replaced
    string Value; //value to replace by

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


//---------------------------------------------------------------------------------------------------------------------
// Class for file buffer
//---------------------------------------------------------------------------------------------------------------------
class tFileBuffer {
private:
    vector<bool> Used;   // flags for lines
    int FirstLine;       // number of the first line in current zone
    int EnfOfZoneLine;   // number of the first line after current zone
    int NextZoneLine;    // number of the line with the begining of next zone after keyword
    int CurrentLine;     // last returned line number
    unsigned int Sum;    // check sum
    string FNameFull;    // file name
    vector<string> Lines;// array of lines
    tMacroSubstitutionVector ListMacroLocal; // macro list for current buffer

    // Переводит файл с bool в формате (+\-)name в формат name 0(1)
    void PatchOldToNewConvert(FILE* pFILE, string &file, const string& fName, bool log);
    void InitMacro(bool log);

public:
    tFileBuffer() { Clear(); }
    tFileBuffer(const string& fName, int CrashIt=IO_CRASH, bool log=true){ 
        Clear(); 
        LoadFile(fName, CrashIt, log); 
    }
    ~tFileBuffer() { Clear(); }

    // общие функции
    inline bool Empty() const { return Lines.empty(); } // проверяет заполнен ли буфер
    inline unsigned int CheckSum(void) { return Sum; } // возвращает сумму всех символов (проверка изменений update)
    inline bool EndOfZone(void) { if(CurrentLine == EnfOfZoneLine || !EnfOfZoneLine) return true; else return false; }
    inline tMacroSubstitutionVector& GetMacroList(void) { return ListMacroLocal; } // возвращает список макроподстановок
    string GetFnameShort(void) const; // возвращает только имя файла без пути
    void GetNextLine(string &line); // возвращает следующую строку (пустые и комментарии не считаются)

    // заполнение буфера данными и очистка
    void AllocText(const string& text); // создает буфер с текстом
    void AllocText(const char* text, size_t length); // создает буфер с текстом
    void AllocText(FILE* file); // создает буфер с текстом
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
};


//----------------------------------------------------------------------------------------------------------------------
// Extandable enum: enum + dictionary with possible external names 
//----------------------------------------------------------------------------------------------------------------------

template <typename T /*base enum type*/ , 
          typename F /*what is registered (usually func pointer)*/ >
class tStupidEnum{
private:
#define STUPID_CODE -1234567 // some code that is supposed to never appear in enums
    int myEnum; // base enum 
    string dict; // base dictionary of text labels for enum int codes
    vector<F> extopts; // array of external options 
    inline tStupidEnum(){} // forbidden default constructor 
public:
    inline const char *GetDictionary(){ return dict.c_str(); }
    inline tStupidEnum(T def, const char *d){ myEnum = (int)def; dict = d; }
    inline int AddExternalOption(const string &s, F o){ 
        if(ProcessEnum(s, dict) != UNDEF_INT) return -1; // check that option name is not in dict
        int code = (int)(STUPID_CODE + (int)extopts.size()); // new code that is supposed to never appear in enums
        string scode = stringprintf("%d", code); 
        ASSERT(ProcessEnum(scode, dict) == UNDEF_INT, "enum code breakdown"); // check that code is not in dict
        dict = dict + " " + s + " " + scode; // adding new option with new nonsense code to dict
        extopts.push_back(o);
        return 0;
    }
    inline int SetExternalOption(const string &s){ // to reassign existing option to external module's option
        int i = ProcessEnum(s, dict);
        if(i == UNDEF_INT) return -1; 
        myEnum = i;
        return 0;
    }
    inline F GetExternalOption() const{ // gets external option if it was set, 0 otherwise 
        int i = (int)myEnum - STUPID_CODE;
        if(i>=0 && i<(int)extopts.size()) return extopts[i];
        return (F) 0;
    }
    inline string GetExternalOptionName() const{ // gets the name of an external option, empty string if not set
        int i = (int)myEnum - STUPID_CODE;
        if(i>=0 && i<(int)extopts.size()) return GetEnumName(myEnum, dict.c_str());
        return string("");
    }
    inline operator int&() { return myEnum; } // cast to int of base enum 
    inline operator T() const { return (T)myEnum; } // cast to base enum
    inline tStupidEnum<T,F>& operator=(const T &object){ myEnum = (int)object; return *this; } // assign base enum 
#undef STUPID_CODE
};

#endif
