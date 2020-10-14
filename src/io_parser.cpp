//======================================================================================================================
//   Parser functions
//======================================================================================================================
#include "parser.h"
#include "const.h"
#include <stdarg.h>
#include <algorithm>
#include <functional>
#include <cctype>
#include <string.h>
#include <math.h>

#define FP_FORMAT "%.15g"
#ifndef ASSERT
#include "lib_base.h" // crash, InformKGB
#include "io_base.h"  // IO_Open - встроенная обёртка от fopen
#endif
#ifdef EXTRAPRECISION_COLESO
#include "extraprecision.h"
#endif

tCmdLine CmdLine;                         // списки параметров командной строки 
static tMacroSubstitutionVector ListMacroGlobal; // список глобальных макроподстановок 

static const string PARAMOFF = "UNDEF";  // key word for param off
static const string PARAMTRUE = "TRUE";  // key word for boolean param`s value
static const string PARAMFALSE= "FALSE";  // key word for boolean param`s value

struct tErrAccumulatorIgnore : public tErrAccumulator {
    tErrAccumulatorIgnore() : tErrAccumulator() { ERRQUIET = ERR_QUIET; }
};
static tErrAccumulatorIgnore EA_IGNORE;

//======================================================================================================================
// Simple error accumulator
//======================================================================================================================

// Crash if errors were found 
void tErrAccumulator::Check_crash() const{ 
    if(ERRCOUNT) crash("%s", ERRMSG.c_str());
}

void IO_error_add_str(tErrAccumulator* ErrAccumulator, string sbuf) {
    if(ErrAccumulator!=NULL) ErrAccumulator->Add(sbuf);
    else crash("%s", sbuf.c_str());
}
void IO_error_add_str(tErrAccumulator& ErrAccumulator, string sbuf) {
    ErrAccumulator.Add(sbuf);
}

//======================================================================================================================
// Basic string processing
//======================================================================================================================

// Check whether a line is empty but a comment
bool isComment(const char* pline){
    int i=0;
    if(pline[0]=='\0') return false;
    while(pline[i]==' '||pline[i]=='\t'||pline[i]==0xA||pline[i]==0xD) i++;
    if(pline[i]=='#' || pline[i]==0) return true;
    else return false;
}

 
// Get next line from file <pf>
// fullLine=0 (default): skip empty or comment lines; fullLine=1: return the line even it is a comment
void GetNextLine(string &sbuf,FILE* pf,bool fullLine){ 
    char cbuf[128]; 
    sbuf.clear();
    do {
        sbuf.clear();
        int eol=0;
        do{
            if(fgets(cbuf,100,pf)==NULL) cbuf[0]='\0';
            sbuf = sbuf + string(cbuf);
            int n = (int)strlen(cbuf);
            if(n==0) break;
            if(cbuf[n-1]==0xA/*LF*/|| cbuf[n-1]==0xD/*CR*/) eol=1;
        }
        while(!eol);

    }while(!fullLine && isComment(sbuf.c_str()));

    if(!fullLine) {
        const char* ppline = sbuf.c_str();
        int i = 0;
        while(ppline[i] >= ' ' || ppline[i] == '\t') i++;
        sbuf = sbuf.substr(0, i);
    }
}
 
// Check whether spaces are present in string
static inline bool CheckSpaces(const char *cbuf){
    int n = (int)strlen(cbuf);
    for(int i=0; i<n; i++) if(cbuf[i]==' ' || cbuf[i]=='\t') return true;
    return false;
}

// Remove spaces at the left and at the right of a string
string& trim_string(string& s, const char *ws){
    if( ws && ws[0] ){
        s.erase(0, s.find_first_not_of(ws));
        return s.erase(s.find_last_not_of(ws) + 1);
    }
    else{
        s.erase(s.begin(), std::find_if(s.begin(), s.end(),
            std::not1(std::ptr_fun<int, int>(std::isspace))));
        s.erase(std::find_if(s.rbegin(), s.rend(),
            std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
        return s;
    }
}

// Remove spaces at the right of a string
string& trim_right_string(string& s, const char *ws){
    if( ws==NULL || ws[0]==0) ws=" "; 
    return s.erase(s.find_last_not_of(ws) + 1);
}

void tolowercase(string& str){ 
    transform(str.begin(), str.end(), str.begin(), ::tolower); 
}

void toLowerCaseTrim(string& str){ 
    trim_string(str, " ");
    tolowercase(str);
}

// Cut the first word from a string
// Returns true if success (string cut), false if error (string unchanged)
bool GetNext(string &sword, string &sline, int lowerCase, string separators){
    std::size_t beg = sline.find_first_not_of(separators);
    if(beg == std::string::npos) return false;

    separators += "#";
    std::size_t end = sline.find_first_of(separators, beg);
    if(end == std::string::npos) {
        sword = sline.substr(beg);
        sline.clear();
    }
    else {
        sword = sline.substr(beg, end - beg);
        end = sline.find_first_not_of(separators, end);
        if(end == std::string::npos) sline.clear();
        else sline = sline.substr(end);
        if(sline[0] == '#') sline.clear();
    }
    if(lowerCase) tolowercase(sword);
    return true;
}

// Cut the first word from a string
// Returns true if success (string cut), false if error (string unchanged)
bool GetNextWord(string &sword, string &sline, int lowerCase){ 
    return GetNext(sword, sline, lowerCase, " \t\n\r");
}

// Cut the first list element from a string
// Returns true if success (string cut), false if error (string unchanged)
bool GetNextElem(string &sword, string &sline, int lowerCase){
    return GetNext(sword, sline, lowerCase, ",\n\r");
}

//======================================================================================================================
// Converting strings to integers and floating point numbers
//======================================================================================================================

// Convert word (array of symbols) to integer
int GetIntFromWord(const char* word, tErrAccumulator* ErrAccumulator){
    int temp;
    int Succeeded = sscanf( word, "%d", &temp);
    if( !Succeeded || Succeeded == EOF ) { // check if something went wrong during the conversion
        IO_error_add("GetIntFromWord: Can't convert %s to integer",word);
        temp=UNDEF_INT;
    }
    return temp;
}

// Convert word (array of symbols) to double
double GetDoubleFromWord(const char* word, tErrAccumulator* ErrAccumulator){
    NativeDouble Result; // number which will contain the result
    int Succeeded = sscanf( word, "%lf", &Result );
    if( !Succeeded || Succeeded == EOF ){ // check if something went wrong during the conversion
        IO_error_add("GetDoubleFromWord: Can't convert %s to double",word);
        MakeNaN(Result);
    }
    return Result;
}

// Check whether a word represents an integer
int IsInteger(const char* word) {
    const char* c = word;
    if(!*c) return 0; // empty string
    while(*c) {
        if((*c=='-' || *c=='+') && c==word) { c++; continue; } // minus or plus in the beginning
        if(*c<'0' || *c>'9') return 0;
        c++;
    }
    return 1;
}

// Check whether a word represents an floating-point number
// Format:
// [{+,-}] <value1> [{E,e} [{+,-}] <value2>]
// where <value1> is a series of digits with not more than one '.' 
// and <value2> is a series of digits
int IsNumber(const char* word) {
    int num_digits = 0, num_dots = 0, exp_part = 0;
    if(!*word) return 0; // empty string is not a number
    for(const char* c = word; *c; c++) {
        if(*c>='0' && *c<='9') { num_digits++; continue; }
        if(*c=='E' || *c=='e') {
            if(exp_part) return 0; // more than one 'E' - error
            if(!num_digits) return 0; // <value1> have no digits - error
            exp_part = 1;
            num_digits = 0;
            continue;
        }
        if(*c=='.') { 
            if(exp_part) return 0; // point after 'E' - error
            if(num_dots) return 0; // two points - error
            num_dots = 1;
            continue;
        }
        if(*c=='+' || *c=='-') {
            if(c==word) continue; // at the begginning - ok
            if(*(c-1)=='E' || *(c-1)=='e') continue; // just after 'E' - ok
            return 0;
        }
        return 0; // wrong symbol
    }
    if(!num_digits) return 0; // <value1> or <value2> (if the latter given) has no digits - error
    return 1;
}


// Convert arithmetic expression (array of symbols) to double
double ParseDoubleFromWord(const char* word, tErrAccumulator* ErrAccumulator) {
    if(IsNumber(word)) return GetDoubleFromWord(word, ErrAccumulator); // use quicklier subroutine for numbers

    tFormula cond;
    if(cond.Parse(word, false /*no log*/, ErrAccumulator)) {
        double val;
        MakeNaN(val);
        return val;
    }
    double val = cond.Evaluate();
    if(IsNaN(val)) IO_error_add("%s is NaN", word);
    return val;
}

// Convert arithmetic expression (array of symbols) to int
int ParseIntFromWord(const char* word, tErrAccumulator* ErrAccumulator) {
    if(IsInteger(word)) return GetIntFromWord(word, ErrAccumulator); // use quicklier subroutine for numbers

    tFormula cond;
    if(cond.Parse(word, false /*no log*/, ErrAccumulator)) return UNDEF_INT;
    double val = cond.Evaluate();
    if(IsNaN(val)) { IO_error_add("%s is NaN", word); return UNDEF_INT; }
    if(val > 0x7FFFFFFE || val < -0x7FFFFFFF) { IO_error_add("%s yields integer overflow", word); return UNDEF_INT; }
    return round_noisette(val);
}

// Get integer parameter <par> with from <line>
// Returns false if line is empty
// Returns true even if word-to-int convertion failed
bool GetIntFromLine(int &par,string &line, tErrAccumulator* ErrAccumulator){
    string word;
    if(!GetNextWord(word,line)) return false;

    par = GetIntFromWord(word.c_str(), ErrAccumulator);
    return true;
}

// Getting two integer parameters <par1>,<par2> from <line>
// Returns false if line has not two words
// Returns true even if word-to-int convertion failed
bool GetIntIntFromLine(int &par1, int &par2, string &line, tErrAccumulator* ErrAccumulator){ 
    string word1, word2;
    if(!GetNextWord(word1,line)) return false;
    if(!GetNextWord(word2,line)) return false;

    par1 = GetIntFromWord(word1.c_str(), ErrAccumulator);
    par2 = GetIntFromWord(word2.c_str(), ErrAccumulator);
    return true;
}

// Get integer and double parameters <par1>,<par2> from <line>
// Returns false if line has not two words
// Returns true even if word-to-int convertion failed
bool GetIntDoubleFromLine(int &par1, double &par2, string &line, tErrAccumulator* ErrAccumulator){ 
    string word1, word2;
    if(!GetNextWord(word1,line)) return false;
    if(!GetNextWord(word2,line)) return false;

    par1 = GetIntFromWord(word1.c_str(), ErrAccumulator);
    par2 = GetDoubleFromWord(word2.c_str(), ErrAccumulator);
    return true;
}

// Getting a list of <maxSize> doubles from <line>
// Returns actual number of elements read
// If an attempt to read the next element was failed, this will be noted in ErrAccumulator
int GetDoubleArrayFromLine(double *par, int maxSize, string &line, tErrAccumulator* ErrAccumulator){
    string word;
    for(int i=0;i<maxSize;i++){
        if(!GetNextWord(word,line)) return i;
        par[i] = GetDoubleFromWord(word.c_str(), ErrAccumulator);
        if( IsNaN(par[i]) ) return i;
    }
    return maxSize;
}

// Getting a list of <maxSize> doubles from <line>
// Returns actual number of elements read
// If an attempt to read the next element was failed, this will be noted in ErrAccumulator
int GetDoubleArrayFromLine(vector<double>& par,const int maxSize,string &line, tErrAccumulator* ErrAccumulator){
    string word;
    for(int i=0;i<maxSize;i++){
        if(!GetNextWord(word,line)) return i;
        double value = GetDoubleFromWord(word.c_str(), ErrAccumulator);
        if( IsNaN(value) ) return i;
        par.push_back(value);
    }
    return par.size();
}

//======================================================================================================================
// Enums processing
//======================================================================================================================

// Getting parameter from line from txt file
string GetEnumName(int val, const char *Dictionary){
    string res = "";
    string pname, pcode; 
    string D = Dictionary;
    
    while(GetNextWord(pname, D)){ //iterating all words in dictionary
        if(!GetNextWord(pcode,D)) return res;
        int code = GetIntFromWord(pcode.c_str());
        if(IS_UNDEF_INT(code)) crash("GetEnumName: dictionary \"%s\" is incorrect", Dictionary);
        if(code == val){ 
            res=pname;
            break; 
        }               
    }
    return res;
}

// Get enum value from a line, admitting literal code only. Returns UNDEF_INT if error
int GetEnumValue(const string& word, string Dictionary) {
    bool isFound = false;
    string pname, pcode;
    if(Dictionary.empty()) crash("Empty dictionary");
    int value = UNDEF_INT;
    while(GetNextWord(pname, Dictionary)){ //iterating all words in dictionary
        // если за строковым выражением не следует число, ошибка
        if(!GetNextWord(pcode, Dictionary)) crash("Error in dictionary (%s)", Dictionary.c_str());
        // если за строковым выражением не следует число, ошибка
        if(!IsInteger(pcode.c_str())) crash("Error in dictionary (%s)", Dictionary.c_str());

        if(CompareWords(pname, word)){ //if word is found
            if(isFound) crash("Error in dictionary (%s): value %s occured twice", Dictionary.c_str(), word.c_str());
            
            //remember the word found,its code,and setting flag to found state,then go on reading 
            isFound = true;
            value = GetIntFromWord(pcode.c_str(), NULL);
        }
    }
    return value;
}

// Returns 1 value <value> presents in <Dictionary>, 0 otherwise
int CheckEnumValue(const int value, string Dictionary) {
    bool isFound = false;
    string pname, pcode;
    if(Dictionary.empty()) crash("Empty dictionary");
    while(GetNextWord(pname, Dictionary)){ //iterating all words in dictionary
        // если за строковым выражением не следует число, ошибка
        if(!GetNextWord(pcode, Dictionary)) crash("Error in dictionary (%s)", Dictionary.c_str());
        // если за строковым выражением не следует число, ошибка
        if(!IsInteger(pcode.c_str())) crash("Error in dictionary (%s)", Dictionary.c_str());

        int tempValue = GetIntFromWord(pcode.c_str(), NULL);
        if(IS_UNDEF_INT(tempValue)) crash("Error in dictionary (%s), wrong code %s", Dictionary.c_str(), pcode.c_str());
        if(tempValue == value) isFound = true; //if word is found
    }
    return isFound;
}

// Get enum value from a word, admitting both digital and literal codes
// Returns UNDEF_INT if the code (digital or literal) is incorrect
int ProcessEnum(const string& word, const string& Dictionary){
    if(IsInteger(word.c_str())) {
        int value = GetIntFromWord(word.c_str());
        if(!CheckEnumValue(value, Dictionary)) return UNDEF_INT;
        else return value;
    }
    else {
        return GetEnumValue(word, Dictionary);
    }
}


//======================================================================================================================
// Formatted output to std::string 
//======================================================================================================================
string stringprintf(const char* fmt, ...) {
#ifdef EXTRAPRECISION // uses function from extraprecision.h
    va_list ap;
    va_start(ap,fmt);
    string S = vstringprintf<qd_real,1>(fmt, ap);
    va_end(ap);
    return S;
}
string _stringprintf(const char* fmt, ...) {
#endif
    int bufsize = 128; // начальный размер буфера
    char *sbuf = new char [bufsize]; //буфер для пихания в чар-строку

    do{//пока успешно не принтанем 
        va_list ap;
        va_start(ap,fmt);
        int nstr = vsnprintf(sbuf,bufsize,fmt,ap);
        va_end(ap);

        if((nstr>=0)&&(nstr<bufsize)){//корректно принтанулось
            string s =  string(sbuf); //результирующий стринг
            delete[] sbuf; 
            return s;            
        }
        if(nstr >= bufsize) bufsize = nstr + 1;
        if(nstr < 0) bufsize *= 2;
        
        delete[] sbuf;
        sbuf = new char [bufsize];

    }while(bufsize < 1000000); //пока строка не разрослась до абсурдных размеров
    delete[] sbuf;
    crash("stringprintf failed! buffer size too big"); 
}

//======================================================================================================================
// tLimitation class methods for initialization of parameters
//======================================================================================================================
void tLimitation::operator =(const tLimitation &aLimitation){
    lowLimit = aLimitation.lowLimit;
    symLowLim = aLimitation.symLowLim;
    symUpLim = aLimitation.symUpLim;
    upLimit = aLimitation.upLimit;
    parType = aLimitation.parType;
}

// Setting the limitations of the parameter
void tLimitation::SetLimits(int _symLowLim,double _lowLim,int _symUpLim,double _upLim){
    if((_symLowLim==LT || _symLowLim==LE || _symLowLim==UNLIM) && (_symUpLim==GT || _symUpLim==GE || _symUpLim==UNLIM)){
        int ib=_symLowLim; _symLowLim=_symUpLim; _symUpLim=ib;
        double db=_lowLim; _lowLim=_upLim; _upLim=db;
    }
    if(_symLowLim!=UNLIM && _symUpLim!=UNLIM && _lowLim>_upLim ) crash("SetLimits: lower limit > upper limit!");
    if(_symLowLim==LT || _symLowLim==LE ) crash("SetLimits: wrong symbol for lower limit!");
    if(_symUpLim==GT || _symUpLim==GE ) crash("SetLimits: wrong symbol for upper limit!");

    parType = PARTYPE_DOUBLE;

    switch( _symLowLim ){
    case UNLIM: symLowLim = UNLIM; break;
    case GT:    symLowLim = GT; lowLimit=_lowLim; break;
    case GE:    symLowLim = GE; lowLimit=_lowLim; break;
    }
    switch( _symUpLim ){
    case UNLIM: symUpLim = UNLIM; break;
    case LT:    symUpLim = LT; upLimit=_upLim; break;
    case LE:    symUpLim = LE; upLimit=_upLim; break;
    }
}

// Setting the limitations of the parameter
void tLimitation::SetLimits(const int _symLowLim,const int _lowLim,const int _symUpLim,const int _upLim){
    SetLimits(_symLowLim,static_cast<double>(_lowLim),_symUpLim,static_cast<double>(_upLim));
    parType = PARTYPE_INT;
}

// Checking the value of parameter
bool tLimitation::CheckValue(const double value) const {
    if(IsNaN(value)) return false; 
    const double tinyloc = 1e-16;
    bool lowGood = false;
    bool upGood = false;
    if( symLowLim==UNLIM ) lowGood = true;
    if( symLowLim==GT && (value>lowLimit) ) lowGood = true;
    if( symLowLim==GE && (value>lowLimit || fabs(value-lowLimit)<tinyloc) ) lowGood = true;
    if( symUpLim==UNLIM ) upGood = true;
    if( symUpLim==LT && (value<upLimit)) upGood = true;
    if( symUpLim==LE && (value<upLimit || fabs(value-upLimit)<tinyloc) ) upGood = true;
    return (lowGood && upGood);
}

bool tLimitation::CheckValue(const int value) const { // Checking the value of parameter
    bool lowGood = false;
    bool upGood = false;
    if( symLowLim==UNLIM ) lowGood = true;
    if( symLowLim==GT && value> round_noisette(lowLimit) ) lowGood = true;
    if( symLowLim==GE && value>=round_noisette(lowLimit) ) lowGood = true;
    if( symUpLim==UNLIM ) upGood = true;
    if( symUpLim==LT && value< round_noisette(upLimit) ) upGood = true;
    if( symUpLim==LE && value<=round_noisette(upLimit) ) upGood = true;
    return (lowGood && upGood);
}

//======================================================================================================================
// tParameter class methods
//======================================================================================================================
tParameter::tParameter(const tParameter &aPar){ //Copy constructor
    ptrDouble = aPar.ptrDouble;
    ptrInt = aPar.ptrInt;
    ptrString = aPar.ptrString;
    ptrStringArray = aPar.ptrStringArray;
    ptrIntArray = aPar.ptrIntArray;
    ptrDoubleArray = aPar.ptrDoubleArray;
    limit = aPar.limit;
   
    parName=aPar.parName; 
    Dictionary=aPar.Dictionary; 
 
    isSet = aPar.isSet;
    crashIt = aPar.crashIt;
    parType = aPar.parType;
}

bool tParameter::CheckParamName(const string& _parName){
    if(CheckSpaces(_parName.c_str()))
        crash("tParameter::initParamBase: param name can't have spaces! [%s]", _parName.c_str());

    if(_parName.empty() == false) {
        // проверка, что имя состоит только из латинских букв, цифр и разрешенных символов
        for(size_t i = 0; i < _parName.length(); i++) { //проверка,что имя состоит только из латинских букв,цифр и '.[]<>
            if(!(_parName[i] == '\'' || _parName[i] == '-' || _parName[i] == '.' || (_parName[i] >= '0' && _parName[i] <= '9')
                || _parName[i] == '<' || _parName[i] == '>' || (_parName[i] >= 'A' && _parName[i] <= 'Z')
                || _parName[i] == '[' || _parName[i] == ']' || _parName[i] == '_'
                || (_parName[i] >= 'a' && _parName[i] <= 'z')))
                crash("tParameter::initParamBase: param name contains illegal characters! [%s]", _parName.c_str());
            // 39 ', 45 -, 46 ., 48-57 0-9, 60 <, 62 >, 65-90 A-Z, 91 [, 93 ], 95 _, 97-122 a-z
        }

        // проверка, что имя начинается с латинской буквы или разрешенного символа
        if(!(_parName[0] == '\'' || _parName[0] == '-' || _parName[0] == '<' || _parName[0] == '>' ||
            (_parName[0] >= 'A' && _parName[0] <= 'Z') || _parName[0] == '[' || _parName[0] == ']' || _parName[0] == '_'
            || (_parName[0] >= 'a' && _parName[0] <= 'z')))
            crash("tParameter::initParamBase: param name starts with illegal character! [%s]", _parName.c_str());
        // 39 ', 45 -, 60 <, 62 >, 65-90 A-Z, 91 [, 93 ], 95 _, 97-122 a-z
    }

    return true;
}

// Setting parameters attributes

void tParameter::initParamBase(const string& _parName, const tParType _parType, const int _crashIt){
    _init();
    if (CheckParamName(_parName)) parName = _parName;
    crashIt = _crashIt;
    parType = _parType;
    if(CheckSpaces(parName.c_str()))
       crash("tParameter::initParamBase: param name can't have spaces! [%s]", parName.c_str());
}

tParameter::tParameter(int &par,const string& _parName,const tParType _parType,const int _crashIt,
                       const int _symLowLim, const int _lowLim,const int _symUpLim,const int _upLim){
    initParamBase(_parName,  _parType, _crashIt);
    ptrInt = &par;
    limit.SetLimits(_symLowLim,_lowLim,_symUpLim,_upLim);
    if(parType!=PARTYPE_INT && parType!=PARTYPE_BOOL) crash("Parser error: incompatible types");
}
tParameter::tParameter(double &par,const string& _parName,const tParType _parType, const int _crashIt,
                       const int _symLowLim, const double _lowLim,const int _symUpLim,const double _upLim){
    initParamBase(_parName,  _parType, _crashIt);
    ptrDouble = &par;
    limit.SetLimits(_symLowLim,_lowLim,_symUpLim,_upLim);
    if(parType!=PARTYPE_DOUBLE) crash("Parser error: incompatible types");
}
tParameter::tParameter(double *par,const string& _parName,const tParType _parType, const int _crashIt){
    initParamBase(_parName,  _parType, _crashIt);
    ptrDouble = par;
    limit.SetLimits(UNLIM,0,UNLIM,0);
    if(parType!=PARTYPE_DOUBLE3) crash("Parser error: incompatible types");
}
tParameter::tParameter(string& par, const string& _parName, const tParType _parType, const int _crashIt) {
    initParamBase(_parName,  _parType, _crashIt);
    ptrString = &par;
    limit.SetLimits(UNLIM,0,UNLIM,0);
    if(parType!=PARTYPE_WORD && parType!=PARTYPE_STRING) crash("Parser error: incompatible types");
}
tParameter::tParameter(vector<string>& par, const string& _parName, const tParType _parType, const int _crashIt) {
    initParamBase(_parName,  _parType, _crashIt);
    ptrStringArray = &par;
    limit.SetLimits(UNLIM,0,UNLIM,0);
    if(parType!=PARTYPE_WORDLIST) crash("Parser error: incompatible types");
}
tParameter::tParameter(vector<int>& par, const string& _parName, const tParType _parType, const int _crashIt) {
    initParamBase(_parName,  _parType, _crashIt);
    ptrIntArray = &par;
    limit.SetLimits(UNLIM,0,UNLIM,0);
    if(parType!=PARTYPE_INTLIST) crash("Parser error: incompatible types");
}
tParameter::tParameter(vector<double>& par, const string& _parName, const tParType _parType, const int _crashIt) {
    initParamBase(_parName,  _parType, _crashIt);
    ptrDoubleArray = &par;
    limit.SetLimits(UNLIM,0,UNLIM,0);
    if(parType!=PARTYPE_DOUBLELIST) crash("Parser error: incompatible types");
}
tParameter::tParameter(int &par,const string& _parName,const tParType _parType,const int _crashIt, const char *dict){
    initParamBase(_parName,  _parType, _crashIt);
    ptrInt = &par;
    limit.SetLimits(UNLIM,0,UNLIM,0);
    if(dict) Dictionary =  string(dict);
    else Dictionary.clear();
    if(parType!=PARTYPE_ENUM) crash("Parser error: incompatible types");
}

bool tParameter::CheckPtr(const string& par) const {
    if(&par == ptrString) return true;
    else return false;
}

// Synonyms for optional parameters
// Attention! Parser doesn't work properly with dd_real and qd_real:
//   first, all the internal procedures are in double precision;
//   second, the lower bytes are nullified when requesting values dd_real or qd_real variables 
//   (if the values will not be read from the file, the lower bytes of the default values will be corrupted)
template<> void tParamManager::Request(double &par, const string& _parName){
    RequestParameter(par, _parName, PARTYPE_DOUBLE, IO_DONTCRASH);
}
#ifdef EXTRAPRECISION // парсер сам работает на qd_real, поэтому с другими типами проблемы
    template<> void tParamManager::Request(NativeDouble&, const string&){ printf("tParamManager::Request: error, only qd_real can be requested in EXTRAPRECISION mode"); exit(0); }
    template<> void tParamManager::Request(dd_real&, const string&){ printf("tParamManager::Request: error, only qd_real can be requested in EXTRAPRECISION mode"); exit(0); }
#else // парсер работает на double, поэтому просто считываем старший член. Младшие приходится занулять, иначе там могут быть большие значения
    #ifdef EXTRAPRECISION_COLESO
        template<> void tParamManager::Request(qd_real &par, const string& _parName){
            par[1] = par[2] = par[3] = 0.0;
            RequestParameter(par[0], _parName, PARTYPE_DOUBLE, IO_DONTCRASH);
        }
        template<> void tParamManager::Request(dd_real &par, const string& _parName){
            par[1] = 0.0;
            RequestParameter(par[0], _parName, PARTYPE_DOUBLE, IO_DONTCRASH);
        }
    #endif
#endif
template<> void tParamManager::Request(vector<double> &par, const string& _parName) {
    RequestParameter(par, _parName, PARTYPE_DOUBLELIST, IO_DONTCRASH);
}
template<> void tParamManager::Request(int& par, const string& _parName){
    RequestParameter(par, _parName, PARTYPE_INT, IO_DONTCRASH);
}
template<> void tParamManager::Request(string& par, const string& _parName){
    RequestParameter(par, _parName, PARTYPE_STRING, IO_DONTCRASH);
}

//Функция для описания параметра типа PARTYPE_INT и PARTYPE_BOOL
void tParamManager::RequestParameter(int &par,const string& _parName,const tParType _parType,const int _crashIt,
    const int _symLowLim,const int _lowLim,const int _symUpLim,const int _upLim){
    if (_parType == PARTYPE_BOOL) {
        if ((_symLowLim != UNLIM && _symLowLim != GE) || _lowLim != 0 || (_symLowLim != UNLIM && _symLowLim != LE) 
            || (_upLim != 0 && _upLim != 1))
            crash("tParamManager::RequestParameter: limitations for type bool are incorrect");
        AddParameter(tParameter(par, _parName, _parType, _crashIt, GE, 0, LE, 1));
    }
    else AddParameter(tParameter(par, _parName, _parType, _crashIt, _symLowLim, _lowLim, _symUpLim, _upLim));
}

//добавляет новый параметр newParameter в вектор tParamsVectorBase
void tParamManager::AddParameter(const tParameter &newParameter) {
    for(size_t i=0; i<size(); i++) {
        if((at(i).GetParName())[0] && CompareWords(at(i).GetParName(), newParameter.GetParName()))
            crash("tParamManager::AddNew: double request of parameter %s!", at(i).GetParName());
    }
    tParamsVectorBase::push_back(newParameter);
}

bool tParameter::GetValue(string pline, bool log, tErrAccumulator* ErrAccumulator) {
    int err = ErrAccumulator ? ErrAccumulator->Count() : 0; // remember current number of errors
    if(isSet) { IO_error_add("Double definition of parameter %s", parName.c_str()); return false; }

    { // проверка на отключение параметра
        string tstring = pline, tword;
        if(GetNextWord(tword, tstring)) {
            if(CompareWords(tword, PARAMOFF)) {
                isSet = false;
                if(log) pprintf("%30s: %s\n", parName.c_str(), PARAMOFF.c_str());
                return true;
            }
        }
    }

    int tempInt;
    double tempDouble;
    string word;
    switch(parType){
    case PARTYPE_INT:
        trim_right_string(pline); // remove spaces at the right
        tempInt = ParseIntFromWord(pline.c_str(), ErrAccumulator);
        if(ErrAccumulator != NULL && ErrAccumulator->Count() > err) break; // an error has occured
        if(!limit.CheckValue(tempInt))
            IO_error_add("Parameter %s=%d is out of limitations!", parName.c_str(), tempInt);
        else
            *ptrInt = tempInt;
        break;
    case PARTYPE_DOUBLE:
        trim_right_string(pline); // remove spaces at the right
        tempDouble = ParseDoubleFromWord(pline.c_str(), ErrAccumulator);
        if(ErrAccumulator != NULL && ErrAccumulator->Count() > err) break; // an error has occured
        if(!limit.CheckValue(tempDouble))
            IO_error_add("Parameter %s=" FP_FORMAT " is out of limitations!", parName.c_str(), tempDouble);
        else
            *ptrDouble = tempDouble;
        break;
    case PARTYPE_DOUBLE3: {
        double vals[3];
        for(int i = 0; i<3; i++) {
            if(!GetNextElem(word, pline)) {
                IO_error_add("Not enough values for parameter \"%s\"", parName.c_str());
                break;
            }
            vals[i] = ParseDoubleFromWord(word.c_str(), ErrAccumulator);
        }

        if(!(ErrAccumulator != NULL && ErrAccumulator->Count() > err)) {
            ptrDouble[0] = vals[0]; ptrDouble[1] = vals[1]; ptrDouble[2] = vals[2];
        }
        break;
    }
    case PARTYPE_STRING: {
        size_t found = pline.find('#');
        if(found == string::npos) found = pline.length();
        while(found>0 && (pline[found - 1] <= 0x20 && pline[found - 1] >= 0)) found--;
        *ptrString = pline.substr(0, found);
        break;
    }
    case PARTYPE_WORD:
        if(!GetNextWord(word, pline, 0)) IO_error_add("Can not read parameter \"%s\"", parName.c_str());
        else {
            trim_string(pline);
            if(!pline.empty()) IO_error_add("Unknown remainder \"%s\" in line after reading the parameter \"%s\"", 
                                             pline.c_str(), parName.c_str());
            else *ptrString = word;
        }
        break;
    case PARTYPE_BOOL:
        if(!GetNextWord(word, pline)) { *ptrInt = 1; break; } // no value = TRUE
        if(CompareWords(word, PARAMTRUE)) { *ptrInt = 1; break; } // explicitly set TRUE
        if(CompareWords(word, PARAMFALSE)) { *ptrInt = 0; break; } // explicitly set FALSE
        tempInt = ParseIntFromWord(word.c_str(), ErrAccumulator);
        if(tempInt != 0 && tempInt != 1) IO_error_add("Parameter %s=%d is out of limitations", parName.c_str(), *ptrInt);
        else *ptrInt = tempInt;
        break;
    case PARTYPE_ENUM:
        if(!GetNextWord(word, pline)) { IO_error_add("Can not read parameter \"%s\"", parName.c_str()); break; }
        tempInt = ProcessEnum(word, Dictionary);
        if(IS_UNDEF_INT(tempInt)) IO_error_add("Wrong value \"%s\" for parameter %s", word.c_str(), parName.c_str());
        else *ptrInt = tempInt;
        break;
    case PARTYPE_INTLIST:
    {
        vector<int> val_list;
        while(GetNextElem(word, pline, false /*no lowercase*/)) {
            int val = ParseIntFromWord(word.c_str(), ErrAccumulator);
            if(ErrAccumulator != NULL && ErrAccumulator->Count() > err) break;
            val_list.push_back(val);
        }
        if(ErrAccumulator != NULL && ErrAccumulator->Count() > err) break;
        // Uncomment the following line to forbid empty list
        //if(val_list.empty()) { IO_error_add("Can not read parameter \"%s\"",parName.c_str()); break; }

        *ptrIntArray = val_list;
    }
    break;
    case PARTYPE_DOUBLELIST:
    {
        vector<double> val_list;
        while(GetNextElem(word, pline, false /*no lowercase*/)) {
            double val = ParseDoubleFromWord(word.c_str(), ErrAccumulator);
            if(ErrAccumulator != NULL && ErrAccumulator->Count() > err) break;
            val_list.push_back(val);
        }
        if(ErrAccumulator != NULL && ErrAccumulator->Count() > err) break;
        *ptrDoubleArray = val_list;
    }
    break;
    case PARTYPE_WORDLIST:
        while(GetNextElem(word, pline, false /*no lowercase*/))
            ptrStringArray->push_back(word);
        break;
    default: crash("tParameter::GetValue: unknown parType %i", parType);
    }

    if(ErrAccumulator != NULL && ErrAccumulator->Count() > err){
        IO_error_add("Reading param %s failed", parName.c_str());
    }
    else {
        isSet = true;
        if(log) PrintParameter();
    }
    return true;
}

void tParameter::PrintParameter() const { // Printing parameter to log file
    pprintf("%30s: ", parName.c_str());
    switch( parType ){
    case PARTYPE_ENUM: 
        pprintf("%d (%s) \n",*ptrInt, GetEnumName(*ptrInt, Dictionary.c_str()).c_str()); break;
    case PARTYPE_INT:
        pprintf("%d\n",*ptrInt); break;
    case PARTYPE_DOUBLE:
        if(IsNaN(*ptrDouble)) pprintf("NaN\n");
        else pprintf(FP_FORMAT "\n",*ptrDouble);
        break;
    case PARTYPE_DOUBLE3:
        pprintf("[" FP_FORMAT ", " FP_FORMAT ", " FP_FORMAT "]\n", ptrDouble[0], ptrDouble[1], ptrDouble[2]); break;
    case PARTYPE_STRING:
    case PARTYPE_WORD:
        pprintf("%s\n",ptrString->c_str()); break;
    case PARTYPE_BOOL:
        pprintf("%s\n",*ptrInt ? "true" : "false"); break;
    case PARTYPE_INTLIST:
        pprintf("[");
        for(size_t i=0; i<ptrIntArray->size(); i++) {
            if(i) pprintf(", ");
            pprintf("%i ", (*ptrIntArray)[i]);
        }
        pprintf("]\n");
        break;
    case PARTYPE_DOUBLELIST:
        pprintf("[");
        for(size_t i=0; i<ptrDoubleArray->size(); i++) {
            if(i) pprintf(", ");
            pprintf(FP_FORMAT, (*ptrDoubleArray)[i]);
        }
        pprintf("]\n");
        break;
    case PARTYPE_WORDLIST:
        pprintf("[");
        for(size_t i=0; i<ptrStringArray->size(); i++) {
            if(i) pprintf(", ");
            pprintf("%s", (*ptrStringArray)[i].c_str());
        }
        pprintf("]\n");
        break;
    case PARTYPE_UNDEFTYPE:  break;
    }
}

//======================================================================================================================
//tFileBuffer class methods
//======================================================================================================================

// создает буффер с текстом
void tFileBuffer::AllocText(const string& text) {
    AllocText(text.c_str(), text.length());
}

// создает буффер с текстом
void tFileBuffer::AllocText(const char* text, size_t length) {
    Clear(); // очистка

    // записываем текст построчно
    string temp;
    size_t nPos;
    for(size_t i = 0; i < length; i++) {
        while(text[i] != 0xA /*LF*/ && text[i] != 0xD /*CR*/ && text[i]!=0) {
            temp += text[i];
            Sum += text[i];
            i++;
        }

        if((nPos = temp.find('#')) != string::npos) temp = temp.substr(0, nPos);
        trim_string(temp);

        if(!temp.empty()) { Lines.push_back(temp); temp.clear(); }
    }

    // настраиваем метки на строки
    EnfOfZoneLine = Lines.size();
    if(EnfOfZoneLine) FirstLine = 0;
    GotoZone();
}

// создает буффер с текстом
void tFileBuffer::AllocText(FILE* file) {
    Clear(); // очистка

    // записываем текст построчно
    string temp;
    size_t nPos;
    char c = (char)fgetc(file);
    while(!feof(file)) {
        while(c != 0xA /*LF*/ && c != 0xD /*CR*/ && !feof(file)) {
            temp += c;
            Sum += c;
            c = (char)fgetc(file);
        }

        while((c == 0xA /*LF*/ || c == 0xD /*CR*/) && !feof(file)) {
            c = (char)fgetc(file);
        }

        if((nPos = temp.find('#')) != string::npos) temp = temp.substr(0, nPos);
        trim_string(temp);

        if(!temp.empty()) { Lines.push_back(temp); temp.clear(); }
    }

    // настраиваем метки на строки
    EnfOfZoneLine = Lines.size();
    if(EnfOfZoneLine) FirstLine = 0;
    GotoZone();
}

void tFileBuffer::InitMacro(bool log) {
    ListMacroLocal.insert(ListMacroLocal.end(), ListMacroGlobal.begin(), ListMacroGlobal.end());
    vector<string> ListCmdMacroLocal;

    // проход по глобальной командной строке - формируем локальные списоки
    ListCmdMacroLocal = CmdLine.GetFileMacros(GetFnameShort());

    // работаем со своим списком макро-подстановок, так как могут быть переопределения в файле
    // проход по буферу - ищем и обрабатываем макро-подстановки
    for(size_t i = 0; i < Lines.size(); i++) {
        // заподозрили макро-подстановку
        // отметка об использовании ставиться автоматически при выделении Used, если строка начиналась с '$',
        // так как в список макроподстановок она попадает всегда
        if(Lines[i][0] == '$') ListMacroLocal.AddOrReplace(Lines[i], log);
    }

    // проход по командной строке - ищем и обрабатываем макро-подстановки
    for(size_t j = 0; j < ListCmdMacroLocal.size(); j++) ListMacroLocal.AddOrReplace(ListCmdMacroLocal.at(j), log);

    ListMacroLocal.Evaluate(log); // вычисление значений всего вектора

    // подставляем значения в текст
    if(ListMacroLocal.size()) {
        for(size_t i = 0; i < Lines.size(); i++) {
            if(ListMacroLocal.Substitution(Lines[i], log) < 0) {
                // произошла ошибка при макроподстановке
                Lines.erase(Lines.begin() + i); i--;
            }
        }
    }

    ListMacroLocal.ErrAccumulator.MoveTo(&ErrAccumulator);
}

// удаление внутренних данных
void tFileBuffer::Clear(void){
    Sum = 0;
    ListMacroLocal.clear();
    Lines.clear();
    if(Used) delete[] Used;
    Used = NULL;
    CurrentLine = FirstLine = EnfOfZoneLine = NextZoneLine = -1;
    FNameFull.clear();
    // убираем префикс с названием файла, использованный в случае возникновения ошибок
    ErrAccumulator.ERRPREF.clear();
}

// делает текущей следующую зону и устанавливает указатель для чтения на ее начало
void tFileBuffer::GotoNextZone(void) {
    if(NextZoneLine >= 0) FirstLine = NextZoneLine;
    else FirstLine = EnfOfZoneLine;
    EnfOfZoneLine = Lines.size();
    NextZoneLine = -1;
    GotoZone();
}

// устанавливает указатель для чтения на начало текущей зоны
void tFileBuffer::GotoZone(void) {
    CurrentLine = FirstLine;
}

// Reading all file to buffer, returns 1 if error
int tFileBuffer::LoadFile(const string& fName, int CrashIt, bool log) {
#ifdef _NOISETTE
    InformKGB
#endif
    if(Lines.size()>0) crash("tFileBuffer: can't read to buffer because it is in use");

    FILE* pFILE = NULL;
#ifdef _NOISETTE // Испольуем встроенную систему ввода-вывода, чтобы вносить исправления в путь
    int FileID = -1;
    if(!fName.empty()) {
        FileID = IO_Open(fName, IO_READ, IO_BINARY, IO_DONTCRASH);
        if(FileID >= 0) pFILE = IO_GetFile(FileID);
    }
#else // Если встроенной системы ввода-вывода нет, используем стандартную функцию
    if(!fName.empty()) pFILE = fopen(fName.c_str(), "rb");
#endif

    if(pFILE == NULL){
        if(CrashIt) /*IO_error_add*/ crash("tFileBuffer::LoadFile: file %s not found\n", fName.c_str());
        if(log) pprintf("tFileBuffer::LoadFile: file %s not found\n", fName.c_str());
        FNameFull = fName;
        return 1;
    }
    if(log) pprintf("tFileBuffer: loading %s\n", fName.c_str());

    { // to remove 01.2021
        string file;
        PatchOldToNewConvert(pFILE, file, fName, log);
        AllocText(file);
    }
    //    AllocText(pFILE);

    FNameFull = fName;
    // устанавливаем префикс для ошибок, чтоб не передавать имя файла
    ErrAccumulator.ERRPREF = string("In file ") + GetFnameShort() + ": ";
    if (!CompareWords(GetFnameShort(),MACROTXT)) InitMacro(log);

#ifdef _NOISETTE
    IO_Close(FileID);
#else
    fclose(pFILE);
#endif

    return 0;
}

// создает массив флагов об использовании строк (если его еще нет)
void tFileBuffer::AllocUsed(void) {
    if(!Used) {
        Used = new bool[Lines.size()];
        for(size_t i = 0; i < Lines.size(); i++) {
            if(Lines[i][0] == '$') Used[i] = true; // макро-подстановка
            else Used[i] = false;
        }
    }
}

// ставит флаг об использовании предыдущей строке
void tFileBuffer::MarkPrevLineAsUsed(void) {
    ASSERT(CurrentLine > 0, "invalid line number");
    if(Used) Used[CurrentLine - 1] = true;
}

// Get next line from buffer <pf>
void tFileBuffer::GetNextLine(string &sbuf){
    if(CurrentLine < 0 || CurrentLine>= EnfOfZoneLine) sbuf.clear();
    else sbuf = Lines.at(CurrentLine++);
}

// Переводит файл с bool в формате (+\-)name в формат name 0(1)
void tFileBuffer::PatchOldToNewConvert(FILE* pFILE, string &file, const string& /*fName*/, bool /*log*/) {
    string line;
    bool plus, minus, comment;
    do {
        plus = false; minus = false; comment = false;
        ::GetNextLine(line, pFILE, true);
        if(line.size() > 0) {
            size_t i = 0;

            while(i < line.size() && line[i] == ' ') i++; // перемотка пробелов
            if(i != line.size()) {

                while(i < line.size() && line[i] == '#') { i++; comment = true; } // перемотка комментариев
                if(i != line.size()) {

                    while(i < line.size() && line[i] == ' ') i++;  // перемотка пробелов
                    if(i != line.size()) {

                        if(line[i] == '+') plus = true;
                        if(line[i] == '-') minus = true;
                        if(plus || minus) {
                            if(i + 1 != line.size()) {
                                if(!(comment == true && (line[i + 1] == '-' || line[i + 1] == '+'))) {
                                    if(line[i] == ' ' || line[i] == '#') 
                                        crash("tFileBuffer: can't convert file because it has single '-' or '+'");

//if(log) pprintf0("=WARNING: in file %s obsolete format has been fixed for the parameter %s", fName.c_str(), line.c_str());
                                    line.erase(i, 1);

                                    while(i < line.size() && (line[i] > ' ' || line[i] < 0)) i++;//перемотка до пробелов
                                    if(plus) line.insert(i, " 1");
                                    if(minus) line.insert(i, " 0");
                                }
                            }
                        }
                    }
                }
            }
        }
        file.append(line);
    } while(!feof(pFILE));
}

//Try to find end of the current zone
int tFileBuffer::MakeZone(const char* KeyWord) {
#ifdef _NOISETTE
    InformKGB
#endif
    string line, word;
    int cCurrentLine = CurrentLine;
    GotoZone();
    while (!EndOfZone()) {
        GetNextLine(line);
        string pline = line;
        if(!GetNextWord(word, pline)) continue;
        if(CompareWords(KeyWord, word)) {
            MarkPrevLineAsUsed();
            EnfOfZoneLine = CurrentLine-1;
            NextZoneLine = CurrentLine;
            CurrentLine = cCurrentLine;
            return 0;
        }
    }
    EnfOfZoneLine = Lines.size();
    NextZoneLine = -1;
    CurrentLine = cCurrentLine;
    return 1;
}

// Copies current zone (part of file) to file buffer FB
void tFileBuffer::CopyActualZoneTo(tFileBuffer& FB) const {
    FB.Clear();
    for (int i = FirstLine; i<EnfOfZoneLine; i++) FB.Lines.push_back(Lines.at(i));
}

// возвращает только имя файла без пути
string tFileBuffer::GetFnameShort(void) const {
    const char* lastSlash = strrchr(FNameFull.c_str(), '/');
    const char* lastBackslash = strrchr(FNameFull.c_str(), '\\');
    size_t strtName;
    if(lastBackslash == lastSlash) strtName = 0;
    else {
        if(lastBackslash>lastSlash) strtName = lastBackslash - FNameFull.c_str() + 1;
        else strtName = lastSlash - FNameFull.c_str() + 1;
    }

    return FNameFull.c_str() + strtName;
}

// выводит неиспользованные строки в текущей зоне
int tFileBuffer::PrintUnknownKeywords(void) {
    if(!Used) return 0;

    GotoZone();
    int printed = 0;
    while (!EndOfZone()) {
        bool disabled = false;
        string line, word, tstring;
        GetNextLine(line);
        if(Used[CurrentLine - 1] || !GetNextWord(word, line, 0)) continue;

        GetNextWord(tstring, line, 0);
        ListMacroLocal.Substitution(tstring, false);
        if(CompareWords(tstring, PARAMOFF)) disabled = true;

        if(!disabled) {
            // check for disabling parameter by "paramname=" option
            tstring = CmdLine.GetFileParam(GetFnameShort(), word);

            if(!tstring.empty()) {
                ListMacroLocal.Substitution(tstring, false);
                if(CompareWords(tstring, PARAMOFF)) disabled = true;
            }
        }

        if(!disabled) {
            pprintf0("=WARNING: UNKNOWN KEYWORD= %s\n", word.c_str());
        }
        printed++;
    }

    CmdLine.PrintUnusedFileArguments(GetFnameShort());
    return printed;
}

//======================================================================================================================
//  Additional methods for buffers
//======================================================================================================================

// Reads list of strings from file buffer FB to list
void parse_words_list(tFileBuffer & fb, vector<string>& list){
    list.clear();
    while(!fb.EndOfZone()) {
        string word, line;
        fb.GetNextLine(line);
        if(GetNextWord(word, line, 0)) list.push_back(word);
    }
}

//======================================================================================================================
//  tCmdLine class methods
//======================================================================================================================

// Добавление всей ком строки (не забывать убирать экзешник при передаче!), т.е. (argc-1,argv+1)
void tCmdLine::AddArguments(int argc, char** argv) {
    for(int iarg = 0; iarg < argc; iarg++) AddArg(argv[iarg]);
    if(ErrAccumulator.Check()){
        IO_error_add("Errors in command line");
        ErrAccumulator.Check_crash();
    }
}

// Добавление одного аргумента
void tCmdLine::AddArgument(char* param) {
    AddArg(param);
    if(ErrAccumulator.Check()){
        IO_error_add("Errors in command line");
        ErrAccumulator.Check_crash();
    }
}

// Добавление одного аргумента (обходим аргументы и растаскиваем на 3 кучи)
void tCmdLine::AddArg(char* param) {
    if(*(param) == '-') { // 1) параметры ком. строки
        char* c = param;
        c++;
        // проверка, что имя параметра не начинается с цифры или точки
        if(*c != 0 && (*c<'0' || *c>'9') && *c != '.') { AddStandardArgument(param); return; }
    }
    if(strchr(param, ':') != NULL) { // 2) оверрайд файлового ввода
        char* c = param;
        while(*c != 0 && *c != ':' && *c != '"' && *c != '=') c++; // проверка, что ':' не находится внутри ""
        if(*c == ':') {
            c++;
            // проверка, что ':' не находится внутри пути (например, d:\...)
            if(*c != '\\' && *c != '/') { AddFileArgument(param); return; }
        }
    }

    AddUnnamedArgument(param); // 3) старорежимные безымянные фиксированные параметры 
}

// Добавить аргумент (оверрайд файлового ввода)
void tCmdLine::AddFileArgument(char* param) {
    string fname, pname, pvalue;
    char* c = param;
    while(*c != 0 && *c != ':') c++;
    if(*c != 0) {
        char buf = *c;
        *c = 0;
        fname = string(param);
        *c = buf;
        c++;
        char* d = c;
        while(*c != 0 && *c != '=') c++;
        if(*c) {
            buf = *c;
            *c = 0;
            pname = string(d);
            *c = buf;
            if(*(c + 1) == 0) pvalue=PARAMOFF;
            else {
                c++;
                d = c;
                pvalue = string(d);
            }
            //else {
            // if (*(c + 1) == '=') 
            //    IO_error_add("tCmdLine::AddFileArgument: incorrect combination of characters '==' in '%s'", param);
            //}
        }
        else {
            pname = string(d);
            pvalue = " ";
        }
        
        trim_string(fname);
        trim_string(pname);

        // проверки имени параметра
        if(pname.empty()) {
            IO_error_add("Missing parameter name with value %s for file %s in command line!", pvalue.c_str(), fname.c_str());
            return;
        }

        if(pname.size() >= 2) {
            if(pname[0] == '$' && pname[pname.size() - 1] == '$') {
                // проверка имени для макро
                if(!tParameter::CheckParamName(pname.substr(1, pname.size()-2))) return;
            }
            else { if(!tParameter::CheckParamName(pname)) return; }
        }
        else { if(!tParameter::CheckParamName(pname)) return; }

        // проверка переопределения параметра
        for(size_t j = 0; j < ListFileArguments.size(); j++) {
            if(CompareWords(ListFileArguments.at(j).Name, pname)) {
                if(CompareWords(ListFileArguments.at(j).Fname, fname) || ListFileArguments.at(j).Fname.empty() || fname.empty()) {
                    IO_error_add("Double definition of parameter %s for file %s in command line!", pname.c_str(), fname.c_str());
                    return;
                }
            }
        }
        ListFileArguments.push_back(tCmdLineArgFile(fname, pname, pvalue));
    }
}

// Добавить аргумент (обычные аргументы)
void tCmdLine::AddStandardArgument(char* param){
    string result(param);
    string pname, pvalue;
    while(result.length() > 0 && result.at(0) == '-') {
        result.erase(0, 1);
    }
    if(result.length() > 0) {
        if(result.find("No-") == 0 || result.find("NO-") == 0 || result.find("no-") == 0 || result.find("nO-") == 0) {
            result.erase(0, 3);
            while(result.length() > 0 && result.at(0) == '-') result.erase(0, 1);
            if(result.find('=') == string::npos) {
                pname = result;
                pvalue = "0";
            }
            else crash("tCmdLine::AddStandardArgument: character '=' and \"no\" in '%s' at the same time", param);
        }
        else {
            size_t fnd=result.find('=');
            if(fnd != string::npos) {
                pname=result.substr(0, fnd); fnd++;
                pvalue = result.substr(fnd, result.size()-fnd);
            }
            else { pname = result; pvalue = " "; }
        }
    }
    else crash("tCmdLine::AddStandardArgument: empty parameter '%s'", param);

    // проверки имени параметра
    trim_string(pname);
    if(pname.empty()) {
        IO_error_add("Missing parameter name with value %s in command line!", pvalue.c_str());
        return;
    }
    if(!tParameter::CheckParamName(pname)) return;

    // проверка переопределения параметра
    for(size_t j = 0; j < ListStandardArguments.size(); j++) {
        if(CompareWords(ListStandardArguments.at(j).Name, pname)) {
            IO_error_add("Double definition of parameter %s in command line!", pname.c_str());
            return;
        }
    }
    ListStandardArguments.push_back(tCmdLineArgStandard(pname, pvalue));
}

// Добавить аргумент (обычные аргументы)
void tCmdLine::AddUnnamedArgument(char* param){
    ListUnnamedArguments.push_back(string(param));
}

// Вернуть следующий неиспользованный безымянник, и выставить ему отметку об использовании
string tCmdLine::GetNextUnnamedArgument() {
    for(size_t i = 0; i < ListUnnamedArguments.size(); i++) {
        if(!ListUnnamedArguments.at(i).Used) {
            ListUnnamedArguments.at(i).SetUsed();
            return ListUnnamedArguments.at(i).Value;
        }
    }
    return string("");
}

// Вернуть стандартный аргумент с именем parName
string tCmdLine::GetStandardArgument(const string &parName) {
    for(size_t i = 0; i < ListStandardArguments.size(); i++) {
        if(CompareWords(parName, ListStandardArguments.at(i).Name)) {
            ListStandardArguments.at(i).SetUsed();
            return ListStandardArguments.at(i).Value;
        }
    }
    return string("");
}

// Печать всех аргументов
void tCmdLine::PrintAllArguments(void) const {
    if(!ListFileArguments.size() && !ListStandardArguments.size() && !ListUnnamedArguments.size()) return; 
    pprintf("tCmdLine: arguments set\n");

#ifdef _NOISETTE 
    pprintf(LINE_SEPARATOR_);
#endif

    for(size_t i = 0; i<ListFileArguments.size(); i++) {
        if(!ListFileArguments.at(i).Fname.empty()) 
            pprintf("%30s: %s (%s)\n",ListFileArguments.at(i).Name.c_str(),
                    ListFileArguments.at(i).Value.c_str(), ListFileArguments.at(i).Fname.c_str());
        else pprintf("%30s: %s\n", ListFileArguments.at(i).Name.c_str(),
                      ListFileArguments.at(i).Value.c_str());
    }
    for(size_t i = 0; i<ListStandardArguments.size(); i++) 
        pprintf("%30s: %s\n", ListStandardArguments.at(i).Name.c_str(),
                 ListStandardArguments.at(i).Value.c_str());
    for(size_t i = 0; i<ListUnnamedArguments.size(); i++)
        pprintf("%30s  %s\n", "", ListUnnamedArguments.at(i).Value.c_str());
}

// Печать неиспользованных аргументов
void tCmdLine::PrintUnusedArguments(void) const {
    for(size_t i = 0; i<ListFileArguments.size(); i++) {
        if(!ListFileArguments.at(i).Used) {
            if(!ListFileArguments.at(i).Fname.empty())
                pprintf0("=WARNING: UNUSED CMDLINE ARGUMENT= %s: %s (%s)\n", ListFileArguments.at(i).Name.c_str(),
                         ListFileArguments.at(i).Value.c_str(), ListFileArguments.at(i).Fname.c_str());
            else pprintf0("=WARNING: UNUSED CMDLINE ARGUMENT= %s: %s\n", ListFileArguments.at(i).Name.c_str(),
                          ListFileArguments.at(i).Value.c_str());
        }
    }
    PrintUnusedStandardArguments();
    PrintUnusedUnnamedArguments();
}

// Печать неиспользованных неименнованных аргументов для файла fileName
void tCmdLine::PrintUnusedUnnamedArguments(void) const {
    for(size_t i = 0; i<ListUnnamedArguments.size(); i++) {
        if(!ListUnnamedArguments.at(i).Used)
            pprintf0("=WARNING: UNUSED CMDLINE ARGUMENT= %s\n", ListUnnamedArguments.at(i).Value.c_str());
    }
}

// Печать неиспользованных стандартных аргументов
void tCmdLine::PrintUnusedStandardArguments(void) const {
    for(size_t i = 0; i < ListStandardArguments.size(); i++) {
        if(!ListStandardArguments.at(i).Used)
            pprintf0("=WARNING: UNUSED CMDLINE ARGUMENT= %s: %s\n", ListStandardArguments.at(i).Name.c_str(),
                     ListStandardArguments.at(i).Value.c_str());
    }
}

// Печать неиспользованных аргументов для файла fileName
void tCmdLine::PrintUnusedFileArguments(const string &fileName) const {
    for(size_t i = 0; i<ListFileArguments.size(); i++) {
        if(CompareWords(fileName, ListFileArguments.at(i).Fname) &&
           !ListFileArguments.at(i).Used) {
            if(!ListFileArguments.at(i).Fname.empty())
                pprintf0("=WARNING: UNUSED CMDLINE ARGUMENT= %s: %s (%s)\n", ListFileArguments.at(i).Name.c_str(),
                         ListFileArguments.at(i).Value.c_str(), ListFileArguments.at(i).Fname.c_str());
            else pprintf0("=WARNING: UNUSED CMDLINE ARGUMENT= name: %s, value: %s\n", ListFileArguments.at(i).Name.c_str(),
                          ListFileArguments.at(i).Value.c_str());
        }
    }
}

// Удаление аргументов из всех списков 
void tCmdLine::ClearAllArguments(void) {
    ClearStandardArguments();
    ClearFileArguments();
    ClearUnnamedArguments();
}

// Удаление аргументов (обычные аргументы ком. строки)
void tCmdLine::ClearStandardArguments() {
    ListStandardArguments.clear();
}

// Удаление аргументов (оверрайд файлового ввода) 
void tCmdLine::ClearFileArguments() {
    ListFileArguments.clear();
}

// Удаление аргументов (старорежимные безымянные фиксированные аргументы)
void tCmdLine::ClearUnnamedArguments() {
    ListUnnamedArguments.clear();
}

// Возвращает макроподстановки (для аргументов с указанным именем файла и с пропущенным именем файла)
vector<string> tCmdLine::GetFileMacros(const string &fileName, bool globalMacro) {
    vector<string> result;
    string buf;
    for(size_t j = 0; j < ListFileArguments.size(); j++){
        if((globalMacro && (ListFileArguments.at(j).Fname.empty() || CompareWords(fileName, ListFileArguments.at(j).Fname))) ||
           (!globalMacro && CompareWords(fileName, ListFileArguments.at(j).Fname))) {
            if(ListFileArguments.at(j).Name[0] == '$') { // заподозрили макро-подстановку
                ListFileArguments.at(j).SetUsed();
                buf = ListFileArguments.at(j).Name;
                buf.append(" ");
                buf.append(ListFileArguments.at(j).Value);
                result.push_back(buf);
            }
        }
    }
    return result;
}

// Возвращает параметр (для аргументов с указанным именем файла и с пропущенным именем файла)
string tCmdLine::GetFileParam(const string &fileName, const string &paramName) {
    for(size_t j = 0; j < ListFileArguments.size(); j++){
        if(ListFileArguments.at(j).Fname.empty() || CompareWords(fileName, ListFileArguments.at(j).Fname)) {
            if(CompareWords(ListFileArguments.at(j).Name, paramName)) {
                ListFileArguments.at(j).SetUsed();
                return ListFileArguments.at(j).Value;
            }
        }
    }
    return string("");
}

//======================================================================================================================
//  tMacroSubstitution class methods
//======================================================================================================================

// Инициализирует макро-подстановку из строки line. Если флаг dollars используется, то имя должно быть обрамлено 
// символами $ с двух сторон, в противном случае - нет
bool tMacroSubstitution::InitMacroSubstitution(const string& line, bool dollars, tErrAccumulator* ErrAccumulator) {
    Value = line;

    // проверяем, что в строке line содержится имя и значение через разделитель
    if(!GetNextWord(Name, Value)) {
        IO_error_add("InitMacroSubstitution: Macro substitution \"%s\" should consist of name and value", line.c_str());
        return false;
    }

    if(!dollars) return true; // если обрамление имени долларами не подразумевается

    // если имя обрамлено $
    // отрезаем начальный и конечный символ $
    if(Name.length() >= 2) {
        if(Name[0] == '$' && Name[Name.size() - 1] == '$') {
            Name = string(Name.c_str(), 1, Name.length() - 2);
            return true;
        }
    }

    IO_error_add("Macro substitution \"%s\" must have a name starting and ending with the character \'$\'", line.c_str());
    return false;
}

//======================================================================================================================
//  tMacroSubstitutionVector class methods
//======================================================================================================================

// Добавляет новую макро-подстановку с текущий вектор из строки line
// log - печать в лог
// fromMacroFile - флаг, определяющий формат имени макро-подстановки (обрамление символом $)
bool tMacroSubstitutionVector::AddNew(const string& line, bool log, bool fromMacroFile) {
    tMacroSubstitution newMacro;

    // пробуем создать новую макро-подстановку из line
    if(newMacro.InitMacroSubstitution(line, !fromMacroFile, &ErrAccumulator) != true) return false;

    // проверяем, что такого имени в списке нет
    for(size_t i = 0; i < size(); i++) {
        if(at(i).Name == newMacro.Name) {
            at(i).Value = newMacro.Value;
            ErrAccumulator.Add(stringprintf("tMacroSubstitutionVector::AddNew: Redefinition of macro substitution \"%s\"\n", line.c_str()));
            return false;
        }
    }

    // добавляем в список
    push_back(newMacro);

    if(log) pprintf("Local macro %s added: %s\n", newMacro.Name.c_str(), newMacro.Value.c_str());
    return true;
}

bool tMacroSubstitutionVector::AddOrReplace(const string& line, bool log) {
    tMacroSubstitution newMacro;
    if(newMacro.InitMacroSubstitution(line,true,&ErrAccumulator) != true) return false;

    // если это замена, то сначала ищем в списке нужную
    for(size_t i = 0; i < size(); i++) {
        if(CompareWords(at(i).Name,newMacro.Name)) {
            at(i).Value = newMacro.Value;
            if(log) pprintf("Local macro %s changed: %s\n", newMacro.Name.c_str(), newMacro.Value.c_str());
            return true;
        }
    }

    return AddNew(line, log, false);
}

bool tMacroSubstitutionVector::Evaluate(bool log) {
    if(!size()) return true;

    enum tEval {
        NO,
        INT,
        DOUBLE
    };
    tEval eval = NO;
    size_t ready = 0, remainder = size();
    static string EVAL = "EVAL ";
    static string EVAL_INT = "EVAL_INT ";

    // Обрезаем пробелы в значении
    for(size_t i = 0; i < size(); i++)  trim_string(at(i).Value);

    do {
        ready = 0; // счетчик готовых макро-подстановок на текущем проходе (те, где нет EVAL)
        for(size_t i = 0; i < size(); i++) {
            string oldValue;
            // Обрезаем пробелы в значении макроса
            if(CompareWords(at(i).Value.substr(0, EVAL.length()), EVAL)) {
                oldValue = at(i).Value;
                at(i).Value = at(i).Value.substr(EVAL.length());
                eval = DOUBLE;
            }
            if(CompareWords(at(i).Value.substr(0, EVAL_INT.length()), EVAL_INT)) {
                oldValue = at(i).Value;
                at(i).Value = at(i).Value.substr(EVAL_INT.length());
                eval = INT;
            }

            // пробуем выполнить подстановку
            if(!oldValue.empty()) {
                int result;
                result=Substitution(at(i).Value,log);
                if(result == 1 || result == 0) { // все возможные подстановки сделаны или нет вложенных подстановок
                    double val = ParseDoubleFromWord(at(i).Value.c_str()); // вычисляем значение выражения
                    ErrAccumulator.Check_crash(); // если в процессе подстановки возникли ошибки
                    // записываем его обратно в строковую переменную
                    if(eval == DOUBLE) at(i).Value = stringprintf("%.15e", val);
                    else at(i).Value = stringprintf("%i", int(val));
                    if (log) pprintf0("Macro evaluated: $%s$ -> %s\n", at(i).Name.c_str(), at(i).Value.c_str());
                    ready++;
                }
            }
            else ready++;
        }
        if(remainder == size() - ready) {
            IO_error_add("tMacroSubstitutionVector::Evaluate: All macros can not be avaluated");
            return false; //за новый проход готовых подстановок не прибавилось
        }
        else remainder = size() - ready;
    } while(ready != size());

    if(eval!=NO && log) pprintf0("All macros evaluated\n");
    return true;
}

// 0 - макросов нет, 1 - выполнена макроподстановка, 2 - выполнены макроподстановки, но есть EVAL, -1 - ошибка
int tMacroSubstitutionVector::Substitution(string &line, bool log) {
    int result = 0, count = 0;
    static int MAX_COUNT = 1000; // счетчик максимального количества подстановок
    size_t posStart = 0;
    string oldline = line;

    while((posStart = line.find('$', posStart)) != string::npos && count < MAX_COUNT) { // ищем начальный доллар, 
                                                                                        // пропуская два доллара подряд
        count++;

        if(posStart == line.find("$$", posStart)) {
            posStart += 2;
            continue;
        }

        size_t posEnd = posStart;
        posEnd = line.find('$', posStart + 1);
        while((posEnd = line.find('$', posEnd)) != string::npos) { // ищем конечный доллар, пропуская два доллара подряд
            if(posEnd != line.find("$$", posEnd)) break;
            else posEnd += 2;
        }
        if(posEnd == string::npos) {
            result = -1;
            IO_error_add("tMacroSubstitutionVector::Substitution: Incorrect using of character '$' in line \"%s\"", 
                         oldline.c_str());
            return result;
        }

        // Ищем макрос по глобальному списку
        size_t i;
        string tempstrig = line.substr(posStart + 1, posEnd - posStart - 1);
        static string EVAL = "EVAL";
        for(i = 0; i < size(); i++) {
            if(CompareWords(at(i).Name, tempstrig)) {
                if(CompareWords(at(i).Value.substr(0, EVAL.length()), EVAL)) {
                    result = 2;
                    return result;
                }

                tempstrig = line.substr(posEnd + 1, line.length());
                line = line.substr(0, posStart);
                line.append(at(i).Value);
                line.append(tempstrig);
                if(result == 0) result = 1;
                break;
            }
        }

        if(i == size()) {
            result = -2;
            IO_error_add("tParamManager::MacroSubstitution: No macro substitution for \"%s\" in line \"%s\"", 
                          tempstrig.c_str(), oldline.c_str());
            return result;
        }
    }

    if(count == MAX_COUNT) {
        result = -3;
        IO_error_add("tParamManager::MacroSubstitution: the number of substitutions exceeded for line \"%s\"",
                     oldline.c_str());
        return result;
    }

    if(result == 1 || result == 0) {
        // макросов нет
        // осталось заменить $$ на $
        while(line.find("$$") != string::npos) {
            string oldValue = line;
            size_t oldPos = oldValue.find("$$");
            line = line.substr(0,oldPos);
            line.append("$");
            line.append(oldValue.substr(oldPos + 2, oldValue.length() - oldPos - 2));
        }

        if(result && log) pprintf("MacroSubstitution: %s -> %s\n", oldline.c_str(), line.c_str());
    }

    return result;
}

void tMacroSubstitutionVector::InitGlobalMacro(string fname, bool CMDLine, bool log) {
    string parseLine;

    tFileBuffer FB(fname, IO_DONTCRASH);

    if(!CompareWords(FB.GetFnameShort(), MACROTXT)) 
        crash("macro substitutions must be described in the file %s", MACROTXT);

    while(!FB.EndOfZone()) {
        FB.GetNextLine(parseLine);
        ListMacroGlobal.AddNew(parseLine, log, true);
    }

    // Обрабатываем командную строку на новые макро-подстановки или переопределение
    if(CMDLine) {
        vector<string> macros = CmdLine.GetFileMacros(FB.GetFnameShort(),true);
        for(size_t j = 0; j <macros.size(); j++){
            ListMacroGlobal.AddOrReplace(macros.at(j), log);
        }
    }

    ListMacroGlobal.Evaluate(log); // вычисление значений для всего вектора
    ListMacroGlobal.ErrAccumulator.Check_crash();
}

 tMacroSubstitutionVector& tMacroSubstitutionVector::GlobalMacro(void) {
     return ListMacroGlobal;
}

//======================================================================================================================
//  tParamManager class methods
//======================================================================================================================

//======================================================================================================================
// Reading params from text file: 
// fileName -- name of the file, 
// log -- печать отчёта в парлог,
// CMDLine -- чтение аргументов командной строки
//======================================================================================================================
void tParamManager::ReadParamsFromFile(const string& fileName, bool log, bool CMDLine){
    int CrashIt = IO_DONTCRASH;
    for(size_t i=0; i < size(); i++)
        if(at(i).GetCrashIt()) CrashIt = IO_CRASH; // set IO_CRASH if any oblibatory parameter is present

    tFileBuffer FB(fileName, CrashIt,log);
    ReadParamsFromBuffer(FB,log,CMDLine);
    if(log) FB.PrintUnknownKeywords();
}

//======================================================================================================================
// Чтение параметров из файлового буфера согласно ранее заполненному перечню
// log -- печать отчёта в парлог
//======================================================================================================================
void tParamManager::ReadParamsFromBuffer(tFileBuffer& FB, bool log, bool CMDLine){
#ifdef _NOISETTE 
    InformKGB
        if(log) pprintf(LINE_SEPARATOR_);
#endif

    // Nothing to read - error
    if(tParamsVectorBase::empty()) crash("params list is empty");

    // Check that all parameters have names
    for(size_t i = 0; i < size(); i++)
        if(!*(at(i).GetParName()))
            crash("can't read unnamed parameter from file buffer");

    FB.AllocUsed();
    // устанавливаем префикс для ошибок, чтоб не передавать имя файла
    ErrAccumulator.ERRPREF = string("In file ") + FB.GetFnameShort() + ": ";

    vector<string> string_values(size());
    FB.GotoZone();
    // проход по файлу - обрабатываем параметры
    while(!FB.EndOfZone()) {
        string parseLine, word;
        FB.GetNextLine(parseLine); // считывание строки для обработки
        if(!GetNextWord(word, parseLine)) continue;

        for(size_t i = 0; i < size(); i++) {
            // если имя параметра в строке не совпадает с именем параметра
            if(!CompareWords(word, at(i).GetParName())) continue; 

            if(string_values[i].empty()) {
                string_values[i] = parseLine;
                if(string_values[i].empty()) string_values[i] = " ";
            }
            else IO_error_add("double definition of parameter %s", word.c_str());
            FB.MarkPrevLineAsUsed();
        }
    }

    // проход по командной строке - определяем или переопределям параметры
    if(CMDLine) {
        for(size_t i = 0; i < size(); i++){
            // parameter that must not be in this file -> skip
            if(at(i).GetAntiCrashIt() || at(i).GetAntiDontCrashIt()) continue;

            string line = CmdLine.GetFileParam(FB.GetFnameShort(), at(i).GetParName());
            if(FB.GetMacroList().Substitution(line, log) < 0) line = ""; // произошла ошибка при макроподстановке
            if(!line.empty()) string_values[i] = line;
        }
    }

    // считываем параметры
    for(size_t i = 0; i < size(); i++){
        if(string_values[i].empty()) continue; // value is not prescribed - skip (checks are below)
        if(at(i).GetAntiCrashIt()) {
            IO_error_add("Obsolete parameter %s must not be there!", at(i).GetParName());
            continue;
        }
        else if(at(i).GetAntiDontCrashIt()) {
            if(log) pprintf0("=WARNING: Obsolete parameter %s must not be there!\n", at(i).GetParName());
            continue;
        }
        // trying to read the value
        at(i).GetValue(string_values[i], log, &ErrAccumulator);
    }

    StatusInformation(log);

    FB.ErrAccumulator.MoveTo(&ErrAccumulator);
    // убираем префикс с названием файла, использованный в случае возникновения ошибок
    ErrAccumulator.ERRPREF.clear();
}

void tParamManager::StatusInformation(bool log) {
    // печать ошибок\предупреждений\информации об обязательных\необязательных\запрещенных параметрах
    bool printed = false;
    for(size_t i = 0; i < size(); i++){
        if(at(i).GetIsSet()) continue; // value is prescribed and read => skip
        if(at(i).GetAntiCrashIt() || at(i).GetAntiDontCrashIt()) { // value must not be set => ok
            // nothing to do
        }
        else if(at(i).GetCrashIt()){ // value must be set but is not set
            IO_error_add("Parameter %s is not set!", at(i).GetParName());
        }
        else { // value is optional and not set. Printing default values
            if(log) {
                if(!printed) {
                    pprintf("Default used for parameters:\n");
                    printed = true;
                }
                at(i).PrintParameter();
            }
        }
    }
}

tParameter& tParamManager::operator[] (const char* parname) {
    for(size_t i = 0; i < size(); i++)
        if(CompareWords(at(i).GetParName(), parname)) return at(i);
    crash("tParamManager: error, can't find parameter %s in the list!", parname);
}

void tParamManager::ReadParamsFromCommandLine(bool log){
    if(tParamsVectorBase::empty()) IO_error_add("ReadParamsFromCommandLine: params list is empty!");
    else {
        // processing unnamed command line parameters  
        // unnamed params must go in a specific fixed order
        // optional unnamed param can't exist if there are obligatory params or named params after it!
        // usage example:
        //    mytool.exe <unnamed param1> <unnamed param2> [unnamed param3] [file:namedparam=value] ...
        bool unnamedOptParam = false; // flag that unnamed optional param appeared

        for(size_t i = 0; i < size(); i++) {
            if(!(at(i).GetParName())[0] && !at(i).GetIsSet()) { // param without name - reading
                if(at(i).GetCrashIt() == 0) unnamedOptParam = true;
                else
                    if(unnamedOptParam)
                        crash("tParamManager::ReadParamsFromCommandLine: cmd line param can't be optional if obligatory param goes after it!");

                string arg = CmdLine.GetNextUnnamedArgument();
                if(arg.empty()) {
                    if(!unnamedOptParam) IO_error_add("Not enough records in command line to read obligatory parameter");
                }
                else {
                    if(ListMacroGlobal.Substitution(arg, log) < 0) {
                        arg = ""; // произошла ошибка при макроподстановке
                        ListMacroGlobal.ErrAccumulator.MoveTo(&ErrAccumulator);
                    }
                    at(i).GetValue(arg.c_str(), log);
                }
            }
        }

        for(size_t i = 0; i < size(); i++) {
            if((at(i).GetParName())[0]) {
                string arg = CmdLine.GetStandardArgument(at(i).GetParName());
                if(ListMacroGlobal.Substitution(arg, log) < 0) {
                    arg = ""; // произошла ошибка при макроподстановке
                    ListMacroGlobal.ErrAccumulator.MoveTo(&ErrAccumulator);
                }
                if(!arg.empty()) at(i).GetValue(arg.c_str(), log);
            }
        }

        StatusInformation(log);

        CmdLine.ErrAccumulator.MoveTo(&ErrAccumulator);
    }
    if(ErrAccumulator.Check()){
        IO_error_add("Errors in command line");
        ErrAccumulator.Check_crash();
    }
}

// Extracts parameters pairs <parName,parValue> to MapString2String
void tParamManager::ExtractParametersPairs(tFileBuffer &fb, MapString2String& par) {
    while(!fb.EndOfZone()) {
        string line, parName;
        fb.GetNextLine(line);
        if(GetNextWord(parName, line, 0)) {
            string parValue, word;
            tolowercase(parName);
            int count = 0;
            while(1) {
                if(!GetNextWord(word, line)) break;
                if(word[0] == '#') break; // Дошли до комментария
                if(count > 0) parValue += " ";
                else count++;
                parValue += word;
            }
            if(!parValue.empty()) par[parName] = parValue;
        }
    }
}

//======================================================================================================================
//Simplified formula parser 
//======================================================================================================================
#define ACT_NUMBER     0 // s_action::number = number  
#define ACT_ARGUMENT   1 // s_action::arg1 = argument's index
#define ACT_FUNCTION   2 // s_action::arg1 = argument's index, arg2 = function code
#define ACT_BINARY     3 // s_action::arg1,arg2 = argument indices, arg3 = operation code
#define ACT_NEGATION   4 // s_action::arg1 = argument's index
#define ACT_UNARYMINUS 5 // s_action::arg1 = argument's index

double func_sin(  double a){ return sin(a); }
double func_cos(  double a){ return cos(a); }
double func_tan(  double a){ return tan(a); }
double func_ctg(  double a){ return 1.0/tan(a); }
double func_asin( double a){ return asin(a); }
double func_acos( double a){ return acos(a); }
double func_atan( double a){ return atan(a); }
double func_actg( double a){ return 0.5*PiNumber-atan(a); }
double func_sh(   double a){ return 0.5*(exp(a)-exp(-a)); }
double func_ch(   double a){ return 0.5*(exp(a)+exp(-a)); }
double func_th(   double a){ return (exp(a)-exp(-a))/(exp(a)+exp(-a));}
double func_log(  double a){ return log(a);}
double func_log10(double a){ return log(a)/log(10.0); }
double func_exp(  double a){ return exp(a);}
double func_sqr(  double a){ return a*a;}
double func_sqrt( double a){ return sqrt(a);}
double func_fabs( double a){ return fabs(a);}

//typedef double(*tFuncDouble2Arg)(double,double);
//double func_min( double a, double b){ return MIN(a,b);}
//double func_max( double a, double b){ return MAX(a,b);}

const int MAX_OPER = 23;
const char operations[MAX_OPER][3] = 
// 0    1    2    3    4    5    6     7    8    9   10   11   12    13   14    15   16    17   18    19   20    21   22 
{ "(", ")", "!", "<", ">", "=", "==", "+", "-", "/","*", "^", "**", "|", "||", "&", "&&", "%", "\\", "<=",">=", "_", "__" };

#define W_RESULT 0
#define W_FUNCTION 1
#define W_BRAKE_OPEN 2
#define W_BRAKE_CLOSE 3
#define W_UNARY 4
#define W_BINARY 5

tFormula::tFormula(int _ntr){
    clear();
    ntr = _ntr;
    argvals.resize(ntr);
    argptrs.resize(ntr);
}

tFormula::tFormula(const tFormula &f) : vector<s_action>(f){
    funcs = f.funcs;
    consts = f.consts; 
    argnames=f.argnames; 
    argvals=f.argvals; 
    argptrs=f.argptrs; 
    argused=f.argused; 
    ntr = f.ntr;
}

void tFormula::clear(){

    ntr = 0; 
    funcs.clear(); 
    consts.clear(); 
    argnames.clear(); 
    argused.clear(); 
    argvals.clear(); 
    argptrs.clear(); 
    vector<s_action>::clear();

    //functions
    funcs.push_back( pair<string, tFuncDouble1Arg>("sin",   func_sin)); 
    funcs.push_back( pair<string, tFuncDouble1Arg>("cos",   func_cos)); 
    funcs.push_back( pair<string, tFuncDouble1Arg>("tan",   func_tan));
    funcs.push_back( pair<string, tFuncDouble1Arg>("tg",    func_tan));
    funcs.push_back( pair<string, tFuncDouble1Arg>("ctg",   func_ctg));
    funcs.push_back( pair<string, tFuncDouble1Arg>("asin",  func_asin));
    funcs.push_back( pair<string, tFuncDouble1Arg>("arcsin",func_asin));
    funcs.push_back( pair<string, tFuncDouble1Arg>("acos",  func_acos));
    funcs.push_back( pair<string, tFuncDouble1Arg>("arccos",func_acos));
    funcs.push_back( pair<string, tFuncDouble1Arg>("atan" , func_atan));
    funcs.push_back( pair<string, tFuncDouble1Arg>("arctan",func_atan));
    funcs.push_back( pair<string, tFuncDouble1Arg>("atg" ,  func_atan));
    funcs.push_back( pair<string, tFuncDouble1Arg>("arctg" ,func_atan));
    funcs.push_back( pair<string, tFuncDouble1Arg>("actg",  func_actg));
    funcs.push_back( pair<string, tFuncDouble1Arg>("arcctg",func_actg));
    funcs.push_back( pair<string, tFuncDouble1Arg>("sh",    func_sh));
    funcs.push_back( pair<string, tFuncDouble1Arg>("ch",    func_ch));
    funcs.push_back( pair<string, tFuncDouble1Arg>("th",    func_th));
    funcs.push_back( pair<string, tFuncDouble1Arg>("ln",    func_log));
    funcs.push_back( pair<string, tFuncDouble1Arg>("log",   func_log));
    funcs.push_back( pair<string, tFuncDouble1Arg>("lg",    func_log10));
    funcs.push_back( pair<string, tFuncDouble1Arg>("exp",   func_exp));
    funcs.push_back( pair<string, tFuncDouble1Arg>("sqr",   func_sqr));
    funcs.push_back( pair<string, tFuncDouble1Arg>("sqrt",  func_sqrt));
    funcs.push_back( pair<string, tFuncDouble1Arg>("abs",   func_fabs));

    //constants
    AddConst("pi", PiNumber);
    AddConst("_e",  exp(1.0));
}

void tFormula::SetNumThreads(int _ntr) {
    argvals.resize(_ntr);
    argptrs.resize(_ntr);
    // если уменьшаем число, то ничего делать не надо
    for(int trn=ntr; trn<_ntr; trn++) {
        argvals[trn].resize(argnames.size(), 0.0);
        argptrs[trn].resize(argnames.size(), NULL);
    }
    ntr = _ntr;
}

inline int tFormula::FindArg(const char *n) const{
    string sn = n;
    for(int i=0; i<(int)argnames.size(); i++) 
        if(CompareWords(sn, argnames[i])) return i;
    return -1;
}

inline int tFormula::FindConst(const char *n) const{
    string sn = n;
    for(int i=0; i<(int)consts.size(); i++) 
        if(CompareWords(sn, consts[i].first)) return i;
    return -1;
}

inline int tFormula::FindFunc(const char *n) const{
    string sn = n;
    for(int i=0; i<(int)funcs.size(); i++) 
        if(CompareWords(sn, funcs[i].first)) return i;
    return -1;
}

int tFormula::AddArg(const char *n, double v, double *p, int DoCheck){
    if(DoCheck)
    if(FindConst(n)>=0 || FindFunc(n)>=0 || FindArg(n)>=0)
        crash("tFormula::AddArg: can't add arg %s, its name is already used\n", n); 

    argnames.push_back(n);
    argused.push_back(0);
    for(int itrn=0; itrn<ntr; itrn++) {
        argvals[itrn].push_back(v);
        argptrs[itrn].push_back(p); 
    }
    return (int) argnames.size() - 1;
}

int tFormula::AddConst(const char *n, double v){
    if(FindConst(n)>=0 || FindFunc(n)>=0 || FindArg(n)>=0)
        crash("tFormula::AddConst: can't add const %s, its name is already used\n", n); 

    consts.push_back(pair<string, double>(n, v));
    return (int) consts.size() - 1;
}

int tFormula::SetArg(const char *n, double  v, int trn){
    if(FindConst(n)>=0 || FindFunc(n)>=0 )
        crash("tFormula::SetArg: can't set arg %s, it isn't an argument\n", n); 

    int i=FindArg(n); 
    if(i>=0){ // аргумент уже присутствует - выставляем значение
        if(trn>=0) {
            argvals[trn][i]=v; argptrs[trn][i]=NULL;
        }
        else for(int itrn=0; itrn<ntr; itrn++) {
            argvals[itrn][i]=v; argptrs[itrn][i]=NULL;
        }
        return i;
    }
    // добавляем новую переменную
    if(trn!=-1) crash("tFormula::SetArg: can't add new variable %s in parallel mode", n);
    return AddArg(n, v, NULL);
}

int tFormula::SetArgPtr(const char *n, double *p, int trn){
    if(FindConst(n)>=0 || FindFunc(n)>=0 )
        crash("tFormula::SetArgPtr: can't set arg %s, it isn't an argument\n", n); 

    int i=FindArg(n); 
    if(i>=0){ // аргумент уже присутствует - выставляем указатель
        if(trn>=0) {
            argvals[trn][i]=0.0; argptrs[trn][i]=p;
        }
        else for(int itrn=0; itrn<ntr; itrn++) {
            argvals[itrn][i]=0.0; argptrs[itrn][i]=p;
        }
        return i;
    }
    // добавляем новую переменную
    if(trn!=-1) crash("tFormula::SetArg: can't add new variable %s in parallel mode", n);
    return AddArg(n, 0.0, p);
}

void tFormula::AddDefaultVars(){ 
    AddArg("index");
    AddArg("x");
    AddArg("y");
    AddArg("z");
    AddArg("label");
}

bool tFormula::CheckArgUsed(const char *n)const{
    int i=FindArg(n); 
    if(i>=0) return argused[i]>0; 
    return 0;
}


// parsing string expression. If error, returns 1 and keeps the formula empty
int tFormula::Parse(const char* str, bool parseLog, tErrAccumulator* ErrAccumulator) {
    if(size() != 0) crash("formula is already in use");
    if(strlen(str)==0) { clear(); return 1; }

    string _str_ = string(" ") + str; // create a string, which will automatically release memory in destruction

    if(parseLog) pprintf("Parsing %s: stage 1\n", str);
    vector< pair<int,int> > P;

    char* c = const_cast<char*>(_str_.c_str()) + 1;
    while(*c && *c!='#') {
        if(*c <= 0x20) { c++; continue; }

        // Cutting symbolic or numeric (number + dot) code
        char* cc = c;
        int dollar = (*cc=='$');
        int letter = (*cc>='a' && *cc<='z') || (*cc>='A' && *cc<='Z') || (*cc == '"');
        int number = (*cc>='0' && *cc<='9') || *cc=='-' || *cc=='.'; // unary minus 
        int parentheses = *cc=='(' || *cc==')';
        int operation = *cc=='!' || *cc=='<' || *cc=='>' || *cc=='='|| *cc=='*'
            || *cc=='+' || *cc=='-' || *cc=='/' || *cc=='^' || *cc=='|' || *cc=='&' || *cc=='%' || *cc=='\\' || *cc=='_';
        if(!(letter || number || operation || dollar || parentheses))
            { IO_error_add("Error: wrong start symbol %c\n", *cc); clear(); return 1; }
        letter*=3; number*=3; operation*=3; dollar*=3; parentheses*=3;

        int quotes_opened = (*cc == '"');
        while((letter | number | operation | dollar | parentheses) & 1) {
            cc++;
            if(quotes_opened) {
                // continue with letters only inside the quotes. Nothing to check
                if(!*cc || *cc=='#') { IO_error_add("Error: unclosed quotes"); clear(); return 1; }
                if(*cc=='"') quotes_opened = 0;
                continue;
            }
            letter &= 2 + ((*cc>='a' && *cc<='z') || (*cc>='A' && *cc<='Z') || (*cc>='0' && *cc<='9'));
            number &= 2 + ((*cc>='0' && *cc<='9') || *cc=='.' || ((*cc=='-' || *cc=='+') && (*(cc-1)=='e' || *(cc-1)=='E')) || *cc=='e' || *cc=='E');
            operation &= 2 + (*cc=='=' || *cc=='|' || *cc=='&' || *cc=='*' || *cc=='_');
            dollar &= 2 + ((*cc>='0' && *cc<='9') || (*cc>='a' && *cc<='z') || (*cc>='A' && *cc<='Z'));
            parentheses &= 2;
        }
        if(number && operation) { // minus - unary or binary
            number = 0; cc = c+1; // унарный минус
        }

        char t = *cc; *cc = 0; // end of line marker
        if(parseLog) pprintf("Found code: %s\n", c);
        int arg = -1, type = -1;

        if(number) {
            type = 0; arg = (int)size();
            s_action A;
            A.cmd = ACT_NUMBER;
            A.number = GetDoubleFromWord(c, ErrAccumulator);
            push_back(A);
        }
        if(letter || dollar) { // argument, $argument or function 
            int cmd=-1;
            //looking for functions
            if(dollar) c++; // skipping $
            else{ //function can't have $ prefix
                for(int i=0; i<(int)funcs.size(); i++){
                    if(CompareWords(c, funcs[i].first)){ 
                        type = W_FUNCTION;
                        arg = cmd = i;
                        break;
                    }
                }
            }
            //if function not found, looking for constants
            if(cmd<0){
                string ccc = string(c);
                trim_string(ccc, "\""); // arguments can be inside the quotes - cut them!
                for(int i=0; i<(int)consts.size(); i++){
                    if(CompareWords(ccc, consts[i].first)){ 
                        type = 0; arg = (int)size();
                        s_action A;
                        cmd = A.cmd = ACT_NUMBER;
                        A.number = consts[i].second;
                        push_back(A);
                    }
                }
            }
            //if not found, looking for arguments
            if(cmd<0){
                string ccc = string(c);
                trim_string(ccc, "\""); // arguments can be inside the quotes - cut them!
                for(int i=0; i<(int)argnames.size(); i++){
                    if(CompareWords(ccc, argnames[i])){ 
                        type = 0; arg = (int)size();
                        s_action A;
                        cmd = A.cmd = ACT_ARGUMENT;
                        A.arg1 = i;
                        argused[i]=1;
                        push_back(A);
                    }
                }
            }
            if(cmd<0){
                IO_error_add("Error: operation or argument %s is undefined", c); clear(); return 1; 
            }
        }
        if(operation || parentheses) {
            int cmd;
            for(cmd=0; cmd<MAX_OPER; cmd++)
                if(CompareWords(c, operations[cmd])) break;
            if(cmd==MAX_OPER) { IO_error_add("Error: operation %s is undefined", c); clear(); return 1; }
            if(cmd==0) type = W_BRAKE_OPEN;
            else if(cmd==1) type = W_BRAKE_CLOSE;
            else if(cmd==2) type = W_UNARY;
            else {
                type = W_BINARY;
                arg = cmd;
            }
        }
        *(c=cc) = t;

        P.push_back(pair<int,int>(type, arg));
    }

    if(parseLog) pprintf("Parsing %s: stage 2\n", str);
    int NP = P.size();
    if(Parse2(P, 0, NP-1, ErrAccumulator)<0) { clear(); return 1; }

    if(parseLog) {
        pprintf("Operation list:\n");
        for(unsigned int i=0; i<size(); i++) {
            pprintf("index = %i  cmd = %i", i, (*this)[i].cmd);
            if((*this)[i].cmd == ACT_NUMBER) pprintf(" number = %e", (*this)[i].number);
            else pprintf(" arg1 = %i arg2 = %i arg3 = %i", (*this)[i].arg1, (*this)[i].arg2, (*this)[i].arg3);
            pprintf("\n");
        }
    }
    return 0;
}

double tFormula::Evaluate(const vector<const double*>& ptrs, const vector<double>& vals) const {
    int N = (int)size();
    if(N==0) { double v; MakeNaN(v); return v; }
#define TCOND_SIZE 16 // для коротких формул избегаем аллокейтов
    double rbuf[TCOND_SIZE];
    
    double* r=NULL;
    
    if(N>=TCOND_SIZE) r = new double[N];
    else{
        for(int i=0; i<TCOND_SIZE; i++) rbuf[i]=0.0;
        r = &rbuf[0];
    }

    for(int i=0; i<N; i++) {
        s_action A = (*this)[i];
        switch(A.cmd) {;
        case ACT_NUMBER: 
            r[i] = A.number; 
            break;
        case ACT_ARGUMENT:
            if(A.arg1<0 || A.arg1>= (int)vals.size() ) 
                crash("tFormula wrong named argument index %d", A.arg1);
            r[i] = (!ptrs.empty() && ptrs[A.arg1] ? *ptrs[A.arg1] : vals[A.arg1]); 
            break;
        case ACT_BINARY: {
            double a = r[A.arg1];
            double b = r[A.arg2];
            switch(A.arg3) {;
            case 3: r[i] = a<b; break;
            case 4: r[i] = a>b; break;
            case 5: case 6: r[i] = fabs(a-b) < 1e-11*(fabs(a)+fabs(b)+1e-11); break;
            case 7: r[i] = a+b; break;
            case 8: r[i] = a-b; break;
            case 9: r[i] = a/b; break;
            case 10: r[i] = a*b; break;
            case 11: case 12: r[i] = pow(a,b); break;
            case 13: case 14: r[i] = round_noisette(a) || round_noisette(b); break;
            case 15: case 16: r[i] = round_noisette(a) && round_noisette(b); break;
            case 17: r[i] = round_noisette(a) / round_noisette(b); break;
            case 18: r[i] = round_noisette(a) % round_noisette(b); break;
            case 19: r[i] = a<=b; break;
            case 20: r[i] = a>=b; break;
            case 21: r[i] = MIN(a,b); break;
            case 22: r[i] = MAX(a,b); break;
            default: crash("Internal error, unknown operation %i", A.arg3);
            }
            //if(IsNaN(r[i])) IO_error_add("Error, NaN detected, operation code = %i, argument = %e", A.arg2, a);
            break; }
        case ACT_FUNCTION: {
            double a = r[A.arg1];
            if(A.arg2<0 || A.arg2>= (int)funcs.size() ) 
                crash("tFormula wrong function index %d", A.arg2);
            // IO_error_add("Internal error, function code %i undefined", A.arg2);
            r[i] = funcs[A.arg2].second(a); 
            //if(IsNaN(r[i])) IO_error_add("Error, NaN detected, function code = %i, argument = %e", A.arg2, a);
            break;}
        case ACT_NEGATION: r[i] = (round_noisette(r[A.arg1])==0); break;
        case ACT_UNARYMINUS: r[i] = -r[A.arg1]; break;
        default: crash("Internal error, unknown action %i", A.cmd);
        }
        if(IsNaN(r[i])) { MakeNaN(r[N-1]); break; } // if NaN is detected, stop computations and return NaN
    }
    double retval = r[N-1];
    if(N>=TCOND_SIZE) delete[] r;
#undef TCOND_SIZE
    return retval;
}

double tFormula::Evaluate(int trn) const { // версия, автоматически подставляющая внутренние массивы
    if(ntr==1 && trn==-1) trn = 0;
    ASSERT(trn>=0 && trn<ntr, "wrong thread index\n"); // это не ошибка ввода, это баг
    return Evaluate(argptrs[trn], argvals[trn]);
}

//specialized backward-compatibility crap
double tFormula::Evaluate(int index, int label, const double* vec, int trn){
    SetArg("index", (double)index, trn);
    if(vec){ SetArg("x", vec[0], trn); SetArg("y", vec[1], trn); SetArg("z", vec[2], trn) ;}
    SetArg("label", (double)label, trn);
    return Evaluate(trn);
}

int tFormula::Test(int index, int label, const double* vec, int trn){
    double res = Evaluate(index, label, vec, trn);
    return fabs(res) > 0.5;
}

// auxiliary subroutine: substitution an expression by it's result
void tFormula::Substitute(vector< pair<int,int> > &P, int Open, int Close, int& Imax, unsigned int newcode) {
    P[Open].first = W_RESULT;
    P[Open].second = newcode;

    int i = Open+1, j = Close+1;
    while(j<=Imax) {
        P[i].first = P[j].first;
        P[i].second = P[j].second;
        i++; j++;
    }
    Imax = i-1;
}

int tFormula::BinaryOperations(vector< pair<int,int> > &P, int Imin, int& Imax, int op_min, int op_max, int op_min2, 
                               int op_max2, bool unMinusFlag, tErrAccumulator* ErrAccumulator) {
    while(1) {
        int Open = -1;
        for(int i=Imin; i<=Imax; i++) {
            if(P[i].first==W_BINARY) {
                if((P[i].second>=op_min && P[i].second<=op_max) || (P[i].second>=op_min2 && P[i].second<=op_max2)) {
                    Open = i;
                    break; 
                }
            }
        }
        if(Open==-1) return 0; // no more operations found

        if(Open==Imin) { IO_error_add("Error: binary operation sign on start"); return -1; }
        if(Open==Imax) { IO_error_add("Error: binary operation sign at the end"); return -1; }
        if(P[Open-1].first!=W_RESULT) { IO_error_add("Error: wrong binary operation argument (left)"); return -1; }
        if(P[Open+1].first!=W_RESULT) { 
            if(unMinusFlag) { Imin = Open+1; continue; } // справа может быть унарный минус - по флагу пропускаем пока
            IO_error_add("Error: wrong binary operation argument (right)");
            return -1;
        }
        s_action A;
        A.cmd = ACT_BINARY;
        A.arg1 = P[Open-1].second;
        A.arg2 = P[Open+1].second;
        A.arg3 = P[Open].second;
        push_back(A);

        Substitute(P, Open-1, Open+1, Imax, (int)size() - 1);
    }
}

// returns number of final action
int tFormula::Parse2(vector< pair<int,int> > &P, int Imin, int Imax, tErrAccumulator* ErrAccumulator) {
// first looking for parentheses
    int Open, Close;
    while(1) {
        Open = -1, Close = -1;
        for(int i=Imin; i<=Imax; i++)
            if(P[i].first==W_BRAKE_OPEN)
            { Open = i; break; }
        if(Open==-1) break; // '(' not found

        int NumOpenBraces = 1;
        for(int i = Open+1; i<=Imax; i++) {
            if(P[i].first==W_BRAKE_OPEN) NumOpenBraces++;
            if(P[i].first==W_BRAKE_CLOSE) NumOpenBraces--;
            if(NumOpenBraces==0) { Close = i; break; }
        }
        if(Close==-1) { IO_error_add("Error: can't find closing brace"); return -1; }

        int result = Parse2(P, Open+1, Close-1, ErrAccumulator);
        if(result<0) return result; // pass error code to the calling function

        // if there is a function before bracket - adding action
        if(Open>0 && P[Open-1].first==W_FUNCTION) {
            s_action A;
            A.cmd = ACT_FUNCTION;
            A.arg1 = result;
            A.arg2 = P[Open-1].second;
            push_back(A);
            result = (int)size() - 1;
            Open--;
        }
        Substitute(P, Open, Close, Imax, result);
    };

// Looking for logic hi-priority unary ops: !
    while(1) {
        Open = -1;
        for(int i=Imin; i<=Imax; i++)
            if(P[i].first==W_UNARY)
                { Open = i; break; }
        if(Open==-1) break; // not found 

        if(Open==Imax) { IO_error_add("Error: binary operation sign on start"); return -1; }
        if(P[Open+1].first!=W_RESULT) { IO_error_add("Error: wrong unary operation argument"); return -1; }
        s_action A;
        A.cmd = ACT_NEGATION;
        A.arg1 = P[Open+1].second;
        push_back(A);

        Substitute(P, Open, Open+1, Imax, (int)size() - 1);
    }

// Looking for binary operations
    // very hi-priority: ^, **
    if(BinaryOperations(P, Imin, Imax, 11, 12, -1, -1, true /*quiet: unMinusFlag on*/, ErrAccumulator)) return -1;
    // hi-priority: *, /
    if(BinaryOperations(P, Imin, Imax, 9, 10, -1, -1, true /*quiet: unMinusFlag on*/, ErrAccumulator)) return -1; 

// Looking for unary minus or unary plus
    while(1) {
        Open = -1;
        for(int i=Imin; i<=Imax; i++) {
            if(P[i].first==W_BINARY && (P[i].second>=7 && P[i].second<=8)) { 
                if(i!=Imin) if(P[i-1].first==W_RESULT) continue; // this is a binary operator
                Open = i;
                break; 
            }
        }
        if(Open==-1) break; // unary operator not found

        if(Open==Imax) { IO_error_add("Error: unary operation sign on the end"); return -1; }
        if(P[Open+1].first!=W_RESULT) { IO_error_add("Error: wrong unary operation argument"); return -1; }
        if(P[Open].second==7) { // plus
            Substitute(P, Open, Open+1, Imax, P[Open+1].second);
        }
        else if(P[Open].second==8) { // minus
            s_action A;
            A.cmd = ACT_UNARYMINUS;
            A.arg1 = P[Open+1].second;
            push_back(A);
            Substitute(P, Open, Open+1, Imax, (int)size() - 1);
        }
    };

// Looking for binary operations
    if(BinaryOperations(P, Imin, Imax, 11, 12, -1, -1, false, ErrAccumulator)) return -1; // very hi-priority: ^, **
    if(BinaryOperations(P, Imin, Imax, 9, 10, -1, -1, false, ErrAccumulator)) return -1;  // hi-priority: *, /
    if(BinaryOperations(P, Imin, Imax, 7, 8, -1, -1, false, ErrAccumulator)) return -1;   // lo-priority: +, -

// Looking for binary operations
    // comparison: >, <, ==, =, >=, <=, _, __
    if(BinaryOperations(P, Imin, Imax, 3, 6, 19, 22, false, ErrAccumulator)) return -1;
    // every remaining binary operations (logic operations)
    if(BinaryOperations(P, Imin, Imax, 0, 999, -1, -1, false, ErrAccumulator)) return -1;

    if(Imax-Imin != 0) {
        IO_error_add("Can't parse expression. Unknown error."); return -1;
    }
    return P[Imin].second;
}
//======================================================================================================================
// Miscellaneous parser subroutines
//======================================================================================================================
void tParamManager::GetUnnamedPathFromCommandLine(std::string &path) {
    size_t i;
    unsigned int count = 0;

    // проверяем все запрошенные параметры
    for(i = 0; i < size(); i++) {
        if(!(at(i).GetParName())[0]) { // находим неименнованные
            if(at(i).CheckPtr(path)) { // проверяем на соответствие 
                string arg = CmdLine.GetNextUnnamedArgument();
                if (arg.empty()) {
                    if (at(i).GetCrashIt())
                        IO_error_add("Not enough records in command line to read directory");
                }
                else {
                    if(ListMacroGlobal.Substitution(arg, false) < 0) {
                        arg = ""; // произошла ошибка при макроподстановке
                        ListMacroGlobal.ErrAccumulator.MoveTo(&ErrAccumulator);
                    }
                    at(i).GetValue(arg.c_str(), false);
                }
                return;
            }
            count++;
        }
    }

    if(i == size()) crash("The requested directory parameter is not specified in ParamManager");
}

//---------------------------------------------------------------------------------------------------------------------
// Выводит предупреждения обо всех параметрах, не заданных в файле. Возвращает число предупреждений
int tParamManager::WarnAboutAllAbsentParameters() {
    int mark = 0;

    // проверяем, что все параметры шайб считаны
    std::vector<tParameter> &_PM = *this; // у tParamManager переопределено [], поэтому кастим к базовому
    for(size_t iparam = 0; iparam < _PM.size(); iparam++) {
        if(!_PM[iparam].GetIsSet()) {
            pprintf("Warning! Parameter is not set. Assuming:\n");
            _PM[iparam].PrintParameter();
            mark++;
        }
    }
    return mark;
}
