#pragma once
#ifndef _COLESO_H
#define _COLESO_H

#include "pointfuncs.h"
#include "es_utils.h"
#include <vector>
#include <string>
#include <math.h>

class tFileBuffer;
#define DESCRIPTION(X) const char* description() const OVERRIDE { return X; }


//======================================================================================================================
// === ACOUSTICS =======================================================================================================
//======================================================================================================================

//======================================================================================================================
template<typename fpv>
struct s_Gaussian2D : tPointFunction_EP<fpv>, tSpaceForm<fpv> {
    DESCRIPTION("Gaussian pulse 2D");
//======================================================================================================================
    fpv FlowVelX, FlowVelY, FlowVelZ; // Background flow

protected:
    tCompoundGaussIntegrator<fpv> GI;

public:
    tFuncType Type() const OVERRIDE { return FUNC_PULSATIONS; } // поскольку наследуем напрямую tPointFunction_EP, указываем тип
    s_Gaussian2D():tSpaceForm<fpv>() { 
        FlowVelX=FlowVelY=FlowVelZ=0.0; 
        tSpaceForm<fpv>::numCoords = 2; // для правильной нормировки амплитуды источника
    }

    const char* filename() const { return "es_gaussian2d.txt"; }
    void ReadParams(tFileBuffer&) OVERRIDE;
    void CalcDirect(fpv t, fpv r, fpv& p, fpv& ur) const;
    //void CalcViaFI(fpv t, fpv r, fpv& p, fpv& ur) const;
    void CalcViaFourier(fpv t, fpv r, fpv& p, fpv& ur) const;
    void CalcViaBesselFourier(fpv t, fpv r, fpv& p, fpv& ur) const;
    void Init(void);
    void PointValue(fpv t, const fpv* coor, fpv* V) const; 
};
//======================================================================================================================


//======================================================================================================================
struct s_Gaussian3D : tPulsFunction, tSpaceForm<double> {
    DESCRIPTION("Gaussian pulse 3D");
//======================================================================================================================
    double FlowVelX, FlowVelY, FlowVelZ; // Скорость фонового потока
    s_Gaussian3D():tSpaceForm<double>() { FlowVelX=FlowVelY=FlowVelZ=0.0; }
    const char* filename() const { return "es_gaussian3d.txt"; }
    void ReadParams(tFileBuffer&) OVERRIDE;
    void CalcValues(double t, double x, double y, double z, double* uex) const;
    void Init(void);
    void PointValue(double t, const double* coor, double* V) const; 
};
//======================================================================================================================


//======================================================================================================================
struct s_EntropyVortex : tPulsFunction {
    DESCRIPTION("entropy and vortex wave");
//======================================================================================================================
    double FlowVelX, FlowVelY; // Скорость фонового потока
    double PerX, PerY; // Период по каждому направлению
    double Xterm, Yterm; // Координаты центра энтропийного и вихревого импульса
    double Bterm, Aterm1, Aterm2; // Полуширина и амплитуды энтропийного и вихревого импульса

    s_EntropyVortex() {
        FlowVelX = FlowVelY = 0.0;
        PerX = PerY = huge;
        Xterm = Yterm = 0.0;
        Bterm = Aterm1 = Aterm2 = 0.0;
    }
    const char* filename() const OVERRIDE { return "es_entropyvortex.txt"; }
    void ReadParams(tFileBuffer&) OVERRIDE;
    void PointValue(double t, const double* coor, double* V) const; 
};
//======================================================================================================================


//======================================================================================================================
struct s_Planar : tPulsFunction { // abstract class - for PlanarGauss and PlanarSinus
    DESCRIPTION("planar wave (sinus/gauss profile)");
//======================================================================================================================
    double Aterm;
    double Xterm, Yterm, Zterm;
    double Dir[3];
// Параметры фонового потока
    double RhoRef, URef[3], PRef, gam;
// For channel: Rmax
    double Rmax;
    s_Planar() {Dir[0] = 1.0; Dir[1] = Dir[2] = 0.0; gam = 1.4; RhoRef = 1.0;
                                             URef[0] = URef[1] = URef[2] = 0.0; PRef = 1.0 / gam; Rmax = huge;}
    const char* filename() const OVERRIDE { return "es_planar.txt"; }
    void ReadParams(tFileBuffer& FB) OVERRIDE;
    void Init() OVERRIDE;
};
//======================================================================================================================

//======================================================================================================================
struct s_PlanarSinus : s_Planar{
    DESCRIPTION("planar wave (sinus profile)");
//====================================================================================================================== 
    double Freq;
    double Phase;
    int AllowNegativePhase;
    s_PlanarSinus() { Freq = 1.0; Phase = 0.0; AllowNegativePhase = 0; }
    const char* filename() const OVERRIDE { return "es_planarsinus.txt"; }
    void ReadParams(tFileBuffer& FB) OVERRIDE;
    void PointValue(double t, const double* coor, double* V) const;
};
//======================================================================================================================

//======================================================================================================================
struct s_PlanarGauss : s_Planar {
    DESCRIPTION("planar wave (Gaussian profile)");
//====================================================================================================================== 
    double Bterm;
    double NormalDistance;  
    s_PlanarGauss() { Bterm = 1.0; NormalDistance = huge; }
    const char* filename() const OVERRIDE { return "es_planargauss.txt"; }
    void ReadParams(tFileBuffer& FB) OVERRIDE;
    void PointValue(double t, const double* coor, double* V) const;
};
//======================================================================================================================

//======================================================================================================================
struct s_Coaxial : tComplexPulsFunction {
    DESCRIPTION("acoustic mode in cylinder or between coaxial cyl-s");
//====================================================================================================================== 
    double FlowVel, SoundVel; // скорость потока и звука в канале
    int AsimuthalMode, RadialMode; // asimuthal and radial modes
    double rmin, rmax, kz; // cylinder radii; axial wave number
    double Ampl; // амплитуда
    double phase; // фаза
    int CoorAxis; // cylinder axis: 0 - X, 1 - Y, 2 - Z

    s_Coaxial() {
        rmin=1.0; rmax=2.0; kz=0.0; MakeNaN(kr); 
        AsimuthalMode = 8; RadialMode=0; Ampl=1.0; phase=0.0; 
        SoundVel=1.0; FlowVel=0.0; 
        CoorAxis = 2;
    }
    const char* filename() const OVERRIDE { return "es_coaxial.txt"; }
    void ReadParams(tFileBuffer& FB) OVERRIDE;
    void Init() OVERRIDE;
    void PointValue(double t, const double* coor, double* V) const;
    void SetEquationConstants(double /*_gam*/, double /*_Rey*/, double /*_Pr*/, double _SoundVel, double _FlowVel) OVERRIDE {
        FlowVel = _FlowVel; SoundVel = _SoundVel;
    }

private:
    double kr; // radial wave number (evaluated at Init())
    static int FindRoot(int nu, double rmin, double rmax, double &kr); // Calculating root of equation using Newton method
};
//======================================================================================================================

//======================================================================================================================
struct s_SinusVisc : tPulsFunction {
    DESCRIPTION("sinus wave with viscousity");
//======================================================================================================================
    double mu, gamma, Prandtl; // viscosity, specific ratio, Prandtl number
    double Ampl; // amplitude of acoustic mode
    double AmplVX, AmplVY, AmplVZ; // amplitudes of vortex mode (no-heat-conductive case only)
    double kx, ky, kz; // wave vector
    double ReOmega, ImOmega; // frequency for no-heat-conductive case. Not need to be specified by user
    double ReCoeff, ImCoeff; // auxiliary data for no-heat-conductive case. Not to be specified by user

    s_SinusVisc() {
        Ampl = 1.0; AmplVX=AmplVY=AmplVZ=0.0; kx=ky=kz=0.0; mu = 0.0;
        gamma=1.4; Prandtl=1e50; MakeNaN(ReOmega); MakeNaN(ImOmega); MakeNaN(ReCoeff); MakeNaN(ImCoeff);
    }
    const char* filename() const OVERRIDE { return "es_sinusvisc.txt"; }
    void SetEquationConstants(double _gam, double _Rey, double _Pr, double /*_SoundVel*/, double /*_FlowVel*/) OVERRIDE {
        mu = 1.0 / _Rey; gamma = _gam; Prandtl = _Pr;
    }

    void ReadParams(tFileBuffer& FB) OVERRIDE;
    void Init(void) OVERRIDE; 
    void PointValue(double t, const double* coor, double* V) const;
};
//======================================================================================================================

//======================================================================================================================
struct s_Polynom : tScalarFunction {
    DESCRIPTION("Chebyshev polynomial planar profile");
//======================================================================================================================
    double SoundSpeed;
    double Aterm;
    double Xmin, Xmax;
    int Degree;

    s_Polynom() { MakeNaN(SoundSpeed); }
    const char* filename() const OVERRIDE { return "es_polynom.txt"; }
    void ReadParams(tFileBuffer& FB) OVERRIDE;
    void PointValue(double t, const double* coor, double* V) const;
};
//======================================================================================================================

//======================================================================================================================
struct s_4peak : tScalarFunction {
    DESCRIPTION("planar wave with 4 peaks of different form");
//======================================================================================================================
    double SoundSpeed;
    double Aterm;
    double Xmin, Xmax;
    double Period;

    s_4peak() { MakeNaN(Period); Aterm=1.0; Xmin=-1.0; Xmax=1.0; SoundSpeed=1.0; }
    const char* filename() const OVERRIDE { return "es_4peak.txt"; }
    void ReadParams(tFileBuffer& FB) OVERRIDE;
    void Init() OVERRIDE { if(IsNaN(Period)) Period = Xmax - Xmin; }
    void PointValue(double t, const double* coor, double* V) const;
};
//======================================================================================================================


//======================================================================================================================
struct s_Cylinder : tPulsFunction {
    DESCRIPTION("Gaussian impulse by cylinder diffraction");
//======================================================================================================================
    double X, Y, Radius; // Cylinder position and radius
    double Xterm, Yterm, Aterm, Bterm; // Gaussian impulse parameters

    int NumK, NumN; // Discretization over radial and angular modes
    double Kmax;

private:
    std::vector<double> Coeff; // Pre-generated coeffitients
    tGaussIntegrator<double> GI;

public:
    s_Cylinder() { 
        X = Y = 0.0; Radius = 0.5;  // geometry
        NumK = NumN = 200; MakeNaN(Kmax); // integration parameters
    }
    const char* filename() const OVERRIDE { return "es_cylinder.txt"; }
    void ReadParams(tFileBuffer& FB) OVERRIDE;
    void Init(void);
    void PrintCoeffMatrix(FILE* f);
    void PointValue(double t, const double* coor, double* V) const;
};
//======================================================================================================================

//======================================================================================================================
template<typename fpv>
struct s_IVP2D : tPointFunction_EP<fpv>, tSpaceForm<fpv> {
    DESCRIPTION("initial-value problem, using Green function");
//======================================================================================================================
    fpv FlowVelX, FlowVelY;
    fpv Hmax; // Source domain width divided by Bterm ( ==1 if cos^2 )
    int mm; // разбиение отрезков интегрирования
    tGaussIntegrator<fpv> GI; // узлы и коэффициенты гауссовой квадратуры

    tFuncType Type() const OVERRIDE { return FUNC_PULSATIONS; } // поскольку наследуем напрямую tPointFunction_EP, указываем тип

    s_IVP2D() { FlowVelX = FlowVelY = 0.0; Hmax=sqrt(9.2*sizeof(fpv)); GI.GR = 4*sizeof(fpv); mm = 3; }
    const char* filename() const OVERRIDE { return "es_ivp2d.txt"; }
    void ReadParams(tFileBuffer& FB) OVERRIDE;
    void Init() OVERRIDE;
    void PointValue(fpv t, const fpv* coor, fpv* V) const OVERRIDE;

protected:
    void CalcIntegralOverPhi(fpv x, fpv r, fpv phimax, fpv* u, int numcontour=1) const;
    void CalcIntegralOverZ(fpv z1, fpv z2, fpv phi, fpv t, fpv x, fpv* u, int numcontour=1) const;
    fpv GetSlope(fpv z, fpv t, fpv x) const;
    fpv FindIntermediatePoint(fpv zl, fpv zr, fpv slope, fpv t, fpv x, int numiters = 50) const;
    void PointValueAA(fpv t, const fpv* coor, fpv* u) const;
};
//======================================================================================================================

//======================================================================================================================
template<typename fpv>
struct s_IVP2DWP : s_IVP2D<fpv> {
    DESCRIPTION("initial-value problem, using Green function (wave potential)");
//======================================================================================================================
    tFuncType Type() const OVERRIDE { return FUNC_SCALAR; }
    void PointValue(fpv t, const fpv* coor, fpv* V) const OVERRIDE;
};
//======================================================================================================================

//======================================================================================================================
template<typename fpv>
struct s_Corner : tPointFunction_EP<fpv>, tSpaceForm<fpv> {
    DESCRIPTION("half-line or corner diffraction");
//======================================================================================================================
    fpv Hmax; // Source domain width divided by Bterm
    int mm; // разбиение отрезков интегрирования
    tGaussIntegrator<fpv> GI; // узлы и коэффициенты гауссовой квадратуры
    fpv angle; // угол (по умолчанию 2*Pi)
    s_Corner() { 
        angle = GetPiNumber2<fpv>(); 
        Hmax=sqrt(9.2*sizeof(fpv)); GI.GR = 4*sizeof(fpv); mm = 3;  
    }

    tFuncType Type() const OVERRIDE { return FUNC_PULSATIONS; } // поскольку наследуем напрямую tPointFunction_EP, указываем тип

    void CalcIntOverPhi_NEW(fpv phi1, fpv phi2, fpv t, fpv R, fpv r, fpv rs, fpv hatPhi, fpv* u) const;
    void CalcIntOverZ_NEW(fpv z1, fpv z2, fpv alpha, fpv t, fpv r, fpv rs, fpv hatPhi, fpv* u) const;
    void MainFunc(fpv t, const fpv* coord, fpv* u, int inv = -1, int clear = 1, fpv* uu = NULL) const;

    const char* filename() const OVERRIDE { return "es_corner.txt"; }
    void ReadParams(tFileBuffer& FB) OVERRIDE;
    void Init() OVERRIDE;
    void PointValue(fpv t, const fpv* coor, fpv* V) const OVERRIDE;
};
//======================================================================================================================

//======================================================================================================================
template<typename fpv>
struct s_CornerWP : s_Corner<fpv> {
    DESCRIPTION("half-line or corner diffraction (wave potential)");
//======================================================================================================================
    tFuncType Type() const OVERRIDE { return FUNC_SCALAR; }
    void PointValue(fpv t, const fpv* coor, fpv* V) const OVERRIDE;
};
//======================================================================================================================


//======================================================================================================================
struct s_CornerPlanar : tPulsFunction {
    DESCRIPTION("planar wave by corner diffraction");
//======================================================================================================================
    double angle; // угол
    double phi0;  // направление, откуда приходит волна
    double X0;    // расстояние до центра волны
    double Bterm; // полуширина волны
    double Aterm; // амплитуда волны
    int NumNodes; // количество узлов для интегрирования

private:
    int m, n; // рациональное представление для угла
    tCompoundGaussIntegrator<double> GI;

public:
    s_CornerPlanar() { angle = Pi2/3.0; phi0 = 0.25*PiNumber; X0=100.0; Bterm=6.0; Aterm = 1.0; 
                       m=n=0; NumNodes=21*sizeof(double); NumNodes *= sizeof(double); } 
    const char* filename() const OVERRIDE { return "es_cornerplanar.txt"; }
    void ReadParams(tFileBuffer& FB) OVERRIDE;
    void Init() OVERRIDE;
    void PointValue(double t, const double* coor, double* V) const;
    void set_reducer(double value) { GI.reducer = value; } // workaround

    // Вспомогательная процедура -- вычисление интегралов
    static double EvalJ0(double a, double b, const tCompoundGaussIntegrator<double>& CGI, int mode = -1);
};
//======================================================================================================================


//======================================================================================================================
// Протяжённый покоящийся точечный источник: 1D и 3D
template<typename fpv>
struct s_Source1D3D : tPointFunction_EP<fpv>, tSpaceForm<fpv> {
    tFuncType Type() const OVERRIDE { return FUNC_PULSATIONS; } // поскольку наследуем напрямую tPointFunction_EP, указываем тип
//======================================================================================================================
    fpv Ampl, Freq, Phase, tmin, tmax; // Параметры, описывающие сигнал во времени
    int SignalType; // 1 - sin, 6 - sin^4

    s_Source1D3D(int _NumCoords=3) { 
        Phase = 0.0; SignalType = 1; tmin = 0.0; tmax = 1e50;  Ampl = Freq = 1.0;
        tSpaceForm<fpv>::numCoords = _NumCoords;
    }
    void ReadParams(tFileBuffer& FB) OVERRIDE;
    fpv TimeForm(fpv t) const;
    fpv TimeFormDeriv(fpv t) const;
    fpv TimeFormDeriv2(fpv t) const;
    fpv TimeFormAntideriv(fpv t) const;
    void PointValue(fpv t, const fpv* coor, fpv* V) const = 0;
    void Init() OVERRIDE { tSpaceForm<fpv>::Init(); GI.Init(); }
protected:
    tGaussIntegrator<fpv> GI;
};
template<typename fpv>
struct s_Source1D : s_Source1D3D<fpv> { 
    DESCRIPTION("1D acoustic source");
    s_Source1D():s_Source1D3D<fpv>(1) {}
    const char* filename() const OVERRIDE { return "es_source1d.txt"; }
    void PointValue1(fpv t, fpv x, fpv* V) const;
    void PointValue(fpv t, const fpv* coor, fpv* V) const;
};
template<typename fpv>
struct s_Source3D : s_Source1D3D<fpv> {
    DESCRIPTION("3D acoustic source");
    s_Source3D():s_Source1D3D<fpv>(3) {}
    const char* filename() const OVERRIDE { return "es_source3d.txt"; }
    void PointValue1(fpv t, fpv r, fpv& p, fpv& urr) const;
    void PointValue(fpv t, const fpv* coor, fpv* V) const;
};
//======================================================================================================================


//======================================================================================================================
// Протяжённый покоящийся точечный источник: 2D
// В отличие от 1D и 3D, начальные условия ненулевые
struct s_Source2D : tPulsFunction, tSpaceForm<double>, tGaussIntegrator<double> {
    DESCRIPTION("2D acoustic source");
//======================================================================================================================
    double Ampl, Freq, Phase; // Параметры, описывающие сигнал во времени
    double Hmax; // Distance from the source center where Gaussian with unit half-width is negligible
    int mm;         // число отрезков в составной формуле Гаусса
    s_Source2D() { numCoords = 2; Phase = 0.0; Ampl = Freq = 1.0; GR = 100; mm=1; Hmax=8.5; }

    void ReadParams(tFileBuffer& FB) OVERRIDE;
    void Init() OVERRIDE { tSpaceForm<double>::Init(); tGaussIntegrator<double>::Init(GR); }
    const char* filename() const OVERRIDE { return "es_source2d.txt"; }
    void PointValue(double t, const double* coor, double* V) const;
};
//======================================================================================================================


//======================================================================================================================
template<typename fpv>
struct s_PointSource : tPointFunction_EP<fpv>, tSpaceForm<fpv> {
    tFuncType Type() const OVERRIDE { return FUNC_PULSATIONS; } // поскольку наследуем напрямую tPointFunction_EP, указываем тип
    DESCRIPTION("point source in a uniform [unsteady] flow");
//======================================================================================================================
    fpv SoundVel, FlowMach; // скорость звука и скорость фонового потока (вдоль оси X)
    fpv Ampl, Freq, Phase, tmin, tmax; // параметры сигнала. Freq - линейная частота
    int SignalType; // 1 : sin, 6 : sin^4

    int UseUnsteady /* =0, если поле стационарное */;
    fpv UnstAmpl, UnstFreq, UnstPhase; // параметры фонового потока, если он нестационарный

    s_PointSource() : tSpaceForm<fpv>() {
        SoundVel = 1.0; FlowMach = 0.0;
        Ampl = 1.0; Freq = 1.0; Phase = 0.0; tmin = 0.0; tmax = 1.0;
        SignalType = 6; // (sin(2*pi*f*t))^4
        UseUnsteady = 0;
        UnstAmpl = 0.0; UnstFreq = 1.0; UnstPhase = 0.0;
    }
    const char* filename() const OVERRIDE { return "es_pointsource.txt"; }
    void ReadParams(tFileBuffer& FB) OVERRIDE;

    fpv FindSourceTime(fpv time, fpv r0x, fpv r0y, fpv* re) const;
    fpv monopole(fpv* r0, fpv time, int diff) const;
    void PointValue1(fpv t, const fpv* coor, fpv* V) const; // Поле от одного источника
    void PointValue(fpv t, const fpv* coor, fpv* V) const OVERRIDE; // front-end
    void Init() OVERRIDE;
};
//======================================================================================================================


//======================================================================================================================
struct s_RotatingDipole : tPulsFunction {
    DESCRIPTION("rotating dipole");
//======================================================================================================================
    double Ampl, Omega, Phase;
    int dir;

    s_RotatingDipole() { Ampl=Omega=1.0; Phase=0.0; dir=2; }
    const char* filename() const OVERRIDE { return "es_rotatingdipole.txt"; }
    void ReadParams(tFileBuffer& FB) OVERRIDE;
    void PointValue(double t, const double* coor, double* V) const;
};
//======================================================================================================================



//======================================================================================================================
struct s_VortexInCylinder : tPulsFunction {
    DESCRIPTION("linear vortex in cylinder");
//======================================================================================================================
    double nu; // параметры уравнения - кинематическая вязкость
    double R; // радиус круга
    int CoorAxis; // направление, нормальное к кругу. Для 2D расчёта = Coor_Z, для RZ расчёта = Coor_X
    int NumModes; // количество мод
    std::string profile; // initial data: vphi(r/R)

private:
    std::vector<double> data; // нули функции Бесселя и амплитуда мод

public:
    s_VortexInCylinder() { CoorAxis = 2; R=1.0; nu=1.0; NumModes=25; profile="x*sin(pi*x)"; }

    const char* filename() const OVERRIDE { return "es_vortexincylinder.txt"; }
    void ReadParams(tFileBuffer& FB) OVERRIDE;
    void Init(void);
    void PointValue(double t, const double* coor, double* V) const;
};
//======================================================================================================================


//======================================================================================================================
template<typename fpv>
struct s_WaveInChannel_PrivateData; // Вспомогательный класс, содержащий комплексные переменные (частота и др.)

template<typename fpv>
struct s_WaveInChannel : tPointFunction_EP<fpv> {
    DESCRIPTION("wave in channel with viscosity");
//======================================================================================================================
    // equation parameters
    fpv nu;       // viscosity
    fpv gamma;    // specific ratio (unused if Pr=infinity)
    fpv Prandtl;  // Prandtl number
    // geometry parameters
    int form;     // 0 - planar channel, 1 - cylindrical channel
    fpv R;        // for cylindrical channel only: cylinder raduis
    fpv Ymin, Ymax; // for planar channel only: channel position
    int CoorAxis; // for cylindrical channel only: axial direction: 0 - X, 1 - Y, 2 - Z
                  // (For 2D computations set CoorAxis=2, when using ZR system set CoorAxis=0)
    // solution parameters
    fpv k;        // axial wave number
    int l;        // for cylindrical channel only: asimuthal mode
    int kmode;    // radial (for cylindrical channel) / transversal (for planar channel) mode 
    fpv ReOmega, ImOmega; // complex frequency. If set, kmode is ignored
    fpv ampl;     // multiplicator
    fpv phase;    // phase (axial shift)

    // additional parameters adjusting iterative process for solving the equation for 'omega'
    fpv _dmumax, _dmumin;  // minimal and maximal value of delta(mu)
    int loglevel; // уровень вывода отчёта об инициализации
    int ControlRootsJump; // enhanced control of switching from one mode to another

    // Data precalculated at Init()
    s_WaveInChannel_PrivateData<fpv>* data; 

    tFuncType Type() const OVERRIDE { return FUNC_PULSATIONS_COMPLEX; }

    // Constructors & destructors
    s_WaveInChannel() { 
        nu=1.0; gamma=1.4; Prandtl=1e50; 
        form=0; R=1.0; Ymin=-0.5; Ymax=0.5; CoorAxis = 2;
        k=0.0; l=8; kmode=1; ampl=1.0; phase=0.0; MakeNaN(ReOmega); MakeNaN(ImOmega);
        MakeNaN(_dmumin); MakeNaN(_dmumax); loglevel=1; ControlRootsJump=0;
        data=NULL;
    }
    s_WaveInChannel<fpv>& operator=(const s_WaveInChannel<fpv> &S);
    s_WaveInChannel(const s_WaveInChannel<fpv>& S) { data=NULL; *this = S; }
    void Free();
    ~s_WaveInChannel() { Free(); }

    // Reading or setting parameters
    const char* filename() const OVERRIDE { return "es_waveinchannel.txt"; }
    void SetEquationConstants(double _gam, double _Rey, double _Pr, double /*_SoundVel*/, double /*_FlowVel*/) OVERRIDE {
        gamma = fpv(_gam); nu = fpv(1.0 / _Rey); Prandtl = fpv(_Pr);
    }
    void ReadParams(tFileBuffer& FB) OVERRIDE;

    // Initialization
    int SolveEquation(void); // returns nonzero if error
    void Init(void) OVERRIDE;

    // Use
    void PointValue(fpv t, const fpv* coor, fpv* V) const;
};
//======================================================================================================================


//======================================================================================================================
struct s_AcousticShock : tPulsFunction { // Проход акустической волны через фронт ударной
    DESCRIPTION("acoustic wave through a shock propagation");
//======================================================================================================================
    double Xterm, Aterm, Bterm; // Acoustic wave parameters
    double MUL[3];         // Values left to the discontinuity
    double MUR[3];         // Values right to the discontinuity
    double DF_Coor;        // Initial discontinuity position
    double gam;

    s_AcousticShock(){
        DF_Coor=0.0; Xterm=-50.0; Aterm=1e-3; Bterm=1.0; gam = 1.4;
        for(int i=0; i<3; i++) { MakeNaN(MUL[i]); MakeNaN(MUR[i]); }
        MakeNaN(VShock);
    }
    const char* filename() const OVERRIDE { return "es_acousticshock.txt"; }
    void ReadParams(tFileBuffer& FB) OVERRIDE;
    void Init() OVERRIDE;
    void PointValue(double t, const double* coor, double* V) const OVERRIDE;
private:
    double VShock;
    void GetConstitutor(double* UR) const;
};
//======================================================================================================================


//======================================================================================================================
// === FULL NONLINEAR SOLUTIONS ========================================================================================
//======================================================================================================================


//======================================================================================================================
struct s_Riemann : tPhysFunction { // Решение задачи Римана о распаде произвольного разрыва
    DESCRIPTION("Riemann problem");
//======================================================================================================================
    double MUL[3];         // Values left to the discontinuity
    double MUR[3];         // Values right to the discontinuity
    double DF_Coor;        // Initial discontinuity position
    double e[3];           // Unit vector normal to front: 0 - X, 1 - Y, 2 - Z
    int SteadyTest;        // Do check that ushock=0 at initialization

    s_Riemann() {
        e[0]=1.0; e[1]=e[2]=0.0; SteadyTest=0; DF_Coor=0.0; 
        for(int i=0; i<3; i++) { MakeNaN(MUL[i]); MakeNaN(MUR[i]); }
        MakeNaN(gam);
    }
    inline s_Riemann(double _gam, const double* L, const double* R) { 
        e[0]=1.0; e[1]=e[2]=0.0; SteadyTest=0; DF_Coor=0.0;
        gam=_gam; for(int i=0; i<3; i++) { MUL[i]=L[i]; MUR[i]=R[i]; }
    }
    int IsShock(double* ushock) const; // Is the solution a shock; if so, returns shock velocity
    int IsUnsteady() const; // Is the solution steady?

    const char* filename() const OVERRIDE { return "es_riemann.txt"; }
    void ReadParams(tFileBuffer& FB) OVERRIDE;
    void Init() OVERRIDE;
    void PointValue(double t, const double* coor, double* V) const OVERRIDE;
};
// Exact solution for Riemann problem (input: conservative vars, output: physical vars)
void RMN_PrimVar(const double* UAL, const double* UAR, const double* normals, double* U0, double lambda, double gam, int Var_NumTurb = 0);
//======================================================================================================================


//======================================================================================================================
struct s_Shock : s_Riemann { // То же, что s_Riemann, только другие параметры считываются
    DESCRIPTION("Riemann problem (shock)");
//======================================================================================================================
    double Mach;

    s_Shock() { MUL[0]=MUR[0]=MUL[2]=MUR[2]=1.0; MUL[1]=MUR[1]=0.0; MakeNaN(Mach); }
    const char* filename() const OVERRIDE { return "es_shock.txt"; }
    void ReadParams(tFileBuffer& FB) OVERRIDE;
    void Init() OVERRIDE;
};
//======================================================================================================================


//======================================================================================================================
struct s_ShockRefl : tPhysFunction { // Задача об отражении ударной волны от твёрдой стенки (1D)
    DESCRIPTION("shock reflection problem");
//======================================================================================================================
    double WallPosition;
    double MUL[3];         // Values left to the discontinuity
    double MUR[3];         // Values right to the discontinuity. MUR[1]=0
    double DF_Coor;        // Initial discontinuity position
    s_ShockRefl() { WallPosition = 0.0; DF_Coor=0.0; for(int i=0; i<3; i++) { MakeNaN(MUL[i]); MakeNaN(MUR[i]); }}

    const char* filename() const OVERRIDE { return "es_shockrefl.txt"; }
    void ReadParams(tFileBuffer& FB) OVERRIDE;
    void Init() OVERRIDE;
    void PointValue(double t, const double* coor, double* V) const;
};
//======================================================================================================================


//======================================================================================================================
struct s_ViscShock : tPhysFunction { // Структура ударной волны
    DESCRIPTION("structure of shock wave");
//======================================================================================================================
    double x0, ushock;     // Shock position at t=0 and velocity
    double MUL[3], MUR[3]; // Infinity values: rho, u, p
    int mode;              // 0 - nu=const, 1 - mu=const
    int SteadyTest;        // Do check that ushock=0 at initialization
    int AutodetectShockPosition; // NOISEtte-only. For steady-state computations: autodetect x0 before calculating ES

    s_ViscShock() {
        mode=-1;
        x0=ushock=0.0; MUL[0]=MUR[0]=MUL[2]=MUR[2]=1.0; MUL[1]=MUR[1]=0.0; 
        SteadyTest=0;
        AutodetectShockPosition=0;
    }
    const char* filename() const OVERRIDE { return "es_viscshock.txt"; }
    void ReadParams(tFileBuffer& FB) OVERRIDE;
    void Init() OVERRIDE;
    void PointValue(double t, const double* coor, double* V) const;
};
//======================================================================================================================


//======================================================================================================================
struct s_CurlFreeCylinder : tPhysFunction {
    DESCRIPTION("irrotational flow around a cylinder");
//======================================================================================================================
    double X, Y, R; // Center coordinates of the cylinder and its radius
    double Gamma; // Circulation
    int PCEmode; // Set rho=1 instead of constant entropy (for pseudo-compressible equations)

    s_CurlFreeCylinder() { X = Y = 0.0; R = 0.5; Gamma = 0.0; PCEmode = 0; }
    const char* filename() const OVERRIDE { return "es_curlfreecylinder.txt"; }
    void ReadParams(tFileBuffer& FB) OVERRIDE;
    void PointValue(double t, const double* coor, double* V) const OVERRIDE;
    void Init(void) OVERRIDE;
};
//======================================================================================================================

//======================================================================================================================
struct s_PotentialSphere : tPhysFunction { // Потенциальное обтекание сферы
    DESCRIPTION("potential flow around a sphere");
//======================================================================================================================
    double X, Y, Z, R; // Координаты и радиус сферы
    int PCEmode; // pseudo-compressible equations

    s_PotentialSphere() { X = Y = Z = 0.0; R = 0.5; PCEmode = 0; }
    const char* filename() const OVERRIDE { return "es_potentialsphere.txt"; }
    void PointValue(double t, const double* coor, double* V) const OVERRIDE;
    void ReadParams(tFileBuffer& FB) OVERRIDE;
    void Init(void) OVERRIDE;
};
//======================================================================================================================


//======================================================================================================================
struct s_Couette : tPhysFunction { // Течение Куэтта между двумя параллельными плоскостями
    DESCRIPTION("Couette flow between two plates");
// Предполагается, что направление скорости по Y, а переменные меняются по X
//======================================================================================================================
    double xL, xR; // Координаты пластин
    double vL, vR; // Скорости пластин
    double tL, tR; // Температура стенок (для изотермических ГУ)
    double pressure; // Давление
    int condL, condR; // Условия на стенках: 0 -- адиабатическое, 1 -- изотермическое
    int ViscType;  // 0 - mu=const; 1 - mu=mu0/sqrt(T)
    #ifdef _NOISETTE
        int AutodetectPressure;
    #endif

    s_Couette() { 
        xL=-0.5; xR=0.5; condL=0; condR=1; gam=1.4; Pr=1.0; vL=0.0; vR=1.0; tL=tR=1.0/gam; ViscType=0; MakeNaN(pressure);
        #ifdef _NOISETTE
            AutodetectPressure=0; 
        #endif
    }
    const char* filename() const OVERRIDE { return "es_couette.txt"; }
    void ReadParams(tFileBuffer& FB) OVERRIDE;
    void Init();
    void PointValue(double t, const double* coor, double* V) const;
};
//======================================================================================================================


//======================================================================================================================
struct s_ViscSphere : tPhysFunction {
    DESCRIPTION("viscous uncompressible flow around sphere");
// Вязкое несжимаемое обтекание сферы при Re << 1, r << R/Re (см. Ландау, Лифшиц, гидродинамика)
//======================================================================================================================
    double R; // радиус сферы
    s_ViscSphere() { R = 0.5; }
    void ReadParams(tFileBuffer& FB) OVERRIDE;
    void PointValue(double t, const double* coor, double* V) const;
};
//======================================================================================================================


//======================================================================================================================
struct tVortexFunction : tPhysFunction {
    double r0[3]; // координаты центра вихря
    double v0[3]; // фоновый поток
    double Mach;  // "сила" вихря (обычно - отношение максимальной скорости к скорости звука на бесконечности)
    double R;     // "радиус" вихря (обычно - расстояние, на котором достигается максимальная скорость)
    int axis;     // ось вихря (0 - X, 1 - Y, 2 - Z)
    double Per[3]; // период по X, Y, Z. Используется только для финитного вихря; для остальных -- huge!

    tVortexFunction() { r0[0]=r0[1]=r0[2]=0.0; Mach=0.1; R=1.0; axis=2; Per[0]=Per[1]=Per[2]=huge; SoundVel=1.0; }
    void ReadParams(tFileBuffer& FB) OVERRIDE;
    void Init() OVERRIDE;
    void PointValue(double t, const double* coor, double* V) const OVERRIDE;
    virtual void Profile(double r, double& rho, double& u_over_r, double& p) const = 0; // собственно профиль вихря
};
//======================================================================================================================

//======================================================================================================================
struct s_RankineVortex : tVortexFunction {
    DESCRIPTION("Rankine vortex");
//======================================================================================================================
    const char* filename() const OVERRIDE { return "es_rankinevortex.txt"; }
    void Init() { SoundVel = 1.0; }
    void Profile(double r, double& rho, double& u_over_r, double& p) const OVERRIDE;
};
//======================================================================================================================

//======================================================================================================================
struct s_GaussianVortex : tVortexFunction {
    DESCRIPTION("vortex with Gaussian circulation profile"); // Oseen-Lamb vortex without dissipation
//======================================================================================================================
protected:
    static const double alpha0;
    static const double Mach2A;

public:
    const char* filename() const OVERRIDE { return "es_gaussianvortex.txt"; }
    void Profile(double r, double& rho, double& u_over_r, double& p) const OVERRIDE;
};
//======================================================================================================================

//======================================================================================================================
struct s_Vortex_BG : tVortexFunction {
    DESCRIPTION("Bosnyakov - Gadjiev vortex");
//======================================================================================================================
    const char* filename() const OVERRIDE { return "es_vortex_bg.txt"; }
    void Init() { SoundVel = 1.0; }
    void Profile(double r, double& rho, double& u_over_r, double& p) const OVERRIDE;
};
//======================================================================================================================

//======================================================================================================================
struct s_FiniteVortex : tVortexFunction {
    DESCRIPTION("finite vortex");
//======================================================================================================================
    int mode; // Режим: 1 - кубический профиль скорости, 0 - профиль степени 2d
    int deg;  // Степень полинома для скорости (для mode=0)

    s_FiniteVortex() : tVortexFunction() { mode=0; deg=3; }
    const char* filename() const OVERRIDE { return "es_finitevortex.txt"; }
    void Profile(double r, double& rho, double& u_over_r, double& p) const OVERRIDE;
    void ReadParams(tFileBuffer& FB) OVERRIDE;
    void Init() { SoundVel = 1.0; }
};
//======================================================================================================================


//======================================================================================================================
struct s_SimpleWave : tPhysFunction {
    DESCRIPTION("simple wave (preprint 2013-53 by Ladonkina et al.)");
    // Attention! The solution is valid until the shock occurs
//======================================================================================================================
    double l, x0, sign_FlowVel;
    s_SimpleWave() { x0 = 0.0; l = 0.2; sign_FlowVel = 1.0; }

    void Init(void) OVERRIDE;
    const char* filename() const OVERRIDE { return "es_simplewave.txt"; }
    void ReadParams(tFileBuffer& FB) OVERRIDE;
    void PointValue(double t, const double* coor, double* V) const;
};
//======================================================================================================================


//======================================================================================================================
struct s_ConcCyl : tPointFunction {
    DESCRIPTION("flow between two cylinders (2D, visc)");
    // ATTENTION! PointValue() returns the solution for velocity and temperature. 
    //            Density and pressure can't be defined uniquely!
//======================================================================================================================
    double X[3]; // a point on the axis
    double RIn, ROut; // cylinder radii
    double OmegaIn, OmegaOut; // angular velocities of cylinders
    double Pr, gam; // Prandtl number and the specific radio
    double Twall; // Wall temperature
    int CoorAxis; // axis direction: 0 - X, 1 - Y, 2 - Z. In NOISEtte, use Coor_Z for 2D and Coor_X for Coor_X

    s_ConcCyl() { X[0]=X[1]=X[2]=0.0; RIn=1.0; ROut=2.0; OmegaIn=OmegaOut=0.0; Pr=1.0; gam=1.4; CoorAxis=2; }
    tFuncType Type() const { return FUNC_TEMPVEL; }
    const char* filename() const OVERRIDE { return "es_conccyl.txt"; }
    void ReadParams(tFileBuffer& FB) OVERRIDE;
    void Init() OVERRIDE;
    void PointValue(double t, const double* coor, double* V) const;
};
//======================================================================================================================


//======================================================================================================================
struct s_Blasius : tPhysFunction {
    DESCRIPTION("solution of Blasius equation");
//======================================================================================================================
    s_Blasius() { FlowVel=1.0; SoundVel = 1e10; Rey=1.0; }
 
    void ReadParams(tFileBuffer& FB) OVERRIDE;
    const char* filename() const OVERRIDE { return "es_blasius.txt"; }
    void Init() OVERRIDE;
    void PointValue(double t, const double* coor, double* V) const OVERRIDE;

private:
    static const double alpha;
    std::vector<double> UV0, UV1;
    double d_eta;
public:
    inline void GetData(std::vector<double>*& pUV0, std::vector<double>*& pUV1, double*& pd_eta) {
        pUV0 = &UV0; pUV1 = &UV1; pd_eta = &d_eta;
    }
};
//======================================================================================================================

//======================================================================================================================
struct s_comprBlasius : tPhysFunction {
    DESCRIPTION("compressible modification of Blasius solution");
// (see Keldysh Intitute Preprint №15 1990)
//======================================================================================================================
    int mu_type;      // type of viscosity ( 0 = power law; 1 = Sutherland law )
    double ps;        // parameter of power law
    double sat;       // parameter of Sutherland law

    // BC parameters
    int bc_type;      // type of BC for h  ( 0 = adiabatic; 1 = isothermic )
    double temp_wall; // temperature value

    // parameters for calculation selfsimilar flow
    double L_sim;     //  mesh size 
    int N_sim;        //  number of intervals 

    // parameters for iterative procedure
    double eps_rel;   // relative accuracy
    double eps_abs;   // absolute accuracy
    int Iter_max;     // max. number of iterations

    s_comprBlasius() : mu_type(0), ps(0.5), sat(0.0), bc_type(0), temp_wall(2.),
        L_sim(5.), N_sim(1000), eps_rel(1e-7), eps_abs(1e-8), Iter_max(100) { FlowVel=SoundVel=1.0; gam=1.4; Pr=1.0; }
 
    void ReadParams(tFileBuffer& FB) OVERRIDE;
    const char* filename() const OVERRIDE { return "es_comprblasius.txt"; }
    void Init() OVERRIDE;
    void PointValue(double t, const double* coor, double* V) const OVERRIDE;

    int DumpTable(const char* fname = NULL) const;

private:
    std::vector<double> table_KSIUVT;
};
//======================================================================================================================

#endif
