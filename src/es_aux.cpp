// *********************************************************************************************************************
// *****                                                                                                           *****
// *****                                  Collection of exact solutions (ColESo)                                   *****
// *****                            Object-free interface (for use in FORTRAN programs)                            *****
// *****                                                                                                           *****
// *********************************************************************************************************************

#include "base_parser.h"
#include "coleso.h"

//======================================================================================================================
// Create object by name
//======================================================================================================================
#define A(N) if(CompareWords(string(FuncName), string(#N))) return (tPointFunction*) new s_##N
tPointFunction* CreateFunction(const char* FuncName) {
    A(Gaussian2D)<double>;
    A(Gaussian3D);
    A(EntropyVortex);
    A(PlanarSinus);
    A(PlanarGauss);
    A(Coaxial);
    A(SinusVisc);
    A(Polynom);
    A(4peak);
    A(Cylinder);
    A(IVP2D)<double>;
    A(Corner)<double>;
    A(CornerPlanar);
    A(Source1D)<double>;
    A(Source3D)<double>;
    A(Source2D);
    A(PointSource)<double>;
    A(RotatingDipole);
    A(VortexInCylinder);
    A(WaveInChannel)<double>;
    A(AcousticShock);
    A(Riemann);
    A(Shock);
    A(ShockRefl);
    A(ViscShock);
    A(CurlFreeCylinder);
    A(PotentialSphere);
    A(Couette);
    A(ViscSphere);
    A(RankineVortex);
    A(GaussianVortex);
    A(Vortex_BG);
    A(FiniteVortex);
    A(SimpleWave);
    A(Blasius);
    A(ConcCyl);
    return NULL;
}

//======================================================================================================================
// Object-free interface (for use in FORTRAN programs)
//======================================================================================================================
static vector<tPointFunction*> VPF;
extern "C" {
    void coleso_add_function(char* FUNCNAME, int* ID) {
        tPointFunction* F = CreateFunction(FUNCNAME);
        if(F==NULL) { 
            if(ID) *ID=-1;
        }
        else {
            if(ID) *ID = (int)VPF.size();
            VPF.push_back(F);
        }
    };
    void coleso_read_file(int ID, char* FILENAME) {
        if(ID<0 || ID>=int(VPF.size())) return;
        VPF[ID]->ReadParamsFromFile(FILENAME);
    }
    void coleso_set_parameter(char* PARAMNAME, char* PARAMVALUE) {
        string S = ":" + string(PARAMNAME) + "=" + string(PARAMVALUE);
        CmdLine.AddArgument((char*)(S.c_str()));
    };
    void coleso_read_set(int ID) {
        if(ID<0 || ID>=int(VPF.size())) return;
        tFileBuffer FB;
        VPF[ID]->ReadParams(FB);
        CmdLine.ClearFileArguments();
    }
    void coleso_init(int ID) {
        if(ID<0 || ID>=int(VPF.size())) return;
        VPF[ID]->Init();
    }
    void coleso_pointvalue(int ID, double T, double* C, double* V) {
        if(ID<0 || ID>=int(VPF.size())) return;
        VPF[ID]->PointValue(T, C, V);
    }
}
