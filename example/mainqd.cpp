// Attention: to build this code, ColESo library should be compiled with 
// EXTRAPRECISION_COLESO flag set in src/base_config.h
#include "coleso.h"
#include "extraprecision.h"

int main(int, char**) {
    s_WaveInChannel<dd_real> S;
    S.CoorAxis = 2;
    S.k = 5.0;
    S.kmode = 1;
    S.l = 1;
    S.nu = 1e-4;
    S.R = 1.0;
    S.gamma = 1.5;
    S.Prandtl = 1.0;
    S.form = 1;
    S._dmumax = 0.02;
    S.Init();

    dd_real C[3]  = {0.9, 0.0, 0.0};
    dd_real time  = 0.0;
    dd_real V[10];

    FILE* out = fopen("data_dd.dat", "wt");
    const dd_real dx = dd_real(1.0) / dd_real(1e4);
    for(C[0]=0.0; C[0]<=1.0; C[0]+=dx) {
        S.PointValue(time, C, V);
        std::string res = C[0].to_string() + " " + V[0].to_string();
        fprintf(out, "%s\n", res.c_str());
    }
    fclose(out);
}
