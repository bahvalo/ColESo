#define _CRT_SECURE_NO_DEPRECATE 1
#define _CRT_SECURE_NO_WARNINGS  1
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

void coleso_add_function(char* FUNCNAME, int* ID);
void coleso_read_file(int ID, char* FILENAME);
void coleso_init(int ID);
void coleso_pointvalue(int ID, double T, double* C, double* V);

#define PiNumber 3.1415926535897932384626433

void LowerCase(char* name) {
    for(; *name; name++) if(*name >= 'A' && *name <= 'Z') name -= 32;
}

#define SOLUTION(X, Y, Z) \
    int ID; \
    printf("%s\n", Z); \
    coleso_add_function((char*)X, &ID);  \
    coleso_read_file(ID, (char*)Y); \
    coleso_init(ID); \
    out = fopen((const char*)Z, "wt"); \
    if(out==NULL) { printf("Error opening file %s for writing\n", Z); exit(0); }

int main(int argc, char** argv) { 
    double T, x, y, z, r, phi, C[3], V[10];
    FILE* out;
    char SolutionName[256] = "";
    if(argc > 1) strncpy(SolutionName, argv[1], 255);
    LowerCase(SolutionName);

    if(!SolutionName[0] || !strcmp(SolutionName, "4peak")) {
        SOLUTION("4peak", "PARAMS/es_4peak.txt", "DATA1D/4peak.dat");
        for(x = -1.0; x<1.0; x+=0.001) {
            T = 0.0; C[0]=x; C[1]=0.0; C[2]=0.0;
            coleso_pointvalue(ID, T, C, V);
            fprintf(out, "%25.15f %25.15f\n", x, V[0]);
        }
        fclose(out);
    }
    if(!SolutionName[0] || !strcmp(SolutionName, "planargauss")) {
        SOLUTION("planargauss", "PARAMS/es_planargauss.txt", "DATA1D/planargauss.dat");
        for(x = 43.0; x<=57.0000001; x+=0.01) {
            T = 25.0; C[0]=x; C[1]=0.0; C[2]=0.0;
            coleso_pointvalue(ID, T, C, V);
            fprintf(out, "%25.15f %25.15f\n", x, V[0]);
        }
        fclose(out);
    }
    if(!SolutionName[0] || !strcmp(SolutionName, "planarsinus")) {
        {
            SOLUTION("planarsinus", "PARAMS/es_planarsinus_1.txt", "DATA1D/planarsinus_1.dat");
            for(x = 0.0; x<=100.000001; x+=0.1) {
                T = 50.0; C[0]=x; C[1]=0.0; C[2]=0.0;
                coleso_pointvalue(ID, T, C, V);
                fprintf(out, "%25.15f %25.15f\n", x, V[0]);
            }
            fclose(out);
        }
        {
            SOLUTION("planarsinus", "PARAMS/es_planarsinus_2.txt", "DATA2D/planarsinus_2.dat");
            for(y = 0.0; y<=1.000000001; y+=0.02) {
                for(x = 0.0; x<=4.000000001; x+=0.02) {
                    T = 10.0; C[0]=x; C[1]=y; C[2]=0.0;
                    coleso_pointvalue(ID, T, C, V);
                    fprintf(out, "%25.15f %25.15f %25.15f\n", x, y, V[0]);
                }
                fprintf(out, "\n");
            }
            fclose(out);
        }
    }
    if(!SolutionName[0] || !strcmp(SolutionName, "entropyvortex")) {
        {
            SOLUTION("entropyvortex", "PARAMS/es_entropyvortex.txt", "DATA2D/entropyvortex_1.dat");
            for(y = 37.0; y<=64.000000001; y+=0.5) {
                for(x = 41; x<=68.000000001; x+=0.5) {
                    T = 5.0; C[0]=x; C[1]=y; C[2]=0.0;
                    coleso_pointvalue(ID, T, C, V);
                    fprintf(out, "%25.15f %25.15f %25.15f\n", x, y, V[1]);
                }
                fprintf(out, "\n");
            }
            fclose(out);
        }
        {
            SOLUTION("entropyvortex", "PARAMS/es_entropyvortex.txt", "DATA2D/entropyvortex_2.dat");
            for(y = 37.0; y<=64.000000001; y+=0.5) {
                for(x = 41; x<=68.000000001; x+=0.5) {
                    T = 5.0; C[0]=x; C[1]=y; C[2]=0.0;
                    coleso_pointvalue(ID, T, C, V);
                    fprintf(out, "%25.15f %25.15f %25.15f\n", x, y, V[0]);
                }
                fprintf(out, "\n");
            }
            fclose(out);
        }
    }
    if(!SolutionName[0] || !strcmp(SolutionName, "source1d")) {
        {
            SOLUTION("source1d", "PARAMS/es_source1d_1.txt", "DATA1D/source1d_1.dat");
            for(x = 0.0; x<=100.000001; x+=0.1) {
                T = 20.0; C[0]=x; C[1]=0.0; C[2]=0.0;
                coleso_pointvalue(ID, T, C, V);
                fprintf(out, "%25.15f %25.15f %25.15f\n", x, V[0], V[1]);
            }
            fclose(out);
        }
        {
            SOLUTION("source1d", "PARAMS/es_source1d_2.txt", "DATA1D/source1d_2.dat");
            for(x = 0.0; x<=100.000001; x+=0.1) {
                T = 20.0; C[0]=x; C[1]=0.0; C[2]=0.0;
                coleso_pointvalue(ID, T, C, V);
                fprintf(out, "%25.15f %25.15f %25.15f\n", x, V[0], V[1]);
            }
            fclose(out);
        }
    }
    if(!SolutionName[0] || !strcmp(SolutionName, "AcousticShock")) {
        double V1[5], V2[5], V3[5];
        SOLUTION("acousticshock", "PARAMS/es_acousticshock.txt", "DATA1D/acousticshock.dat");
        for(x = 0.0; x<=460.000001; x+=1.0) {
            C[0]=x; C[1]=0.0; C[2]=0.0;
            coleso_pointvalue(ID,   0.0, C, V1);
            coleso_pointvalue(ID,  35.0, C, V2);
            coleso_pointvalue(ID, 105.0, C, V3);
            fprintf(out, "%25.15f %25.15f %25.15f %25.15f\n", x, V1[0], V2[0], V3[0]);
        }
        fclose(out);
    }
    if(!SolutionName[0] || !strcmp(SolutionName, "gaussian2d")) {
        SOLUTION("gaussian2d", "PARAMS/es_gaussian2d.txt", "DATA2D/gaussian2d.dat");
        for(y = 0.; y<=34.000000001; y+=0.5) {
            for(x = 0.; x<=40.000000001; x+=0.5) {
                T = 30.0; C[0]=x; C[1]=y; C[2]=0.0;
                coleso_pointvalue(ID, T, C, V);
                fprintf(out, "%25.15f %25.15f %25.15f\n", x, y, V[0]);
            }
            fprintf(out, "\n");
        }
        fclose(out);
    }
    if(!SolutionName[0] || !strcmp(SolutionName, "source2d")) {
        SOLUTION("source2d", "PARAMS/es_source2d.txt", "DATA2D/source2d.dat");
        for(y = 0.; y<=100.000000001; y+=1.0) {
            for(x = 0.; x<=100.000000001; x+=1.0) {
                T = 12.0; C[0]=x; C[1]=y; C[2]=0.0;
                coleso_pointvalue(ID, T, C, V);
                fprintf(out, "%25.15f %25.15f %25.15f\n", x, y, V[0]);
            }
            fprintf(out, "\n");
        }
        fclose(out);
    }
    if(!SolutionName[0] || !strcmp(SolutionName, "gaussian3d")) {
        SOLUTION("gaussian3d", "PARAMS/es_gaussian3d.txt", "DATA2D/gaussian3d.dat");
        for(y = 0.; y<=34.000000001; y+=0.5) {
            for(x = 0.; x<=40.000000001; x+=0.5) {
                T = 30.0; C[0]=x; C[1]=y; C[2]=0.0;
                coleso_pointvalue(ID, T, C, V);
                fprintf(out, "%25.15f %25.15f %25.15f\n", x, y, V[0]);
            }
            fprintf(out, "\n");
        }
        fclose(out);
    }
    if(!SolutionName[0] || !strcmp(SolutionName, "source3d")) {
        SOLUTION("source3d", "PARAMS/es_source3d.txt", "DATA2D/source3d.dat");
        for(r = 0.; r<=50.000000001; r+=1.5) {
            for(x = -25.; x<=100; x+=1.5) {
                T = 25.0; C[0]=x; C[1]=r; C[2]=0.0;
                coleso_pointvalue(ID, T, C, V);
                fprintf(out, "%25.15f %25.15f %25.15f\n", x, r, V[1]);
            }
            fprintf(out, "\n");
        }
        fclose(out);
    }
    if(!SolutionName[0] || !strcmp(SolutionName, "pointsource")) {
        SOLUTION("pointsource", "PARAMS/es_pointsource.txt", "DATA2D/pointsource.dat");
        for(r = 0.; r<=150.000000001; r+=2.0) {
            for(x = 0.; x<=300; x+=2.) {
                T = 50.0; C[0]=x; C[1]=r; C[2]=0.0;
                coleso_pointvalue(ID, T, C, V);
                fprintf(out, "%25.15f %25.15f %25.15f\n", x, r, V[4]);
            }
            fprintf(out, "\n");
        }
        fclose(out);
    }
    if(!SolutionName[0] || !strcmp(SolutionName, "rotatingdipole")) {
        SOLUTION("rotatingdipole", "PARAMS/es_rotatingdipole.txt", "DATA2D/rotatingdipole.dat");
        for(y = -1.005; y<=1.000000001; y+=0.01) {
            for(x = -1.005; x<=1.00000001; x+=0.01) {
                T = 0.0; C[0]=x; C[1]=y; C[2]=0.0;
                coleso_pointvalue(ID, T, C, V);
                fprintf(out, "%25.15f %25.15f %25.15f\n", x, y, V[4]);
            }
            fprintf(out, "\n");
        }
        fclose(out);
    }
    if(!SolutionName[0] || !strcmp(SolutionName, "coaxial")) {
        {
            SOLUTION("coaxial", "PARAMS/es_coaxial_1.txt", "DATA2D/coaxial_1.dat");
            for(r = 0.3; r<=1.000000001; r+=0.01) {
                for(phi = -PiNumber; phi<=PiNumber+0.00000001; phi+=0.01*PiNumber) {
                    T = 1.0; C[0]=r*cos(phi); C[1]=r*sin(phi); C[2]=0.0;
                    coleso_pointvalue(ID, T, C, V);
                    fprintf(out, "%25.15f %25.15f %25.15f\n", C[0], C[1], V[0]);
                }
                fprintf(out, "\n");
            }
            fclose(out);
        }
        {
            SOLUTION("coaxial", "PARAMS/es_coaxial_2.txt", "DATA2D/coaxial_2.dat");
            for(r = 10.; r<=110.000000001; r+=1.) {
                for(z = 0.; z<=50.000000001; z+=1.) {
                    T = 17.0; C[0]=r; C[1]=0.0; C[2]=z;
                    coleso_pointvalue(ID, T, C, V);
                    fprintf(out, "%25.15f %25.15f %25.15f\n", r, z, V[0]);
                }
                fprintf(out, "\n");
            }
            fclose(out);
        }
    }
    if(!SolutionName[0] || !strcmp(SolutionName, "corner")) {
        {
            SOLUTION("corner", "PARAMS/es_corner_1.txt", "DATA2D/corner_1.dat");
            for(y = -42.; y<=66.; y+=2.) {
                for(x = -85.0; x<=110.00000001; x+=2.) {
                    T = 69.0; C[0]=x; C[1]=y; C[2]=0.0;
                    coleso_pointvalue(ID, T, C, V);
                    fprintf(out, "%25.15f %25.15f %25.15f\n", x, y, V[0]);
                }
                fprintf(out, "\n");
            }
            fclose(out);
        }
        {
            SOLUTION("corner", "PARAMS/es_corner_2.txt", "DATA2D/corner_2.dat");
            for(r = 0.; r<=1.000000001; r+=0.01) {
                double phimax;
                phimax = 2.*PiNumber/3.;
                for(phi = 0.; phi<=phimax+0.000001; phi+=0.01*phimax) {
                    T = 0.9; C[0]=r*cos(phi); C[1]=r*sin(phi); C[2]=0.0;
                    coleso_pointvalue(ID, T, C, V);
                    fprintf(out, "%25.15f %25.15f %25.15f\n", C[0], C[1], V[0]);
                }
                fprintf(out, "\n");
            }
            fclose(out);
        }
    }
    if(!SolutionName[0] || !strcmp(SolutionName, "cornerplanar")) {
        {
            SOLUTION("cornerplanar", "PARAMS/es_cornerplanar_1.txt", "DATA2D/cornerplanar_1.dat");
            for(y = -13.; y<=16.; y+=0.5) {
                for(x = -25.0; x<=20.00000001; x+=0.5) {
                    T = 60.0; C[0]=x; C[1]=y; C[2]=0.0;
                    coleso_pointvalue(ID, T, C, V);
                    fprintf(out, "%25.15f %25.15f %25.15f\n", x, y, V[0]);
                }
                fprintf(out, "\n");
            }
            fclose(out);
        }
        {
            SOLUTION("cornerplanar", "PARAMS/es_cornerplanar_2.txt", "DATA2D/cornerplanar_2.dat");
            for(r = 0.; r<=1.000000001; r+=0.01) {
                double phimax;
                phimax = 2.*PiNumber/3.;
                for(phi = 0.; phi<=phimax+0.000001; phi+=0.01*phimax) {
                    T = 0.4; C[0]=r*cos(phi); C[1]=r*sin(phi); C[2]=0.0;
                    coleso_pointvalue(ID, T, C, V);
                    fprintf(out, "%25.15f %25.15f %25.15f\n", C[0], C[1], V[0]);
                }
                fprintf(out, "\n");
            }
            fclose(out);
        }
    }
    if(!SolutionName[0] || !strcmp(SolutionName, "cylinder")) {
        SOLUTION("cylinder", "PARAMS/es_cylinder.txt", "DATA2D/cylinder.dat");
        for(r = 0.5; r<=10; r+=0.1) {
            double phimax;
            phimax = PiNumber;
            for(phi = 0.; phi<=phimax+0.000001; phi+=0.005*phimax) {
                T = 8.; C[0]=r*cos(phi); C[1]=r*sin(phi); C[2]=0.0;
                coleso_pointvalue(ID, T, C, V);
                fprintf(out, "%25.15f %25.15f %25.15f\n", C[0], C[1], V[0]);
            }
            fprintf(out, "\n");
        }
        fclose(out);
    }
    if(!SolutionName[0] || !strcmp(SolutionName, "sinusvisc")) {
        SOLUTION("sinusvisc", "PARAMS/es_sinusvisc.txt", "DATA2D/sinusvisc.dat");
        for(y = 0.0; y<=1.000000001; y+=0.01) {
            for(x = 0.0; x<=1.00000001; x+=0.01) {
                T = 0.0; C[0]=x; C[1]=y; C[2]=0.0;
                coleso_pointvalue(ID, T, C, V);
                fprintf(out, "%25.15f %25.15f %25.15f\n", x, y, V[0]);
            }
            fprintf(out, "\n");
        }
        fclose(out);
    }
    if(!SolutionName[0] || !strcmp(SolutionName, "vortexincylinder")) {
        double V1[5], V2[5];
        SOLUTION("vortexincylinder", "PARAMS/es_vortexincylinder.txt", "DATA1D/vortexincylinder.dat");
        for(x = 0.0; x<=1.000001; x+=0.01) {
            C[0]=x; C[1]=0.0; C[2]=0.0;
            coleso_pointvalue(ID,  0.0, C, V1);
            coleso_pointvalue(ID,  0.1, C, V2);
            fprintf(out, "%25.15f %25.15f %25.15f\n", x, V1[2], V2[2]);
        }
        fclose(out);
    }
    if(!SolutionName[0] || !strcmp(SolutionName, "waveinchannel")) {
        {
            SOLUTION("waveinchannel", "PARAMS/es_waveinchannel_1.txt", "DATA2D/waveinchannel_1.dat");
            for(y = 0.0; y<=0.5; y+=(0.02*(0.5-y)+0.001)) {
                for(x = -0.25; x<=0.2500000001; x+=0.01) {
                    T = 0.5; C[0]=x; C[1]=y; C[2]=0.0;
                    coleso_pointvalue(ID, T, C, V);
                    fprintf(out, "%25.15f %25.15f %25.15f\n", x, y, V[1]);
                }
                fprintf(out, "\n");
            }
            fclose(out);
        }
        {
            SOLUTION("waveinchannel", "PARAMS/es_waveinchannel_2.txt", "DATA2D/waveinchannel_2.dat");
            for(y = 0.0; y<=0.5; y+=0.01) {
                for(x = -0.25; x<=0.2500000001; x+=0.01) {
                    T = 1.0; C[0]=x; C[1]=y; C[2]=0.0;
                    coleso_pointvalue(ID, T, C, V);
                    fprintf(out, "%25.15f %25.15f %25.15f\n", x, y, V[1]);
                }
                fprintf(out, "\n");
            }
            fclose(out);
        }
        {
            SOLUTION("waveinchannel", "PARAMS/es_waveinchannel_3.txt", "DATA1D/waveinchannel_3.dat");
            for(r = 0.0; r<=1.0; r+=0.001) {
                T = 0.1; C[0]=0.0; C[1]=r; C[2]=0.0;
                coleso_pointvalue(ID, T, C, V);
                fprintf(out, "%25.15f %25.15f %25.15f %25.15f\n", r, V[0], V[2], V[4]);
            }
            fclose(out);
        }
        {
            SOLUTION("waveinchannel", "PARAMS/es_waveinchannel_4.txt", "DATA2D/waveinchannel_4.dat");
            for(r = 0.0; r<=1.000000001; r+=0.01*(1-r)+1e-5) {
                for(phi = 0; phi<=PiNumber/4.+0.0000001; phi+=0.01) {
                    T = 0.0; C[0]=r*cos(phi); C[1]=r*sin(phi); C[2]=0.0;
                    coleso_pointvalue(ID, T, C, V);
                    fprintf(out, "%25.15f %25.15f %25.15f\n", C[0], C[1], V[0]);
                }
                fprintf(out, "\n");
            }
            fclose(out);
        }
    }    
    if(!SolutionName[0] || !strcmp(SolutionName, "rankinevortex") || !strcmp(SolutionName, "gaussianvortex") || !strcmp(SolutionName, "finitevortex")) {
        int ID_F, ID_R, ID_G;
        double VF[5], VG[5], VR[5];
        printf("%s\n", "vortexes.dat");
        coleso_add_function((char*)"rankinevortex", &ID_R);
        coleso_add_function((char*)"gaussianvortex", &ID_G);
        coleso_add_function((char*)"finitevortex", &ID_F);
        coleso_read_file(ID_R, (char*)"PARAMS/es_rankinevortex.txt");
        coleso_read_file(ID_F, (char*)"PARAMS/es_finitevortex.txt");
        coleso_read_file(ID_G, (char*)"PARAMS/es_gaussianvortex.txt");
        coleso_init(ID_R);
        coleso_init(ID_F);
        coleso_init(ID_G);
        out = fopen((const char*)"DATA1D/vortexes.dat", "wt");
        if(out==NULL) { printf("Error opening file %s for writing\n", "DATA1D/vortexes.dat"); exit(0); }
        for(r = 0.0; r<=8.0; r+=0.1) {
            T = 0.0; C[0]=r; C[1]=0.0; C[2]=0.0;
            coleso_pointvalue(ID_R, T, C, VR);
            coleso_pointvalue(ID_G, T, C, VG);
            coleso_pointvalue(ID_F, T, C, VF);
            fprintf(out, "%25.15f %25.15f %25.15f %25.15f %25.15f %25.15f %25.15f\n", r, VF[2], VG[2], VR[2], VF[4], VG[4], VR[4]);
        }
        fclose(out);
    }
    if(!SolutionName[0] || !strcmp(SolutionName, "riemann")) {
        SOLUTION("riemann", "PARAMS/es_riemann.txt", "DATA1D/riemann.dat");
        for(x = 175.0; x<=310.000001; x+=1.0) {
            T = 4.0; C[0]=x; C[1]=0.0; C[2]=0.0;
            coleso_pointvalue(ID, T, C, V);
            fprintf(out, "%25.15f %25.15f\n", x, V[0]);
        }
        fclose(out);
    }
    if(!SolutionName[0] || !strcmp(SolutionName, "simplewave")) {
        double V1[5], V2[5];
        SOLUTION("simplewave", "PARAMS/es_simplewave.txt", "DATA1D/simplewave.dat");
        for(x = -0.7; x<=0.5000001; x+=0.01) {
            C[0]=x; C[1]=0.0; C[2]=0.0;
            coleso_pointvalue(ID,  0.0, C, V1);
            coleso_pointvalue(ID,  0.09, C, V2);
            fprintf(out, "%25.15f %25.15f %25.15f %25.15f %25.15f\n", x, V1[0], V1[4], V2[0], V2[4]);
        }
        fclose(out);
    }
    if(!SolutionName[0] || !strcmp(SolutionName, "viscshock")) {
        {
            SOLUTION("viscshock", "PARAMS/es_viscshock_1.txt", "DATA1D/viscshock_1.dat");
            for(x = -20.0; x<=20.0; x+=0.1) {
                T = 0.0; C[0]=x; C[1]=0.0; C[2]=0.0;
                coleso_pointvalue(ID, T, C, V);
                fprintf(out, "%25.15f %25.15f %25.15f %25.15f\n", x, V[0], V[1], V[4]);
            }
            fclose(out);
        }
        {
            SOLUTION("viscshock", "PARAMS/es_viscshock_2.txt", "DATA1D/viscshock_2.dat");
            for(x = -20.0; x<=20.0; x+=0.1) {
                T = 0.0; C[0]=x; C[1]=0.0; C[2]=0.0;
                coleso_pointvalue(ID, T, C, V);
                fprintf(out, "%25.15f %25.15f %25.15f %25.15f\n", x, V[0], V[1], V[4]);
            }
            fclose(out);
        }
        {
            SOLUTION("viscshock", "PARAMS/es_viscshock_3.txt", "DATA1D/viscshock_3.dat");
            for(x = -20.0; x<=20.0; x+=0.1) {
                T = 0.0; C[0]=x; C[1]=0.0; C[2]=0.0;
                coleso_pointvalue(ID, T, C, V);
                fprintf(out, "%25.15f %25.15f %25.15f %25.15f\n", x, V[0], V[1], V[4]);
            }
            fclose(out);
        }
    }
    if(!SolutionName[0] || !strcmp(SolutionName, "conccyl")) {
        SOLUTION("conccyl", "PARAMS/es_conccyl.txt", "DATA2D/conccyl.dat");
        for(phi = PiNumber/4.0; phi<=7*PiNumber/12.+0.000000001; phi+=0.02*PiNumber/3) {
            for(r = 1.0; r<=2.00000001; r+=0.01) {
                T = 0.0; C[0]=r*cos(phi); C[1]=r*sin(phi); C[2]=0.0;
                coleso_pointvalue(ID, T, C, V);
                fprintf(out, "%25.15f %25.15f %25.15f\n", C[0], C[1], V[4]);
            }
            fprintf(out, "\n");
        }
        fclose(out);
    }
    if(!SolutionName[0] || !strcmp(SolutionName, "CurlFreecylinder")) {
        SOLUTION("curlfreecylinder", "PARAMS/es_curlfreecylinder.txt", "DATA2D/curlfreecylinder.dat");
        for(phi = 0; phi<=PiNumber+0.0000001; phi+=PiNumber/50) {
            for(r = 0.5; r<=5.00000001; r+=0.01) {
                T = 0.0; C[0]=r*cos(phi); C[1]=r*sin(phi); C[2]=0.0;
                coleso_pointvalue(ID, T, C, V);
                fprintf(out, "%25.15f %25.15f %25.15f\n", C[0], C[1], V[4]);
            }
            fprintf(out, "\n");
        }
        fclose(out);
    }
    if(!SolutionName[0] || !strcmp(SolutionName, "potentialsphere")) {
        SOLUTION("potentialsphere", "PARAMS/es_potentialsphere.txt", "DATA2D/potentialsphere.dat");
        for(phi = 0; phi<=PiNumber+0.0000001; phi+=PiNumber/50) {
            for(r = 0.5; r<=5.00000001; r+=0.025) {
                T = 0.0; C[0]=r*cos(phi); C[1]=r*sin(phi); C[2]=0.0;
                coleso_pointvalue(ID, T, C, V);
                fprintf(out, "%25.15f %25.15f %25.15f\n", C[0], C[1], V[4]);
            }
            fprintf(out, "\n");
        }
        fclose(out);
    }
    if(!SolutionName[0] || !strcmp(SolutionName, "viscsphere")) {
        SOLUTION("viscsphere", "PARAMS/es_viscsphere.txt", "DATA2D/viscsphere.dat");
        for(phi = 0; phi<=PiNumber+0.0000001; phi+=PiNumber/50) {
            for(r = 0.5; r<=5.00000001; r+=0.025) {
                T = 0.0; C[0]=r*cos(phi); C[1]=r*sin(phi); C[2]=0.0;
                coleso_pointvalue(ID, T, C, V);
                fprintf(out, "%25.15f %25.15f %25.15f\n", C[0], C[1], V[4]);
            }
            fprintf(out, "\n");
        }
        fclose(out);
    }
    return 0;
}
