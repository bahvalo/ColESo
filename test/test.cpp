#include "coleso.h"
#include "base_parser.h"
using namespace std;

void ColESo_Tests(int argc, char** argv);

#define SOLUTION(X, Y, Z) \
    X S; \
    S.ReadParamsFromFile(#Y); \
    S.Init(); \
    pprintf0("%s\n", Z); \
    FILE* out = fopen((const char*)Z, "wt"); \
    if(out==NULL) crash("Error opening file %s for writing\n", Z);

int main(int argc, char** argv) {
    double V[10];
    string SolutionName = "";
    if(argc > 1) SolutionName = argv[1];
    if(CompareWords(SolutionName, "test")) {
        ColESo_Tests(argc-1, argv+1);
        return 0;
    }

    if(SolutionName.empty() || CompareWords(SolutionName, "4peak")) {
        SOLUTION(s_4peak, PARAMS/es_4peak.txt, "DATA1D/4peak.dat");
        for(double x = -1.0; x<1.0; x+=0.001) {
            double T = 0.0, C[3] = {x, 0.0, 0.0};
            S.PointValue(T, C, V);
            fprintf(out, "%25.15f %25.15f\n", x, V[0]);
        }
        fclose(out);
    }
    if(SolutionName.empty() || CompareWords(SolutionName, "planargauss")) {
        SOLUTION(s_PlanarGauss, PARAMS/es_planargauss.txt, "DATA1D/planargauss.dat");
        for(double x = 43.0; x<=57.0000001; x+=0.01) {
            double T = 25.0, C[3] = {x, 0.0, 0.0};
            S.PointValue(T, C, V);
            fprintf(out, "%25.15f %25.15f\n", x, V[0]);
        }
        fclose(out);
    }
    if(SolutionName.empty() || CompareWords(SolutionName, "planarsinus")) {
        {
            SOLUTION(s_PlanarSinus, PARAMS/es_planarsinus_1.txt, "DATA1D/planarsinus_1.dat");
            for(double x = 0.0; x<=100.000001; x+=0.1) {
                double T = 50.0, C[3] = {x, 0.0, 0.0};
                S.PointValue(T, C, V);
                fprintf(out, "%25.15f %25.15f\n", x, V[0]);
            }
            fclose(out);
        }
        {
            SOLUTION(s_PlanarSinus, PARAMS/es_planarsinus_2.txt, "DATA2D/planarsinus_2.dat");
            for(double y = 0.0; y<=1.000000001; y+=0.02) {
                for(double x = 0.0; x<=4.000000001; x+=0.02) {
                    double T = 10.0, C[3] = {x, y, 0.0};
                    S.PointValue(T, C, V);
                    fprintf(out, "%25.15f %25.15f %25.15f\n", x, y, V[0]);
                }
                fprintf(out, "\n");
            }
            fclose(out);
        }
    }
    if(SolutionName.empty() || CompareWords(SolutionName, "entropyvortex")) {
        {
            SOLUTION(s_EntropyVortex, PARAMS/es_entropyvortex.txt, "DATA2D/entropyvortex_1.dat");
            for(double y = 37.0; y<=64.000000001; y+=0.5) {
                for(double x = 41; x<=68.000000001; x+=0.5) {
                    double T = 5.0, C[3] = {x, y, 0.0};
                    S.PointValue(T, C, V);
                    fprintf(out, "%25.15f %25.15f %25.15f\n", x, y, V[1]);
                }
                fprintf(out, "\n");
            }
            fclose(out);
        }
        {
            SOLUTION(s_EntropyVortex, PARAMS/es_entropyvortex.txt, "DATA2D/entropyvortex_2.dat");
            for(double y = 37.0; y<=64.000000001; y+=0.5) {
                for(double x = 41; x<=68.000000001; x+=0.5) {
                    double T = 5.0, C[3] = {x, y, 0.0};
                    S.PointValue(T, C, V);
                    fprintf(out, "%25.15f %25.15f %25.15f\n", x, y, V[0]);
                }
                fprintf(out, "\n");
            }
            fclose(out);
        }
    }
    if(SolutionName.empty() || CompareWords(SolutionName, "source1d")) {
        {
            SOLUTION(s_Source1D<double>, PARAMS/es_source1d_1.txt, "DATA1D/source1d_1.dat");
            for(double x = 0.0; x<=100.000001; x+=0.1) {
                double T = 20.0, C[3] = {x, 0.0, 0.0};
                S.PointValue(T, C, V);
                fprintf(out, "%25.15f %25.15f %25.15f\n", x, V[0], V[1]);
            }
            fclose(out);
        }
        {
            SOLUTION(s_Source1D<double>, PARAMS/es_source1d_2.txt, "DATA1D/source1d_2.dat");
            for(double x = 0.0; x<=100.000001; x+=0.1) {
                double T = 20.0, C[3] = {x, 0.0, 0.0};
                S.PointValue(T, C, V);
                fprintf(out, "%25.15f %25.15f %25.15f\n", x, V[0], V[1]);
            }
            fclose(out);
        }
    }
    if(SolutionName.empty() || CompareWords(SolutionName, "AcousticShock")) {
        SOLUTION(s_AcousticShock, PARAMS/es_acousticshock.txt, "DATA1D/acousticshock.dat");
        double V1[5], V2[5], V3[5];
        for(double x = 0.0; x<=460.000001; x+=1.0) {
            double C[3] = {x, 0.0, 0.0};
            S.PointValue(  0.0, C, V1);
            S.PointValue( 35.0, C, V2);
            S.PointValue(105.0, C, V3);
            fprintf(out, "%25.15f %25.15f %25.15f %25.15f\n", x, V1[0], V2[0], V3[0]);
        }
        fclose(out);
    }
    if(SolutionName.empty() || CompareWords(SolutionName, "gaussian2d")) {
        SOLUTION(s_Gaussian2D<double>, PARAMS/es_gaussian2d.txt, "DATA2D/gaussian2d.dat");
        for(double y = 0.; y<=34.000000001; y+=0.5) {
            for(double x = 0.; x<=40.000000001; x+=0.5) {
                double T = 30.0, C[3] = {x, y, 0.0};
                S.PointValue(T, C, V);
                fprintf(out, "%25.15f %25.15f %25.15f\n", x, y, V[0]);
            }
            fprintf(out, "\n");
        }
        fclose(out);
    }
    if(SolutionName.empty() || CompareWords(SolutionName, "source2d")) {
        SOLUTION(s_Source2D, PARAMS/es_source2d.txt, "DATA2D/source2d.dat");
        for(double y = 0.; y<=100.000000001; y+=1.0) {
            for(double x = 0.; x<=100.000000001; x+=1.0) {
                double T = 12.0, C[3] = {x, y, 0.0};
                S.PointValue(T, C, V);
                fprintf(out, "%25.15f %25.15f %25.15f\n", x, y, V[0]);
            }
            fprintf(out, "\n");
        }
        fclose(out);
    }
    if(SolutionName.empty() || CompareWords(SolutionName, "gaussian3d")) {
        SOLUTION(s_Gaussian3D, PARAMS/es_gaussian3d.txt, "DATA2D/gaussian3d.dat");
        for(double y = 0.; y<=34.000000001; y+=0.5) {
            for(double x = 0.; x<=40.000000001; x+=0.5) {
                double T = 30.0, C[3] = {x, y, 0.0};
                S.PointValue(T, C, V);
                fprintf(out, "%25.15f %25.15f %25.15f\n", x, y, V[0]);
            }
            fprintf(out, "\n");
        }
        fclose(out);
    }
    if(SolutionName.empty() || CompareWords(SolutionName, "source3d")) {
        SOLUTION(s_Source3D<double>, PARAMS/es_source3d.txt, "DATA2D/source3d.dat");
        for(double r = 0.; r<=50.000000001; r+=1.5) {
            for(double x = -25.; x<=100; x+=1.5) {
                double T = 25.0, C[3] = {x, r, 0.0};
                S.PointValue(T, C, V);
                fprintf(out, "%25.15f %25.15f %25.15f\n", x, r, V[1]);
            }
            fprintf(out, "\n");
        }
        fclose(out);
    }
    if(SolutionName.empty() || CompareWords(SolutionName, "pointsource")) {
        SOLUTION(s_PointSource<double>, PARAMS/es_pointsource.txt, "DATA2D/pointsource.dat");
        for(double r = 0.; r<=150.000000001; r+=2.0) {
            for(double x = 0.; x<=300; x+=2.) {
                double T = 50.0, C[3] = {x, r, 0.0};
                S.PointValue(T, C, V);
                fprintf(out, "%25.15f %25.15f %25.15f\n", x, r, V[4]);
            }
            fprintf(out, "\n");
        }
        fclose(out);
    }
    if(SolutionName.empty() || CompareWords(SolutionName, "rotatingdipole")) {
        SOLUTION(s_RotatingDipole, PARAMS/es_rotatingdipole.txt, "DATA2D/rotatingdipole.dat");
        for(double y = -1.005; y<=1.000000001; y+=0.01) {
            for(double x = -1.005; x<=1.00000001; x+=0.01) {
                double T = 0.0, C[3] = {x, y, 0.0};
                S.PointValue(T, C, V);
                fprintf(out, "%25.15f %25.15f %25.15f\n", x, y, V[4]);
            }
            fprintf(out, "\n");
        }
        fclose(out);
    }
    if(SolutionName.empty() || CompareWords(SolutionName, "coaxial")) {
        {
            SOLUTION(s_Coaxial, PARAMS/es_coaxial_1.txt, "DATA2D/coaxial_1.dat");
            for(double r = 0.3; r<=1.000000001; r+=0.01) {
                for(double phi = -PiNumber; phi<=PiNumber+0.00000001; phi+=0.01*PiNumber) {
                    double T = 1.0, C[3] = {r*cos(phi), r*sin(phi), 0.0};
                    S.PointValue(T, C, V);
                    fprintf(out, "%25.15f %25.15f %25.15f\n", C[0], C[1], V[0]);
                }
                fprintf(out, "\n");
            }
            fclose(out);
        }
        {
            SOLUTION(s_Coaxial, PARAMS/es_coaxial_2.txt, "DATA2D/coaxial_2.dat");
            for(double r = 10.; r<=110.000000001; r+=1.) {
                for(double z = 0.; z<=50.000000001; z+=1.) {
                    double T = 17.0, C[3] = {r, 0.0, z};
                    S.PointValue(T, C, V);
                    fprintf(out, "%25.15f %25.15f %25.15f\n", r, z, V[0]);
                }
                fprintf(out, "\n");
            }
            fclose(out);
        }
    }
    if(SolutionName.empty() || CompareWords(SolutionName, "corner")) {
        {
            SOLUTION(s_Corner<double>, PARAMS/es_corner_1.txt, "DATA2D/corner_1.dat");
            for(double y = -42.; y<=66.; y+=2.) {
                for(double x = -85.0; x<=110.00000001; x+=2.) {
                    double T = 69.0, C[3] = {x, y, 0.0};
                    S.PointValue(T, C, V);
                    fprintf(out, "%25.15f %25.15f %25.15f\n", x, y, V[0]);
                }
                fprintf(out, "\n");
            }
            fclose(out);
        }
        {
            SOLUTION(s_Corner<double>, PARAMS/es_corner_2.txt, "DATA2D/corner_2.dat");
            for(double r = 0.; r<=1.000000001; r+=0.01) {
                double phimax = 2.*PiNumber/3.;
                for(double phi = 0.; phi<=phimax+0.000001; phi+=0.01*phimax) {
                    double T = 0.9, C[3] = {r*cos(phi), r*sin(phi), 0.0};
                    S.PointValue(T, C, V);
                    fprintf(out, "%25.15f %25.15f %25.15f\n", C[0], C[1], V[0]);
                }
                fprintf(out, "\n");
            }
            fclose(out);
        }
    }
    if(SolutionName.empty() || CompareWords(SolutionName, "cornerplanar")) {
        {
            SOLUTION(s_CornerPlanar<double>, PARAMS/es_cornerplanar_1.txt, "DATA2D/cornerplanar_1.dat");
            for(double y = -13.; y<=16.; y+=0.5) {
                for(double x = -25.0; x<=20.00000001; x+=0.5) {
                    double T = 60.0, C[3] = {x, y, 0.0};
                    S.PointValue(T, C, V);
                    fprintf(out, "%25.15f %25.15f %25.15f\n", x, y, V[0]);
                }
                fprintf(out, "\n");
            }
            fclose(out);
        }
        {
            SOLUTION(s_CornerPlanar<double>, PARAMS/es_cornerplanar_2.txt, "DATA2D/cornerplanar_2.dat");
            for(double r = 0.; r<=1.000000001; r+=0.01) {
                double phimax = 2.*PiNumber/3.;
                for(double phi = 0.; phi<=phimax+0.000001; phi+=0.01*phimax) {
                    double T = 0.4, C[3] = {r*cos(phi), r*sin(phi), 0.0};
                    S.PointValue(T, C, V);
                    fprintf(out, "%25.15f %25.15f %25.15f\n", C[0], C[1], V[0]);
                }
                fprintf(out, "\n");
            }
            fclose(out);
        }
    }
    if(SolutionName.empty() || CompareWords(SolutionName, "cylinder")) {
        SOLUTION(s_Cylinder, PARAMS/es_cylinder.txt, "DATA2D/cylinder.dat");
        for(double r = 0.5; r<=10; r+=0.1) {
            double phimax = PiNumber;
            for(double phi = 0.; phi<=phimax+0.000001; phi+=0.005*phimax) {
                double T = 8., C[3] = {r*cos(phi), r*sin(phi), 0.0};
                S.PointValue(T, C, V);
                fprintf(out, "%25.15f %25.15f %25.15f\n", C[0], C[1], V[0]);
            }
            fprintf(out, "\n");
        }
        fclose(out);
    }
    if(SolutionName.empty() || CompareWords(SolutionName, "sinusvisc")) {
        SOLUTION(s_SinusVisc, PARAMS/es_sinusvisc.txt, "DATA2D/sinusvisc.dat");
        for(double y = 0.0; y<=1.000000001; y+=0.01) {
            for(double x = 0.0; x<=1.00000001; x+=0.01) {
                double T = 0.0, C[3] = {x, y, 0.0};
                S.PointValue(T, C, V);
                fprintf(out, "%25.15f %25.15f %25.15f\n", x, y, V[0]);
            }
            fprintf(out, "\n");
        }
        fclose(out);
    }
    if(SolutionName.empty() || CompareWords(SolutionName, "vortexincylinder")) {
        SOLUTION(s_VortexInCylinder, PARAMS/es_vortexincylinder.txt, "DATA1D/vortexincylinder.dat");
        double V1[5], V2[5];
        for(double x = 0.0; x<=1.000001; x+=0.01) {
            double C[3] = {x, 0.0, 0.0};
            S.PointValue( 0.0, C, V1);
            S.PointValue( 0.1, C, V2);
            fprintf(out, "%25.15f %25.15f %25.15f\n", x, V1[2], V2[2]);
        }
        fclose(out);
    }
    if(SolutionName.empty() || CompareWords(SolutionName, "waveinchannel")) {
        {
            SOLUTION(s_WaveInChannel<double>, PARAMS/es_waveinchannel_1.txt, "DATA2D/waveinchannel_1.dat");
            for(double y = 0.0; y<=0.5; y+=(0.02*(0.5-y)+0.001)) {
                for(double x = -0.25; x<=0.2500000001; x+=0.01) {
                    double T = 0.5, C[3] = {x, y, 0.0};
                    S.PointValue(T, C, V);
                    fprintf(out, "%25.15f %25.15f %25.15f\n", x, y, V[1]);
                }
                fprintf(out, "\n");
            }
            fclose(out);
        }
        {
            SOLUTION(s_WaveInChannel<double>, PARAMS/es_waveinchannel_2.txt, "DATA2D/waveinchannel_2.dat");
            for(double y = 0.0; y<=0.5; y+=0.01) {
                for(double x = -0.25; x<=0.2500000001; x+=0.01) {
                    double T = 1.0, C[3] = {x, y, 0.0};
                    S.PointValue(T, C, V);
                    fprintf(out, "%25.15f %25.15f %25.15f\n", x, y, V[1]);
                }
                fprintf(out, "\n");
            }
            fclose(out);
        }
        {
            SOLUTION(s_WaveInChannel<double>, PARAMS/es_waveinchannel_3.txt, "DATA1D/waveinchannel_3.dat");
            for(double r = 0.0; r<=1.0; r+=0.001) {
                double T = 0.1, C[3] = {0, r, 0.0};
                S.PointValue(T, C, V);
                fprintf(out, "%25.15f %25.15f %25.15f %25.15f\n", r, V[0], V[2], V[4]);
            }
            fclose(out);
        }
        {
            SOLUTION(s_WaveInChannel<double>, PARAMS/es_waveinchannel_4.txt, "DATA2D/waveinchannel_4.dat");
            for(double r = 0.0; r<=1.000000001; r+=0.01*(1-r)+1e-5) {
                for(double phi = 0; phi<=PiNumber/4.+0.0000001; phi+=0.01) {
                    double T = 0.0, C[3] = {r*cos(phi), r*sin(phi), 0.0};
                    S.PointValue(T, C, V);
                    fprintf(out, "%25.15f %25.15f %25.15f\n", C[0], C[1], V[0]);
                }
                fprintf(out, "\n");
            }
            fclose(out);
        }
    }
    if(SolutionName.empty() || CompareWords(SolutionName, "rankinevortex") || CompareWords(SolutionName, "gaussianvortex") || CompareWords(SolutionName, "finitevortex")) {
        s_RankineVortex SR;
        s_GaussianVortex SG;
        s_FiniteVortex SF;
        SR.ReadParamsFromFile("PARAMS/es_rankinevortex.txt");
        SG.ReadParamsFromFile("PARAMS/es_gaussianvortex.txt");
        SF.ReadParamsFromFile("PARAMS/es_finitevortex.txt");
        SR.Init();
        SG.Init();
        SF.Init();
        FILE* out = fopen("DATA1D/vortexes.dat", "wt");
        if(out==NULL) crash("Error opening file %s for writing\n", "DATA1D/vortexes.dat");
        printf("DATA1D/vortexes.dat\n");
        double VR[5], VG[5], VF[5];
        for(double r = 0.0; r<=8.0; r+=0.1) {
            double T = 0.0, C[3] = {r, 0.0, 0.0};
            SR.PointValue(T, C, VR);
            SG.PointValue(T, C, VG);
            SF.PointValue(T, C, VF);
            fprintf(out, "%25.15f %25.15f %25.15f %25.15f %25.15f %25.15f %25.15f\n", r, VF[2], VG[2], VR[2], VF[4], VG[4], VR[4]);
        }
        fclose(out);
    }
    if(SolutionName.empty() || CompareWords(SolutionName, "riemann")) {
        SOLUTION(s_Riemann, PARAMS/es_riemann.txt, "DATA1D/riemann.dat");
        for(double x = 175.0; x<=310.000001; x+=1.0) {
            double T = 4.0, C[3] = {x, 0.0, 0.0};
            S.PointValue(T, C, V);
            fprintf(out, "%25.15f %25.15f\n", x, V[0]);
        }
        fclose(out);
    }
    if(SolutionName.empty() || CompareWords(SolutionName, "simplewave")) {
        SOLUTION(s_SimpleWave, PARAMS/es_simplewave.txt, "DATA1D/simplewave.dat");
        double V1[5], V2[5];
        for(double x = -0.7; x<=0.5000001; x+=0.01) {
            double C[3] = {x, 0.0, 0.0};
            S.PointValue( 0.0, C, V1);
            S.PointValue( 0.09, C, V2);
            fprintf(out, "%25.15f %25.15f %25.15f %25.15f %25.15f\n", x, V1[0], V1[4], V2[0], V2[4]);
        }
        fclose(out);
    }
    if(SolutionName.empty() || CompareWords(SolutionName, "viscshock")) {
        {
            SOLUTION(s_ViscShock, PARAMS/es_viscshock_1.txt, "DATA1D/viscshock_1.dat");
            for(double x = -20.0; x<=20.0; x+=0.1) {
                double T = 0.0, C[3] = {x, 0.0, 0.0};
                S.PointValue(T, C, V);
                fprintf(out, "%25.15f %25.15f %25.15f %25.15f\n", x, V[0], V[1], V[4]);
            }
            fclose(out);
        }
        {
            SOLUTION(s_ViscShock, PARAMS/es_viscshock_2.txt, "DATA1D/viscshock_2.dat");
            for(double x = -20.0; x<=20.0; x+=0.1) {
                double T = 0.0, C[3] = {x, 0.0, 0.0};
                S.PointValue(T, C, V);
                fprintf(out, "%25.15f %25.15f %25.15f %25.15f\n", x, V[0], V[1], V[4]);
            }
            fclose(out);
        }
        {
            SOLUTION(s_ViscShock, PARAMS/es_viscshock_3.txt, "DATA1D/viscshock_3.dat");
            for(double x = -20.0; x<=20.0; x+=0.1) {
                double T = 0.0, C[3] = {x, 0.0, 0.0};
                S.PointValue(T, C, V);
                fprintf(out, "%25.15f %25.15f %25.15f %25.15f\n", x, V[0], V[1], V[4]);
            }
            fclose(out);
        }
    }
    if(SolutionName.empty() || CompareWords(SolutionName, "conccyl")) {
        SOLUTION(s_ConcCyl, PARAMS/es_conccyl.txt, "DATA2D/conccyl.dat");
        for(double phi = PiNumber/4.0; phi<=7*PiNumber/12.+0.000000001; phi+=0.02*PiNumber/3) {
            for(double r = 1.0; r<=2.00000001; r+=0.01) {
                double T = 0.0, C[3] = {r*cos(phi), r*sin(phi), 0.0};
                S.PointValue(T, C, V);
                fprintf(out, "%25.15f %25.15f %25.15f\n", C[0], C[1], V[4]);
            }
            fprintf(out, "\n");
        }
        fclose(out);
    }
    if(SolutionName.empty() || CompareWords(SolutionName, "curlfreecylinder")) {
        SOLUTION(s_CurlFreeCylinder, PARAMS/es_curlfreecylinder.txt, "DATA2D/curlfreecylinder.dat");
        for(double phi = 0; phi<=PiNumber+0.0000001; phi+=PiNumber/50) {
            for(double r = 0.5; r<=5.00000001; r+=0.01) {
                double T = 0.0, C[3] = {r*cos(phi), r*sin(phi), 0.0};
                S.PointValue(T, C, V);
                fprintf(out, "%25.15f %25.15f %25.15f\n", C[0], C[1], V[4]);
            }
            fprintf(out, "\n");
        }
        fclose(out);
    }
    if(SolutionName.empty() || CompareWords(SolutionName, "potentialsphere")) {
        SOLUTION(s_PotentialSphere, PARAMS/es_potentialsphere.txt, "DATA2D/potentialsphere.dat");
        for(double phi = 0; phi<=PiNumber+0.0000001; phi+=PiNumber/50) {
            for(double r = 0.5; r<=5.00000001; r+=0.025) {
                double T = 0.0, C[3] = {r*cos(phi), r*sin(phi), 0.0};
                S.PointValue(T, C, V);
                fprintf(out, "%25.15f %25.15f %25.15f\n", C[0], C[1], V[4]);
            }
            fprintf(out, "\n");
        }
        fclose(out);
    }
    if(SolutionName.empty() || CompareWords(SolutionName, "viscsphere")) {
        SOLUTION(s_ViscSphere, PARAMS/es_viscsphere.txt, "DATA2D/viscsphere.dat");
        for(double phi = 0; phi<=PiNumber+0.0000001; phi+=PiNumber/50) {
            for(double r = 0.5; r<=5.00000001; r+=0.025) {
                double T = 0.0, C[3] = {r*cos(phi), r*sin(phi), 0.0};
                S.PointValue(T, C, V);
                fprintf(out, "%25.15f %25.15f %25.15f\n", C[0], C[1], V[4]);
            }
            fprintf(out, "\n");
        }
        fclose(out);
    }
    return 0;
}
