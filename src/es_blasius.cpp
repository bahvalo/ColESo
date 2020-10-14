// *********************************************************************************************************************
// *****                                                                                                           *****
// *****                                  Collection of exact solutions (ColESo)                                   *****
// *****                      Blasius profile and compessible modification of Blasius profile                      *****
// *****                                                                                                           *****
// *********************************************************************************************************************

#include "parser.h"
#include "coleso.h"
#ifdef _NOISETTE
#include "lib_base.h"
#endif

// *********************************************************************************************************************
// *****                                                                                                           *****
// *****                                      Blasius boundary layer                                               *****
// *****                                                                                                           *****
// *********************************************************************************************************************

//======================================================================================================================
void s_Blasius::ReadParams(tFileBuffer& FB) {
    tParamManager PM;
    PM.Request(Rey, "re");
    PM.Request(gam, "gam");
    PM.Request(FlowVel, "FlowVel");
    PM.Request(SoundVel, "SoundVel");
    PM.ReadParamsFromBuffer(FB);
}
//======================================================================================================================

//======================================================================================================================
// Calculating the solution of Blasius problem: 
// 2 f''' + f'' f = 0 with boundary conditions f(0) = f'(0) = 0, f(+infty) = 1.
// In fact, we know that f''(0) = 0.3320573362151946 and solve the Cauchy problem
// F' = G(F), F = (f, f', f''), G(f) = (f', f'', -f f''/2)
// Output: UV0 = f', UV1 = 1/2*(eta*f' - f), where eta is the independent variable
//======================================================================================================================
const double s_Blasius::alpha = 0.3320573362151946;
void s_Blasius::Init() {
    // Parameters of the numerical solution. Chosen in a way to get the solution with double precision
    const double L_sim = 12.5;  // size of region
    const int N_sim = 4096;     // number of intervals
    const double h = L_sim / N_sim; // mesh step
    d_eta = h;

    UV0.resize(N_sim + 1);
    UV1.resize(N_sim + 1);
    double F[3] = {0., 0., alpha}; // initial data
    UV0[0] = UV1[0] = 0.0;

    for (int i = 1; i <= N_sim; ++i){
        double k1[3] = {h*F[1], h*F[2], -0.5*h*F[0]*F[2]};
        double k2[3] = {h*(F[1]+0.5*k1[1]), h*(F[2]+0.5*k1[2]), -0.5*h*((F[0]+0.5*k1[0])*(F[2]+0.5*k1[2]))};
        double k3[3] = {h*(F[1]+0.5*k2[1]), h*(F[2]+0.5*k2[2]), -0.5*h*((F[0]+0.5*k2[0])*(F[2]+0.5*k2[2]))};
        double k4[3] = {h*(F[1]+    k3[1]), h*(F[2]+    k3[2]), -0.5*h*((F[0]+    k3[0])*(F[2]+    k3[2]))};
        for(int j=0; j<3; j++) F[j] += C1_6*(k1[j]+k4[j]) + C1_3*(k2[j]+k3[j]);

        UV0[i] = F[1];
        UV1[i] = 0.5*(i*h*F[1] - F[0]);
    }
}
//======================================================================================================================


//======================================================================================================================
// Interpolation of function with values V[0], V[1], V[2], V[3], V[4] (at points 0, 1, 2, 3, 4) to point x
// using 4-th order Lagrange polynomial
//======================================================================================================================
static double FivePointInterpolation(const double* V, double x) {
    double sum = 0.0;
    for(int i=0; i<5; i++) {
        double m = 1.0;
        for(int j=0; j<5; j++) {
            if(j==i) continue;
            m *= (x - j) / (i - j);
        }
        sum += V[i] * m;
    }
    return sum;
}
//======================================================================================================================

//======================================================================================================================
// Blasius profile: obtaining the solution in the given point using interpolation
//======================================================================================================================
void s_Blasius::PointValue(double /*t*/, const double* coor, double* V) const{
    if(UV0.size()<5 || UV0.size()!=UV1.size()) crash("s_Blasius: Init not done");
    V[Var_R] = 1.0;
    V[Var_U] = V[Var_V] = V[Var_W] = 0.0;
    V[Var_P] = SoundVel*SoundVel / gam;

    double x = coor[0], y = fabs(coor[1]);
    if (y < 1e-50) return;

    const double eta = (x < 1e-10) ? huge : sqrt(Rey) * y / sqrt(x);
    const int i = int(eta / d_eta);
    if(i < 3) { // Use 2 terms of Taylor series
        double eta3 = eta*eta*eta;
        double f  = 0.5 * alpha*eta*eta * (1.0 - alpha*eta3 / 120.0);
        double df = alpha*eta* (1.0 - alpha*eta3 / 48.0);
        V[Var_U] = df;
        V[Var_V] = 0.5*(eta*df - f);
    }
    else if(i > int(UV0.size() - 4)) { // Free stream values
        V[Var_U] = 1.;
        V[Var_V] = UV1[UV0.size() - 1];
    }
    else { // Interpolation
        double dd = eta / d_eta - i; // fractional part, 0 <= dd <= 1
        // Interpolation over ((i-2) - (i-1) - (i) - (i+1) - (i+2)) and ((i-1) - (i) - (i+1) - (i+2) - (i+3))
        V[Var_U] = FivePointInterpolation(&(UV0[i-2]),2.0+dd)*(1.0-dd) + FivePointInterpolation(&(UV0[i-1]),1.0+dd)*dd;
        V[Var_V] = FivePointInterpolation(&(UV1[i-2]),2.0+dd)*(1.0-dd) + FivePointInterpolation(&(UV1[i-1]),1.0+dd)*dd;
    }

    V[Var_U] *= FlowVel;
    V[Var_V] *= FlowVel / sqrt(Rey * x);
    if(coor[1]<0) V[Var_V] = - V[Var_V];
}
//======================================================================================================================



// *********************************************************************************************************************
// *****                                                                                                           *****
// *****                              Compressible Blasius boundary layer                                          *****
// *****                                                                                                           *****
// *********************************************************************************************************************


//======================================================================================================================
void s_comprBlasius::ReadParams(tFileBuffer& FB) {
    tParamManager PM;
    PM.Request(Rey, "re"); // Reynolds number
    //PM.Request(FlowVel, "FlowVel"); // Flow speed = 1.0
    PM.Request(SoundVel, "SoundVel"); // Sound speed (1 / Mach)
    PM.Request(Pr, "pr"); // Prandtl number
    PM.Request(gam, "gam"); //  Specific ratio

    // Viscosity
    PM.Request(mu_type, "mu_type");  // type of viscosity ( 0 = power law; 1 = Sutherland law )
    PM.Request(ps, "ps");            //  parameter of power law(mu = T^ps)
    PM.Request(sat, "sat");          // parameter of Sutherland law

    // BC parameters
    PM.Request(bc_type, "bc_type");  // type of BC for heat ( 0 = adiabatic; 1 = isothermic )
    PM.Request(temp_wall, "temp_wall");  //  temperature value  (if bc_type == 1 )  

    // parameters for calculation selfsimilar flow
    PM.Request(L_sim, "L_sim");       //  mesh size 
    PM.Request(N_sim, "N_sim");       //  number of intervals 
    // parameters for iterative procedure
    PM.Request(eps_rel, "eps_rel");   // relative accuracy
    PM.Request(eps_abs, "eps_abs");   // absolute accuracy
    PM.Request(Iter_max, "Iter_max"); // max. number of iterations

    PM.ReadParamsFromBuffer(FB);
}
//======================================================================================================================

//======================================================================================================================
double PowerLaw(double x, double y) {     //  h^(ps-1)
    return pow(x, y);
}

double Sutherland(double x, double y) {     // h*sqrt(h)*(1+sat)/(h+sat)
    return x*sqrt(x)*(1.0 + y) / (x + y);
}

typedef double(*FuncDoubleArgs)(double, double);

//======================================================================================================================
void s_comprBlasius::Init() {
    if(IsNaN(SoundVel)) crash("SoundVel not set");
    const double ma = 1.0 / SoundVel;
    table_KSIUVT.resize((N_sim + 1)*4);

    FuncDoubleArgs ViscousLaw = NULL;
    double pn, buf;

    const int n2 = 2 * N_sim + 1;

    //  iteration  i
    vector<double> z(n2);
    vector<double> h(n2);
    vector<double> fi(n2);
    //  iteration  i+1     
    vector<double> z1(n2);
    vector<double> h1(n2);
    vector<double> fi1(n2);

    //  прогоночные коэффициенты
    vector<double> al(n2);
    vector<double> bet(n2);

    double gg = gam - 1.;
    double ama = 0.25 * gg * Pr * ma * ma;

    if (mu_type == 0) {
        ViscousLaw = PowerLaw;   // h^(n-1)
        pn = ps - 1.;
    }
    else {
        ViscousLaw = Sutherland;   // h^1.5  (1 + sat) / ( h + sat )
        pn = sat;  //(=Ts)
    }

    // ====== Cycle: two meshes
    // kp=0    base (coarse) mesh (h);
    // kp=1    fine mesh (0.5 h)

    for (int kp = 0, n = N_sim + 1; kp < 2; kp++, n += N_sim) {
        double hz = L_sim / (n - 1);
        double hh = hz * hz;
        double ss = sqrt(L_sim);

        // ------ Initial distribution
        if (bc_type == 0) temp_wall = 1.0 + 0.5 * gg * sqrt(Pr) * ma * ma;
        for (int i = 0; i < n; i++) {
            h[i] = (1.0 - temp_wall) * (double)(i)* hz / L_sim + temp_wall;   // h    h(0)=hw    h(L_sim)=1
            z[i] = 2.0 * sqrt((double)(i)* hz) / ss;       // z = fi'
            fi[i] = 0.0;       // fi
        }

        //   ====== Cycle: nonlinear iterations
        //  h1,fi1,z1 - предыдущая итерация
        //  h,fi,z    - новая итерация   

        bool nonliniter = true;
        int k = 0;
        while (nonliniter) {
            k++;
            if(k > Iter_max) crash("s_comprBlasius::Init error, maximal number of iterations exceeded");

            // ------ Computation of function FI' = z (Cauchy problem: trapezoid rule)
            fi1[0] = 0.0;
            for (int i = 1; i < n; i++) {
                fi1[i] = fi1[i - 1] + 0.5 * hz * (z[i] + z[i - 1]);
            }

            // ------ Computation of function Z=2*U (Two-point boundary value problem: Samarsky's monotonic scheme and sweep method)
            //                s+1    s
            //         DUMAX=| Z  -  Z|
            double hs, hl, hp;   double r, kap, a, b, c, d;

            al[1] = 0.0;      //bet[1] = 0.0;
            for (int i = 1; i < n - 1; i++) {
                hs = ViscousLaw(h[i], pn);
                hl = 0.5 * (ViscousLaw(h[i - 1], pn) + hs);
                hp = 0.5 * (ViscousLaw(h[i + 1], pn) + hs);
                r = 0.5 * hz * fi1[i] / hs;
                kap = 1.0 / (1.0 + r);
                a = hl * kap / hh;        // a*z_{i-1} + c*z_i + b*z_{i+1} = 0
                b = hp * (kap + hz * fi1[i] / hs) / hh;
                c = a + b;
                d = c - al[i] * a;
                al[i + 1] = b / d;   // bet[i + 1] = a * bet[i] / d;
            }

            //------ обратный ход прогонки
            z1[n - 1] = 2.0;
            double dumax = fabs(z1[n - 1] - z[n - 1]);
            for (int i = n - 2; i >= 0; i--) {
                z1[i] = al[i + 1] * z1[i + 1];   // +bet[i + 1];
                dumax = MAX(dumax, fabs(z1[i] - z[i]));
            }

            // Computation of function H (Two-point boundary value problem: Samarsky's monotonic scheme and sweep method)
            //                s+1    s
            //         DHMAX=| H  -  H|

            double hs0, hs1, f;
            if (bc_type == 0) { // adiabatic BC
                al[1] = 1.0;
                hs0 = ViscousLaw(h[0], pn);
                hs1 = ViscousLaw(h[1], pn);
                bet[1] = 0.5 * ama * hs0 * (z1[1] * z1[1]) / hs1;
            }
            else {
                al[1] = 0.0;
                bet[1] = temp_wall;
            }

            for (int i = 1; i < n - 1; i++) {
                hs = ViscousLaw(h[i], pn);
                hl = 0.5 * (ViscousLaw(h[i - 1], pn) + hs);
                hp = 0.5 * (ViscousLaw(h[i + 1], pn) + hs);
                r = 0.5 * hz * Pr * fi1[i] / hs;
                kap = 1.0 / (1.0 + r);
                a = hl * kap / hh;
                b = hp * (kap + hz * Pr * fi1[i] / hs) / hh;
                c = a + b;
                buf = 5.0 * (z1[i + 1] - z1[i - 1]) / (10.0 * hz);
                f = ama * hs * buf * buf;

                d = c - al[i] * a;
                al[i + 1] = b / d;
                bet[i + 1] = (a * bet[i] + f) / d;
            }
            //  backward substitution
            h1[n - 1] = 1.0;
            double dhmax = fabs(h1[n - 1] - h[n - 1]);
            for (int i = n - 2; i >= 0; i--) {
                h1[i] = al[i + 1] * h1[i + 1] + bet[i + 1];
                dhmax = MAX(dhmax, fabs(h1[i] - h[i]));
            }

            // ====== Condition of nonlinear iteration finish
            nonliniter = false;
            for (int i = 0; i < n; i++) {
                double err1 = eps_rel * fabs(z[i]) + eps_abs;
                double err2 = eps_rel * fabs(h[i]) + eps_abs;
                if (dumax > err1 || dhmax > err2) {
                    nonliniter = true;
                    break;
                }
            }

            if (nonliniter) {  // ------ update: s+1 iteration 
                for (int i = 0; i < n; i++) {
                    fi[i] = fi1[i];
                    h[i] = h1[i];
                    z[i] = z1[i];
                }
            }

        } // end of nonlinear iterations
          // ------ Go to the mesh in 2 times smaller step
        if (kp == 0) {
            for (int i = 0; i < n; i++) {
                double* table_i = &(table_KSIUVT[4*i]);
                table_i[1] = z1[i];   // U
                table_i[2] = fi1[i];  // V
                table_i[3] = h1[i];    // T
            }
        }

    } // End mesh cycle

    // ====== The Runge method of increasing the order of accuracy
    for (int i = 0; i < N_sim + 1; i++) {
        double* table_i = &(table_KSIUVT[4*i]);
        table_i[1] = C4_3 * z1[2 * i] - C1_3 * table_i[1];    // U
        table_i[2] = C4_3 * fi1[2 * i] - C1_3 * table_i[2];    // V
        table_i[3] = C4_3 * h1[2 * i] - C1_3 * table_i[3];    // T
    }

    // ======= calculate table of physical variables
    double hz = L_sim / N_sim;
    for (int i = 0; i < N_sim + 1; i++) {
        double* table_i = &(table_KSIUVT[4*i]);
        table_i[0] = 0.0;  // KSI
        fi1[i] = log(table_i[3]);
        table_i[1] = 0.5 * table_i[1];   // u = 0.5 * fi' = 0.5 * z
        table_i[2] = i * hz * table_i[1] * (1.0 - fi1[i]) - 0.5 * table_i[2];  //  integral(0-eta) [log h] d(eta)
    }

    fi[0] = 0.0;
    for (int i = 1; i < N_sim + 1; i++) {
        double* table_i  = &(table_KSIUVT[4*i]);
        double* table_i1 = &(table_KSIUVT[4*(i-1)]);
        table_i[0] = table_i1[0] + 0.5 * (table_i1[3] + table_i[3]) * hz * 2.;   // integral(0-eta) [h] d(eta)
        fi[i] = fi[i - 1] + 0.5 * (fi1[i - 1] + fi1[i]) * hz;    // integral(0-eta) [log h] d(eta)
        table_i[2] = (table_i[2] + table_i[1] * fi[i]) * table_i[3];
    }
}   
//======================================================================================================================

//======================================================================================================================
void s_comprBlasius::PointValue(double /*t*/, const double* coor, double* V) const {
    V[Var_P] = SoundVel*SoundVel / gam;
    V[Var_W] = 0.0;

    double x = coor[0], y = fabs(coor[1]);
    const double eps = 1e-7;
    if (x < eps) {   
        V[Var_R] = 1. / table_KSIUVT[3]; //  rho
        V[Var_U] = 1.; //  U
        V[Var_V] = 0.; //  V
        return;
    }

    const double ksi_current = sqrt(Rey) * y / sqrt(x);

    int ii = 0; 
    const double* table_ii = &(table_KSIUVT[0]);
    const double* table_ii1 = table_ii;
    for (int i = 0; i < N_sim; i++) {
        if (ksi_current >= table_ii[0] && ksi_current <= table_ii1[0]) break;
        table_ii = table_ii1;   
        table_ii1 = table_ii + 4;        
        ii = i + 2;       
     }

    if (fabs(ksi_current - table_ii[0]) < eps) {
        V[Var_R] = 1. / table_ii[3];   //  rho
        V[Var_U] = table_ii[1];
        V[Var_V] = table_ii[2];      
    }
    else if (fabs(ksi_current - table_ii1[0]) < eps) {
        V[Var_R] = 1. / table_ii1[3];
        V[Var_U] = table_ii1[1];
        V[Var_V] = table_ii1[2];    
    }
    else if (ii < N_sim) {
        double delta_ksi = table_ii1[0] - table_ii[0];
        double coeff_left = (table_ii1[0] - ksi_current) / delta_ksi;
        double coeff_right = 1. - coeff_left;     
        double rho_inv = table_ii[3] * coeff_left + table_ii1[3] * coeff_right;
        V[Var_R] = 1. / rho_inv;
        V[Var_U] = table_ii[1] * coeff_left + table_ii1[1] * coeff_right;
        V[Var_V] = table_ii[2] * coeff_left + table_ii1[2] * coeff_right;
    }
    else {
        V[Var_R] = 1. / table_ii[3];
        V[Var_U] = table_ii[1];
        V[Var_V] = table_ii[2];
    }
    V[Var_V] /= sqrt(Rey * x);
    if(coor[1]<0) V[Var_V] *= -1.0;
}
//======================================================================================================================

//======================================================================================================================
// Dump of internal data for debug purpose
//======================================================================================================================
int s_comprBlasius::DumpTable(const char* fname) const {
    if(fname==NULL) fname = "s_comprBlasius_data.txt";
    FILE* F = fopen(fname, "wt");
    if(F==NULL) return 1;

    const double hz = L_sim / N_sim;
    fprintf(F, "%15s %15s %15s %15s %15s\n", "ETA", "KSI", "U", "T", "V");
    for (int i = 0; i < N_sim + 1; i++) {
        const double* table_i = &(table_KSIUVT[i*4]);
        fprintf(F, "% 15.5f % 15.5f % 15.5f % 15.5f % 15.5f\n", 
            i * hz, table_i[0], table_i[1], table_i[3], table_i[2]);
    }
    fclose(F);
    return 0;
}
//======================================================================================================================
