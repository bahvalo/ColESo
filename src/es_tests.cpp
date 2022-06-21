// *********************************************************************************************************************
// *****                                                                                                           *****
// *****                                  Collection of exact solutions (ColESo)                                   *****
// *****                                             Self-test module                                              *****
// *****                                                                                                           *****
// *********************************************************************************************************************

#include "base_parser.h"
#include "coleso.h"
#include "spaceform.h"
#ifdef EXTRAPRECISION_COLESO
#include "extraprecision.h"
#endif

//======================================================================================================================
// Applying operator of linearized or full Euler equations (depending on solution type) to the exact solution 
// Input:  S       -- an exact solution to verify
//         T, coor -- point in time in space where derivatives should be calculated
//         numcoor -- number of coordinates in use (not to waste time for 'z' derivatives in 2D)
//         h       -- step for numerical differentiation
//         MF      -- for linearized equations, background field (rho, u, v, w, p).
//                    If MF==NULL, default values are used: 1, 0, 0, 0, 1/gamma
//         gamma   -- specific ratio (in use only for full Euler equations or if MF==NULL)
//         method  -- code of a formula for derivatives computation
// Output: V       -- computed values (5 variables)
//======================================================================================================================
template<typename fpv>
void Check_EE_LEE(const tPointFunction_EP<fpv>& S, fpv T, const fpv* coor, fpv* V, int numcoor, fpv h, 
              const fpv* MF, fpv gamma, int method) {
    static const fpv c_1_3  = fpv(1.0) / fpv(3.0);
    static const fpv c_4_3  = fpv(4.0) / fpv(3.0);
    static const fpv c_20_3 = fpv(20.0) / fpv(3.0);
    static const fpv c_3_5  = fpv(3.0) / fpv(5.0);
    static const fpv c_6_5  = fpv(6.0) / fpv(5.0);
    static const fpv c_1_6  = fpv(1.0) / fpv(6.0);
    static const fpv c_1_10 = fpv(1.0) / fpv(10.0);
    static const fpv c_25_12 = fpv(25.0) / fpv(12.0);
    static const fpv c_49_20 = fpv(49.0) / fpv(20.0);

    if(numcoor<1 || numcoor>3) crash("Check_EE_LEE error: wrong numcoor = %i", numcoor);
    if(method!=2 && method!=4 && method!=6) crash("Check_EE_LEE error: wrong method = %i", method);

    static const int N = tPointFunction::NumVarsMax;
    int n = S.NumVars();
    if(S.Type() == FUNC_PULSATIONS_COMPLEX) n = 5; // we'll look on the real part only
    if(n!=5) crash("Check_EE_LEE error: wrong NumVars = %i", n);

    int LEEmode = 0;
    if(S.Type()==FUNC_PULSATIONS || S.Type() == FUNC_PULSATIONS_COMPLEX) LEEmode = 1;
    else if(S.Type()==FUNC_PHYSICAL) LEEmode = 0;
    else crash("Check_EE_LEE error: wrong solution type");

    const fpv bufB[4] = {coor[0], coor[1], coor[2], T};
    const fpv invh = fpv(1.0)/h;
    fpv derivs[4][N];
    for(int icoor=0; icoor<4; icoor++) {
        // Do not compute the derivatives that are known to be zero
        if(icoor!=3 /*time*/ && icoor>=numcoor) {
            for(int ivar=0; ivar<n; ivar++) derivs[icoor][ivar] = 0.0;
            continue;
        }

        // If time is too small to use central derivatives, then we will use forward derivatives
        if(icoor==3 && T < h*(0.5*method + 0.5)) { 
            fpv VV[7][N];
            for(int i=0; i<=method; i++)
                S.PointValue(T+i*h, coor, VV[i]);
            for(int ivar=0; ivar<n; ivar++) {
                #define W(z) VV[z][ivar]
                switch(method) {
                case 2:
                    derivs[icoor][ivar] = -0.5*W(2) + 2.*W(1) - 1.5*W(0);
                    break;
                case 4:
                    derivs[icoor][ivar] = -0.25*W(4) + c_4_3*W(3) - 3.*W(2) + 4.*W(1) - c_25_12*W(0);
                    break;
                case 6:
                    derivs[icoor][ivar] = -c_1_6*W(6) + c_6_5*W(5) - 3.75*W(4) + c_20_3*W(3) - 7.5*W(2) + 6.*W(1) - c_49_20*W(0);
                    break;
                default: crash("Check_EE_LEE error: wrong method = %i", method);
                }
                derivs[icoor][ivar] *= invh;
                #undef W
            }
        }
        // Main branch: using central derivatives of the order 'method'
        else {
            fpv buf[4] = {coor[0], coor[1], coor[2], T};
            fpv& t = buf[3];
            fpv* c = buf;

            // Compute the values in a neighborhood of the point given
            fpv fm3[N], fm2[N], fm1[N], f1[N], f2[N], f3[N];
            fpv d1[N], d2[N], d3[N]; // 1-st order derivatives using 2-nd order central differences with step 2h, 4h, 6h
            buf[icoor] = bufB[icoor] - h;  S.PointValue(t, c, fm1);
            buf[icoor] = bufB[icoor] + h;  S.PointValue(t, c, f1);
            for(int ivar=0; ivar<n; ivar++) d1[ivar] = (f1[ivar] - fm1[ivar]) / (fpv(2.)*h);
            if(method>=4) {
                buf[icoor] = bufB[icoor] - 2.*h;  S.PointValue(t, c, fm2);
                buf[icoor] = bufB[icoor] + 2.*h;  S.PointValue(t, c, f2);
                for(int ivar=0; ivar<n; ivar++) d2[ivar] = (f2[ivar] - fm2[ivar]) / (fpv(4.)*h);
            }
            if(method>=6) {
                buf[icoor] = bufB[icoor] - 3.*h;  S.PointValue(t, c, fm3);
                buf[icoor] = bufB[icoor] + 3.*h;  S.PointValue(t, c, f3);
                for(int ivar=0; ivar<n; ivar++) d3[ivar] = (f3[ivar] - fm3[ivar]) / (fpv(6.)*h);
            }

            // Calculation of the first derivatives with high (2-nd, 4-th or 6-th) order using 2-nd order CD on different meshes
            for(int ivar=0; ivar<n; ivar++) {
                switch(method) {
                case 2:
                    derivs[icoor][ivar] = d1[ivar];
                    break;
                case 4:
                    derivs[icoor][ivar] = c_4_3*d1[ivar] - c_1_3*d2[ivar];
                    break;
                case 6:
                    derivs[icoor][ivar] = (fpv(1.5))*d1[ivar] - c_3_5*d2[ivar] + c_1_10*d3[ivar];
                    break;
                default: crash("Check_EE_LEE error: wrong method = %i", method);
                }
            }
        }
    }

    if(LEEmode) {
        const fpv rho_ref = MF==NULL ? fpv(1.0) : MF[Var_R];
        const fpv c_ref2 = MF==NULL ? fpv(1.0) : gamma*MF[Var_P]/MF[Var_R];

        // Now calculating the LEE operator applied to the solution given using the derivatives computed above
        const fpv div_u = derivs[0][Var_U] + derivs[1][Var_V] + derivs[2][Var_W];
        V[Var_R] = derivs[3][Var_R] + rho_ref*div_u;
        V[Var_U] = derivs[3][Var_U] + derivs[0][Var_P] / rho_ref;
        V[Var_V] = derivs[3][Var_V] + derivs[1][Var_P] / rho_ref;
        V[Var_W] = derivs[3][Var_W] + derivs[2][Var_P] / rho_ref;
        V[Var_P] = derivs[3][Var_P] + rho_ref*c_ref2*div_u;
        if(MF!=NULL) { // Drift
            for(int ivar=0; ivar<5; ivar++) 
                V[ivar] += derivs[0][ivar]*MF[Var_U] + derivs[1][ivar]*MF[Var_V] + derivs[2][ivar]*MF[Var_W];
        }
    }
    else {
        fpv U[5];
        S.PointValue(T, coor, U);
        // d(rho)/dt + d(rho ux)/dx + d(rho uy)/dy + d(rho uz)/dz
        const fpv div_u = derivs[0][Var_U] + derivs[1][Var_V] + derivs[2][Var_W];
        V[Var_R] = derivs[3][Var_R] + U[Var_R]*div_u
            + U[Var_U]*derivs[0][Var_R] + U[Var_V]*derivs[1][Var_R] + U[Var_W]*derivs[2][Var_R];
        // d(rho ux)/dt + d(rho ux^2 + p)/dx + d(rho ux uy)/dy + d(rho ux uz)/dz
        V[Var_U] = U[Var_R]*derivs[3][Var_U] + U[Var_U]*derivs[3][Var_R]
            + derivs[0][Var_R]*U[Var_U]*U[Var_U] + 2.*U[Var_R]*U[Var_U]*derivs[0][Var_U] + derivs[0][Var_P]
            + derivs[1][Var_R]*U[Var_U]*U[Var_V] + U[Var_R]*(U[Var_V]*derivs[1][Var_U] + U[Var_U]*derivs[1][Var_V])
            + derivs[2][Var_R]*U[Var_U]*U[Var_W] + U[Var_R]*(U[Var_W]*derivs[2][Var_U] + U[Var_U]*derivs[2][Var_W]);
        V[Var_V] = U[Var_R]*derivs[3][Var_V] + U[Var_V]*derivs[3][Var_R]
            + derivs[0][Var_R]*U[Var_U]*U[Var_V] + U[Var_R]*(U[Var_V]*derivs[0][Var_U] + U[Var_U]*derivs[0][Var_V])
            + derivs[1][Var_R]*U[Var_V]*U[Var_V] + 2.*U[Var_R]*U[Var_V]*derivs[1][Var_V] + derivs[1][Var_P]
            + derivs[2][Var_R]*U[Var_V]*U[Var_W] + U[Var_R]*(U[Var_W]*derivs[2][Var_V] + U[Var_V]*derivs[2][Var_W]);
        V[Var_W] = U[Var_R]*derivs[3][Var_W] + U[Var_W]*derivs[3][Var_R]
            + derivs[0][Var_R]*U[Var_U]*U[Var_W] + U[Var_R]*(U[Var_W]*derivs[0][Var_U] + U[Var_U]*derivs[0][Var_W]) 
            + derivs[1][Var_R]*U[Var_W]*U[Var_V] + U[Var_R]*(U[Var_V]*derivs[1][Var_W] + U[Var_W]*derivs[1][Var_V])
            + derivs[2][Var_R]*U[Var_W]*U[Var_W] + 2.*U[Var_R]*U[Var_W]*derivs[2][Var_W] + derivs[2][Var_P];
        // dE/dt + d((E+p) ux)/dx + d((E+p) uy)/dy + d((E+p) uz)/dz
        const fpv E = 0.5*U[Var_R]*VELDOT(U,U) + U[Var_P]/(gamma - 1.0);
        fpv gradE[4];
        for(int icoor=0; icoor<4; icoor++) {
            gradE[icoor] = 0.5*derivs[icoor][Var_R]*VELDOT(U,U) 
                + U[Var_R]*VELDOT(U, derivs[icoor]) 
                + derivs[icoor][Var_P] / (gamma - 1.0);
        }
        V[Var_E] = gradE[3] + (E+U[Var_P])*div_u
            + U[Var_U]*(gradE[0]+derivs[0][Var_P]) 
            + U[Var_V]*(gradE[1]+derivs[1][Var_P]) 
            + U[Var_W]*(gradE[2]+derivs[2][Var_P]);
    }
}


//======================================================================================================================
// Applying the operator of linearized Navier-Stokes equations on the uniform background field to a function
// Input:  S       -- an exact solution to verify
//         t, coor -- point in time in space where derivatives should be calculated
//         numcoor -- number of coordinates in use (not to waste time for 'z' derivatives in 2D)
//         h       -- step for numerical differentiation
//         MF      -- background field (rho, u, v, w, p). If MF==NULL, default values are used: 1, 0, 0, 0, 1/gamma
//         nu, Pr, gamma -- equation parameters
//         method  -- code of a formula for derivatives computation
// Output: V       -- computed values (5 variables)
//======================================================================================================================
template<typename fpv>
void CheckLNSE(const tPointFunction_EP<fpv>& S, fpv t, const fpv* coor, fpv* V, int numcoor, fpv h, 
               const fpv* MF, fpv nu, fpv Pr, fpv gamma, int method) {
    // Calculate the convective part and time derivative first
    Check_EE_LEE(S, t, coor, V, numcoor, h, MF, gamma, method);

    static const int N = tPointFunction::NumVarsMax;
    int n = S.NumVars();
    if(S.Type() == FUNC_PULSATIONS_COMPLEX) n = 5; // we'll look on the real part only
    if(n!=5) crash("CheckLNSE error: wrong NumVars = %i", n);

    tFixBlock<fpv,5> f[7][7][7]; // stored values. Positions: {-3, -2, -1, 0, 1, 2, 3} at each direction
    for(int i=0; i<7; i++) for(int j=0; j<7; j++) for(int k=0; k<7; k++) for(int ivar=0; ivar<5; ivar++) MakeNaN(f[i][j][k][ivar]);

    // Compute the values in a neighborhood of the point given
    for(int i=0; i<7; i++) for(int j=0; j<7; j++) for(int k=0; k<7; k++) {
        if(numcoor!=3 && k!=3) { f[i][j][k] = fpv(0.0); continue; }
        if(numcoor==1 && j!=3) { f[i][j][k] = fpv(0.0); continue; }
        if(i!=3 && j!=3 && k!=3) continue;
        if(method==2 && (i==0 || i==1 || i==5 || i==6 || j==0 || j==1 || j==5 || j==6 || k==0 || k==1 || k==5 || k==6)) continue;
        if(method==4 && (i==0 || i==6 || j==0 || j==6 || k==0 || k==6)) continue;

        fpv buf[4] = {coor[0]+(i-3)*h, coor[1]+(j-3)*h, coor[2]+(k-3)*h};
        fpv val[N];
        S.PointValue(t, buf, val);
        f[i][j][k] = val;
    }

    // Calculation of direct second derivatives
    static const fpv c_1_3  = fpv(1.0) / fpv(3.0);
    static const fpv c_4_3  = fpv(4.0) / fpv(3.0);
    static const fpv c_3_5  = fpv(3.0) / fpv(5.0);
    static const fpv c_1_10 = fpv(1.0) / fpv(10.0);
    static const fpv c_1_12 = fpv(1.0) / fpv(12.0);
    static const fpv c_3_20 = fpv(3.0) / fpv(20.0);
    static const fpv c_1_90 = fpv(1.0) / fpv(90.0);
    static const fpv c_49_18 = fpv(49.0) / fpv(18.0);

    tFixBlock<fpv,5> D2V[3][3]; // second order derivatives (symmetric matrix 3x3 of N-vectors)
    switch(method) {
    case 2:
        D2V[0][0] = f[4][3][3]-fpv(2.0)*f[3][3][3]+f[2][3][3];
        D2V[1][1] = f[3][4][3]-fpv(2.0)*f[3][3][3]+f[3][2][3];
        D2V[2][2] = f[3][3][4]-fpv(2.0)*f[3][3][3]+f[3][3][2];
        break;
    case 4:
        D2V[0][0] = -fpv(2.5)*f[3][3][3] + c_4_3*(f[2][3][3]+f[4][3][3]) - c_1_12*(f[1][3][3]+f[5][3][3]);
        D2V[1][1] = -fpv(2.5)*f[3][3][3] + c_4_3*(f[3][2][3]+f[3][4][3]) - c_1_12*(f[3][1][3]+f[3][5][3]);
        D2V[2][2] = -fpv(2.5)*f[3][3][3] + c_4_3*(f[3][3][2]+f[3][3][4]) - c_1_12*(f[3][3][1]+f[3][3][5]);
        break;
    case 6:
        D2V[0][0] = -c_49_18*f[3][3][3] + fpv(1.5)*(f[2][3][3]+f[4][3][3]) - c_3_20*(f[1][3][3]+f[5][3][3]) + c_1_90*(f[0][3][3]+f[6][3][3]);
        D2V[1][1] = -c_49_18*f[3][3][3] + fpv(1.5)*(f[3][2][3]+f[3][4][3]) - c_3_20*(f[3][1][3]+f[3][5][3]) + c_1_90*(f[3][0][3]+f[3][6][3]);
        D2V[2][2] = -c_49_18*f[3][3][3] + fpv(1.5)*(f[3][3][2]+f[3][3][4]) - c_3_20*(f[3][3][1]+f[3][3][5]) + c_1_90*(f[3][3][0]+f[3][3][6]);
        break;
    default: crash("CheckLNSE error: wrong method = %i", method);
    }

    // Calculation of cross derivatives
    for(int icoor=1; icoor<3; icoor++) for(int jcoor=0; jcoor<icoor; jcoor++) {
        // Calculating of cross derivatives using 2-nd order formula on different stencils
        int pp[3]={3,3,3}, pm[3]={3,3,3}, mp[3]={3,3,3}, mm[3]={3,3,3}; 
        tFixBlock<fpv,5> D1, D2, D3;
        {
            pp[icoor]++; pp[jcoor]++;
            pm[icoor]++; pm[jcoor]--;
            mp[icoor]--; mp[jcoor]++;
            mm[icoor]--; mm[jcoor]--;
            D1 = f[pp[0]][pp[1]][pp[2]] - f[pm[0]][pm[1]][pm[2]] - f[mp[0]][mp[1]][mp[2]] + f[mm[0]][mm[1]][mm[2]];
            D1 *= fpv(0.25);
        }
        if(method>=4) {
            pp[icoor]++; pp[jcoor]++;
            pm[icoor]++; pm[jcoor]--;
            mp[icoor]--; mp[jcoor]++;
            mm[icoor]--; mm[jcoor]--;
            D2 = f[pp[0]][pp[1]][pp[2]] - f[pm[0]][pm[1]][pm[2]] - f[mp[0]][mp[1]][mp[2]] + f[mm[0]][mm[1]][mm[2]];
            D2 *= fpv(0.0625);
        }
        if(method>=6) {
            pp[icoor]++; pp[jcoor]++;
            pm[icoor]++; pm[jcoor]--;
            mp[icoor]--; mp[jcoor]++;
            mm[icoor]--; mm[jcoor]--;
            D3 = f[pp[0]][pp[1]][pp[2]] - f[pm[0]][pm[1]][pm[2]] - f[mp[0]][mp[1]][mp[2]] + f[mm[0]][mm[1]][mm[2]];
            D3 *= fpv(1.)/fpv(36.);
        }

        // Combining to get high-order approximation
        switch(method) {
        case 2:
            D2V[icoor][jcoor] = D1;
            break;
        case 4:
            D2V[icoor][jcoor] = c_4_3*D1 - c_1_3*D2;
            break;
        case 6:
            D2V[icoor][jcoor] = (fpv(1.5))*D1 - c_3_5*D2 + c_1_10*D3;
            break;
        default: crash("CheckLNSE error: wrong method = %i", method);
        }
        D2V[jcoor][icoor] = D2V[icoor][jcoor];
    }

    const fpv invh2 = 1.0 / (h*h);
    for(int i=0; i<3; i++) for(int j=0; j<3; j++) D2V[i][j] *= invh2;
    tFixBlock<fpv,5> L = D2V[0][0] + D2V[1][1] + D2V[2][2]; // Laplassian

    // Now calculating the viscous & heat operator applied to the solution given using the derivatives computed above
    V[Var_U] -= nu*L[1] + c_1_3*nu*(D2V[0][0][1]+D2V[0][1][2]+D2V[0][2][3]);
    V[Var_V] -= nu*L[2] + c_1_3*nu*(D2V[1][0][1]+D2V[1][1][2]+D2V[1][2][3]);
    V[Var_W] -= nu*L[3] + c_1_3*nu*(D2V[2][0][1]+D2V[2][1][2]+D2V[2][2][3]);
    V[Var_P] -= (nu / Pr) * (gamma*L[4] - L[0]);
}


//======================================================================================================================
// Instantiating functions of the template class
//======================================================================================================================
#define INSTANTIATE(fpv) \
    template void Check_EE_LEE(const tPointFunction_EP<fpv>& S, fpv T, const fpv* coor, fpv* V, int numcoor, fpv h, \
              const fpv* MF, fpv gamma, int method); \
    template void CheckLNSE(const tPointFunction_EP<fpv>& S, fpv t, const fpv* coor, fpv* V, int numcoor, fpv h, \
              const fpv* MF, fpv nu, fpv Pr, fpv gamma, int method);


INSTANTIATE(NativeDouble)
#ifdef EXTRAPRECISION_COLESO
INSTANTIATE(dd_real)
INSTANTIATE(qd_real)
#endif

template<typename fpv>
inline int CheckErr(const fpv* V, int N, fpv max_err) {
    for(int i=0; i<N; i++) if(IsNaN(V[i])) return 2;
    for(int i=0; i<N; i++) if(fabs(V[i]) > max_err) return 1;
    return 0;
}

// *********************************************************************************************************************
// *****                                                                                                           *****
// *****                                       Tests for specific solutions:                                       *****
// *****                             check that the solution really satisfies the equation                         *****
// *****                                                                                                           *****
// *********************************************************************************************************************

template<typename fpv>
void VerifyLEE_Gaussian2D(int NTests) {
    fpv warn_eps = pow_to_const_int<sizeof(fpv)/8>(1e-8);
    fpv c[3] = {0.,0.,0.}, t;
    
    for(int i=-NTests; i; i++) {
        s_Gaussian2D<fpv> S;
        S.Aterm = 1.0;
        S.Bterm = 1.0;
        S.FlowVelX = S.FlowVelY = S.FlowVelZ = 0.0;
        S.r0[0] = S.r0[1] = S.r0[2] = 0.0;
        S.Form = S.FORM_GAUSSIAN;
        S.NormalizeForm = 0;
        S.Init();

        t    = exp((fpv(rand()) / RAND_MAX) * 6.0);
        for(int ir=0; ir<3; ir++) {
            fpv r = S.Bterm;
            if(ir==1) r = exp((fpv(rand()) / RAND_MAX) * 6.0);
            if(ir==2) r = S.Bterm * pow((fpv(rand()) / RAND_MAX), 8.);
            if(i==-1) { r = 0.0; t = 300.0; if(ir) break; }
            fpv phi = (fpv(rand()) / RAND_MAX) * Pi2;
            c[0] = r * cos(phi);
            c[1] = r * sin(phi);

            fpv eps = pow_to_const_int<sizeof(fpv)/8>(1e-3);

            fpv LEEerr[5];
            Check_EE_LEE<fpv>((const tPointFunction_EP<fpv> &)(tPointFunction&)S, t, c, LEEerr, 2, eps, NULL, 1.4, 6);
            if(CheckErr(LEEerr, 5, warn_eps)) {
                Check_EE_LEE<fpv>((const tPointFunction_EP<fpv> &)(tPointFunction&)S, t, c, LEEerr, 2, 0.1*eps, NULL, 1.4, 6);
            }
            if(CheckErr(LEEerr, 5, warn_eps)) {
                pprintf0("WARNING! (t = %f, r = %f)\n", double(t), double(r));
                pprintf0("err: % 6.2e % 6.2e % 6.2e\n", double(LEEerr[0]), double(LEEerr[1]), double(LEEerr[2]));
            }
        }
    }
}

template<typename fpv>
void VerifyLEE_IVP2D(int NTests) {
    fpv warn_eps = pow_to_const_int<sizeof(fpv)/8>(1e-8);
    fpv c[3] = {0.,0.,0.}, t;
    
    for(int i=-NTests; i; i++) {
        s_IVP2D<fpv> S;
        S.FlowVelX = 0.1; S.FlowVelY = 0.1;
        S.Aterm = 1.0;
        S.Bterm = 1.0;
        S.r0[0] = S.r0[1] = S.r0[2] = 0.0;
        S.Form = S.FORM_COS2;
        S.NormalizeForm = 1;
        //S.GI.GR = 256;
        S.Init();

        const fpv gamma = 1.4;
        const fpv MF[5] = {1.0, S.FlowVelX, S.FlowVelY, 0., 1.0/gamma};

        t    = exp((fpv(rand()) / RAND_MAX) * 6.0);
        for(int ir=0; ir<3; ir++) {
            fpv r = S.Bterm;
            if(ir==1) r = exp((fpv(rand()) / RAND_MAX) * 6.0);
            if(ir==2) r = S.Bterm * pow((fpv(rand()) / RAND_MAX), 8.);
            if(i==-1) { r = 7.271832339860; t = 14.972247327990; if(ir) break; }
            fpv phi = (fpv(rand()) / RAND_MAX) * Pi2;
            c[0] = r * cos(phi);
            c[1] = r * sin(phi);

            fpv eps = pow_to_const_int<sizeof(fpv)/8>(1e-3);

            fpv LEEerr[5];
            Check_EE_LEE<fpv>((const tPointFunction_EP<fpv> &)(tPointFunction&)S, t, c, LEEerr, 2, eps, MF, gamma, 6);
            if(CheckErr(LEEerr, 5, warn_eps)) {
                Check_EE_LEE<fpv>((const tPointFunction_EP<fpv> &)(tPointFunction&)S, t, c, LEEerr, 2, 0.1*eps, MF, gamma, 6);
            }
            if(CheckErr(LEEerr, 5, warn_eps)) {
                pprintf0("WARNING! (t = %17.12f, r = %17.12f, phi = %10.6f)\n", double(t), double(r), double(phi));
                pprintf0("err: % 6.2e % 6.2e % 6.2e % 6.2e % 6.2e\n", double(LEEerr[0]), double(LEEerr[1]), double(LEEerr[2]), double(LEEerr[3]), double(LEEerr[4]));
            }
        }
    }
}

//template<typename fpv>
void VerifyLEE_CornerPlanar(int NTests) {
    typedef double fpv;
    fpv c[3] = {0.,0.,0.}, t;
    
    for(int i=-NTests; i; i++) {
//        pprintf0("i = %i\n", -i);
        s_CornerPlanar S;
        int half_line_case = i&1;
        if(i==-1) half_line_case = 1; // or 0
        if(half_line_case) {
            S.angle = Pi2;
            S.phi0 = 0.25*PiNumber;
            S.X0 = 50.0;
            S.Bterm = 2.0;
            S.Aterm = 1.0;
        }
        S.set_reducer(0.1);
        S.Init();

        const fpv gamma = 1.4; // unused

        t    = exp((double(rand()) / RAND_MAX) * 6.0);
        for(int ir=0; ir<3; ir++) {
            double r = fabs(S.X0 - t);
            if(ir==1) r = exp((double(rand()) / RAND_MAX) * 6.0);
            if(ir==2) r = t + S.Bterm * pow((double(rand()) / RAND_MAX), 8.);
            double phi = (double(rand()) / RAND_MAX) * S.angle;
            if(i==-1) { r = 397.926243684650; t = 1.846064024587; phi = 2.285317621112; if(ir) break; }
            c[0] = r * cos(phi);
            c[1] = r * sin(phi);

            fpv eps = pow_to_const_int<sizeof(fpv)/8>(1e-3);

            // Check that out point is not too near to the corner boundary
            if(S.angle < 3.0) {
                double dist = MIN(c[1], r * sin(S.angle - phi));
                if(dist < 3.01*eps) continue; // multiplicator 3 is due to 6-th order derivatives
            }
            if(S.angle > 5.0) {
                double dist = r; 
                if(c[0] > 0.0) dist = MIN(r, fabs(c[1]));
                if(dist < 3.01*eps) continue; // multiplicator 3 is due to 6-th order derivatives
            }
            fpv warn_eps = pow_to_const_int<sizeof(fpv)/8>(1e-8) / MIN(r, 1.0);

            fpv LEEerr[5];
            Check_EE_LEE<fpv>((const tPointFunction_EP<fpv> &)(tPointFunction&)(tPointFunction&)(S), t, c, LEEerr, 2, eps, NULL, gamma, 6);
            if(fabs(LEEerr[0]) > warn_eps || fabs(LEEerr[1]) > warn_eps || fabs(LEEerr[2]) > warn_eps) {
                Check_EE_LEE<fpv>((const tPointFunction_EP<fpv> &)(tPointFunction&)S, t, c, LEEerr, 2, 0.1*eps, NULL, gamma, 6);
            }
            if(CheckErr(LEEerr, 5, warn_eps)) {
                pprintf0("WARNING! (t = %17.12f, r = %17.12f, \n          phi = %17.12f, angle = %17.12f)\n", 
                    double(t), double(r), double(phi), double(S.angle));
                pprintf0("err: % 6.2e % 6.2e % 6.2e\n", double(LEEerr[0]), double(LEEerr[1]), double(LEEerr[2]));
            }
        }
    }
}


template<typename fpv>
void VerifyLEE_Corner(int NTests) {
    fpv c[3] = {0.,0.,0.}, t;
    
    for(int i=-NTests; i; i++) {
        //pprintf0("i = %i\n", -i);
        s_Corner<fpv> S;
        int half_line_case = i&1;
        if(i==-1) half_line_case = 0;
        if(half_line_case) {
            S.angle = Pi2;
        }
        else {
            S.angle = Pi2 / 3.0;
        }
        S.Aterm = 1.0;
        S.r0[0] = 50.0/sqrt(2.0);
        S.r0[1] = 50.0/sqrt(2.0);
        S.Bterm = 2.0;
        S.Init();

        const fpv gamma = fpv(7.)/fpv(5.); // unused

        t = (double(rand()) / RAND_MAX);
        if(t > 0.5) t = exp((2.*t-1.0) * 6.0);
        else t = exp(-20*t);

        for(int ir=0; ir<3; ir++) {
            double r = fabs(sqrt(SQR(S.r0[0])+SQR(S.r0[1])) - t);
            if(ir==1) r = exp((double(rand()) / RAND_MAX) * 6.0);
            if(ir==2) r = (double)(t + S.Bterm * pow((double(rand()) / RAND_MAX), 8.));
            double phi = (double(rand()) / RAND_MAX) * S.angle;
            if(i==-1) { r = 220.251500183477; t = 270.251500183477; phi = 1.430864246412; if(ir) break; }
            c[0] = r * cos(phi);
            c[1] = r * sin(phi);

            fpv eps = pow_to_const_int<sizeof(fpv)/8>(1e-3);

            // Check that out point is not too near to the corner boundary
            if(S.angle < 3.0) {
                double dist = MIN(c[1], r * sin(S.angle - phi));
                if(dist < 3.01*eps) continue; // multiplicator 3 is due to 6-th order derivatives
            }
            if(S.angle > 5.0) {
                double dist = r;
                if(c[0] > 0.0) dist = MIN(r, fabs(double(c[1])));
                if(dist < 3.01*eps) continue; // multiplicator 3 is due to 6-th order derivatives
            }
            fpv warn_eps = pow_to_const_int<sizeof(fpv)/8>(1e-8) / MIN(r, 1.0);

            fpv LEEerr[5];
            Check_EE_LEE<fpv>((const tPointFunction_EP<fpv> &)(tPointFunction&)S, t, c, LEEerr, 2, eps, NULL, gamma, 6);
            if(fabs(LEEerr[0]) > warn_eps || fabs(LEEerr[1]) > warn_eps || fabs(LEEerr[2]) > warn_eps) {
                Check_EE_LEE<fpv>((const tPointFunction_EP<fpv> &)(tPointFunction&)S, t, c, LEEerr, 2, 0.1*eps, NULL, gamma, 6);
            }
            if(CheckErr(LEEerr, 5, warn_eps)) {
                pprintf0("WARNING! (t = %17.12f, r = %17.12f, \n          phi = %17.12f, angle = %17.12f)\n", 
                    double(t), double(r), double(phi), double(S.angle));
                pprintf0("err: % 6.2e % 6.2e % 6.2e\n", double(LEEerr[0]), double(LEEerr[1]), double(LEEerr[2]));
            }
        }
    }
}


//template<typename fpv>
void VerifyLEE_Gaussian3D(int NTests) {
    typedef double fpv;
    fpv warn_eps = pow_to_const_int<sizeof(fpv)/8>(1e-8);
    fpv c[3] = {0.,0.,0.}, t;
    
    for(int i=-NTests; i; i++) {
        s_Gaussian3D S;
        S.Aterm = 1.0;
        S.Bterm = 1.0;
        S.FlowVelX = S.FlowVelY = S.FlowVelZ = 0.0;
        S.r0[0] = S.r0[1] = S.r0[2] = 0.0;
        S.Form = S.FORM_GAUSSIAN;
        S.NormalizeForm = 0;
        S.Init();

        t    = exp((double(rand()) / RAND_MAX) * 6.0);
        for(int ir=0; ir<2; ir++) {
            double r = S.Bterm;
            if(ir==1) r = exp((double(rand()) / RAND_MAX) * 4.0);
            else      r = S.Bterm * pow((double(rand()) / RAND_MAX), 8.);
            double phi = (double(rand()) / RAND_MAX) * Pi2;
            c[0] = r * cos(phi);
            c[1] = r * sin(phi);

            double eps = sizeof(fpv)==8 ? 1e-2 : 3e-5;

            fpv LEEerr[5];
            Check_EE_LEE<fpv>((const tPointFunction_EP<fpv> &)(tPointFunction&)S, t, c, LEEerr, 3, eps, NULL, 1.4, 6);
            if(fabs(LEEerr[0]) > warn_eps || fabs(LEEerr[1]) > warn_eps || fabs(LEEerr[2]) > warn_eps) {
                Check_EE_LEE<fpv>((const tPointFunction_EP<fpv> &)(tPointFunction&)S, t, c, LEEerr, 3, 0.1*eps, NULL, 1.4, 6);
            }
            if(CheckErr(LEEerr, 5, warn_eps)) {
                pprintf0("WARNING! (t = %f, r = %f)\n", double(t), double(r));
                pprintf0("err: % 6.2e % 6.2e % 6.2e\n", double(LEEerr[0]), double(LEEerr[1]), double(LEEerr[2]));
            }
        }
    }
}


template<typename fpv>
void VerifyLEE_Source1D(int NTests) {
    fpv warn_eps = pow_to_const_int<sizeof(fpv)/8>(1e-8);
    fpv c[3] = {0.,0.,0.}, t = 0.;
    
    for(int i=-NTests; i; i++) {
        s_Source1D<fpv> S;
        S.Form = tSpaceForm<fpv>::FORM_GAUSSIAN;
        S.Aterm = 1.0;
        S.Bterm = 5.0;
        S.NormalizeForm = 0;
        S.r0[0] = S.r0[1] = S.r0[2] = 0.0;
        S.Ampl = 1.0;
        S.Freq = fpv(1.0)/fpv(10.0);
        S.Phase = 0.0;
        S.tmin = 0.0; S.tmax = 1e50;
        S.SignalType = 1;
        S.Init();

        if(i!=-1) t    = exp((double(rand()) / RAND_MAX) * 4.0);
        for(int ir=0; ir<((i==-1)?1:2); ir++) {
            fpv r, phi=0.0;
            if(i==-1) {
                t = 33, r = 50;
            }
            else {
                r = S.Bterm;
                if(ir==1) r = exp((double(rand()) / RAND_MAX) * 4.0);
            }
            c[0] = r * cos(phi);
            c[1] = r * sin(phi);

            fpv eps = sizeof(fpv)==8 ? 1e-2 : 3e-5;

            fpv LEEerr[5];
            Check_EE_LEE<fpv>((const tPointFunction_EP<fpv> &)(tPointFunction&)S, t, c, LEEerr, 1, eps, NULL, 1.4, 6);
            fpv source = S.SpaceForm(c) * S.TimeForm(t);
            LEEerr[Var_R] -= source;
            LEEerr[Var_P] -= source;
            if(CheckErr(LEEerr, 5, warn_eps)) {
                pprintf0("WARNING! (t = %25.15f, r = %25.15f, \n    phi = %25.15f)\n", double(t), double(r), double(phi));
                pprintf0("err: % 6.2e % 6.2e % 6.2e\n", double(LEEerr[0]), double(LEEerr[1]), double(LEEerr[2]));
            }
        }
    }
}

void VerifyLEE_Source2D(int NTests) {
    typedef double fpv;
    fpv warn_eps = pow_to_const_int<sizeof(fpv)/8>(2e-7);
    fpv c[3] = {0.,0.,0.}, t = 0.;
    
    for(int i=-NTests; i; i++) {
        s_Source2D S;
        S.Form = tSpaceForm<fpv>::FORM_COS2;
        S.Aterm = 1.0;
        S.Bterm = 5.0;
        S.NormalizeForm = 0;
        S.r0[0] = S.r0[1] = S.r0[2] = 0.0;
        S.Ampl = 1.0;
        S.Freq = fpv(1.0)/fpv(10.0);
        S.Phase = 0.0;
        S.GR = 200;
        S.mm = 1;
        S.Init();

        if(i!=-1) t    = exp((fpv(rand()) / RAND_MAX) * 4.0);
        for(int ir=0; ir<((i==-1)?1:2); ir++) {
            fpv r, phi=0.0;
            if(i==-1) {
                t = 1.785, r = 5.00115;
            }
            else {
                r = S.Bterm;
                if(ir==1) r = exp((fpv(rand()) / RAND_MAX) * 4.0);
                else      r = S.Bterm * pow((fpv(rand()) / RAND_MAX), 8);
            }
            c[0] = r * cos(phi);
            c[1] = r * sin(phi);

            fpv eps = pow_to_const_int<sizeof(fpv)/8>(5e-3);

            fpv LEEerr[5];
            Check_EE_LEE<fpv>((const tPointFunction_EP<fpv> &)(tPointFunction&)S, t, c, LEEerr, 3, eps, NULL, 1.4, 6);
            fpv source = S.SpaceForm(c) * sin(Pi2*S.Freq*t);
            LEEerr[Var_R] -= source;
            LEEerr[Var_P] -= source;
            if(CheckErr(LEEerr, 5, warn_eps)) {
                pprintf0("WARNING! (t = %25.15f, r = %25.15f, \n    phi = %25.15f)\n", double(t), double(r), double(phi));
                pprintf0("err: % 6.2e % 6.2e % 6.2e\n", double(LEEerr[0]), double(LEEerr[1]), double(LEEerr[2]));
            }
        }
    }
}

template<typename fpv> 
void VerifyLEE_Source3D(int NTests) {
    fpv warn_eps = pow_to_const_int<sizeof(fpv)/8>(1e-8);
    fpv c[3] = {0.,0.,0.}, t = 0.;
    
    for(int i=-NTests; i; i++) {
        s_Source3D<fpv> S;
        S.Form = tSpaceForm<fpv>::FORM_GAUSSIAN;
        S.Aterm = 1.0;
        S.Bterm = 5.0;
        S.NormalizeForm = 0;
        S.r0[0] = S.r0[1] = S.r0[2] = 0.0;
        S.Ampl = 1.0;
        S.Freq = fpv(1.0)/fpv(10.0);
        S.Phase = 0.0;
        S.tmin = 0.0; S.tmax = 1e50;
        S.SignalType = 1;
        S.Init();

        if(i!=-1) t    = exp((fpv(rand()) / RAND_MAX) * 4.0);
        for(int ir=0; ir<((i==-1)?1:2); ir++) {
            fpv r, phi=0.0;
            if(i==-1) {
                t = 19.8, r = 31.055;
            }
            else {
                r = S.Bterm;
                if(ir==1) r = exp((fpv(rand()) / RAND_MAX) * 4.0);
                else      r = S.Bterm * pow((fpv(rand()) / RAND_MAX), 8);
            }
            c[0] = r * cos(phi);
            c[1] = r * sin(phi);

            fpv eps = pow_to_const_int<sizeof(fpv)/8>(5e-3);

            fpv LEEerr[5];
            Check_EE_LEE<fpv>((const tPointFunction_EP<fpv> &)(tPointFunction&)S, t, c, LEEerr, 3, eps, NULL, 1.4, 6);
            fpv source = S.SpaceForm(c) * S.TimeForm(t);
            LEEerr[Var_R] -= source;
            LEEerr[Var_P] -= source;
            if(CheckErr(LEEerr, 5, warn_eps)) {
                pprintf0("WARNING! (t = %25.15f, r = %25.15f, \n    phi = %25.15f)\n", double(t), double(r), double(phi));
                pprintf0("err: % 6.2e % 6.2e % 6.2e\n", double(LEEerr[0]), double(LEEerr[1]), double(LEEerr[2]));
            }
        }
    }
}

void VerifyLEE_RotatingDipole(int NTests) {
    typedef double fpv;
    fpv warn_eps = pow_to_const_int<sizeof(fpv)/8>(1e-8);
    fpv c[3] = {0.,0.,0.}, t = 0.;
    
    for(int i=-NTests; i; i++) {
        s_RotatingDipole S;
        S.Init();

        t = exp((double(rand()) / RAND_MAX) * 4.0);
        fpv r = 1.0 + (double(rand()) / RAND_MAX) * 10.0;
        fpv phi = Pi2 * double(rand()) / RAND_MAX;
        fpv theta = -0.5*PiNumber + PiNumber * double(rand()) / RAND_MAX;

        c[0] = r * cos(theta) * cos(phi);
        c[1] = r * cos(theta) * sin(phi);
        c[2] = r * sin(theta);

        fpv eps = pow_to_const_int<sizeof(fpv)/8>(5e-3);

        fpv LEEerr[5];
        Check_EE_LEE<fpv>((const tPointFunction_EP<fpv> &)(tPointFunction&)S, t, c, LEEerr, 3, eps, NULL, 1.4, 6);
        if(CheckErr(LEEerr, 5, warn_eps)) {
            pprintf0("WARNING! (t = %25.15f, r = %25.15f, \n    phi = %25.15f)\n", double(t), double(r), double(phi));
            pprintf0("err: % 6.2e % 6.2e % 6.2e\n", double(LEEerr[0]), double(LEEerr[1]), double(LEEerr[2]));
        }
    }
}

template<typename fpv>
void VerifyLEE_PointSource(int NTests) {
    fpv warn_eps = 1e7 /* no ideas why are the arithmetic errors so big */ * pow_to_const_int<sizeof(fpv)/8>(1e-12);
    fpv c[3] = {0.,0.,0.}, t = 0.;
    
    for(int i=-NTests; i; i++) {
        s_PointSource<fpv> S;
        S.FlowMach = (fpv(rand()) / RAND_MAX) * 0.8;
        S.SoundVel = 0.1 + (fpv(rand()) / RAND_MAX) * 10;
        S.SignalType = 6;
        S.tmin = -10000.0;
        S.tmax = 10000.0;

        t = exp((fpv(rand()) / RAND_MAX) * 4.0);
        fpv r = 1.0 + (fpv(rand()) / RAND_MAX) * 10.0;
        fpv phi = Pi2 * fpv(rand()) / RAND_MAX;
        fpv theta = -0.5*PiNumber + PiNumber * fpv(rand()) / RAND_MAX;

        if(i==-1) {
            t = 44.026136598983335, r =         6.813776055177465, phi =  2.881287344757850, theta =        -0.222194310577543;
            S.FlowMach = 0.888024536881619;
            S.SoundVel = 0.122888882106998;
        }

        S.Init();

        c[0] = r * cos(theta) * cos(phi);
        c[1] = r * cos(theta) * sin(phi);
        c[2] = r * sin(theta);

        fpv gamma = 1.1 + double(rand()) / RAND_MAX;
        fpv MF[5] = {fpv(1.0), S.FlowMach*S.SoundVel, fpv(0.0), fpv(0.0), SQR(S.SoundVel) / gamma};

        fpv LEEerr[5];
        fpv eps = pow_to_const_int<sizeof(fpv)/8>((double)(1e-3 * (1.0 - S.FlowMach) * MIN(fpv(1.0),S.SoundVel) / S.Freq));
        Check_EE_LEE<fpv>((const tPointFunction_EP<fpv> &)(tPointFunction&)S, t, c, LEEerr, 3, eps, MF, gamma, 6);
        if(CheckErr(LEEerr, 5, warn_eps)) {
            eps *= 10;
            Check_EE_LEE<fpv>((const tPointFunction_EP<fpv> &)(tPointFunction&)S, t, c, LEEerr, 3, eps, MF, gamma, 6);
        }
        if(CheckErr(LEEerr, 5, warn_eps)) {
            eps *= 0.01;
            Check_EE_LEE<fpv>((const tPointFunction_EP<fpv> &)(tPointFunction&)S, t, c, LEEerr, 3, eps, MF, gamma, 6);
        }
        if(CheckErr(LEEerr, 5, warn_eps)) {
            pprintf0("WARNING! (t = %25.15f, r = %25.15f, \n    phi = %25.15f, theta = %25.15f)\n    FlowMach = %25.15f, SoundVel = %25.15f\n", 
                double(t), double(r), double(phi), double(theta), double(S.FlowMach), double(S.SoundVel));
            pprintf0("err: % 6.2e % 6.2e % 6.2e % 6.2e % 6.2e\n", double(LEEerr[0]), double(LEEerr[1]), double(LEEerr[2]), double(LEEerr[3]), double(LEEerr[4]));
        }
    }
}

void VerifyLEE_Coaxial(int NTests) {
    typedef double fpv;
    fpv warn_eps = 2e-6; // high tolerance due to the inaccuracies in Bessel functions computation
    fpv c[3] = {0.,0.,0.}, t = 0.;
    
    for(int i=-NTests; i; i++) {
        s_Coaxial S;
        S.FlowVel = double(rand()) / RAND_MAX;
        S.SoundVel = 0.2 + double(rand()) / RAND_MAX;
        S.AzimuthalMode = 3; S.RadialMode = 4;
        S.rmin = 0.2; S.rmax = 2.;
        S.kz = GetPiNumber<fpv>() * (0.5 + double(rand()) / RAND_MAX);
        S.Ampl = 1.0; S.phase = 0.4; 
        S.CoorAxis = 2;
        S.Init();

        fpv eps = pow_to_const_int<sizeof(fpv)/8>(5e-3);

        t = exp((double(rand()) / RAND_MAX) * 4.0);
        fpv r = S.rmin + 5.*eps + (double(rand()) / RAND_MAX) * (S.rmax - S.rmin - 10.*eps);
        fpv phi = Pi2 * double(rand()) / RAND_MAX;
        fpv z = double(rand()) / RAND_MAX;

        if(i==-1) r = 0.28;

        c[0] = r * cos(phi);
        c[1] = r * sin(phi);
        c[2] = z;

        fpv gamma = 1.1 + double(rand()) / RAND_MAX;
        fpv MF[5] = {1.0, 0.0, 0.0, S.FlowVel, SQR(S.SoundVel) / gamma};

        fpv LEEerr[5];
        Check_EE_LEE<fpv>((const tPointFunction_EP<fpv> &)(tPointFunction&)S, t, c, LEEerr, 3, eps, MF, gamma, 6);
        if(CheckErr(LEEerr, 5, warn_eps)) {
            pprintf0("WARNING! (t = %f, r = %f)\n", double(t), double(r));
            pprintf0("err: % 6.2e % 6.2e % 6.2e % 6.2e % 6.2e\n", double(LEEerr[0]), double(LEEerr[1]), double(LEEerr[2]), double(LEEerr[3]), double(LEEerr[4]));
        }
    }
}

void VerifyLNSE_SinusVisc(int NTests) {
    typedef double fpv;
    fpv warn_eps = pow_to_const_int<sizeof(fpv)/8>(1e-10);
    
    const fpv gamma = 1.4; // unused
    for(int i=-NTests; i; i++) {
        s_SinusVisc S;
        S.mu = 0.01;
        S.Ampl = double(rand()) / RAND_MAX;
        S.AmplVX = double(rand()) / RAND_MAX;
        S.AmplVY = double(rand()) / RAND_MAX;
        S.AmplVZ = double(rand()) / RAND_MAX;
        S.kx = S.ky = S.kz = Pi2 / 3.0;
        S.Init();

        fpv t = double(rand()) / RAND_MAX;
        fpv c[3];
        c[0] = double(rand()) / RAND_MAX;
        c[1] = double(rand()) / RAND_MAX;
        c[2] = double(rand()) / RAND_MAX;

        fpv eps = pow_to_const_int<sizeof(fpv)/8>(3e-3);

        fpv LEEerr[5];
        CheckLNSE<fpv>((const tPointFunction_EP<fpv> &)(tPointFunction&)S, t, c, LEEerr, 3, eps, NULL, S.mu, 1e50 /*no heat*/, gamma, 6);
        if(CheckErr(LEEerr, 5, warn_eps)) {
            pprintf0("WARNING! \n");
            pprintf0("err: % 6.2e % 6.2e % 6.2e % 6.2e % 6.2e\n", double(LEEerr[0]), double(LEEerr[1]), double(LEEerr[2]), double(LEEerr[3]), double(LEEerr[4]));
        }
    }
}

template<typename fpv>
void VerifyLNSE_WaveInChannel(int mode) {
    s_WaveInChannel<fpv> S;
    S.CoorAxis = 2;
    S.k = 0.0;
    S.kmode = 1;
    S.l = 8;
    S.nu = 1e-5;
    S.R = 1.0;
    S.gamma = 1.4;
    S.Prandtl = 1.0;
    S.form = 1;
    S._dmumax = 0.02;
    S.Init();

    fpv C[3]  = {0.999, 0.0, 0.0};
    fpv time  = 0.0;

    fpv V[10];
    if(mode<0) {
        fpv e = 1e-2;
        if(e > 0.25*(S.R-C[0])) e = 0.25*(S.R-C[0]);
        for(int k=0; k<10; k++) {
            CheckLNSE((const tPointFunction_EP<fpv>&)S, time, C, V, 3, e, (const fpv*)NULL, S.nu, S.Prandtl, S.gamma, 6);
            printf("%.7e %.7e %.7e %.7e %.7e\n", (double)fabs(V[0]), (double)fabs(V[1]), (double)fabs(V[2]), (double)fabs(V[3]), (double)fabs(V[4]));
            e /= fpv(10.0);
        }

        FILE* out = fopen("QA_log.txt", "wt");
        fprintf(out, "err_rho %g\n", double(fabs(V[0])));
        fprintf(out, "err_u   %g\n", double(fabs(V[1])));
        fprintf(out, "err_v   %g\n", double(fabs(V[2])));
        fprintf(out, "err_w   %g\n", double(fabs(V[3])));
        fprintf(out, "err_p   %g\n", double(fabs(V[4])));
        fclose(out);
    }
    else {
        fpv e = pow_to_const_int<sizeof(fpv)/8>(3e-3);
        if(e > 0.25*(S.R-C[0])) e = 0.25*(S.R-C[0]);
        fpv eps = pow_to_const_int<sizeof(fpv)/8>(1e-4);
        CheckLNSE((const tPointFunction_EP<fpv>&)S, time, C, V, 3, e, (const fpv*)NULL, S.nu, S.Prandtl, S.gamma, 6);
        double err = double(fabs(V[0]) + fabs(V[1]) + fabs(V[2]) + fabs(V[3]) + fabs(V[4]));
        if(err > eps) {
            pprintf0("WARNING! \n");
            pprintf0("err: % 6.2e % 6.2e % 6.2e % 6.2e % 6.2e\n", double(V[0]), double(V[1]), double(V[2]), double(V[3]), double(V[4]));
        }
    }
}

void VerifyEE_Vortexes(int NTests) {
    typedef double fpv;
    fpv c[3] = {0.,0.,0.}, t = 0.;

    for(int i=-NTests; i; i++) {
        fpv warn_eps = 1e-8;
        tVortexFunction* VF;
        switch(i&3) {
        case 0: {
            VF = new s_FiniteVortex;
            s_FiniteVortex& S = *((s_FiniteVortex*)VF);
            S.mode = 0;
            S.deg = 3;
            warn_eps = 1e-7;
            break;
                }
        case 1:
            VF = new s_RankineVortex;
            break;
        case 2:
            VF = new s_GaussianVortex;
            break;
        case 3:
            VF = new s_Vortex_BG;
            break;
        default: crash("ERROR");
        }

        tVortexFunction& S = *VF;
        S.r0[0] = S.r0[1] = S.r0[2] = 0.0;
        S.v0[0] = 0.1; S.v0[1] = 0.2; S.v0[2] = 0.3;
        S.Mach = 0.6543;
        S.R = 1. + (double(rand()) / RAND_MAX) * 10;
        S.axis = 1;
        S.gam = 1.1 + double(rand()) / RAND_MAX;

        fpv eps = pow_to_const_int<sizeof(fpv)/8>(5e-3);

        fpv r = ((double(rand()) / RAND_MAX) + 0.5) * S.R;
        if((i&3)==1 && fabs(r - S.R) < eps * 4.0) r = 0.5 * S.R; // Rankine vortex is not smooth at r=R
        fpv phi = Pi2 * double(rand()) / RAND_MAX;
        fpv z = double(rand()) / RAND_MAX;

        if(i==-3) { 
            r = 4.04;
            S.R = 4.21;
            phi = 2.456745205712603;
            S.gam = 1.4;
        }

        S.Init();

        c[0] = r * sin(phi);
        c[1] = z;
        c[2] = r * cos(phi);

        fpv EEerr[5];
        Check_EE_LEE<fpv>((const tPointFunction_EP<fpv> &)(tPointFunction&)S, t, c, EEerr, 3, eps, NULL, S.gam, 6);
        if(CheckErr(EEerr, 5, warn_eps)) {
            pprintf0("WARNING! (t = %f, r = %20.15f, R = %20.15f, \n   phi = %20.15f, gam = %20.15f vortex = %i)\n", double(t), double(r), double(S.R), double(phi), double(S.gam), i&3);
            pprintf0("err: % 6.2e % 6.2e % 6.2e % 6.2e % 6.2e\n", double(EEerr[0]), double(EEerr[1]), double(EEerr[2]), double(EEerr[3]), double(EEerr[4]));
        }

        delete VF;
    }
}

void VerifyEE_SimpleWave(int NTests) {
    typedef double fpv;
    fpv c[3] = {0.,0.,0.}, t = 0.;

    for(int i=-NTests; i; i++) {
        fpv warn_eps = 1e-8;
        s_SimpleWave S;
        S.x0 = 0.0; // perturbation center
        S.l = 0.2; 
        S.gam = 5./3.; // specific ratio
        S.SoundVel = sqrt(10.)/3.;
        S.FlowVel = -sqrt(10.);
        S.a = pow(2.,C1_3) - 1.0;
        S.log = 0;

        fpv eps = pow_to_const_int<sizeof(fpv)/8>(5e-5);

        fpv r = ((double(rand()) / RAND_MAX) - 0.5) * S.l * 2.0;

        S.Init();

        c[0] = r;
        c[1] = c[2] = 0.0;
        t = 0.001;

        fpv EEerr[5];
        Check_EE_LEE<fpv>((const tPointFunction_EP<fpv> &)(tPointFunction&)S, t, c, EEerr, 3, eps, NULL, S.gam, 6);
        if(CheckErr(EEerr, 5, warn_eps)) {
            pprintf0("WARNING! (t = %f, r = %20.15f, l = %20.15f)\n", double(t), double(r), double(S.l));
            pprintf0("err: % 6.2e % 6.2e % 6.2e % 6.2e % 6.2e\n", double(EEerr[0]), double(EEerr[1]), double(EEerr[2]), double(EEerr[3]), double(EEerr[4]));
        }
    }
}

// *********************************************************************************************************************
// *****                                                                                                           *****
// *****                                       Interface to the test module                                        *****
// *****                                                                                                           *****
// *********************************************************************************************************************

void ColESo_Tests(int argc, char** argv) { 
    if(argc > 1) {
        string Case;
        int NTests = -1; // infinite
        tParamManager PM;
        CmdLine.AddArguments(argc-1, argv+1);
        PM.RequestParameter(Case, "", PARTYPE_WORD, IO_DONTCRASH);
        if(argc >= 3) PM.RequestParameter(NTests, "", PARTYPE_INT, IO_DONTCRASH);
        PM.ReadParamsFromCommandLine(false /*log*/);
        CmdLine.PrintUnusedArguments();
        CmdLine.ClearAllArguments();

        // Linear cases
        if(CompareWords(Case, "Gaussian2D"))        { VerifyLEE_Gaussian2D<double>(NTests); return; }
        if(CompareWords(Case, "Source1D"))          { VerifyLEE_Source1D<double>(NTests); return; }
        if(CompareWords(Case, "Source2D"))          { VerifyLEE_Source2D(NTests); return; }
        if(CompareWords(Case, "Source3D"))          { VerifyLEE_Source3D<double>(NTests); return; }
        if(CompareWords(Case, "Gaussian3D"))        { VerifyLEE_Gaussian3D(NTests); return; }
        if(CompareWords(Case, "RotatingDipole"))    { VerifyLEE_RotatingDipole(NTests); return; }
        if(CompareWords(Case, "IVP2D"))             { VerifyLEE_IVP2D<double>(NTests); return; }
        if(CompareWords(Case, "Corner"))            { VerifyLEE_Corner<double>(NTests); return; }
        if(CompareWords(Case, "CornerPlanar"))      { VerifyLEE_CornerPlanar(NTests); return; }
        if(CompareWords(Case, "SinusVisc"))         { VerifyLNSE_SinusVisc(NTests); return; }
        if(CompareWords(Case, "Coaxial"))           { VerifyLEE_Coaxial(NTests); return; }
        if(CompareWords(Case, "PointSource"))       { VerifyLEE_PointSource<double>(NTests); return; }
        if(CompareWords(Case, "WaveInChannel"))     { VerifyLNSE_WaveInChannel<double>(1); return; }

        #ifdef EXTRAPRECISION_COLESO
        if(CompareWords(Case, "Gaussian2D-DD"))     { VerifyLEE_Gaussian2D<dd_real>(NTests); return; }
        if(CompareWords(Case, "Source1D-DD"))       { VerifyLEE_Source1D<dd_real>(NTests); return; }
        if(CompareWords(Case, "Source3D-DD"))       { VerifyLEE_Source3D<dd_real>(NTests); return; }
        if(CompareWords(Case, "IVP2D-DD"))          { VerifyLEE_IVP2D<dd_real>(NTests); return; }
        if(CompareWords(Case, "PointSource-DD"))    { VerifyLEE_PointSource<dd_real>(NTests); return; }
        if(CompareWords(Case, "PointSource-QD"))    { VerifyLEE_PointSource<qd_real>(NTests); return; }
        if(CompareWords(Case, "Corner-DD"))         { VerifyLEE_Corner<dd_real>(NTests); return; }
        if(CompareWords(Case, "Corner-QD"))         { VerifyLEE_Corner<qd_real>(NTests); return; }
        if(CompareWords(Case, "WaveInChannel-DD"))  { VerifyLNSE_WaveInChannel<dd_real>(1); return; }
        if(CompareWords(Case, "WaveInChannel-QD"))  { VerifyLNSE_WaveInChannel<qd_real>(1); return; }
        #endif

        // WaveInChannel: printing residual on a succession of meshes
        if(CompareWords(Case, "WaveInChannel-succ"))     { VerifyLNSE_WaveInChannel<double>(-1 /*special mode*/); return; }
        #ifdef EXTRAPRECISION_COLESO
        if(CompareWords(Case, "WaveInChannel-succ-DD"))  { VerifyLNSE_WaveInChannel<dd_real>(-1); return; }
        if(CompareWords(Case, "WaveInChannel-succ-QD"))  { VerifyLNSE_WaveInChannel<qd_real>(-1); return; }
        #endif

        // Nonlinear cases
        if(CompareWords(Case, "Vortexes"))          { VerifyEE_Vortexes(NTests); return; }
        if(CompareWords(Case, "SimpleWave"))        { VerifyEE_SimpleWave(NTests); return; }
        crash("Wrong case %s\n", Case.c_str());
    }

    #define LAUNCH(CASE, NUMTESTS) \
        { char C[] = #CASE, N[] = #NUMTESTS; printf("Launching %s (%s tests)\n", C, N); char* buf[3] = {NULL, C, N}; ColESo_Tests(3, buf); }
    printf("Launching all tests. Type ColESo <testname> <num_tests> for a specific test\n");
    // Linear cases
    LAUNCH(Gaussian2D, 100);
    LAUNCH(IVP2D, 100);
    LAUNCH(Gaussian3D, 100);
    LAUNCH(Source1D, 100);
    LAUNCH(Source2D, 40);
    LAUNCH(Source3D, 40);
    LAUNCH(RotatingDipole, 100);
    LAUNCH(Corner, 50);
    LAUNCH(CornerPlanar, 50);
    LAUNCH(SinusVisc, 1000);
    LAUNCH(Coaxial, 1000);
    LAUNCH(PointSource, 1000);
    LAUNCH(WaveInChannel, 1);

    // Nonlinear cases
    LAUNCH(Vortexes, 100000);
    LAUNCH(SimpleWave, 5000);
}

