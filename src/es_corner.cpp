// *********************************************************************************************************************
// *****                                                                                                           *****
// *****                                  Collection of exact solutions (ColESo)                                   *****
// *****                            Wave propagation from 2D initial axisymmetric pulse;                           *****
// *****              Gaussian pulse and planar wave with Gaussian profile diffraction inside a sector             *****
// *****                                                                                                           *****
// *********************************************************************************************************************

#include "pointfuncs.h"
#include "base_parser.h"
#include "coleso.h"
#include "geom_primitive.h"
#ifdef EXTRAPRECISION_COLESO
#include "extraprecision.h"
#endif
#include "es_specfunc.h"

//======================================================================================================================
// Common subroutine: solution computation via Taylor expansion near t=0
// uex = (rho', u', v', wave potential, p')
//======================================================================================================================
template<typename fpv>
void PointValueTaylor(fpv T, const fpv* coord, fpv* uex, const tSpaceForm<fpv>& S) {
    // временно: только для гауссиана
    if(S.Form!=tSpaceForm<fpv>::FORM_GAUSSIAN) crash("PointValueTaylor: Form!=FORM_GAUSSIAN");

    fpv ll = SQR(coord[0] - S.r0[0]) + SQR(coord[1] - S.r0[1]);
    fpv l = sqrt(ll);
    fpv invBTerm = 1.0/fpv(S.Bterm);
    fpv alpha = - GetLn2<fpv>() * invBTerm * invBTerm;
    fpv f0 = exp(alpha*ll); // значение
    fpv f1_l = 2.0*alpha*f0; // 1-я производная, делённая на l
    //fpv f1 = f1_l / l; // 1-я производная
    fpv f2 = 2.0*alpha*(1.0 + 2.0*alpha*l)*f0; // 2-я производная
    //fpv f3 = 8.0*alpha*alpha*(1.0 - 2.0*alpha*l)*f0; // 3-я производная

    uex[Var_R] = uex[Var_P] = f0 + 0.5*T*T*(f2 + f1_l);
    uex[Var_W] = f0*T + T*T*T*(f2 + f1_l) / fpv(6.0); // минус волновой потенциал
    uex[Var_U] = - T * f1_l * (coord[0] - S.r0[0]);
    uex[Var_V] = - T * f1_l * (coord[1] - S.r0[1]);
}
//======================================================================================================================


// *********************************************************************************************************************
// *****                                                                                                           *****
// *****                Wave propagation from initial pulse - via integrating of Green function                    *****
// *****                                                                                                           *****
// *********************************************************************************************************************

//======================================================================================================================
template<typename fpv>
void s_IVP2D<fpv>::ReadParams(tFileBuffer& FB) {
    tParamManager PM;
    PM.Request(tSpaceForm<fpv>::r0[0], "Xterm");        // pulse center
    PM.Request(tSpaceForm<fpv>::r0[1], "Yterm");        // pulse center
    PM.Request(tSpaceForm<fpv>::Aterm, "Aterm");        // pulse amplitude
    PM.Request(tSpaceForm<fpv>::Bterm, "Bterm");        // pulse radius
    PM.Request(tSpaceForm<fpv>::Form, "Form", "");      // form of the pulse
    PM.Request(FlowVelX, "FlowVelX");                   // background flow velocity
    PM.Request(FlowVelY, "FlowVelY");                   // background flow velocity
    // Approximation parameters
    PM.Request(Hmax, "Hmax");                           // maximal radius of the pulse (it is about 0 if r>Hmax*Bterm)
    PM.Request(GI.GR, "GR");                            // number of nodes of Gaussian quadratures
    PM.Request(mm, "mm");                               // number of segments in compound Gaussian quadratures
    PM.ReadParamsFromBuffer(FB);
}
//======================================================================================================================

//======================================================================================================================
template<typename fpv>
void s_IVP2D<fpv>::Init() {
    if(tSpaceForm<fpv>::Form!=tSpaceForm<fpv>::FORM_GAUSSIAN) Hmax=1.0;
    if(Hmax <= 0.0) crash("Error: Hmax <= 0.0");

    tSpaceForm<fpv>::Init();
    GI.Init(GI.GR, GI_LEGENDRE); // init of Gauss-Legendre quadrature formula on [0,1]
}
//======================================================================================================================


//======================================================================================================================
// Purpose: calculate phi_max
// phi_max is half the central angle of a circumference of radius 'L'
//   subtended by its intersection with a disk of radius 'rad'
// The distance between the circumference and the disk is 'r0'
// If there is no intersection, returns 0.0 or -1.0
//======================================================================================================================
template<typename fpv>
inline fpv FindPhiMax(fpv rad, fpv r0, fpv L) {
    fpv num = rad*rad + r0*r0 - L*L;
    fpv denom = 2. * rad * r0;

    // we need to compute arccos(num/denom)
    if(denom < 1e-50) { // treating the case of zero denominator just for a case
        if(num > 0.0) return -1.0; // no intersection
        else return GetPiNumber<fpv>(); // full intersection
    }

    if(fabs(num) < 0.9*denom) {
        fpv cosphi = num / denom;
        return acos(cosphi);
    }
    else { // in computing arccos(z) with |z|~1 we would get a big error, so we work with sine
        const fpv& r = rad;
        const fpv& x = r0;
        fpv aux = (r+x-L)*(r+x+L)*(L+r-x)*(L+x-r);
        // The condition (aux < 0) is equivalent to (abs(num) > denom)
        // If (r ~ x >> L), we can't check (abs(num) > denom) directly due to round-up errors
        if(aux < 0.0) {
            if(num > 0.0) return -1.0; // no intersection
            else return GetPiNumber<fpv>(); // full intersection
        }
        fpv sinphi = sqrt(aux) / denom;
        fpv as = asin(sinphi);
        if(num < 0) return GetPiNumber<fpv>() - as;
        else return as;
    }
}
//======================================================================================================================


//======================================================================================================================
// Purpose: calculate d(phi_max)/dr
// phi_max is half the central angle of a circumference of radius 'r'
//   subtended by its intersection with a disk of radius 'L'
// The distance between the circumference and the disk is 'x'
// Assume: x>=L, i. e. the center of the circumference is outside the disk
// Used:
// a) for a free space problem, 'x' is the distance between the observer and the pulse center, 'L' is the pulse radius
// b) for a diffraction problem, 'x' is the distance between the observer and the corner vertex,
//      'L' is the diffraction zone radius
//======================================================================================================================
template<typename fpv>
inline fpv Find_dPhiMax(fpv r, fpv x, fpv L) {
    if(x<L) crash("Find_dPhiMax error: only for the case x >= L");

    if(r <= x-L) return huge;
    if(r >= x+L) return -huge;

    // need to compute cosphi = (r*r + x*x - L*L) / (2.*r*x);
    // from the conditions (x>=L) and (r <= x+L) we have (x >= r/2)
    // thus the denominator is equal to zero iff (r = 0), and in this case there problem has no sense
    if(r < 1e-30) return 0.0;

    fpv sinphi = sqrt((r+x-L)*(r+x+L)*(r+L-x)*(L+x-r)) / (2.*r*x);
    fpv dcosphidr = (r*r - x*x + L*L) / (2.*x*r*r);
    if(sinphi < 1e-30) return r<= x ? huge : -huge;
    return -dcosphidr / sinphi;
}
//======================================================================================================================

//======================================================================================================================
// Auxiliary subroutine for s_IVP2D::PointValueAA: calculation of the integral by angle
// Return value: intergal value; its derivative by r; its derivative by 'x'
//======================================================================================================================
template<typename fpv>
void s_IVP2D<fpv>::CalcIntegralOverPhi(fpv x, fpv r, fpv phimax, fpv* sum, int numcontur) const {
    const fpv beff = tSpaceForm<fpv>::Bterm * Hmax;
    sum[0] = sum[1] = sum[2] = 0.0;

    if(phimax < 0.) phimax = FindPhiMax(r, x, beff); // Integral limits by \phi
    if(phimax<=0.) return;
    const fpv dphi = phimax / numcontur;
    for(int j=0; j<numcontur; j++) {
        for(int i=0; i<GI.GR; i++) {
            fpv phi = (j+GI.GN[i])* dphi;
            fpv ll = r*r + x*x - 2.*r*x*cos(phi);
            fpv val[3] = {0.0, 0.0, 0.0};
            val[0] = tSpaceForm<fpv>::SpaceForm_rr(ll, &(val[1])); // Pulse value at a point given
            val[2] = val[1];
            val[1] *= 2.*(r - x*cos(phi));
            val[2] *= 2.*(x - r*cos(phi));
            for(int k=0; k<3; k++)
                sum[k] += 2.0 * val[k] * GI.GC[i] * dphi; // 2. is due to symmetry
        }
    }
}
//======================================================================================================================

//======================================================================================================================
// Auxiliary subroutine for s_IVP2D::PointValueAA: calculation of the integral by 'z'
// Return value:
//     sum[0] is the integral value;
//     sum[1] is its derivative by t (i. e. derivative by 'r' multiplied by (1-z^2))
//     sum[2] is its derivative by x
//======================================================================================================================
template<typename fpv>
void s_IVP2D<fpv>::CalcIntegralOverZ(fpv z1, fpv z2, fpv phi, fpv t, fpv x, fpv* sum, int numcontur) const {
    sum[0] = sum[1] = sum[2] = 0.0;
    const fpv dz = (z2 - z1) / numcontur;
    for(int j=0; j<numcontur; j++) {
        for(int i=0; i<GI.GR; i++) {
            fpv z = z1 + (j+GI.GN[i])* dz;
            fpv r = t * (1.0 - z*z);
            fpv ll = r*r + x*x - 2.*r*x*cos(phi);
            fpv val[3] = {0.0, 0.0, 0.0};
            val[0] = tSpaceForm<fpv>::SpaceForm_rr(ll, &(val[1])); // Pulse value at a point given
            val[2] = val[1];
            val[1] *= 2*(r - x*cos(phi))*(1-z*z);
            val[2] *= 2*(x - r*cos(phi));
            for(int k=0; k<3; k++) {
                val[k] *= (1 - z*z) / sqrt(2.0 - z*z);
                sum[k] += val[k] * dz * GI.GC[i];
            }
        }
    }
}
//======================================================================================================================


//======================================================================================================================
// Calculation of dphi_max/dz for a given z in the case x>=L
//======================================================================================================================
template<typename fpv>
fpv s_IVP2D<fpv>::GetSlope(fpv z, fpv t, fpv x) const {
    fpv L = tSpaceForm<fpv>::Bterm * Hmax; // Source radius
    fpv r = t*(1-z*z);
    if(x<L) crash("s_IVP2D::GetSlope error: only for the case x >= L");
    double ddr = Find_dPhiMax<fpv>(r, x, L);
    return (-2.*t*z) * ddr;
}

//======================================================================================================================
// Find the point z \in (zl, zr) with a given 'slope'=dphi_max/dz using dichotomy
//======================================================================================================================
template<typename fpv>
fpv s_IVP2D<fpv>::FindIntermediatePoint(fpv zl, fpv zr, fpv slope, fpv t, fpv x, int numiters) const {
    fpv vl = GetSlope(zl, t, x) - slope;
    fpv vr = GetSlope(zr, t, x) - slope;
    // If slopes at the both ends of the interval have the same sign, dichotony is not possible
    // This may occur, for instance, if the function is about zero for all z, and due to the arithmetic errors we get
    // different signs. Thus we make return and not crash
    if((vl>0.0) == (vr>0.0)) return -huge;
    for(int i=0; i<numiters; i++) {
        fpv zmid = 0.5*(zl + zr);
        fpv V = GetSlope(zmid, t, x) - slope;
        if((V>0.0)==(vl>0.0)) { zl = zmid; vl = V; }
        else { zr = zmid; vr = V; }
    }
    if(fabs(zl-zr) > 1e-8) crash("s_IVP2D::FindIntermediatePoint: process unconverged");
    return 0.5*(zl+zr);
}


//======================================================================================================================
// Calculation of wave potential for the problem of localized pulse distribution or its diffraction in corner
// Task: compute
// (-1/t) \int\limits_{|r|<t} f(|r+x|) / (2\pi\sqrt{t^2-r^2}) d^2r
// Result is 3 elements: integral value and its derivative in x and t
//======================================================================================================================
template<typename fpv>
void s_IVP2D<fpv>::PointValueAA(fpv t, const fpv* coord, fpv* u) const {
    if(GI.GN==NULL) crash("s_IVP2D::PointValueAA: init not done");
    fpv L = tSpaceForm<fpv>::Bterm * Hmax; // Source radius
    fpv x = sqrt(SQR(coord[0])+SQR(coord[1])); // Distance from source center to observer

    const int printfile = 0; // Debug output of phi_max(z)
    if(printfile){
        FILE* f = fopen("phimax.dat", "wt");
        for(double z=0.0; z<=1.0; z+=0.001) {
            double r = t*(1-z*z);
            double arg = x<tiny ? SIGN((r*r+x*x-L*L))*huge : double((r*r+x*x-L*L)/(2*r*x));
            double val = (arg<-1.0) ? PiNumber : (arg>1.0 ? 0.0 : acos(arg));
            fprintf(f, "%e %e\n", z, val);
        }
        fclose(f);
    }

    if(x > t+L) { // Signal from source has not reached the observer yet
        u[0] = u[1] = u[2] = 0.0;
        return;
    }

    if(x>L) { // Consider the case x >= L, i. e. observer is outside the source
        // Determine domain to integrate over
        const fpv frac_phimax = 0.75;
        fpv phimax_max;

        if(t*t > x*x - L*L) { // 3-zone configuration or 2-zone configuration with a local maximum
            phimax_max = asin(L/x); // maximal value of phimax
        }
        else { // 2-zone configuration without a local maximum
            phimax_max = acos((t*t+x*x-L*L)/(2*t*x)); // value of phimax at zero
        }
        fpv phimax_maxL = frac_phimax * phimax_max;
        fpv phimax_maxR = phimax_maxL;

        // Solving quadratic equation to get the inner part of the integration domain,
        // i. e. the part where the outer integral will be in z and the inner one in phi
        fpv zLL,zRR;
        {
            fpv b2 = x*cos(phimax_maxL);
            fpv d4 = b2*b2 - (x*x - L*L);
            if(d4 < 0) crash("s_IVP2D::PointValueAA: internal error, negative determinant");
            fpv rLL = b2 - sqrt(d4);
            fpv rRR = b2 + sqrt(d4); // 3-zone configuration
            if(t < x+L) rRR = t; // 2-zone configuration
            // converting to the values of z
            //fpv zL  = rL  > t ? 0.0 : sqrt(1 - rL/t);
            //fpv zR  = rR  > t ? 0.0 : sqrt(1 - rR/t);
            zLL = rLL > t ? fpv(0.0) : sqrt(1 - rLL/t);
            zRR = rRR > t ? fpv(0.0) : sqrt(1 - rRR/t);
        }

        fpv r_coor_max = sqrt(x*x - L*L); // coordinate of the point with phi=phi_max
        fpv z_coor_max = 1 - r_coor_max / t;
        if(z_coor_max >= 0.0) z_coor_max = sqrt(z_coor_max);
        else z_coor_max = 0.0;
        // looking for a point with dphi_max/dz = -1.0
        fpv zLLnew = FindIntermediatePoint(z_coor_max, 1.0, -1.0, t, x);
        if(zLLnew >= 0.0) zLL = zLLnew;

        if(t > x+L) {
            fpv myslope = 2.0;
            if(1) { //if(t < 1.1*(x+L)) { // parameter tuning
                fpv zR  = sqrt(1 - (x+L)/t);
                myslope *= phimax_max / (z_coor_max - zR);
            }
            fpv zRRnew = FindIntermediatePoint(0.0, z_coor_max, myslope, t, x);
            if(zRRnew >= 0.0) zRR = zRRnew;
        }

        { // recalculating phimax_maxL
            fpv r = t*(1-zLL*zLL);
            phimax_maxL = acos((r*r+x*x-L*L)/(2*r*x));
            r = t*(1-zRR*zRR);
            phimax_maxR = acos((r*r+x*x-L*L)/(2*r*x));
        }

        // workaround: calculate all with the outer integral in z
        //zLL = zL; zRR = zR;

        fpv summ[3] = {0.0, 0.0, 0.0};

        // считаем внутренний интеграл: вначале разбиваем отрезок, потом применяем квадратуру Гаусса
        fpv dz = (zLL - zRR) / fpv(mm);
        for(int j=0; j<mm; j++) {
            for(int i=0; i<GI.GR; i++) {
                fpv z = zRR + (j + GI.GN[i]) * dz;
                fpv r = t * (1-z*z); // обратно пересчитываем в r
                fpv val[3];
                CalcIntegralOverPhi(x, r, -1.0, val, mm);
                val[1] *= (1 - z*z); // производную по r переводим в производную по t
                for(int k=0; k<3; k++) {
                    val[k] *= 0.5; // убираем симметрию
                    val[k] *= (1 - z*z) / sqrt(2 - z*z);
                    summ[k] += dz * GI.GC[i] * val[k];
                }
            }
        }

        // считаем внешние интегралы: внешний интеграл по phi, внутренний по r
        fpv dphi = phimax_maxR / fpv(mm);
        fpv summ1[3] = {0.0, 0.0, 0.0}, summ2[3] = {0.0, 0.0, 0.0};
        for(int j=0; j<mm; j++) {
            for(int i=0; i<GI.GR; i++) {
                fpv phi = (j + GI.GN[i]) * dphi;
                // решаем квадратное уравнение для определения пределов интегрирования по r
                fpv b2 = x*cos(phi);
                fpv d4 = b2*b2 - (x*x - L*L);
                if(d4 < 0) crash("s_IVP2D::PointValueAA: internal error, negative determinant (2)");
                const fpv rrr = b2 + sqrt(d4);
                // пересчитываем в z (имеем zR <= zrr <= zRR <= zLL <= zll <= zL)
                fpv zrr = rrr > t ? fpv(0.0) : sqrt(1 - rrr/t);
                // интегрируем по z от zrr до zRR
                if(zRR - zrr < tiny) continue;
                fpv val[3];
                CalcIntegralOverZ(zrr, zRR, phi, t, x, val, mm);
                for(int k=0; k<3; k++)
                    summ1[k] += dphi * GI.GC[i] * val[k];
            }
        }

        dphi = phimax_maxL / fpv(mm);
        for(int j=0; j<mm; j++) {
            for(int i=0; i<GI.GR; i++) {
                fpv phi = (j + GI.GN[i]) * dphi;
                // решаем квадратное уравнение для определения пределов интегрирования по r
                fpv b2 = x*cos(phi);
                fpv d4 = b2*b2 - (x*x - L*L);
                if(d4 < 0) crash("s_IVP2D::PointValueAA: internal error, negative determinant (2)");
                const fpv rll = b2 - sqrt(d4);
                // пересчитываем в z (имеем zR <= zrr <= zRR <= zLL <= zll <= zL)
                fpv zll = rll > t ? fpv(0.0) : sqrt(1 - rll/t);
                // интегрируем по z от zrr до zRR
                if(zll - zLL < tiny) continue;
                fpv val[3];
                CalcIntegralOverZ(zLL, zll, phi, t, x, val, mm);
                for(int k=0; k<3; k++)
                    summ2[k] += dphi * GI.GC[i] * val[k];
            }
        }

        //pprintf("rL = %e  rLL = %e rRR = %e rR = %e\n", rL, rLL, rRR, rR);
        //pprintf("zR = %e  zRR = %e zLL = %e zL = %e\n", zR, zRR, zLL, zL);
        for(int k=0; k<3; k++)
            u[k] = -(summ[k]+summ1[k]+summ2[k]) * fpv(2.0) / GetPiNumber<fpv>();
        //pprintf("summs = %e, %e, %e\n", -(summ) * 2.0*t / PiNumber, -(summ1) * 2.0*t / PiNumber, -(summ2) * 2.0*t / PiNumber);
    }
    else { // рассматриваем случай x < L, то есть приёмник находится внутри источника
        fpv summ1[3] = {0.0, 0.0, 0.0}, summ2[3] = {0.0, 0.0, 0.0};
        // внешний интеграл по phi, внутренний по z
        const fpv phimax_max = PiNumber;
        fpv phimin = 0.0;
        if(t < L + x && t > L - x) {
            if(t < tiny || x < tiny) phimin = GetPiNumber<fpv>();
            else phimin = 0.5*(GetPiNumber<fpv>() + acos((t*t + x*x - L*L)/(2.*t*x)));
        }
        int extra_division = (t > L+x) && (t < 1.5*(L+x));
        for(int jj=0; jj<2; jj++) { // дополнительное разбиение
            fpv _phimin = phimin;
            fpv _phimax = phimax_max;
            if(extra_division) {
                if(jj) _phimin = 5.0*sqrt(1 - (x+L)/t);
                else   _phimax = 5.0*sqrt(1 - (x+L)/t);
            }
            else {
                if(jj) break;
            }

            fpv dphi = (_phimax - _phimin) / fpv(mm);
            for(int j=0; j<mm; j++) {
                for(int i=0; i<GI.GR; i++) {
                    fpv phi = _phimin + (j + GI.GN[i]) * dphi;
                    // решаем квадратное уравнение для определения пределов интегрирования по r
                    fpv b2 = x*cos(phi);
                    fpv d4 = b2*b2 - (x*x - L*L);
                    if(d4 < 0) crash("s_IVP2D::PointValueAA: internal error, negative determinant (3)");
                    const fpv rR = b2 + sqrt(d4); // второй корень отрицательный
                    // пересчитываем в z (имеем zR <= zrr <= zRR <= zLL <= zll <= zL)
                    fpv zR = rR> t ? fpv(0.0) : sqrt(1 - rR/t);
                    fpv zL = 1.0;
                    // интегрируем по z от zR до 1
                    fpv val[3] = {0.0, 0.0, 0.0};
                    if(zL - zR < tiny) continue;
                    CalcIntegralOverZ(zR, zL, phi, t, x, val, mm);
                    for(int k=0; k<3; k++)
                        summ1[k] += dphi * GI.GC[i] * val[k];
                }
            }
        }
        // внешний интеграл по z, внутренний по phi
        if(t < L+x && t > L-x) {
            // интегрирование по z ведётся от 0 до 1, но с выделением точки, в которой нарушается гладкость
            // решаем квадратное уравнение для определения пределов интегрирования по r
            fpv b2 = x*cos(phimin);
            fpv d4 = b2*b2 - (x*x - L*L);
            if(d4 < 0) crash("s_IVP2D::PointValueAA: internal error, negative determinant (4)");
            fpv rR = b2 + sqrt(d4); // второй корень отрицательный
            if(rR > t) {
                if(rR > t*1.0000001) crash("s_IVP2D::PointValueAA: internal error, rR > t");
                rR = t;
            }
            fpv zR = sqrt(1 - rR/t); // пересчитываем в z

            // первая часть интервала
            fpv dz = zR / fpv(mm);
            for(int j=0; j<mm; j++) {
                for(int i=0; i<GI.GR; i++) {
                    fpv z = (j + GI.GN[i]) * dz;
                    fpv r = t * (1-z*z); // обратно пересчитываем в r
                    fpv val[3] = {0.0, 0.0, 0.0};
                    CalcIntegralOverPhi(x, r, -1.0, val, mm);
                    val[1] *= (1 - z*z); // производную по r переводим в производную по t
                    for(int k=0; k<3; k++) {
                        val[k] *= 0.5; // 0.5 - обратно убираем симметрию
                        val[k] *= (1 - z*z) / sqrt(2 - z*z);
                        summ2[k] += dz * GI.GC[i] * val[k];
                    }
                }
            }
            // вторая часть интервала (прямоугольная)
            dz = (1.0 - zR) / fpv(mm);
            for(int j=0; j<mm; j++) {
                for(int i=0; i<GI.GR; i++) {
                    fpv z = zR + (j + GI.GN[i]) * dz;
                    fpv r = t * (1-z*z); // обратно пересчитываем в r
                    fpv val[3] = {0.0, 0.0, 0.0};
                    CalcIntegralOverPhi(x, r, phimin, val, mm);
                    val[1] *= (1 - z*z); // производную по r переводим в производную по t
                    for(int k=0; k<3; k++) {
                        val[k] *= 0.5; // 0.5 - обратно убираем симметрию
                        val[k] *= (1 - z*z) / sqrt(2 - z*z);
                        summ2[k] += dz * GI.GC[i] * val[k];
                    }
                }
            }
        }
        for(int k=0; k<3; k++)
            u[k] = -(summ1[k] + summ2[k]) * 2.0 / GetPiNumber<fpv>();
    }
    for(int k=0; k<3; k++)
        u[k] *= -1; // выводим минус волновой потенциал и его производные
}
//======================================================================================================================


//======================================================================================================================
// Calculation of the wave potential for the following problem:
// evolution of a localized pulse of an arbitrary form
//======================================================================================================================
template<typename fpv>
void s_IVP2DWP<fpv>::PointValue(fpv t, const fpv* coord, fpv* uex) const {
    fpv u[3];
    s_IVP2D<fpv>::PointValueAA(t, coord, u);
    uex[Var_R] = u[0] * t;
}


//======================================================================================================================
// Task: canclulate
// \int\limits_{|r|<t} f(|r-r0|) / \sqrt{t^2-r^2} d^2r
// Here origin is the observer point,
// r0 is the vector to the source center
//======================================================================================================================
template<typename fpv>
void s_IVP2D<fpv>::PointValue(fpv t, const fpv* coord, fpv* uex) const {
    for(int ivar=0; ivar<5; ivar++) uex[ivar] = 0.0;

    for(int iPerX = -tSpaceForm<fpv>::MaxPer[0]; iPerX <= tSpaceForm<fpv>::MaxPer[0]; iPerX++)
    for(int iPerY = -tSpaceForm<fpv>::MaxPer[1]; iPerY <= tSpaceForm<fpv>::MaxPer[1]; iPerY++) {
        if(tSpaceForm<fpv>::Checkerboard) if((iPerX+iPerY)&1) continue;
        if(fabs(tSpaceForm<fpv>::PerX) > 0.5*huge && iPerX) continue;
        if(fabs(tSpaceForm<fpv>::PerY) > 0.5*huge && iPerY) continue;
        // Adjusting coordinates taking into account the background flow and the lattice of the initial pulses
        fpv coor[2];
        coor[0] = coord[0] - tSpaceForm<fpv>::r0[0] - FlowVelX*t + tSpaceForm<fpv>::PerX * iPerX;
        coor[1] = coord[1] - tSpaceForm<fpv>::r0[1] - FlowVelY*t + tSpaceForm<fpv>::PerY * iPerY;
        fpv phi1 = GetAngle(coor[0], coor[1]);

        fpv dp, dur;
        if(1) {
            fpv summs[3];
            PointValueAA(t, coor, summs); // -W/t; derivative by t of (-W/t), derivative by x of (-W/t)

            dp = summs[0] + summs[1] * t;
            dur = - summs[2] * t;
        }
        else {
            // debug code: numerical differentiation of the wave potential
            fpv eps = 1e-4;

            fpv summs[3], summ1, summ2;
            PointValueAA(t, coor, summs); summ1 = summs[0] * t;
            PointValueAA(t+eps, coor, summs); summ2 = summs[0] * (t+eps);
            dp = (summ2 - summ1) / eps;
            coor[0] += eps * cos(phi1);
            coor[1] += eps * sin(phi1);
            PointValueAA(t, coor, summs); summ2 = summs[0] * t;
            dur = - (summ2 - summ1) / eps;
        }
        uex[Var_R] += dp;
        uex[Var_P] += dp;
        uex[Var_U] += dur * cos(phi1);
        uex[Var_V] += dur * sin(phi1);
    }
}
//======================================================================================================================



// *********************************************************************************************************************
// *****                                                                                                           *****
// *****                        Planar wave diffraction inside a sector with angle 2*Pi/n                          *****
// *****                                                                                                           *****
// *********************************************************************************************************************

//======================================================================================================================
template<typename fpv>
void s_CornerPlanar<fpv>::ReadParams(tFileBuffer& FB) {
    tParamManager PM;
    PM.Request(angle, "angle"); // angle (can be only 2*Pi/N with natural N). If N is even, then is no diffraction
    PM.Request(phi0, "phi0");   // direction the incoming wave comes from
    PM.Request(X0, "X0");       // distance between the Gaussian center at t=0 and the corner vertex
    PM.Request(Bterm, "Bterm"); // half-width of the incoming Gaussian
    PM.Request(Aterm, "Aterm"); // amplitude of the incoming Gaussian
    PM.ReadParamsFromBuffer(FB);
}
//======================================================================================================================

//======================================================================================================================
template<typename fpv>
void s_CornerPlanar<fpv>::Init() {
    if(angle <= 1e-10 || angle > Pi2 + 1e-10)
        crash("s_CornerPlanar error: 2*Pi/angle should be a natural number");

    // Identifying of the rational number
    double _n = Pi2 / angle;
    n = round_noisette(_n);
    if(fabs(_n - double(n)) > 1e-10)
        crash("s_CornerPlanar error: 2*Pi/angle should be a natural number");

    // If angle = Pi/n, there is no diffraction
    if(n&1) m = 2;             // 2*Pi/n
    else { n = n / 2; m = 1; } // Pi/n

    // Check that the incoming wave is inside the sector
    if(phi0 < 0 || phi0 > angle)
        crash("s_CornerPlanar error: phi0 is outside the corner");

    // Init of Gaussian quadrature rule
    H = sqrt(-2.*log(0.5*NativeDouble(get_eps<fpv>())));
    GJI.Init(int(0.58*H*H)+1, GI_JACOBI1);

    // Just constants
    pi = GetPiNumber<fpv>();
    pi2 = GetPiNumber2<fpv>();
    const_sqrt_2ln2 = sqrt(2.*log(fpv(2.0)));
}
//======================================================================================================================

//======================================================================================================================
// Calculation of  E(y) = int_{0}^{infty} exp(-(x-y)^2/2) x^{-1/2} dx
//======================================================================================================================
template<typename fpv>
fpv calc_E(fpv y, const tGaussIntegrator<fpv>& GJI) {
    NativeDouble H = sqrt(-2.*log(0.5*NativeDouble(get_eps<fpv>())));
    if(H*NativeDouble(y) > H*H+1.) {
        // x=0 is outside the essential support of the Gaussian
        // rewrite E(y) = int_{-infty}^{infty} exp(-x^2/2) f(x) dx, f(x) = Heaviside(x+y)/sqrt(x+y)
        // and use the integration formula with a uniform step
        int n = int(H*H/Pi2)+1;
        fpv h = fpv(H/n);
        fpv h2_2 = 0.5*h*h;
        fpv E = 1./sqrt(y);
        for(int j=1; j<=n; j++) {
            E += exp(-(j*j)*h2_2) * (1./sqrt(y+j*h) + 1./sqrt(y-j*h));
        }
        return E * h;
    }
    else {
        // x=0 is inside the essential support of the Gaussian
        // use the Gauss -- Jacobi quadrature rule
        fpv b = y + H;
        if(b<=0.0) return 0.0;
        fpv E = 0.0;
        for(int i=0; i<GJI.GR; i++) {
            fpv tmp = b*GJI.GN[i] - y;
            E += GJI.GC[i] * exp(-0.5*tmp*tmp);
        }
        return E * sqrt(b);
    }
}

//======================================================================================================================
// Calculating J0 = int f0(eta) d eta, from 0 to infinity, where
// f0(eta) = exp(-(alpha*eta-beta)^2 / 2) / sqrt(eta) / (eta + 1) / Pi
//======================================================================================================================
template<typename fpv>
fpv calc_J0(fpv alpha, fpv beta, const tGaussIntegrator<fpv>& GJI) {
    if(alpha < get_eps<fpv>()*get_eps<fpv>()) {
        return exp(-0.5*beta*beta);
    }

    NativeDouble H = sqrt(-2.*log(0.5*NativeDouble(get_eps<fpv>())));
    if(H*beta > H*H+1.0) { // beta > H+1/H
        if(alpha < 1e-100) return 0.0;
        const fpv inv_alpha = 1.0 / alpha;
        // use the uniform step
        fpv J0 = 0.0;
        int n = int(H*H/Pi2);
        fpv h = fpv(H/n);
        for(int i=-n; i<n; i++) {
            fpv x = i*h;
            fpv eta = (x+beta)*inv_alpha;
            J0 += exp(-0.5*x*x) / (sqrt(eta)*(eta+1));
        }
        return J0*h*inv_alpha / GetPiNumber<fpv>();
    }
    else {
        // use the Gauss -- Jacobi formula
        fpv b = (beta+H) / alpha;
        if(b<=0.0) return fpv(0.0);
        fpv J0 = 0.0;
        fpv JJ = 0.0;
        for(int i=0; i<GJI.GR; i++) {
            fpv eta = b*GJI.GN[i];
            fpv tmp = eta*alpha - beta;
            eta = GJI.GC[i] / (eta+1.0);
            J0 += eta * exp(-0.5*tmp*tmp);
            JJ += eta;
        }
        fpv sqrt_b = sqrt(b);
        J0 *= sqrt_b;
        // 0.2071 ~= 0.5*(sqrt(2.0)-1.0)
        if(!(alpha+beta>H || alpha>0.2071*(beta+H))) {
            JJ *= sqrt_b;
            JJ -= 2.0*atan(sqrt_b);
            JJ *= exp(-0.5*(alpha+beta)*(alpha+beta));
            J0 -= JJ;
        }
        return J0 / GetPiNumber<fpv>();
    }
}

//======================================================================================================================
template<typename fpv>
void s_CornerPlanar<fpv>::PointValue(fpv t, const fpv* coord, fpv* uex) const {
    if(!n) crash("s_CornerPlanar error: init not done");
    for(int ivar=Var_R; ivar<=Var_P; ivar++) uex[ivar] = 0.0;

    // Special case: wave comes parallel to a corner ray
    // Allow to have a reflection of the sector and reflected solution there
    if(fabs(phi0) < 1e-50 && !(m==2 && n==1) && coord[1]<0.0) {
        fpv coord_new[3] = {coord[0], -coord[1], 0.0};
        PointValue(t, coord_new, uex);
        uex[Var_V] = -uex[Var_V];
        return;
    }

    // Direction to the observer
    fpv r = sqrt(SQR(coord[0]) + SQR(coord[1]));
    fpv ref_Phi = GetAngle(coord[0], coord[1]);
    if(ref_Phi > angle+tinyflt) return; // Error: point is outside the corner

    const fpv alpha0 = const_sqrt_2ln2 / Bterm;
    for(int iwave = 0; iwave < 2*n; iwave++) { // Loop over waves (i. e. over reflections of the incoming wave)
        // Detect the wave direction on the Riemann surface
        fpv phi_j = phi0 + (iwave>>1) * 2.0*angle;
        if(iwave&1) phi_j = pi2*m - phi_j;

        fpv psi = ref_Phi - phi_j; // (-4*pi, 2*pi)
        if(psi < 0.0) psi += pi2*m; // (0, 4*pi)
        int gpsi = (psi<pi || psi > 3.*pi);

        // 1. Incident (or reflected) wave
        fpv form = 0.0;
        if(gpsi) {
            fpv phase = t - X0 + r * cos(psi);
            form = exp(- 0.5*SQR(alpha0*phase));
        }

        // 2. Diffraction term #1
        if(m==2) { // if m==1, then the angle of the corner is pi/n => no diffraction
            fpv a = t - r - X0;
            fpv b = cos(0.5*psi); b = 2*r*b*b;
            fpv j0 = calc_J0(alpha0*b, alpha0*a, GJI);
            form += j0 * (0.5-gpsi);
        }

        uex[Var_P] += form;
        uex[Var_U] -= form * cos(phi_j);
        uex[Var_V] -= form * sin(phi_j);

        // 3. Diffraction term #2 (~ 1/sqrt(r)) -- only for the half-line case
        if(m==2 && n==1 && r>1e-200) {
            fpv delta = alpha0*(t - r - X0);
            fpv E = calc_E(delta, GJI) / sqrt(2.*alpha0*pi*pi*r);
            fpv dur = E * cos(0.5*psi);
            fpv duphi = - E * sin(0.5*psi);
            uex[Var_U] -= dur * cos(ref_Phi) - duphi * sin(ref_Phi);
            uex[Var_V] -= duphi * cos(ref_Phi) + dur * sin(ref_Phi);
        }
    }
    uex[Var_U] *= Aterm;
    uex[Var_V] *= Aterm;
    uex[Var_P] *= Aterm;
    uex[Var_R] = uex[Var_P];
}
//======================================================================================================================



// *********************************************************************************************************************
// *****                                                                                                           *****
// *****                        Gaussian pulse diffraction inside a sector with angle 2*Pi/n                       *****
// *****                                                                                                           *****
// *********************************************************************************************************************


//======================================================================================================================
template<typename fpv>
void s_Corner<fpv>::ReadParams(tFileBuffer& FB) {
    tParamManager PM;
    PM.Request(tSpaceForm<fpv>::r0[0], "Xterm");        // pulse center
    PM.Request(tSpaceForm<fpv>::r0[1], "Yterm");        // pulse center
    PM.Request(tSpaceForm<fpv>::Aterm, "Aterm");        // pulse amplitude
    PM.Request(tSpaceForm<fpv>::Bterm, "Bterm");        // pulse radius
    PM.Request(angle, "angle"); // corner angle. Must be 2*PI/n with a natural n
    // Approximation parameters
    PM.Request(Hmax, "Hmax");                           // maximal radius of the pulse (it is about 0 if r>Hmax*Bterm)
    PM.Request(GI.GR, "GR");                            // number of nodes of Gaussian quadratures
    PM.Request(mm, "mm");                               // number of segments in compound Gaussian quadratures
    PM.ReadParamsFromBuffer(FB);
}
//======================================================================================================================

//======================================================================================================================
template<typename fpv>
void s_Corner<fpv>::Init() {
    tSpaceForm<fpv>::Init();
    GI.Init(GI.GR, GI_LEGENDRE); // init of Gauss-Legendre quadrature formula on [0,1]

    if(tSpaceForm<fpv>::Form != tSpaceForm<fpv>::FORM_GAUSSIAN) crash("s_Corner::Init() error: only Gaussian form is supported");
    if(Hmax <= 0.0) crash("Error: Hmax <= 0.0");

    int n = round_noisette(double(Pi2 / angle));
    if(fabs(angle*n - GetPiNumber2<fpv>()) > 1e-10)
        crash("s_Corner error: 2*Pi/angle should be a natural number");
    angle = GetPiNumber2<fpv>() / fpv(n); // пересчитываем angle на случай, если число было задано с погрешностью

    fpv dsh = sqrt(SQR(MIN(tSpaceForm<fpv>::r0[0],fpv(0.0))) + SQR(tSpaceForm<fpv>::r0[1]));
    if(dsh < tSpaceForm<fpv>::Bterm * Hmax) {
        fpv HmaxNew = dsh / tSpaceForm<fpv>::Bterm;
        if(HmaxNew < 3.0) crash("Error: impulse is too near to the halfline (dist/Bterm = %f)", double(HmaxNew));
        pprintf0("Warning: impulse is too near to the halfline:\n  Hmax given = %f, fact = dist/Bterm = %f\n", double(Hmax), double(HmaxNew));
        Hmax = HmaxNew;
    }
}
//======================================================================================================================


//======================================================================================================================
// Intersection of intervals (a,b) and (c,d) on a circle
// Assumptions: 0 < b-a < 2*Pi, 0 < d-c < 2*Pi
// Returns the intervals in the intersection (0 or 1 or 2) and their boundaries (up to 4 numbers)
//======================================================================================================================
template <class fpv>
int IntersectArcs(fpv a, fpv b, fpv c, fpv d, fpv* limits) {
    // if the second interval starts before the first one, swap the intervals
    if(c < a) { SWAP(a, c); SWAP(b, d); }
    // reduce the beginning of the second interval to (a, 2*Pi + a)
    const fpv _2pi = GetPiNumber2<fpv>(); // 2*Pi
    while(c >= _2pi + a) { c -= _2pi; d -= _2pi; }

    if(b < c) {
        d -= _2pi;
        if(d < a) return 0;
        *(limits++) = a; *(limits++) = MIN(b, d);
        return 1;
    }
    else {
        if(d < b) { *(limits++) = c; *(limits++) = d; return 1; }
        *(limits++) = c; *(limits++) = b;
        if(d > _2pi+a) { *(limits++) = _2pi+a; *(limits++) = d; return 2; }
        else return 1;
    }
}
//======================================================================================================================


//======================================================================================================================
// Auxiliary subroutine
// u[0] -- value of the integral
// u[1] -- value of the integral of the time derivative
// u[2] -- value of the integral of the derivative in r
// u[3] -- corrected value of the integral of the derivative in phi
//======================================================================================================================
template<typename fpv>
void s_Corner<fpv>::CalcIntOverZ_NEW(fpv z1, fpv z2, fpv alpha, fpv t, fpv r, fpv rs, fpv hatPhi, fpv* u) const {
    u[0] = u[1] = u[2] = u[3] = 0.0;
    fpv cosalpha = cos(alpha), sinalpha = sin(alpha);
    const fpv coshatPhi = cos(hatPhi), sinhatPhi = sin(hatPhi);
    const fpv xs = r + rs * coshatPhi, ys = rs * sinhatPhi;

    // пытаемся ограничить область интегрирования
    fpv limits[2] = {z1, z2};
    if(z2 < z1) crash("z2 < z1");

    if(1) {
        fpv L = tSpaceForm<fpv>::Bterm*this->Hmax;
        fpv dist2 = xs*xs+ys*ys;
        const fpv b2 = -cos(alpha - GetAngle(xs, ys)) * sqrt(dist2);
        fpv d4 = b2*b2 - (dist2 - L*L);
        if(d4<0.0) return; // луч при данном phi не пересекает носитель начальных данных
        d4 = sqrt(d4);
        limits[1] = MAX(t*(1.0-z2*z2), -b2-d4);
        limits[0] = MIN(t*(1.0-z1*z1), -b2+d4);
        // область интегрирования по r от limits[1] до limits[0]
        if(limits[0] < limits[1]) return; // интервалы не пересекаются
        limits[0] = sqrt(1.0 - limits[0]/t);
        limits[1] = sqrt(1.0 - limits[1]/t);
    }

    const fpv dz = (limits[1]-limits[0]) / this->mm;
    for(int jjj=0; jjj<this->mm; jjj++) {
        for(int iii=0; iii<this->GI.GR; iii++) {
            fpv z = limits[0] + (jjj+this->GI.GN[iii])* dz;
            fpv R = 1.0 - z*z;

            fpv deltax = t*R*cosalpha - xs;
            fpv deltay = t*R*sinalpha - ys;
            fpv ll = SQR(deltax) + SQR(deltay);

            fpv val[4], dF;
            val[0] = tSpaceForm<fpv>::SpaceForm_rr(ll, &dF); // Амплитуда источника в заданной точке
            val[1] = 2.*dF * R * (deltax * cosalpha + deltay * sinalpha);
            val[2] = - 2.*dF * deltax;
            val[3] = 2.*dF * deltay; // multiplication by sign is in calling subroutine
            fpv mult = R / (GetPiNumber<fpv>() * sqrt(1.0 + R)) * dz * this->GI.GC[iii];
            for(int k=0; k<4; k++)
                u[k] += val[k] * mult;
        }
    }
}


//======================================================================================================================
// Auxiliary subroutine
// u[0] -- value of the integral
// u[1] -- value of the integral of the time derivative
// u[2] -- value of the integral of the derivative in r
// u[3] -- corrected value of the integral of the derivative in phi
// Input:
// phi1, phi2 -- integration limits
// t, R, r, rs, hatPhi -- parameters (see math. documentation)
//======================================================================================================================
template<typename fpv>
void s_Corner<fpv>::CalcIntOverPhi_NEW(fpv phi1, fpv phi2, fpv t, fpv R, fpv r, fpv rs, fpv hatPhi, fpv* u) const {
    u[0] = u[1] = u[2] = u[3] = 0.0;

    const fpv coshatPhi = cos(hatPhi), sinhatPhi = sin(hatPhi);
    const fpv xs = r + rs * coshatPhi, ys = rs * sinhatPhi;
    const fpv PHI1 = GetAngle(xs, ys);

    // пытаемся ограничить область интегрирования
    fpv limits[4] = {phi1, phi2};
    int numlimits = 1;
    if(1) {
        fpv L = tSpaceForm<fpv>::Bterm*this->Hmax;
        fpv dist2 = xs*xs+ys*ys;
        fpv num = -(L*L - dist2 - t*t*R*R);
        fpv denom = 2.*t*R*sqrt(dist2);
        fpv cosmax;
        if(denom < fabs(num)) cosmax = SIGN(num);
        else cosmax = num / denom;
        fpv acosine = acos(cosmax);
        numlimits = IntersectArcs(PHI1 - acosine, PHI1 + acosine, phi1, phi2, limits);
    }
    for(int il=0; il<numlimits; il++) {
        const fpv dalpha = (limits[il*2+1]-limits[il*2]) / this->mm;
        for(int jjj=0; jjj<this->mm; jjj++) for(int iii=0; iii<this->GI.GR; iii++) {
            fpv alpha = limits[il*2] + (jjj+this->GI.GN[iii])* dalpha;
            fpv cosalpha = cos(alpha), sinalpha = sin(alpha);

            fpv deltax = t*R*cosalpha - xs;
            fpv deltay = t*R*sinalpha - ys;
            fpv ll = SQR(deltax) + SQR(deltay);

            fpv val[4] = {0.0, 0.0, 0.0, 0.0}, dF;
            val[0] = tSpaceForm<fpv>::SpaceForm_rr(ll, &dF); // Impulse form in a point given
            val[1] = 2.*dF * R * (deltax * cosalpha + deltay * sinalpha);
            val[2] = -2.*dF * deltax;
            val[3] = 2.*dF * deltay; // multiplication by sign is in calling subroutine
            for(int k=0; k<4; k++)
                u[k] += val[k] * this->GI.GC[iii] * dalpha;
        }
    }
}

//======================================================================================================================
// Auxiliary subroutine
// Returns 1 if no intersection
// mode: 0 - intersection with circumference (for boundary integral), 1 - intersection with disk (for 2D integral)
//======================================================================================================================
template<typename fpv>
int BoundPhiZone(fpv r0, fpv rs, fpv L, fpv hatPhi, fpv& psimin, fpv& psimax, int mode) {
    if(mode==1 && rs < L) return 0; // full intersection
    fpv a = SQR(r0) + SQR(rs);
    fpv b = 2.0 * r0 * rs;
    fpv amb = SQR(r0 - rs); // a-b
    fpv L2 = L*L; // L^2
    fpv buf;
    if(mode==1 && r0*r0 > rs*rs - L2) {
        buf = asin(L / rs); // if rs < L then there was a return statement above
    }
    else {
        buf = (a - L2) / (b + get_min_value<fpv>());
        if(buf > 1.0) return 1; // no intersection
        if(buf <= -1.0) return 0; // full intersection

        if(fabs(buf) < 0.5) buf = acos(buf); // for ||buf|-1| << 1 we have troubles with machine precision here
        else { // so we make a regularization
            buf = (b+a-L2)*(L2-amb);
            if(buf < 0.0) buf = 0.0;
            buf = sqrt(buf) / (b + get_min_value<fpv>());
            if(buf > 1.0) buf = 1.0;
            buf = asin(buf);
            if(a < L2) buf = -buf;
        }
    }

    fpv limits[4];
    int numlimits = IntersectArcs<fpv>(0.0, GetPiNumber<fpv>(), hatPhi - buf, hatPhi + buf, limits);
    if(numlimits==0) return 1;
    if(numlimits==2) crash("Internal error: numlimits==2");
    psimin = limits[0]; psimax = limits[1];

    const fpv _Pi2 = GetPiNumber2<fpv>();
    while(psimin<0.0) { psimin += _Pi2; psimax += _Pi2; } // otherwise a sign of sin(0.5*psi) will be wrong
    while(psimin>_Pi2) { psimin -= _Pi2; psimax -= _Pi2; } // otherwise a sign of sin(0.5*psi) will be wrong
    return 0;
}
//======================================================================================================================


//======================================================================================================================
// Main subroutine to calculate the solution for Gaussian pulse on half-line diffraction
// Input:
//   t, coor -- time and coordinates
//   inv     -- 0: cancluation of the signal from a real source, 1: from an imaginary source, -1: from all sources
//   clear   -- 0: write to the output array, 1: increment to the output array
// Output:
//   u[0] -- Q = W/t, where W is the wave potential
//   u[1] -- dQ/dt
//   u[2] -- dQ/d(coord[0])
//   u[3] -- dQ/d(coord[1])
// If uu!=NULL, then for each u[i] returns its four terms (16 values total)
//======================================================================================================================
template<typename fpv>
void s_Corner<fpv>::MainFunc(fpv t, const fpv* coor, fpv* u, int inv, int clear, fpv* _uu) const {
    const fpv tiny_loc = sizeof(fpv)==8 ? tiny : 1e-30;
    const fpv _Pi2 = GetPiNumber2<fpv>();
    const fpv _PiNumber = GetPiNumber<fpv>();

    if(clear) {
        u[0] = u[1] = u[2] = u[3] = 0.0;
        if(_uu != NULL) for(int i=0; i<16; i++) _uu[i] = 0.0;
    }

    // Определяем рациональное представление угла: 2*pi/n
    int n = round_noisette(double(_Pi2 / angle));
    if(inv==-1) {
        for(int i=0; i<2*n; i++)
            MainFunc(t, coor, u, i, 0, _uu);
        if(fabs(coor[0])+fabs(coor[1]) < tiny_loc) u[2] = u[3] = 0.0;
        return;
    }
    if(inv<0 || inv>=2*n) crash("s_Corner<fpv>::MainFunc error: wrong source number %i", inv);

    fpv uu[16]; // основной массив для результатов. Старший индекс = номер производной, младший = номер составляющей
    for(int i=0; i<16; i++) uu[i] = 0.0;

    // Расстояние от наблюдателя до вершины полупрямой
    const fpv dist_observer_vertex = sqrt(coor[0]*coor[0]+coor[1]*coor[1]);
    // Направление из вершины полупрямой на наблюдателя
    fpv phi_observer = GetAngle(coor[0], coor[1]); // 0 <= phi_observer < 2*Pi. В тексте: \phi
    // Корректируем положение наблюдателя на случай, если он из-за какой-нибудь погрешности оказался чуть-чуть вне угла
    if(n>1 && phi_observer > angle) {
        if(phi_observer > 1.5*PiNumber) phi_observer = 0.0;
        else phi_observer = angle;
    }
    fpv cosphi = 1.0, sinphi = 0.0;
    if(dist_observer_vertex > get_min_value<fpv>()) {
        cosphi = coor[0] / dist_observer_vertex;
        sinphi = coor[1] / dist_observer_vertex;
    }

    // Положение рассматриваемого (действительного или мнимого) начального импульса
    // -- переводим центр импульса в полярные координаты
    const fpv dist_src_vertex = sqrt(SQR(tSpaceForm<fpv>::r0[0]) + SQR(tSpaceForm<fpv>::r0[1])); // In text: rs
    fpv phis = GetAngle(tSpaceForm<fpv>::r0[0], tSpaceForm<fpv>::r0[1]);
    // -- поворачиваем и отражаем. Получаем угол в диапазоне (0, 4*Pi)
    phis += fpv(inv >> 1) * fpv(2.0*angle);
    if(inv&1) phis = 4.0*_PiNumber - phis;

    // Разность углов между направлением на наблюдателя и направлением на центр мнимого источника. В тексте: \Phi
    fpv Phi = phi_observer - phis; // (-4*Pi, 2*Pi)
    if(Phi < 0.0) Phi += 4.*_PiNumber; // (0, 4*Pi)
    // Определяем точку Pi*J, лежащую в диапазоне (Phi-Pi/2, Phi+Pi/2). В тексте: J
    // По J определяем геометрические множители. В тексте: g_+, g_-
    int J = 4, gp = 1, gm = 1;
    if(Phi < 0.5*_PiNumber) J = 0;
    else if(Phi < 1.5*_PiNumber) { J = 1; gp = 0; }
    else if(Phi < 2.5*_PiNumber) { J = 2; gp = gm = 0; }
    else if(Phi < 3.5*_PiNumber) { J = 3; gm = 0; }

    // Радиус дифракционной зоны, т. е. множества точек,
    // расположенный в которых импульс к моменту времени t достигнет наблюдателя через вершину полупрямой
    fpv difr_zone_radius = t - dist_observer_vertex;
    // Ширина источника (т. е. расстояние от его центра, на котором наличием импульса можно пренебречь
    const fpv L = tSpaceForm<fpv>::Bterm*this->Hmax;

    // -- возвращаемся в декартовые координаты
    const fpv SXterm = dist_src_vertex * cos(phis);
    const fpv SYterm = dist_src_vertex * sin(phis);
    // Расстояние от наблюдателя до центра начального импульса
    const fpv dist = sqrt(SQR(coor[0] - SXterm) + SQR(coor[1] - SYterm));
    if(t+L <= dist) return;

    // считаем первый интеграл - для подобластей phi>0 и phi<0
    for(int ilight = 0; ilight<=1; ilight++) {
        if((ilight ? gm : gp)==0) continue; // dark zone - skipping
        fpv sum1[4] = {0.0, 0.0, 0.0, 0.0};

        fpv rmin = MAX(fpv(0.0), dist-L);
        fpv rmax = MIN(t, dist+L);
        fpv zmin = sqrt(1.0 - rmax/t);
        fpv zmax = sqrt(1.0 - rmin/t);

        const fpv PhiHatSign = (ilight ^ (J&1)) ? -1.0 : 1.0;
        const fpv hatPhi = _PiNumber - PhiHatSign * Phi; // In text: \hat{\Phi}

        const fpv dz = (zmax - zmin)/ fpv(mm);
        for(int ii=0; ii<mm; ii++) for(int iii=0; iii<GI.GR; iii++) {
            fpv z = (ii+GI.GN[iii]) * dz + zmin;
            fpv R = 1 - z*z; // In text: R
            fpv mult = R / (sqrt(1+R)*_PiNumber) * GI.GC[iii]*dz;

            fpv summLoc[4]; // wave potential, its derivative on t, on r and on phi
            CalcIntOverPhi_NEW(0, _PiNumber, t, R, dist_observer_vertex, dist_src_vertex, hatPhi, summLoc);
            summLoc[3] *= -PhiHatSign;
            for(int k=0; k<4; k++) sum1[k] += summLoc[k] * mult;
        }
        uu[0+ilight]  += sum1[0]; // wave potential over t
        uu[4+ilight]  += sum1[1]; // time derivative of (wave potential over t)
        // converting derivatives over r and phi to the derivatives over source Cartesian coordinates
        uu[8+ilight]  += cosphi * sum1[2] - sinphi * sum1[3];
        uu[12+ilight] += sinphi * sum1[2] + cosphi * sum1[3];
    }

    // если дифракционная зона отсутствует, всё уже сделано
    if(t <= dist_observer_vertex) goto fin;

    if(t < 2.25*dist_observer_vertex) { // t < 2.25*dov
        // считаем второй интеграл - для подобластей phi>0 и phi<0
        fpv sumPlus2[4] = {0.0, 0.0, 0.0, 0.0};
        fpv sumMinus2[4] = {0.0, 0.0, 0.0, 0.0};

        // интеграл по области r > dov, где r -- расстояние от наблюдателя
        // считаем интеграл вначале по phi, потом по z
        // определяем границу области r > dov
        fpv zM = sqrt(1.0 - dist_observer_vertex / t);

        {
            // делаем поправку на конечность носителя источника
            fpv zmin = 0.0, zmax = zM;
            if(dist+L < t) zmin = sqrt(1. - (dist+L)/t);
            if(dist > L) {
                fpv z = sqrt(1. - (dist-L)/t);
                zmax = MIN(zM, z);
            }
            if(zmax < zmin) zmax = zmin;

            fpv dz = (zmax - zmin) / fpv(this->mm);
            for(int j=0; j<this->mm; j++) for(int i=0; i<this->GI.GR; i++) {
                fpv z = zmin + (j + this->GI.GN[i]) * dz;
                fpv R = 1-z*z;

                if(R >= 1.0 || t*R <= dist_observer_vertex*2 - t) continue;
                fpv phimax = FindPhiMax(dist_observer_vertex, t*R, difr_zone_radius);
                if(phimax < 0.0) continue;

                fpv mult = 2. / (_Pi2 * sqrt(1.+R)) * R * dz * this->GI.GC[i];
                for(int ilight=0; ilight<=1; ilight++) {
                    fpv summLoc[4];
                    const fpv PhiHatSign = (ilight ^ (J&1)) ? -1.0 : 1.0;
                    const fpv hatPhi = _PiNumber - PhiHatSign * Phi; // In text: \hat{\Phi}
                    CalcIntOverPhi_NEW(0.0, phimax, t, R, dist_observer_vertex, dist_src_vertex, hatPhi, summLoc);
                    summLoc[3] *= -PhiHatSign;
                    for(int k=0; k<4; k++)
                        (ilight ? sumMinus2 : sumPlus2)[k] += summLoc[k] * mult;
                }
            }
        }

        fpv sumPM2a[2][4] = {{0.,0.,0.,0.}, {0.,0.,0.,0.}};
        for(int ilight=0; ilight<2; ilight++) {
            // интеграл по области r < dov, где r -- расстояние от наблюдателя
            // пусть x -- направление от наблюдателя на вершину полупрямой, y -- перпендикулярное ему
            // считаем интеграл вначале по x, потом по y

            // определяем предел интегрирования по y
            if(dist_observer_vertex <= 0.5*difr_zone_radius) crash("s_Corner::MainFunc: dist_observer_vertex <= 0.5*difr_zone_radius");
            fpv ymax_dz = (difr_zone_radius / dist_observer_vertex) * sqrt(SQR(dist_observer_vertex) - 0.25*SQR(difr_zone_radius));

            // корректируем с учётом финитности начальных данных
            fpv ymin = 0.0, ymax = ymax_dz;

            const fpv PhiHatSign = (ilight ^ (J&1)) ? -1.0 : 1.0;
            const fpv hatPhi = _PiNumber - PhiHatSign * Phi; // In text: \hat{\Phi}

            const fpv coshatPhi = cos(hatPhi), sinhatPhi = sin(hatPhi);
            const fpv r = dist_observer_vertex, rs = dist_src_vertex;
            const fpv xs = r + rs * coshatPhi, ys = rs * sinhatPhi;

            if(1) {
                ymin = MAX(ymin, ys - L);
                ymax = MIN(ymax, ys + L);
                if(ymax < ymin) continue;
            }

            fpv dy = (ymax - ymin) / fpv(this->mm);
            for(int jy=0; jy<this->mm; jy++) for(int iy=0; iy<GI.GR; iy++) {
                fpv y = ymin + (jy + GI.GN[iy]) * dy;

                // определяем пределы интегрирования по x
                fpv xmin = dist_observer_vertex - sqrt(SQR(t - dist_observer_vertex) - y*y);
                fpv xmax = sqrt(SQR(dist_observer_vertex) - y*y);

                // корректируем с учётом финитности начальных данных
                {
                    fpv Lloc = sqrt(L*L - SQR(y - ys));
                    xmin = MAX(xmin, xs - Lloc);
                    xmax = MIN(xmax, xs + Lloc);
                    if(xmax < xmin) continue;
                }

                fpv dx = (xmax - xmin) / fpv(this->mm);
                fpv _sum[4] = {0.,0.,0.,0.};
                for(int jx=0; jx<this->mm; jx++) for(int ix=0; ix<this->GI.GR; ix++) {
                    fpv x = xmin + (jx + this->GI.GN[ix]) * dx;

                    fpv deltax = x - xs, deltay = y - ys;
                    fpv ll = SQR(deltax) + SQR(deltay);

                    fpv dF;
                    fpv F = tSpaceForm<fpv>::SpaceForm_rr(ll, &dF);
                    fpv rr = x*x + y*y;
                    fpv sq = sqrt(t*t - rr);

                    fpv val[4];
                    fpv R = sqrt(rr) / t;
                    fpv cosalpha = x / (sqrt(rr) + tiny_loc);
                    fpv sinalpha = y / (sqrt(rr) + tiny_loc);
                    val[0] = F; // Impulse form in a point given
                    val[1] = 2.*dF * R * (deltax * cosalpha + deltay * sinalpha);
                    val[2] = -2.*dF * deltax;
                    val[3] = 2.*dF * deltay;
                    val[3] *= -PhiHatSign;

                    for(int k=0; k<4; k++)
                        _sum[k] += val[k] / (t * sq) * this->GI.GC[ix] * dx;
                }

                for(int k=0; k<4; k++)
                    sumPM2a[ilight][k] += _sum[k] * this->GI.GC[iy] * dy;
            }
        }

        for(int k=0; k<4; k++) {
            uu[k*4+2] = (sumPlus2[k]  + sumPM2a[0][k] / _Pi2) * (0.5 - gp);
            uu[k*4+3] = (sumMinus2[k] + sumPM2a[1][k] / _Pi2) * (0.5 - gm);
        }
        // converting derivatives over r and phi to the derivatives over source Cartesian coordinates
        for(int ilight=2; ilight<4; ilight++) {
            fpv sum1[4] = {0.,0.,uu[8+ilight],uu[12+ilight]};
            uu[8+ilight]  = cosphi * sum1[2] - sinphi * sum1[3];
            uu[12+ilight] = sinphi * sum1[2] + cosphi * sum1[3];
        }
    }
    else { // t >= 2.25*dist_observer_vertex
#if 0
        // определяем половину угла, под которым носитель начальных данных виден из наблюдателя
        fpv _phi_ = GetPiNumber<fpv>(); // pi, если наблюдатель внутри носителя
        if(dist > L) _phi_ = asin(L / dist); // от 0 до pi/2
        for(int ilight=0; ilight<2; ilight++) {
            fpv sum2[4] = {0.,0.,0.,0.};

            const fpv PhiHatSign = (ilight ^ (J&1)) ? -1.0 : 1.0;
            const fpv hatPhi = _PiNumber - PhiHatSign * Phi; // In text: \hat{\Phi}

            // считаем внешние интегралы: внешний интеграл по phi, внутренний по r
            fpv alphamin = 0.0, alphamax = _PiNumber;
            fpv limits[4] = {alphamin, alphamax, 0., 0.};
            int numlimits = 1;
            {
                const fpv coshatPhi = cos(hatPhi), sinhatPhi = sin(hatPhi);
                const fpv r = dist_observer_vertex, rs = dist_src_vertex;
                const fpv xs = r + rs * coshatPhi, ys = rs * sinhatPhi;
                const fpv PHI1 = GetAngle(xs, ys);
                numlimits = IntersectArcs<fpv>(alphamin, alphamax, PHI1-_phi_, PHI1+_phi_, limits);
            }
            const int nmm = s_IVP2D<fpv>::mm;
            //const int nmm = 40;
            for(int ilim=0; ilim<numlimits; ilim++) {
                fpv dalpha = (limits[ilim*2+1] - limits[ilim*2]) / fpv(nmm);
                for(int j=0; j<nmm; j++) for(int i=0; i<this->GI.GR; i++) {
                    fpv alpha = limits[ilim*2] + (j + this->GI.GN[i]) * dalpha;
                    fpv cosalpha = cos(alpha);
                    fpv sinalpha = sin(alpha);

                    // решаем квадратное уравнение для определения пределов интегрирования по r
                    fpv r_t = dist_observer_vertex / t;
                    fpv b2 = r_t*cosalpha; // b/2
                    fpv d4 = SQR(1.0 - r_t) - SQR(r_t*sinalpha);
                    if(d4 < 0) crash("s_Corner::MainFunc: internal error, negative determinant (22)");
                    fpv Rmax = b2 + sqrt(d4);
                    if(Rmax > 1) crash("s_Corner::MainFunc: internal error, Rmax>1");
                    fpv zll = sqrt(1.0 - Rmax);

                    // интегрируем по z от zll до 1
                    fpv val[4];
                    CalcIntOverZ_NEW(zll, 1.0, alpha, t, dist_observer_vertex, dist_src_vertex, hatPhi, val);
                    val[3] *= -PhiHatSign;
                    for(int k=0; k<4; k++) sum2[k] += val[k] * dalpha * this->GI.GC[i];
                }
            }

            for(int k=0; k<4; k++) sum2[k] *= 0.5 - (ilight ? gm : gp);

            uu[2+ilight]  += sum2[0]; // wave potential over t
            uu[6+ilight]  += sum2[1]; // time derivative of (wave potential over t)
            // converting derivatives over r and phi to the derivatives over source Cartesian coordinates
            uu[10+ilight] += cosphi * sum2[2] - sinphi * sum2[3];
            uu[14+ilight] += sinphi * sum2[2] + cosphi * sum2[3];
        }
#else // NEW VARIANT
        const fpv t_over_r_minus_1 = (dist_observer_vertex<t*1e-70) ? 1e70 : t / dist_observer_vertex - 1.0;
        for(int ilight=0; ilight<2; ilight++) {
            if(t_over_r_minus_1 < get_min_value<fpv>()) continue; // do nothing
            fpv sum2[4] = {0.,0.,0.,0.};

            const fpv PhiHatSign = (ilight ^ (J&1)) ? -1.0 : 1.0;
                  fpv hatPhi = _PiNumber - PhiHatSign * Phi; // In text: \hat{\Phi}
            const fpv coshatPhi = cos(hatPhi), sinhatPhi = sin(hatPhi);
            while(hatPhi<0.0)  hatPhi += _Pi2;
            while(hatPhi>_Pi2) hatPhi -= _Pi2;


            fpv psimin = 0.0, psimax = _PiNumber;
            // пытаемся ограничить область интегрирования
            if(BoundPhiZone(t - dist_observer_vertex, dist_src_vertex, L, hatPhi, psimin, psimax, 1)) continue;

            const int nmm = this->mm;
            const int nmm2 = this->mm;
            const fpv dpsi = (psimax - psimin) / fpv(nmm2);
            for(int j=0; j<nmm2; j++) for(int i=0; i<this->GI.GR; i++) {
                fpv psi = psimin + (j + this->GI.GN[i]) * dpsi;
                fpv cospsi = cos(psi);
                fpv sinpsi = sin(psi);
                fpv sinpsi2 = sin(0.5*psi);
                const fpv a = - dist_observer_vertex * cospsi;
                const fpv b = sqrt(t*t - SQR(dist_observer_vertex*sinpsi));

                fpv ximin = atan(cospsi / sqrt(t_over_r_minus_1*(t_over_r_minus_1+2.0)));
                fpv ximax = atan((t_over_r_minus_1 + cospsi) / (2.0*sinpsi2*sqrt(t_over_r_minus_1)));

                do { // we don't need precision here
                    if(b < get_min_value<fpv>()) break;
                    fpv cospp = cospsi * coshatPhi + sinpsi * sinhatPhi; // cos(psi - hat{Phi})
                    fpv sinpp = sinpsi * coshatPhi - cospsi * sinhatPhi; // sin(psi - hat{Phi})
                    //fpv rmin = dist_src_vertex-L;
                    //fpv rmax = dist_src_vertex+L;
                    fpv aux = L*L - SQR(dist_src_vertex*sinpp);
                    if(aux < 0.0) { ximax = -huge; break; } // no intersection
                    aux = sqrt(aux);
                    fpv rmin = dist_src_vertex*cospp - aux;
                    fpv rmax = dist_src_vertex*cospp + aux;
                    rmin = MAX(fpv(0.0), rmin);
                    rmax = MIN(t - dist_observer_vertex, rmax);
                    fpv sin_ximin_new = (rmin - a) / b;
                    fpv sin_ximax_new = (rmax - a) / b;
                    if(sin_ximin_new < -1.0) sin_ximin_new = -1.0;
                    if(sin_ximin_new > 1.0) break; // it should not occur
                    if(sin_ximax_new < -1.0)break; // it should not occur
                    if(sin_ximax_new > 1.0) sin_ximax_new = 1.0;
                    sin_ximin_new = asin(sin_ximin_new);
                    sin_ximax_new = asin(sin_ximax_new);
                    ximin = MAX(ximin, sin_ximin_new);
                    ximax = MIN(ximax, sin_ximax_new);
                } while(false);
                if(ximax < ximin) continue;

                fpv sum[4] = {0.,0.,0.,0.};
                fpv dxi = (ximax - ximin) / nmm;
                for(int jjj=0; jjj<nmm; jjj++) for(int iii=0; iii<this->GI.GR; iii++) {
                    fpv xi = ximin + (jjj+this->GI.GN[iii])* dxi;
                    fpv R0 = a + b * sin(xi);

                    fpv deltax = R0*cospsi - dist_src_vertex*coshatPhi;
                    fpv deltay = R0*sinpsi - dist_src_vertex*sinhatPhi;
                    fpv ll = SQR(deltax) + SQR(deltay);

                    fpv val[4], dF;
                    val[0] = tSpaceForm<fpv>::SpaceForm_rr(ll, &dF); // Амплитуда источника в заданной точке
                    val[1] = 2.*dF * (deltax * (R0*cospsi + dist_observer_vertex) + deltay * R0*sinpsi) / t;
                    val[2] = - 2.*dF * deltax;
                    val[3] = 2.*dF * deltay; // multiplication by sign is in calling subroutine
                    fpv mult = R0 * dxi * this->GI.GC[iii];
                    for(int k=0; k<4; k++)
                        sum[k] += val[k] * mult;
                }

                sum[3] *= -PhiHatSign;
                for(int k=0; k<4; k++) sum2[k] += sum[k] * dpsi * this->GI.GC[i];
            }

            for(int k=0; k<4; k++) sum2[k] *= 0.5 - (ilight ? gm : gp);
            for(int k=0; k<4; k++) sum2[k] /= (_Pi2*t);

            uu[2+ilight]  += sum2[0]; // wave potential over t
            uu[6+ilight]  += sum2[1]; // time derivative of (wave potential over t)
            // converting derivatives over r and phi to the derivatives over Cartesian coordinates
            uu[10+ilight] += cosphi * sum2[2] - sinphi * sum2[3];
            uu[14+ilight] += sinphi * sum2[2] + cosphi * sum2[3];
        }
#endif
    }

    if(n==1) { // Boundary integrals: half-line case only
        for(int ilight=0; ilight<2; ilight++) {
            const fpv PhiHatSign = (ilight ^ (J&1)) ? -1.0 : 1.0;
            const fpv hatPhi = _PiNumber - PhiHatSign * Phi; // In text: \hat{\Phi}

            int nmm_local = 5*this->mm; // no ideas why here the number of intervals should be greater than in other integral approximations
            fpv psimin = 0.0, psimax = _PiNumber;
            // пытаемся ограничить область интегрирования
            if(BoundPhiZone(t - dist_observer_vertex, dist_src_vertex, L, hatPhi, psimin, psimax, 0)) continue;

            //fpv a = SQR(t - dist_observer_vertex) + SQR(dist_src_vertex);
            fpv b = 2.0 * (t - dist_observer_vertex) * dist_src_vertex;
            fpv amb = SQR(t - dist_observer_vertex - dist_src_vertex); // a-b

            fpv bint = 0.0, sum = 0.0;
            fpv dpsi = (psimax - psimin) / fpv(nmm_local);
            for(int j=0; j<nmm_local; j++) for(int i=0; i<this->GI.GR; i++) {
                fpv psi = psimin + (j + this->GI.GN[i]) * dpsi;
                // fpv ll = a - b * cos(psi - hatPhi); -- lead to troubles with machine precision when a ~ b >> 1
                fpv ll = amb + 2. * b * SQR(sin(0.5*(psi - hatPhi)));
                fpv f = tSpaceForm<fpv>::SpaceForm_rr(ll);
                bint += f * sin(0.5*psi) * this->GI.GC[i] * dpsi;
                sum  += f * cos(0.5*psi) * this->GI.GC[i] * dpsi;
            }
            bint *= sqrt(t - dist_observer_vertex) / (_Pi2 * t);
            bint *= 0.5 - (ilight ? gm : gp);
            sum *=  PhiHatSign * sqrt(difr_zone_radius) / (_Pi2 * t);
            sum *=  0.5 - (ilight ? gm : gp);

            uu[6+ilight]  += sqrt(dist_observer_vertex) / t * bint; // time derivative of (wave potential over t)
            uu[10+ilight] += (-cosphi * bint - sinphi * sum) / MAX(sqrt(dist_observer_vertex), tiny_loc);
            uu[14+ilight] += (-sinphi * bint + cosphi * sum) / MAX(sqrt(dist_observer_vertex), tiny_loc);
        }
    }

fin:
    for(int i=0; i<16; i++) if(IsNaN(uu[i])) crash("s_Corner: NaN detected");
    if(_uu!=NULL) for(int i=0; i<16; i++) _uu[i] += uu[i];
    u[0] += uu[ 0] + uu[ 1] + uu[ 2] + uu[ 3];
    u[1] += uu[ 4] + uu[ 5] + uu[ 6] + uu[ 7];
    u[2] += uu[ 8] + uu[ 9] + uu[10] + uu[11];
    u[3] += uu[12] + uu[13] + uu[14] + uu[15];

#if 0 // debug print
    printf0("t = %22.12e; c = %22.12e %22.12e\n", t, coor[0], coor[1]);
    printf0("u = \n");
    printf0("%e %e %e %e\n", double(uu[0]), double(uu[1]), double(uu[2]), double(uu[3]));
    printf0("%e %e %e %e\n", double(uu[4]), double(uu[5]), double(uu[6]), double(uu[7]));
    printf0("%e %e %e %e\n", double(uu[8]), double(uu[9]), double(uu[10]), double(uu[11]));
    printf0("%e %e %e %e\n", double(uu[12]), double(uu[13]), double(uu[14]), double(uu[15]));
#endif
}

//======================================================================================================================
// Sommerfeld problem of the diffraction on a half-line
// Calculation of the wave potential
//======================================================================================================================
template<typename fpv>
void s_CornerWP<fpv>::PointValue(fpv t, const fpv* coor, fpv* V) const {
    // Если время слишком маленькое, то получаем решение разложением в ряд Тейлора вблизи нуля
    // Так как время маленькое, то наличие полупрямой можно не учитывать
    if(t < tSpaceForm<fpv>::Bterm * (sizeof(fpv)==8 ? 1e-6 : 1e-22)) {
        fpv u[5];
        PointValueTaylor(t, coor, u, (const tSpaceForm<fpv>&)(*this));
        V[0] = u[Var_W]; // PointValueTaylor возвращает волновой потенциал в позицию Var_W
        return;
    }

    fpv u[4];
    s_Corner<fpv>::MainFunc(t, coor, u, -1 /*two sides*/, 1 /*clear*/);
    V[0] = u[0] * t;
}
//======================================================================================================================


//======================================================================================================================
// Sommerfeld problem of the diffraction on a half-line
// Front-end
//======================================================================================================================
template<typename fpv>
void s_Corner<fpv>::PointValue(fpv T, const fpv* coor, fpv* uex) const {
    // If time is too small, then the solution is by the Taylor series at t=0
    // Since time is small, we don't take into account the presence of the half-line
    static const fpv timemin = pow_to_const_int<sizeof(fpv)/8>(1e-6);
    if(T < tSpaceForm<fpv>::Bterm * timemin) {
        PointValueTaylor(T, coor, uex, (tSpaceForm<fpv>&)(*this));
        uex[Var_W] = 0.0; // PointValueTaylor returns the wave potential at this position, removing
        return;
    }

    #if 1
        fpv summs[4];
        MainFunc(T, coor, summs); // -W/t, derivative of (-W/t) in t, derivatives of (-W/t) in x and in y

        uex[Var_R] = uex[Var_P] = summs[0] + summs[1] * T;
        uex[Var_U] = - summs[2] * T;
        uex[Var_V] = - summs[3] * T;
        uex[Var_W] = 0.0;
    #else
        // debug mode: numerical differentiation of the wave potential
        const fpv eps = sizeof(fpv)==8 ? 1e-6 : 1e-30;

        fpv summs[4], summ1, summ2;
        MainFunc(T, coor, summs); summ1 = summs[0] * T;
        MainFunc(T+eps, coor, summs); summ2 = summs[0] * (T+eps);
        uex[Var_R] = uex[Var_P] = (summ2 - summ1) / eps;

        fpv CU[2] = {coor[0]+eps, coor[1]};
        MainFunc(T, CU, summs); summ2 = summs[0] * T;
        uex[Var_U] = - (summ2 - summ1) / eps;

        fpv CV[2] = {coor[0], coor[1]+eps};
        MainFunc(T, CV, summs); summ2 = summs[0] * T;
        uex[Var_V] = - (summ2 - summ1) / eps;
        uex[Var_W] = 0.0;
    #endif
}
//======================================================================================================================

// Instantiation of template classes
template struct s_IVP2D<double>;
template struct s_Corner<double>;
template struct s_CornerPlanar<double>;
template struct s_IVP2DWP<double>;
template struct s_CornerWP<double>;
#ifdef EXTRAPRECISION_COLESO
template struct s_IVP2D<qd_real>;
template struct s_Corner<qd_real>;
template struct s_CornerPlanar<qd_real>;
template struct s_IVP2DWP<qd_real>;
template struct s_CornerWP<qd_real>;
template struct s_IVP2D<dd_real>;
template struct s_Corner<dd_real>;
template struct s_CornerPlanar<dd_real>;
template struct s_IVP2DWP<dd_real>;
template struct s_CornerWP<dd_real>;
#endif

// Verification subroutine
#ifdef EXTRAPRECISION_COLESO
#define real_type1 double
#define real_type2 dd_real
#define mult 1.05 // multiplication step
#define NMAX 1000
#define threshold 1.5e-15
void CheckCornerPlanar() {
    s_CornerPlanar<real_type1> S1;
    s_CornerPlanar<real_type2> S2;
    S1.Init();
    S2.Init();
    const tGaussIntegrator<real_type1>& I1 = S1.get_GJI();
    const tGaussIntegrator<real_type2>& I2 = S2.get_GJI();

    if(1) {
        real_type2 absdelta;
        for(int sign=0; sign<2; sign++) for(int n=-NMAX; n<=NMAX; n++) {
            real_type2 beta = pow(mult, n);
            if(sign) beta *= -1.0;
            real_type1 E1 = calc_E<real_type1>(real_type1(beta), I1);
            real_type2 E2 = calc_E<real_type2>(beta, I2);
            double err = fabs(double(E2-E1));
            if(err>threshold)
                printf("% e  % 25.15e % 25.15e\n", double(beta), double(E2), err);
        }
    }
    if(1) {
        real_type2 alpha, absbeta;
        for(int sign=0; sign<2; sign++) for(int n=-NMAX; n<=NMAX; n++) {
            real_type2 beta = pow(mult, n);
            if(sign) beta *= -1.0;
            for(int na=-NMAX; na<=NMAX; na++) {
                real_type2 alpha = pow(mult, na);
                real_type1 E1 = calc_J0<real_type1>(real_type1(alpha), real_type1(beta), I1);
                real_type2 E2 = calc_J0<real_type2>(alpha, beta, I2);
                double err = fabs(double(E2-E1));
                if(err>threshold)
                    printf("% e % e  % 25.15e % 25.15e\n", double(alpha), double(beta), double(E2), err);
            }
        }
    }
}
#endif
