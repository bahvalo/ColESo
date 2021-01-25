// *********************************************************************************************************************
// *****                                                                                                           *****
// *****                                  Collection of exact solutions (ColESo)                                   *****
// *****                             Point source in a flow (moving point source) - 3D                             *****
// *****                                  Source of a finite size in 1D, 2D, 3D                                    *****
// *****                                               Rotating dipole                                             *****
// *****                                                                                                           *****
// *********************************************************************************************************************

#include "coleso.h"
#include "parser.h"
#include "geom_primitive.h"
#ifdef _NOISETTE
#include "lib_base.h"
#endif
#ifdef EXTRAPRECISION_COLESO
#include "extraprecision.h" 
#endif
#include "es_specfunc.h"
#include <cstdlib>

//======================================================================================================================
// Считывание параметров
//======================================================================================================================
template<typename fpv>
void s_PointSource<fpv>::ReadParams(tFileBuffer& FB){
    tSpaceForm<fpv>::Read(FB, IO_DONTCRASH);
    tParamManager PM; 
    PM.Request(Ampl, "Ampl");
    PM.Request(Freq, "Freq");
    PM.Request(Phase, "Phase");
    PM.Request(tmin, "tmin");
    PM.Request(tmax, "tmax");
    PM.Request(SignalType, "SignalType");

    PM.Request(SoundVel, "SoundVel");
    PM.Request(FlowMach, "FlowMach");
    PM.Request(UnstAmpl, "UnstAmpl");
    PM.Request(UnstFreq, "UnstFreq");
    PM.Request(UnstPhase, "UnstPhase");
    PM.ReadParamsFromBuffer(FB);

    if(PM["UnstAmpl"].GetIsSet() || PM["UnstFreq"].GetIsSet() || PM["UnstPhase"].GetIsSet()) UseUnsteady = 1;
}

//======================================================================================================================
// Проверка параметров
//======================================================================================================================
template<typename fpv>
void s_PointSource<fpv>::Init() {
    if(SignalType!=1 && SignalType!=6) crash("s_PointSource::Init error: wrong SignalType");
    // tSpaceForm<double>::Init(); -- не вызываем! Наоборот, нам нужна L1-норма от формы импульса
    if(!this->NormalizeForm) { this->Aterm *= tSpaceForm<fpv>::AmpliduteLinf2L1(this->Bterm); this->NormalizeForm=1; }
}

//======================================================================================================================
// Определение времени запаздывания и радиуса эмиссии
// Вход: время и координаты наблюдателя относительно источника
//======================================================================================================================
template<typename fpv>
fpv s_PointSource<fpv>::FindSourceTime(fpv time, fpv r0x, fpv r0y, fpv* re) const {
    // Для покоящегося источника в равномерном потоке
    if(!UseUnsteady) {
        fpv X2 = r0x*r0x+r0y*r0y;
        fpv SQ = sqrt(SQR(r0x*FlowMach) + (1.-FlowMach*FlowMach)*X2);
        fpv tau = (r0x*FlowMach - SQ) / (1.-FlowMach*FlowMach); // < 0
        re[0] = r0x + tau*FlowMach;
        re[1] = r0y;
        return time + tau;
    }

    // Для источника, двигающегося с переменной скоростью
    // Начальное приближение
    const s_PointSource* S = this;
    fpv Ufreq = S->UnstFreq / SoundVel;
    fpv Uampl = S->UnstAmpl / SoundVel;
    fpv tau = sqrt(SQR(r0x) + SQR(r0y)) / (FlowMach - Uampl); // > 0

    fpv dtauold = huge;
    fpv L = huge, Lold = huge;
    for(int iter=0; iter<1000; iter++) {
        fpv M = FlowMach + Uampl*sin(Pi2*Ufreq*(time-tau) + S->UnstPhase);
        fpv F = FlowMach*tau + Uampl/(Pi2*Ufreq)*(cos(Pi2*Ufreq*(time-tau) + S->UnstPhase) - cos(Pi2*Ufreq*time + S->UnstPhase));

        re[0] = r0x - F;
        re[1] = r0y;
        fpv re2 = SQR(re[0]) + SQR(re[1]);
        Lold = L;
        L = tau*tau - re2;
        // printf("L = %e\n", L);
        if(fabs(L) < 1e-10) break; // TODO check it!
        fpv dLdTau = 2.0*(tau + re[0]*M);
        fpv dtau = - L / dLdTau;
        if(fabs(L) > fabs(Lold)) crash("FindSourceTime error: Newton process unconverges! q = %f", double(dtau/dtauold));
        dtauold = dtau;
        tau += (iter<50 ? 0.2 : 1.0) * dtau;
        if(iter==200) crash("FindSourceTime error: Newton is under stagnation, dtau = %e", double(fabs(dtau)));
    }
    if(tau < 0) crash("FindSourceTime error: negative Tau found: %f", double(tau));
    return time-tau;
}
//======================================================================================================================


//======================================================================================================================
// Решение волнового уравнения для задачи с движущимся по закону x(t) точечным источником
// Источник в волновом уравнении: функция SourceFunction(t)*Dirac(r - x(t))
// Скорость звука равна 1
// Требуется найти значение переменной в заданный момент времени в заданной точке
// Вход:  time - момент наблюдения
//        R    - двумерный вектор (z, sqrt(x^2+y^2)) из источника к наблюдателю в момент наблюдения
// Input: diff - 0 for function, 1 for X derivative, 2 for radial deriv. in cylindrical coords,
//             - 3 for time deriv. at steady system (with constant R),
//             - 4 for time deriv. at flux system (with constant R-Mt)
//======================================================================================================================
template<typename fpv>
fpv s_PointSource<fpv>::monopole(fpv* R, fpv time, int diff) const {
    const fpv _Pi2 = GetPiNumber2<fpv>();
    fpv freq  = Freq / SoundVel;
 
    // Вычислим момент эмиссии и положение источника относительно наблюдателя в этот момент
    fpv r[2];
    fpv t_star = FindSourceTime(time, R[0], R[1], r);
    fpv absr = sqrt(r[0]*r[0]+r[1]*r[1]);
    {
        // На всякий случай проверяем, что уравнение решено правильно
        fpv TT = time - absr;
        if(fabs(TT-t_star) > 1e-8*(fabs(TT)+fabs(t_star)+1.0))
            crash("s_PointSource::monopole: internal error: TT = %25.15e, t=%25.15e", double(TT), double(t_star));
    }

    // Вычислим функцию источника и её производную в момент эмиссии t* = time - |r*|
    fpv sf = 0.0, dsf = 0.0;

    if(t_star >= tmin && t_star <= tmax) {
        switch(SignalType) {;
        case 1:
            sf = sin(_Pi2*freq*t_star+Phase);
            dsf = _Pi2*freq* cos(_Pi2*freq*t_star+Phase);
            break;
        case 6:
            sf = sin(_Pi2*freq*t_star+Phase);
            sf = sf*sf*sf*sf;
            dsf = _Pi2*freq*(sin(2.*(_Pi2*freq*t_star+Phase)) - 0.5*sin(4.0*(_Pi2*freq*t_star+Phase)));
            break;
        default:
            crash("s_PointSource::monopole: SignalType %i is undefined", SignalType);
        }
    }

    // Считаем M = - dx/dt в момент эмиссии. Ещё считаем dM/dt 
    fpv M = FlowMach, dMdt = 0.0;
    if(UseUnsteady) {
        fpv Uampl = UnstAmpl / SoundVel;
        M += fpv(Uampl*sin(Pi2*freq*t_star + UnstPhase));
        dMdt = fpv(Uampl*Pi2*freq*cos(Pi2*freq*t_star+ UnstPhase));
    }

    fpv rxi = absr + r[0]*M; // знаменатель в формуле
    fpv xi = rxi / absr; // допплеровский множитель

    fpv solut = 0.0;
    switch (diff) {;
    case 0: // волновой потенциал
        solut = sf / rxi;
        break;
    case 1: // Z derivative in cylindrical coords
        solut = - dsf * r[0] / (rxi*rxi) -
            sf / (rxi*rxi) * (M + (1.-M*M) * r[0] / rxi) + sf*r[0]*r[0]/(rxi*rxi*rxi)*dMdt;
        break;
    case 2: // radial deriv. in cylindrical coords
        solut = - dsf * r[1] / (rxi*rxi) - 
            sf / (rxi*rxi) * (    (1.-M*M) * r[1] / rxi) + sf*r[1]*r[0]/(rxi*rxi*rxi)*dMdt;
        break;
    case 3: // time deriv. at steady system (with constant R)
        solut = dsf / rxi;
        break;
    case 4: // for time deriv. at flow system (with constant R-Mt)
        solut = dsf/(xi*rxi) - sf/(rxi*rxi*rxi)*(M*r[0] + M*M*absr) - sf/(xi*rxi*rxi)*dMdt*r[0];
        break;
    default:
        crash("s_PointSource::monopole: error DIFF code");
    }

    if(IsNaN(solut)) {
        pprintf("R = %e %e\n", double(R[0]), double(R[1]));
        pprintf("r  = %e %e\n",  double(r[0]),  double(r[1]));
        pprintf("M  = %e, dMdt = %e\n", double(M), double(dMdt));
        pprintf("time = %e\n", double(time));
        crash("s_PointSource::monopole: NaN detected! See log for details");
    }
    return solut * Ampl * tSpaceForm<fpv>::Aterm / (4.0*GetPiNumber<fpv>());
}
//======================================================================================================================


//======================================================================================================================
// Точное решение для задачи о распространении звука от точечного источника
// Обезразмеривание произвольное
//======================================================================================================================
template<typename fpv>
void s_PointSource<fpv>::PointValue1(fpv t, const fpv* coor, fpv* uex) const {
    for(int ivar=Var_R; ivar<=Var_P; ivar++) uex[ivar] = 0.0;

    // Переобезразмериваем время на скорость звука
    fpv time = t * SoundVel;

    fpv r[2] = {coor[0], sqrt(SQR(coor[1]) + SQR(coor[2]))};
    if(fabs(r[0]) < tiny && fabs(r[1]) < tiny) return; // В точке источника решение не определено. Возвращаем ноль

    // Считаем пульсации плотности и цилиндрических компонент скоростей
    uex[Var_R] = uex[Var_P] = monopole(r, time, 4); // Rho'
    uex[Var_U] = - monopole(r, time, 1); // Uz'
    uex[Var_V] = - monopole(r, time, 2); // Ur'

    // Переходим от цилиндрических компонент к декартовым
    fpv solut = uex[Var_V];
    if(r[1] < get_eps<fpv>())
        uex[Var_V] = uex[Var_W] = 0.0;
    else {
        uex[Var_V] = solut * coor[1] / r[1];
        uex[Var_W] = solut * coor[2] / r[1];
    }

    // Учитываем исходное обезразмеривание
    uex[Var_R] /= SoundVel;
    uex[Var_P] *= SoundVel;
}
//======================================================================================================================

//======================================================================================================================
template<typename fpv>
void s_PointSource<fpv>::PointValue(fpv t, const fpv* Coords, fpv* uex) const {
//======================================================================================================================
    fpv buf[Var_NumMax];
    for(int ivar=0; ivar<Var_N; ivar++) uex[ivar] = 0.0;

    for(int i=0; i<tSpaceForm<fpv>::NAngularPeriods; i++) {
        fpv r[3] = {tSpaceForm<fpv>::r0[0], tSpaceForm<fpv>::r0[1], tSpaceForm<fpv>::r0[2]}; // положение источника
        if(i) RotateVector2D<fpv>(r, i*GetPiNumber2<fpv>()/tSpaceForm<fpv>::NAngularPeriods, r); // поворачиваем положение источника, если периодика по углу
        for(int icoor=0; icoor<3; icoor++) r[icoor] = Coords[icoor] - r[icoor]; // вектор от источника к приёмнику
        PointValue1(t, r, buf);
        for(int ivar=0; ivar<Var_N; ivar++) uex[ivar] += buf[ivar];
    }
}
//======================================================================================================================



// *********************************************************************************************************************
// *****                                                                                                           *****
// *****                            Steady source of finite width (1D, 2D, 3D)                                     *****
// *****                                                                                                           *****
// *********************************************************************************************************************

//======================================================================================================================
template<typename fpv>
void s_Source1D3D<fpv>::ReadParams(tFileBuffer& FB) {
    tParamManager PM;
    // Зависимость источника от координаты
    tSpaceForm<fpv>::Read(FB, true /* periodics is permitted */);
    // Зависимость источника от времени
    PM.Request(Ampl, "Ampl");
    PM.Request(Freq, "Freq");
    PM.Request(Phase, "Phase");
    PM.Request(tmin, "tmin");
    PM.Request(tmax, "tmax");
    PM.Request(SignalType, "SignalType");
    PM.ReadParamsFromBuffer(FB);
}
//======================================================================================================================


//======================================================================================================================
// Зависимость сигнала от времени
//======================================================================================================================
template<typename fpv>
fpv s_Source1D3D<fpv>::TimeForm(fpv time) const {
    if(time < tmin || time > tmax) return 0.0;
    
    switch (SignalType){
    case 0:
        return Ampl;
    case 1:
        return Ampl*sin(GetPiNumber2<fpv>()*Freq*time+Phase);
    case 6: {
        fpv s = sin(GetPiNumber2<fpv>()*Freq*time+Phase);
        return Ampl*s*s*s*s;
    }
    default: crash("s_Source1D3D::TimeForm: SignalType %i is undefined", SignalType);
    }
}
//======================================================================================================================

//======================================================================================================================
// Зависимость сигнала от времени. Вычисление производной
//======================================================================================================================
template<typename fpv>
fpv s_Source1D3D<fpv>::TimeFormDeriv(fpv time) const {
    if(time < tmin || time > tmax) return 0.0;
    fpv omega = GetPiNumber2<fpv>()*Freq;
    fpv phase = omega*time+Phase;
    
    switch (SignalType){
    case 0:
        return 0.0;
    case 1:
        return Ampl*cos(phase)*omega;
    case 6: {
        fpv s = sin(phase);
        return 4.*Ampl*s*s*s*cos(phase)*omega;
    }
    default: crash("s_Source1D3D::TimeFormDeriv: SignalType %i is undefined", SignalType);
    }
}
//======================================================================================================================

//======================================================================================================================
// Зависимость сигнала от времени. Вычисление второй производной
//======================================================================================================================
template<typename fpv>
fpv s_Source1D3D<fpv>::TimeFormDeriv2(fpv time) const {
    if(time < tmin || time > tmax) return 0.0;
    fpv omega = GetPiNumber2<fpv>()*Freq;
    fpv phase = omega*time+Phase;
    
    switch (SignalType){
    case 0:
        return 0.0;
    case 1:
        return -Ampl*sin(phase)*omega*omega;
    case 6: // (sin(x))^4 = cos(4x)/8 + cos(2x)/2 + 3/8
        return 2.*Ampl*(cos(4.*phase) - cos(2.*phase))*omega*omega;
    default: crash("s_Source1D3D::TimeFormDeriv: SignalType %i is undefined", SignalType);
    }
}
//======================================================================================================================

//======================================================================================================================
// Зависимость сигнала от времени. Вычисление первообразной
//======================================================================================================================
template<typename fpv>
fpv s_Source1D3D<fpv>::TimeFormAntideriv(fpv time) const {
    if(time < tmin) return 0.0;
    if(time > tmax) time = tmax; // первообразная стагнирует, когда функция нулевая
    
    fpv omega = GetPiNumber2<fpv>()*Freq;
    fpv phase = omega*time+Phase;
    switch (SignalType){
    case 0:
        return Ampl*(time - tmin);
    case 1:
        return Ampl*(cos(Phase) - cos(phase))/omega;
    case 6: {
        return Ampl * (0.375*time - (sin(2*phase)-sin(2*Phase)) / (4.*omega) + (sin(4*phase)-sin(4*Phase)) / (32.*omega));
    }
    default: crash("s_Source1D3D::TimeFormAntideriv: SignalType %i is undefined", SignalType);
    }
}
//======================================================================================================================


//======================================================================================================================
// Сигнал для одного источника на расстоянии x от него
//======================================================================================================================
template<typename fpv>
void s_Source1D<fpv>::PointValue1(fpv t, fpv x, fpv* V) const {
    const tGaussIntegrator<fpv>& _GI = s_Source1D3D<fpv>::GI;
    if(_GI.GN==NULL) crash("s_Source1D: Init not done");
    static fpv cnst_radius = sqrt(-2.*log(get_eps<fpv>())/GetLn2<fpv>()); // sqrt(-2*log2(eps)) ~ 10.21
    fpv beff = tSpaceForm<fpv>::Bterm;
    if(tSpaceForm<fpv>::Form==0) beff *= cnst_radius; // Радиус источника: для гауссиана выбирается радиус обрезки
    const fpv q = MAX(tSpaceForm<fpv>::InvBTerm, 4.0*Pi2*s_Source1D3D<fpv>::Freq);

    fpv sum[2] = {0.0, 0.0};
    for(int iii=0; iii<2; iii++) { // iii=0: intergation over y<x; iii=1: integration over y>x
        fpv yL = -beff, yR = beff;
        if(iii==0) {
            yR = MIN(yR, x);
            yL = MAX(yL, (s_Source1D3D<fpv>::tmin-(t-x)));
            yR = MIN(yR, (s_Source1D3D<fpv>::tmax-(t-x)));
        }
        else {
            yL = MAX(yL, x);
            yL = MAX(yL, (t+x-s_Source1D3D<fpv>::tmax));
            yR = MIN(yR, (t+x-s_Source1D3D<fpv>::tmin));
        }
        if(yR <= yL) continue;

        fpv dy0 = yR-yL; // длина отрезка интегрирования
        int numk = int(0.68*q*dy0) + 1;
        fpv dy = dy0 / numk;
        for(int k=0; k<numk; k++) for(int i=0; i<_GI.GR; i++) {
            fpv y = yL + (k + _GI.GN[i]) * dy;
            fpv f = tSpaceForm<fpv>::SpaceForm(y);
            fpv g = this->TimeForm(t - fabs(y - x));
            sum[iii] += f * g * _GI.GC[i];
        }
        sum[iii] *= 0.5 * dy;
    }

    V[Var_R] = V[Var_P] = sum[0] + sum[1];
    V[Var_U] = sum[0] - sum[1];
    V[Var_V] = V[Var_W] = 0.0;
}
//======================================================================================================================

template<typename fpv>
static int _compare_fpv(const void* a, const void* b){
    if( *((fpv*) a) <  *((fpv*) b) ) return -1;
    if( *((fpv*) a) == *((fpv*) b) ) return 0;
    return 1;
}

//======================================================================================================================
// Сигнал для одного источника: вычисляем p и ur/r
//======================================================================================================================
template<typename fpv>
void s_Source3D<fpv>::PointValue1(fpv t, fpv r, fpv& _p, fpv& _urr) const {
    static fpv cnst_radius = sqrt(-2.*log(get_eps<fpv>())/GetLn2<fpv>()); // sqrt(-2*log2(eps)) ~ 10.21
    const tGaussIntegrator<fpv>& _GI = s_Source1D3D<fpv>::GI;
    if(_GI.GN==NULL) crash("s_Source3D: Init not done");
    fpv beff = tSpaceForm<fpv>::Bterm;
    if(tSpaceForm<fpv>::Form==0) beff *= cnst_radius; // Радиус источника: для гауссиана выбирается радиус обрезки

    const fpv q = MAX(tSpaceForm<fpv>::InvBTerm, 4.0*Pi2*s_Source1D3D<fpv>::Freq);
    // Пока в _urr вычисляем радиальную компоненту скорости, не делённую на r. Поделим потом.
    // Это не приводит к потере точности, потому что для вычисления декартовых компонент идёт обратно умножение на r
    _p = _urr = 0.0; 

    const fpv invr = 1.0/MAX(r, fpv(1e-50));
    // Вычисляем слагаемое, не содержащее Delta(...)
    {
        // Определяем все возможные точки разрывов или излома подынтегральных функций
        // Левый предел = 0, правый = MIN(beff, 2*r); промежуточные точки \pm r \pm (t-tmin), \pm r \pm (t-tmax)
        fpv ximax = MIN(beff*invr, fpv(2.0));
        fpv X[11] = {0, 1.0, ximax, 
            -1.0-invr*(t-s_Source1D3D<fpv>::tmax), -1.0-invr*(t-s_Source1D3D<fpv>::tmin), 
            -1.0+invr*(t-s_Source1D3D<fpv>::tmax), -1.0+invr*(t-s_Source1D3D<fpv>::tmin), 
             1.0-invr*(t-s_Source1D3D<fpv>::tmax),  1.0-invr*(t-s_Source1D3D<fpv>::tmin),  
             1.0+invr*(t-s_Source1D3D<fpv>::tmax),  1.0+invr*(t-s_Source1D3D<fpv>::tmin)};
        std::qsort(X, 11, sizeof(fpv), _compare_fpv<fpv>);

        for(int iinterval=0; iinterval<10; iinterval++) {
            const fpv xiL = X[iinterval];
            const fpv xiR = X[iinterval+1];
            if(xiL < 0.0 || xiR > ximax) continue;
            if(xiR <= xiL) continue;

            fpv p = 0.0, ur = 0.0;
            const fpv dy0 = xiR - xiL;
            int numk = int(0.68*q*r*dy0) + 1;
            fpv dy = dy0 / numk;
            for(int k=0; k<numk; k++) for(int i=0; i<_GI.GR; i++) {
                fpv xi = xiL + (k + _GI.GN[i]) * dy;
                fpv x = xi * r;
                fpv f   = tSpaceForm<fpv>::SpaceForm(x);
                fpv df  = tSpaceForm<fpv>::dSpaceForm(x);
                fpv g   = this->TimeForm(t - fabs(r - x));
                fpv G   = this->TimeFormAntideriv(t - fabs(r - x));
                p  += 0.5 * r * xi * f * g * _GI.GC[i];
                ur += 0.5 * (- xi*f*G - x*xi*df*G + x*fabs(1.0-xi)*f*g) * _GI.GC[i];
            }
            _p   += p  * dy;
            _urr += ur * dy;
        }
    }
    // Одиночное слагаемое
    _urr += tSpaceForm<fpv>::SpaceForm(2.0*r) * this->TimeFormAntideriv(t - r);

    // Вычисляем слагаемое, содержащее Delta(...) или его производную по r
    {
        // Определяем все возможные точки разрывов или излома подынтегральных функций
        const fpv xr = beff+r;
        // Левый предел = r, правый = xr; промежуточные точки \pm (t-tmin), \pm (t-tmax)
        fpv X[4] = {r, xr, fabs(t-s_Source1D3D<fpv>::tmin), fabs(t-s_Source1D3D<fpv>::tmax)};
        std::qsort(X, 4, sizeof(fpv), _compare_fpv<fpv>);

        for(int iinterval=0; iinterval<3; iinterval++) {
            const fpv xL = X[iinterval];
            const fpv xR = X[iinterval+1];
            if(xL < r || xR > xr) continue;
            if(xR <= xL) continue;

            fpv p = 0.0, ur = 0.0;
            const fpv dy0 = xR - xL;
            int numk = int(0.68*q*dy0) + 1;
            fpv dy = dy0 / numk;
            for(int k=0; k<numk; k++) for(int i=0; i<_GI.GR; i++) {
                fpv x = xL + (k + _GI.GN[i]) * dy;
                fpv g   = this->TimeForm(t - x);
                fpv G   = this->TimeFormAntideriv(t - x);
                fpv Delta = 0.0, dDelta = 0.0;
                if(r > 0.1*tSpaceForm<fpv>::Bterm) { // Считаем обычным образом
                    Delta = ((x+r)*tSpaceForm<fpv>::SpaceForm(x+r) - (x-r)*tSpaceForm<fpv>::SpaceForm(x-r)) / (2.0 * r);
                    dDelta = - Delta / r + ((x+r)*tSpaceForm<fpv>::dSpaceForm(x+r) + tSpaceForm<fpv>::SpaceForm(x+r) + 
                                            (x-r)*tSpaceForm<fpv>::dSpaceForm(x-r) + tSpaceForm<fpv>::SpaceForm(x-r)) / (2.0 * r);
                }
                else {
                    // Нужно корректно определить предены интегрирования, чтобы не вылететь за beff
                    fpv ximax = (beff - x) * invr;
                    if(ximax > 1.0) ximax = 1.0;
                    if(ximax < -1.0) continue; // Delta = dDelta = 0.0;
                    for(int j=0; j<_GI.GR; j++) {
                        fpv xi = _GI.GN[j]*(ximax+1.0)-1.0;
                        fpv y = x + r * xi;
                        fpv f   = tSpaceForm<fpv>::SpaceForm(y);
                        fpv df  = tSpaceForm<fpv>::dSpaceForm(y);
                        fpv ddf = tSpaceForm<fpv>::ddSpaceForm(y);
                        Delta += (y*df+f) * _GI.GC[j];
                        dDelta += (y*ddf+2*df) * xi * _GI.GC[j];
                    }
                    Delta *= 0.5*(ximax+1.0);
                    dDelta *= 0.5*(ximax+1.0);
                }
                p  += Delta * g * _GI.GC[i];
                ur -= dDelta * G * _GI.GC[i];
            }
            _p   += p  * dy;
            _urr += ur * dy;
        }
    }
    
    _urr *= invr;
}
//======================================================================================================================


//======================================================================================================================
// Сигнал от решётки источников
//======================================================================================================================
template<typename fpv>
void s_Source1D<fpv>::PointValue(fpv t, const fpv* coor, fpv* uex) const {
    fpv buf[Var_NumMax];
    for(int ivar=0; ivar<Var_N; ivar++) uex[ivar] = 0.0;

    for(int i=-s_Source1D3D<fpv>::MaxPer[0]; i<=s_Source1D3D<fpv>::MaxPer[0]; i++) {
        fpv xs = (i==0) ? s_Source1D3D<fpv>::r0[0] : s_Source1D3D<fpv>::r0[0] + s_Source1D3D<fpv>::PerX * i;
        PointValue1(t, coor[0] - xs, buf);
        for(int ivar=0; ivar<Var_N; ivar++) uex[ivar] += buf[ivar];
    }
}
//======================================================================================================================


//======================================================================================================================
// Сигнал от решётки источников
//======================================================================================================================
template<typename fpv>
void s_Source3D<fpv>::PointValue(fpv t, const fpv* coor, fpv* uex) const {
    fpv rs[3] = {tSpaceForm<fpv>::r0[0], tSpaceForm<fpv>::r0[1], tSpaceForm<fpv>::r0[2]};

    for(int ivar=0; ivar<Var_N; ivar++) uex[ivar] = 0.0;
    for(int i=0; i<tSpaceForm<fpv>::NAngularPeriods; i++) {
        if(i) RotateVector2D(rs, GetPiNumber2<fpv>()/tSpaceForm<fpv>::NAngularPeriods, rs);
        fpv rr[3] = {coor[0]-rs[0], coor[1]-rs[1], coor[2]-rs[2]};
        fpv r = sqrt(VDOT(rr,rr));
        fpv p, urr;
        PointValue1(t, r, p, urr);
        uex[Var_R] += p;
        uex[Var_U] += urr * rr[0];
        uex[Var_V] += urr * rr[1];
        uex[Var_W] += urr * rr[2];
        uex[Var_P] += p;
    }
}
//======================================================================================================================



//======================================================================================================================
void s_Source2D::ReadParams(tFileBuffer& FB) {
    tParamManager PM; 
    tSpaceForm<double>::Read(FB, false /* no periodics allowed */);
    PM.Request(Ampl, "Ampl");
    PM.Request(Freq, "Freq");
    PM.Request(Phase, "Phase");
    PM.Request(GR, "GR"); // параметры интегрирования
    PM.Request(mm, "mm"); // параметры интегрирования
    PM.Request(Hmax, "Hmax");
    PM.ReadParamsFromBuffer(FB);
}
//======================================================================================================================

//======================================================================================================================
void s_Source2D::PointValue(double t, const double* coor, double* V) const {
    const double beff = Bterm * (Form==FORM_GAUSSIAN ? Hmax : 1.0); // Радиус источника: для гауссиана выбирается радиус обрезки
    const double omega = Pi2 * Freq;
    const double r = sqrt(SQR(coor[0]-r0[0]) + SQR(coor[1]-r0[1]));

    const double x = omega * r;
    const double J0x = BesselJ0(x);
    const double N0x = BesselN0(x);
    const double J1x = BesselJ1(x);
    const double N1x = BesselN1(x);

    double re = 0.0, im = 0.0; // действительная и мнимая часть интеграла
    double red = 0.0, imd = 0.0; // действительная и мнимая часть интеграла от продифференцированной функции
    // цикл по двум интервалам: от 0 до r и от r до бесконечности
    for(int iint=0; iint<2; iint++) {
        double limL = iint ? r : 0.0, limR = iint ? huge : r;
        if(limL > beff) continue;
        if(limR > beff) limR = beff;

        const double dr00 = (limR - limL) / mm;
        if(dr00 < tiny) continue;
        for(int ii=0; ii<mm; ii++) for(int jj=0; jj<GR; jj++) {
            const double r00 = limL + (ii + GN[jj]) * dr00;
            const double f = tSpaceForm<double>::SpaceForm(r00);
            const double x0 = omega * r00;
            const double J0x0 = BesselJ0(x0);
            const double N0x0 = BesselN0(x0);
            const double J1x0 = BesselJ1(x0);
            const double N1x0 = BesselN1(x0);
            double denom = J1x0*N0x0 - J0x0*N1x0;
            if(fabs(denom)<tiny || IsNaN(denom)) crash("s_Source2D::PointValue: Wronskian for Bessel equation failed");
            const double mult = f * GC[jj] * dr00 / denom;
            re += mult * J0x0*J0x;
            im -= mult * (iint ? N0x0*J0x : N0x*J0x0);
            red -= mult * J0x0*J1x;
            imd += mult * (iint ? N0x0*J1x : N1x*J0x0);
        }
    }

    const double phase = omega*t+Phase;
    const double c = cos(phase), s = sin(phase);
    double P = re * s + im * c;
    double U = red * c - imd * s;
    
    V[Var_R] = V[Var_P] = P;
    if (r < tiny) V[Var_U] = V[Var_V] = 0;
    else {
        V[Var_U] = U * (coor[0] - r0[0]) / r;
        V[Var_V] = U * (coor[1] - r0[1]) / r;
    }
    V[Var_W] = 0;
}
//======================================================================================================================



//======================================================================================================================
void s_RotatingDipole::PointValue(double t, const double* Coords, double* uex) const {
//======================================================================================================================
    for(int i=0; i<5; i++) uex[i]=0.0;

    double x = Coords[(dir+1)%3], y = Coords[(dir+2)%3], z = Coords[dir];
    double R = sqrt(x*x + y*y + z*z);
    if(R<tiny) return;

    double psi = Omega*(t - R) - Phase;
    double c = cos(psi);
    double s = sin(psi);
    double a = -s*x+c*y;
    double b = s*y+c*x;

    // Плотность и давление
    uex[Var_R] = uex[Var_P] = Omega*(a - b*Omega*R)/(R*R*R);

    // Скорости
    double aux = ((3.0 - Omega*Omega*R*R)*b + 3.0*Omega*R*a)/(R*R*R*R*R);
    uex[Var_U + (dir+1)%3] = x*aux - (c - Omega*R*s)/(R*R*R);
    uex[Var_U + (dir+2)%3] = y*aux - (s + Omega*R*c)/(R*R*R);
    uex[Var_U + dir]       = z*aux;
}
//======================================================================================================================

//======================================================================================================================
void s_RotatingDipole::ReadParams(tFileBuffer& FB) {
//======================================================================================================================
    tParamManager PM;
    PM.RequestParameter(Ampl, "Ampl", PARTYPE_DOUBLE, IO_DONTCRASH);
    PM.RequestParameter(Omega, "Omega", PARTYPE_DOUBLE, IO_DONTCRASH);
    PM.RequestParameter(Phase, "Phase", PARTYPE_DOUBLE, IO_DONTCRASH);
    PM.RequestParameter(dir, "dir", PARTYPE_INT, IO_DONTCRASH, GE, 0, LE, 2);
    PM.ReadParamsFromBuffer(FB);
}
//======================================================================================================================



#define INSTANTIATE(T) \
    template struct s_Source1D3D<T>; \
    template struct s_Source1D<T>; \
    template struct s_Source3D<T>; \
    template struct s_PointSource<T>;

INSTANTIATE(NativeDouble)
#ifdef EXTRAPRECISION_COLESO
INSTANTIATE(dd_real)
INSTANTIATE(qd_real)
#endif
