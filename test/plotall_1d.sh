#!/bin/bash

cd ./DATA1D
../plot1d.sh "plot '4peak.dat' using 1:2 w l ls 100 ti 'rho_{exact}'" "x" "rho_{exact}" 4peak.png
../plot1d.sh "plot 'acousticshock.dat' using 1:2 w l ls 101 ti 't=0', 'acousticshock.dat' using 1:3 w l ls 104 ti 't=35', 'acousticshock.dat' using 1:4 w l ls 102 ti 't=105'" "x" "rho\'" acousticshock.png
../plot1d.sh "plot 'planargauss.dat' using 1:2 w l ls 100 ti 'rho''" "x" "rho\'" planargauss.png
../plot1d.sh "plot 'planarsinus_1.dat' using 1:2 w l ls 100 ti 'rho''" "x" "rho\'" planarsinus_1.png
../plot1d.sh "plot 'riemann.dat' using 1:2 w l ls 100 ti 'rho''" "x" "rho\'" riemann.png
../plot1d.sh "plot 'simplewave.dat' using 1:2 w l ls 100 ti 'rho(t=0)', 'simplewave.dat' using 1:3 w l ls 102 ti 'p(t=0)', 'simplewave.dat' using 1:4 w l ls 104 ti 'rho(t=0.1)', 'simplewave.dat' using 1:5 w l ls 101 ti 'p(t=0.1)'" "x" "rho, p" simplewave.png
../plot1d.sh "plot 'source1d_1.dat' using 1:2 w l ls 104 ti 'rho''', 'source1d_1.dat' using 1:3 w l ls 102 ti 'u''" "x" "rho\', u\'" source1d_1.png
../plot1d.sh "plot 'source1d_2.dat' using 1:2 w l ls 104 ti 'rho''', 'source1d_2.dat' using 1:3 w l ls 102 ti 'u''" "x" "rho\', u\'" source1d_2.png
../plot1d.sh "plot 'viscshock_1.dat' using 1:2 w l ls 104 ti 'rho', 'viscshock_1.dat' using 1:3 w l ls 102 ti 'u', 'viscshock_1.dat' using 1:4 w l ls 101 ti 'p'" "x" "rho, u, p" viscshock_1.png
../plot1d.sh "plot 'viscshock_2.dat' using 1:2 w l ls 104 ti 'rho', 'viscshock_2.dat' using 1:3 w l ls 102 ti 'u', 'viscshock_2.dat' using 1:4 w l ls 101 ti 'p'" "x" "rho, u, p" viscshock_2.png
../plot1d.sh "plot 'viscshock_3.dat' using 1:2 w l ls 104 ti 'rho', 'viscshock_3.dat' using 1:3 w l ls 102 ti 'u', 'viscshock_3.dat' using 1:4 w l ls 101 ti 'p'" "x" "rho, u, p" viscshock_3.png
../plot1d.sh "plot 'vortexes.dat' using 1:2 w l ls 104 ti 'FiniteVortex', 'vortexes.dat' using 1:3 w l ls 102 ti 'GaussianVortex', 'vortexes.dat' using 1:4 w l ls 101 ti 'RankineVortex'" "r" "u_{phi}" vortexes_u.png
../plot1d.sh "plot 'vortexes.dat' using 1:5 w l ls 104 ti 'FiniteVortex', 'vortexes.dat' using 1:6 w l ls 102 ti 'GaussianVortex', 'vortexes.dat' using 1:7 w l ls 101 ti 'RankineVortex'" "r" "p" vortexes_p.png
../plot1d.sh "plot 'vortexincylinder.dat' using 1:3 w l ls 104 ti 'u_{phi}', 'vortexincylinder.dat' using 1:2 w l ls 102 ti 'initial data'" "r" "u_{phi}" vortexincylinder.png
../plot1d.sh "plot 'waveinchannel_3.dat' using 1:2 w l ls 101 ti 'rho''', 'waveinchannel_3.dat' using 1:3 w l ls 102 ti 'v''', 'waveinchannel_3.dat' using 1:4 w l ls 104 ti 'p'''" "r" "rho\', v\', p\'" waveinchannel_3.png
