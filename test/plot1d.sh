#!/bin/bash

# _gnuplot <plotcmd> <xtitle> <ytitle> <outfile>

unset LD_LIBRARY_PATH

plotcmd=$1
xtitle=$2
ytitle=$3
outfile=$4

set termtype_cmd="set terminal dumb";
set outfile_ext="";

termtype_cmd_png="set terminal png size 1400,1000; set output \"$outfile\"";

gnuplot <<- EOF
#cat <<- EOF
        $termtype_cmd_png
        set encoding utf8
        unset style line
        set style function lines;
        set style line 100 lt 1 lc rgb "black" lw 2
        set style line 101 lt 1 lc rgb "blue" lw 2
        set style line 102 lt 1 lc rgb "#00AA00" lw 2  # dark-green
        set style line 103 lt 1 lc rgb "#00AAAA" lw 2  # dark-cyan
        set style line 104 lt 1 lc rgb "#CC0000" lw 2  # dark-red

        set key box linestyle 10
        set key samplen 2.5
        set key spacing 2.5
        set key top right
        
        # set style line 5 lt 0
        # set xtics 0.1
        # set ytics $cl_tick_size
        # set mytics 5
        # set mxtics 5
        # set grid xtics ytics mxtics mytics back ls 5
        set xlabel "x" font "Helvetica,16"

#        set xrange [${x_min}:${x_max}]
#        set yrange [${y_min}:${y_max}]
        set xlabel "${xtitle}" font "Helvetica,16"
        set ylabel "${ytitle}" font "Helvetica,16"
#        plot '${datafile}' using 1:2 w l ls 11
        ${plotcmd}
EOF
