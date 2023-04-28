reset

# png
set terminal pngcairo dashed size 1000,800 enhanced font "Verdana,30"
#set terminal pngcairo dashed size 2500,1750 enhanced font "Verdana,50"
# dashed option enables dashed linestyle in pngcairo

# Custom line styles

LW1= 4.0
LW2 = 4.0
PS = 4

set style line 1 lt 1 dt 1 lw LW1 lc rgb '#0072bd' # blue
set style line 2 lt 1 dt 9 lw LW2 lc rgb '#d95319' # orange
set style line 3 lt 1 lw LW1 lc rgb '#edb120' # yellow
set style line 4 lt 1 lw LW1 lc rgb '#7e2f8e' # purple
set style line 5 lt 1 lw LW1 lc rgb '#77ac30' # green
set style line 6 lt 1 lw LW1 lc rgb '#4dbeee' # light-blue
set style line 7 lt 1 lw LW1 lc rgb '#a2142f' # red

set style line 102 lc rgb '#808080' lt 0 lw 3
#set grid back ls 102

# Set the border using the linestyle 80 that we defined
# 3 = 1 + 2 (1 = plot the bottom line and 2 = plot the left line)
# back means the border should be behind anything else drawn
# set style line 80 lt 0 lw 3 lc rgb "#808080"
# set border 3 back ls 80

# MACROS
# x- and ytics for each row resp. column
NOXTICS = "set xtics ('' 1.0,'' 1.4,'' 1.8); \
          unset xlabel; \
          set mxtics 4 "
XTICS = "set xtics 1.0,0.4,2; \
          set xlabel 'r_{ij} (Å)'; \
          set mxtics 4 "
NOYTICS = "set ytics ( '' -30, '' 0, '' 30 ); \
           set mytics 2 ; \
           unset ylabel"
YTICS = "set ytics auto"

NOKEY = "unset key"
KEY = "set key top left"

# Margins for each row resp. column
# ---- 0.95
# |  |
# ---- 0.70
# |  |
# ---- 0.45
# |  |
# ---- 0.20

TMARGIN = "set tmargin at screen 0.95; set bmargin at screen 0.70"
MMARGIN = "set tmargin at screen 0.70; set bmargin at screen 0.45"
BMARGIN = "set tmargin at screen 0.45; set bmargin at screen 0.20"
LMARGIN = "set lmargin at screen 0.15; set rmargin at screen 0.55"
RMARGIN = "set lmargin at screen 0.55; set rmargin at screen 0.95"

# Placement of the a,b,c,d labels in the graphs
POS = "at graph 0.75, 0.80 font 'helvetica, 30'"
POS2 = "at graph 0.70, 0.55 font 'helvetica, 25'"
POS3 = "at graph 0.70, 0.30 font 'helvetica, 25'"

# Enable the use of macros
set macros

set output "fit_plots.png"

#YBOT = 0
#YTOP = 10
#set yrange [YBOT : YTOP]
XBOT = 1.1
XTOP = 1.9
set xrange [XBOT : XTOP]

# Fit this function
f(x)=A*exp(-(x-mu)**2/(2*sigma**2))

# Start multiplot
set multiplot layout 3,1 rowsfirst

# --- GRAPH a
@TMARGIN; @LMARGIN
@NOXTICS; @NOYTICS
@NOKEY
fit f(x) 'r12_acc.dat' using 1:2 via A,mu,sigma
set label 1 'r_{12}' @POS
set label 2 sprintf("m = %3.2f", mu) @POS2
set label 3 sprintf("s = %3.2f", sigma) @POS3
set arrow from mu,0 to mu,f(mu) nohead lw 3
plot "r12_acc.dat" u 1:2, f(x)
unset arrow

# --- GRAPH b
@MMARGIN; @LMARGIN
@NOXTICS; @NOYTICS
set ylabel "1/\chi^2"
@NOKEY
fit f(x) 'r23_acc.dat' using 1:2 via A,mu,sigma
set label 1 'r_{23}' @POS
set label 2 sprintf("m = %3.2f", mu) @POS2
set label 3 sprintf("s = %3.2f", sigma) @POS3
set arrow from mu,0 to mu,f(mu) nohead lw 3
plot "r23_acc.dat" u 1:2, f(x)
unset arrow

# --- GRAPH c
@BMARGIN; @LMARGIN
@XTICS; @NOYTICS
@NOKEY
fit f(x) 'r34_acc.dat' using 1:2 via A,mu,sigma
set label 1 'r_{34}' @POS
set label 2 sprintf("m = %3.2f", mu) @POS2
set label 3 sprintf("s = %3.2f", sigma) @POS3
set arrow from mu,0 to mu,f(mu) nohead lw 3
plot "r34_acc.dat" u 1:2, f(x)
unset arrow

# --- GRAPH d
@TMARGIN; @RMARGIN
@NOXTICS; @NOYTICS
@NOKEY
fit f(x) 'r45_acc.dat' using 1:2 via A,mu,sigma
set label 1 'r_{45}' @POS
set label 2 sprintf("m = %3.2f", mu) @POS2
set label 3 sprintf("s = %3.2f", sigma) @POS3
set arrow from mu,0 to mu,f(mu) nohead lw 3
plot "r45_acc.dat" u 1:2, f(x)
unset arrow

# --- GRAPH e
@MMARGIN; @RMARGIN
@NOXTICS; @NOYTICS
@NOKEY
fit f(x) 'r56_acc.dat' using 1:2 via A,mu,sigma
set label 1 'r_{56}' @POS
set label 2 sprintf("m = %3.2f", mu) @POS2
set label 3 sprintf("s = %3.2f", sigma) @POS3
set arrow from mu,0 to mu,f(mu) nohead lw 3
plot "r56_acc.dat" u 1:2, f(x)
unset arrow

# --- GRAPH f
@BMARGIN; @RMARGIN
@XTICS; @NOYTICS
@NOKEY
fit f(x) 'r16_acc.dat' using 1:2 via A,mu,sigma
set label 1 'r_{16}' @POS
set label 2 sprintf("m = %3.2f", mu) @POS2
set label 3 sprintf("s = %3.2f", sigma) @POS3
set arrow from mu,0 to mu,f(mu) nohead lw 3
plot "r16_acc.dat" u 1:2, f(x)
unset arrow

unset multiplot
### End multiplot
