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
set grid back ls 102

# Set the border using the linestyle 80 that we defined
# 3 = 1 + 2 (1 = plot the bottom line and 2 = plot the left line)
# back means the border should be behind anything else drawn
# set style line 80 lt 0 lw 3 lc rgb "#808080"
# set border 3 back ls 80

# MACROS
# x- and ytics for each row resp. column
NOXTICS = "set xtics ('' 0.5,'' 1,'' 1.5,'' 2,'' 2.5,'' 3); \
          unset xlabel; \
          set mxtics 2 "
XTICS = "set xtics 0.5,0.5,3; \
          set xlabel 'r_{16} (Å)'; \
          set mxtics 2 "
NOYTICS = "set ytics ( '-30' -30, '0' 0, '30' 30 ); \
           set mytics 2 ; \
           unset ylabel"
YTICS = "set ytics 0,1.0,3; \
           set mytics 2"

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
LMARGIN = "set lmargin at screen 0.15; set rmargin at screen 0.95"
#RMARGIN = "set lmargin at screen 0.55; set rmargin at screen 0.95"

# Placement of the a,b,c,d labels in the graphs
POS = "at graph 0.85, 0.80 font 'helvetica, 40'"

# Enable the use of macros
set macros

set output "PLOT_R16_R12_R23_R34.png"

YBOT = 0.5
YTOP = 2.9
set yrange [YBOT : YTOP]
set xrange [0.5 : 3]

# Start multiplot
set multiplot layout 3,1 rowsfirst
# --- GRAPH a
@TMARGIN; @LMARGIN
@NOXTICS; @YTICS
set ylabel "r_{12} (Å)"
@NOKEY
set label 1 '' @POS
#set logscale x 10
#set logscale y 10
plot "< paste r16.dat r12.dat" u 2:4  w p ls 1
# --- GRAPH b
@MMARGIN; @LMARGIN
@NOXTICS; @YTICS
set ylabel "r_{23} (Å)"
@NOKEY
plot "< paste r16.dat r23.dat" u 2:4  w p ls 1
set label 1 '' @POS
# --- GRAPH c
@BMARGIN; @LMARGIN
@XTICS; @YTICS
set ylabel "r_{34} (Å)"
@NOKEY
set label 1 '' @POS
plot "< paste r16.dat r34.dat" u 2:4  w p ls 1
unset multiplot
### End multiplot
