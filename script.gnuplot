
# "time", "s", "v", "a", "j"

set term wxt size 1000, 1000

set style line 1 lc rgb "black"
set style line 2 lc rgb "blue"
set style line 3 lc rgb "violet"
set style line 4 lc rgb "green"

set multiplot layout 4,1

unset title

# Read action timestamps from second datablock
set table $actions
plot datafile index 1 with table
unset table

# Time related data
set offsets graph 0, 0, 0, 0
set size noratio
set lmargin 10
set xlabel "time [s]"

# Mark actions (red vertical lines)
set for [w in $actions[*]] arrow from first w, graph 0.0 to first w, graph 1.0 linecolor rgb 'red' nohead

# Plot data
set ylabel "position [m]"
plot datafile index 0 using 1:2 with lines ls 1 title "s"

set ylabel "velocity [m/s]"
plot datafile index 0 using 1:3 with lines ls 1 title "v"

set ylabel "acceleration [m/s²]"
plot datafile index 0 using 1:4 with lines ls 1 title "a"

set ylabel "jerk [m/s³]"
#unset yrange
#stats datafile index 0 using 1:5 nooutput
#datamargin = (STATS_max_y - STATS_min_y)*0.05
#set yrange [STATS_min_y - datamargin:STATS_max_y + datamargin]
plot datafile index 0 using 1:5 with lines ls 1 title "j"

unset multiplot
pause mouse close
