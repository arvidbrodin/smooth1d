This code is a work in progress, meant to become a jerk-limited jogging planner for linuxcnc.

Besides being functionally incomplete, it needs all assertions replaced with proper error handling. It is also licensed under GPLv3 at the moment since the code in poly.rs is ported from the GNU Science Library which is released under GPLv3+.

(LinuxCNC is GPLv2 which is incompatible with GPLv3.)

## Building

The code is written in Rust so you need a Rust compiler installed on your system. This involves running the rustup script - search the net for instructions. Then run 'cargo build' to build the code. There are no dependencies on other software packages.

## Testing

Run 'cargo test' to execute the tests in path.rs. Each test runs through the trajectory using a dt of 1 ms, checking jerk and acceleration limits at each point. Most tests also check specific values at key points (such as the position at end of move) and makes sure the move finishes.

### Viewing trajectory plots

You can view the trajectory of each test with the gnuplot script supplied. E.g. to view the trajectory produced by the jlim_continued_move test:

gnuplot -e "datafile='jlim_continued_move.data'" script.gnuplot

Each red vertical line in the plot marks a replot() command.
