Build the solution, inside folder Q3
chmod 777 build.sh
./build.sh

The output png files include-
"cmap3a.png" "3dplota.png" which store the output of part a, serial code with delta = 0.1
"cmap3c_redblack.png" "3cplot3d_redblack.png"  and "cmap3c_diagonal.png" "3cplota_diagonal.png"
which has plots verifying that the parallel code and serial code give same output (part c).
"cmap_accuratesol.png" "3cplota_accuratesol.png" which has a highly accurate solution to the poisson's equation.
"timings_d.png"  which has the timings asked for in part (d)
"timings_c.png"  which has the timings asked for in part (c)
The timing have been plotted seperately for serial, red-black and d
