g++ -O3 -fopenmp Poisson.cpp main.cpp utils.cpp 
rm -f out.bin out.txt prof.log outtxt
./a.out 1  0.1 s >> log.log
python3 Analyze.py 21 21 "cmap3a.png" "3dplota.png"  "Slice of solution">> log.log


./a.out 1  0.1 rb >> log.log
python3 Analyze.py 21 21 "cmap3c_redblack.png" "3cplot3d_redblack.png" "Slice of solution"
./a.out 1  0.1 d >> log.log
python3 Analyze.py 21 21 "cmap3c_diagonal.png" "3cplota_diagonal.png" "Slice of solution"
rm -f prof.log 
for mode in s rb d; do
    for error in 0.1 0.01 0.005; do
        ./a.out 8 $error $mode >> log.log
    done
done
python3 timingsc.py "timings_c.png" 
rm -f prof.log

for mode in rb d; do
    for threads in 2 4 8 16; do
        ./a.out $threads 0.005 $mode >> log.log
    done
done
python3 Analyze.py 401 401 "cmap_accuratesol.png" "3cplota_accuratesol.png" "Slice of solution"
python3 timingsd.py "timings_d.png" 

rm -f out.bin out.txt prof.log outtxt a.out 