g++ main.cpp
./a.out q
python3 plot.py "Serial QUICK Scheme"
rm out.txt
./a.out u
python3 plot.py "Serial upwind Scheme"
rm out.txt

mpicc -g -Wall upwind.c -lm
mpiexec -n 2 ./a.out
python3 plot.py "Parallel Upwind scheme, 2 threads"
rm out.txt
mpicc -g -Wall upwind.c -lm
mpiexec -n 4 ./a.out
python3 plot.py "Parallel Upwind scheme, 4 threads"
rm out.txt

mpicc -g -Wall QUICK.c -lm
mpiexec -n 2 ./a.out
python3 plot.py "Parallel QUICK scheme, 2 threads"
rm out.txt
mpicc -g -Wall QUICK.c -lm
mpiexec -n 4 ./a.out
python3 plot.py "Parallel QUICK scheme, 4 threads"
rm out.txt