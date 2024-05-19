rm -f timing.txt
mpicc -O3 mpi_redblack.c -lm
for threads in 2 4 8 10 ; do
    mpiexec -n $threads ./a.out 0.005 a.binn >> timing.txt
done
mpicc -O3 mpi_jb.c -lm
for threads in 2 4 8 10 ; do
    mpiexec -n $threads ./a.out 0.005 a.bin >> timing.txt
done