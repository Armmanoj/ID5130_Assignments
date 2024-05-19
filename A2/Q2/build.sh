g++ serial_jb.cpp 
./a.out 0.01 >> timing.txt 
rm -f timing.txt
mpicc -O3 mpi_redblack.c -lm
for threads in 2 4 8; do
    mpiexec -n $threads ./a.out 0.01 a.bin >> timing.txt
    python3 Analyze.py 201 201 rb$threads.png out.bin a.bin
done

mpicc -O3 mpi_jb.c -lm
for threads in 2 4 8; do
    mpiexec -n $threads ./a.out 0.01 a.bin >> timing.txt
    python3 Analyze.py 201 201 jb$threads.png out.bin a.bin
done