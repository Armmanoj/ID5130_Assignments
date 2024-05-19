g++ -fopenmp -O3 main.cpp Util.cpp Q2.cpp
rm -f debug.log out.txt profile.log
./a.out 1  25
python3 Analysis.py "Q2a.png"
rm profile.log
./a.out 2 100
python3 Analysis.py "Q2bn_100.png"
rm profile.log
./a.out 2 1000
./a.out 4 1000
./a.out 8 1000

#for i in $(seq 2 13); do 
#   ./a.out $i 100000
#done

python3 timing.py "Q2btiming.png"
rm -f debug.log out.txt profile.log a.out