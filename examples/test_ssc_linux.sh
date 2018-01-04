lsb_release -a
gcc -L../build_linux -Wl,-rpath=../build_linux -Wall -o test ver.c ../build_linux/ssc.so -ldl
./test
rm ./test