main_iplib : ip_lib.c bmp.c main_iplib.c
	gcc main_iplib.c ip_lib.c bmp.c -lm -g -Wall -o main_iplib

ip_lib.o : ip_lib.c
	gcc -c ip_lib.c -o ip_lib.o

bmp.o : bmp.c
	gcc -c bmp.c -o bmp.o

main_iplib.o : main_iplib.c
	gcc -c main_iplib.c -o main_iplib.o

clean :
	rm *.o main_iplib
