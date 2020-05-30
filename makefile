main_iplib : ip_lib.c bmp.c main_iplib.c ip_lib.h bmp.h
	gcc main_iplib.c ip_lib.c bmp.c -lm -g -Wall -o main_iplib

ip_lib.o : ip_lib.c ip_lib.h
	gcc -c --ansi --pedantic -ggdb ip_lib.c -o ip_lib.o

bmp.o : bmp.c bmp.h
	gcc -c bmp.c -o bmp.o

clean :
	rm *.o main_iplib
