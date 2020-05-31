CC = gcc
CFLAGS = -Wall -lm -g -O1

main_iplib: main_iplib.o bmp.o ip_lib.o
	$(CC) $(CFLAGS) main_iplib.o ip_lib.o bmp.o -lm -o main_iplib

main_iplib.o:  main_iplib.c ip_lib.h bmp.h
	$(CC) -Wall -g3 -c main_iplib.c -o main_iplib.o -lm

ip_lib.o : ip_lib.c ip_lib.h bmp.h
	$(CC) $(CFLAGS) -c ip_lib.c -o ip_lib.o -lm

bmp.o : bmp.c bmp.h
	$(CC) -Wall -g3 -c bmp.c -o bmp.o -lm

clean :
	rm *.o main_iplib
