CC = gcc
CFLAGS = -Wall --ansi --pedantic -g3 -O3 -fsanitize=address -fsanitize=undefined -std=gnu89 -Wextra

main_iplib: main_iplib.o bmp.o ip_lib.o
	#gcc main_iplib.c ip_lib.c bmp.c -lm -g -Wall -o main_iplib
	$(CC) $(CFLAGS) main_iplib.o ip_lib.o bmp.o -lm -o main_iplib

main_iplib.o:  main_iplib.c
	$(CC)  -Wall -g3 -O3 -c main_iplib.c -o main_iplib.o -lm

ip_lib.o : ip_lib.c ip_lib.h
	$(CC) $(CFLAGS) -c ip_lib.c -o ip_lib.o -lm
	#gcc -Wall --ansi --pedantic -lm -g3 -O3 -fsanitize=address -fsanitize=undefined -std=gnu89 -Wextra ip_lib.c -o ip_lib.o

bmp.o : bmp.c bmp.h
	$(CC) -Wall -g3 -O3 -c bmp.c -o bmp.o -lm
	#gcc -c -Wall --ansi --pedantic -lm -g3 -O3 -fsanitize=address -fsanitize=undefined -std=gnu89 -Wextra bmp.c -o bmp.o
clean :
	rm *.o main_iplib
