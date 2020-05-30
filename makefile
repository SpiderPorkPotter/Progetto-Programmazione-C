CC = gcc
CFLAGS = -Wall --ansi --pedantic -lm -g3 -O3 -fsanitize=address -fsanitize=undefined -std=gnu89 -Wextra

main_iplib: main_iplib.c
	#gcc main_iplib.c ip_lib.c bmp.c -lm -g -Wall -o main_iplib
	$(CC) $(CFLAGS) main_iplib.c ip_lib.c bmp.c -o main_iplib

ip_lib.o : ip_lib.c ip_lib.h
	$(CC) $(CFLAGS) ip_lib.c -o ip_lib.o
	#gcc -Wall --ansi --pedantic -lm -g3 -O3 -fsanitize=address -fsanitize=undefined -std=gnu89 -Wextra ip_lib.c -o ip_lib.o

bmp.o : bmp.c bmp.h
	$(CC) $(CFLAGS) bmp.c -o bmp.o
	#gcc -c -Wall --ansi --pedantic -lm -g3 -O3 -fsanitize=address -fsanitize=undefined -std=gnu89 -Wextra bmp.c -o bmp.o
clean :
	rm *.o main_iplib
