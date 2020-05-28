all: ip_lib.o bmp.o
	gcc main_iplib.c ip_lib.c bmp.c -Wall --ansi --pedantic -lm -g3 -O3 -fsanitize=address -fsanitize=undefined -std=gnu89 -Wextra -g -Wall -O1 -o test

ip_lib.o :
	gcc -c -Wall --ansi --pedantic -lm -g3 -O3 -fsanitize=address -fsanitize=undefined -std=gnu89 -Wextra ip_lib.c

bmp.o :
	gcc -c -Wall --ansi --pedantic -lm -g3 -O3 -fsanitize=address -fsanitize=undefined -std=gnu89 -Wextra bmp.c

clean:
	rm ip_lib.o
	rm bmp.o
