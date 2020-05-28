obj = vm.o execute.o machinecode.o print.o checkregisters.o
ex = vm
flgs = -g3 -fsanitize=address -fsanitize=undefined -stf=gnu89 -pedantic-errors -Wall -Wextra
libs = -lm

$(ex) : $(obj)
	gcc $(flgs) -o $(ex) $(obj) $(libs)

clean :
	rm $(ex) $(obj)

.c.o :
	gcc -c $(flgs) $<

.h.c :
	touch $<
