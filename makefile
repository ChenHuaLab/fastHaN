c1 = gcc -ggdb -w -O2 -mavx -mavx2 -mfma
c2 = gcc -ggdb -w -O2 -mavx -mavx2 -mfma -D_use_bz2
libs1 = -lm -lpthread -lz
libs2 = -lm -lpthread -lz -lbz2

c=$(c1)
libs=$(libs1)

cmd_rm=rm
ifeq ($(OS), Windows_NT)
	libs+= -static
	cmd_rm=del
else
	c+= -D_linux
	cmd_rm=rm -f
endif

objects = myfunction.o myutil.o tcs.o mjn.o network.o
soft=fastHapNetwork

$(soft) : $(objects)
	$(c) -o $(soft) $(objects) $(libs)

myfunction.o : myfunction.h myfunction.c
	$(c) -c myfunction.c
myutil.o : myutil.h myutil.c
	$(c) -c myutil.c
tcs.o : tcs.h tcs.c
	$(c) -c tcs.c
mjn.o : mjn.h mjn.c
	$(c) -c mjn.c
network.o : network.h network.c
	$(c) -c network.c

clean :
	$(cmd_rm) $(soft) *.exe *.o
