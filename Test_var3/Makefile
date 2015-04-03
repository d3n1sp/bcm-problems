f       = f
c       = cpp
o       = o
a       = a
r       = exe

.SUFFIXES: .$f .$c .$o

CC      = mpicxx -c -O3 -I../Bcmlib -I../qd -I../Bslib
CFCF    =
LINK    = mpicxx -L/usr/lib -Wl,-rpath,\$$ORIGIN

RUN     = lame2d-par_comput
RUNOBJS = Test.$o
MYLIBS  = bcmlib.$a qd.$a bslib.$a lapacksm.$a /usr/lib/libmetis.$a

$(RUN): $(RUNOBJS) $(MYLIBS)
	$(LINK) -o $(RUN) $(RUNOBJS) $(MYLIBS) -lpthread -lrt -lm -lgfortran 

.$f.$o:
	$(CF) $<

.$c.$o:
	$(CC) $<

clean:
	-rm *.o
	-rm $(RUN)
	make -C ../Bcmlib clean
	make -C ../qd clean
	make -C ../Bslib clean
	make -C ../lapacksm clean

run:
	./$(RUN)

mpirun:
	mpirun -np 4 ./$(RUN)

