f       = f
c       = cpp
o       = o
a       = a
r       = exe

.SUFFIXES: .$f .$c .$o

CC      = mpicxx -c -O3 -I../Bcmlib -I../qd -I../Bslib
CFCF    =
LINK    = mpicxx -L/usr/lib -Wl,-rpath,\$$ORIGIN
mpirun := mpirun -np 4

RUN     = lame2d-par_comput
RUNOBJS = Test.$o
MYLIBS  = bcmlib.$a qd.$a bslib.$a lapacksm.$a /usr/lib/libmetis.$a

$(RUN): $(RUNOBJS) $(MYLIBS)
	$(LINK) -o $(RUN) $(RUNOBJS) $(MYLIBS) -lpthread -lrt -lm -lgfortran

$(MYLIBS):
	make -C ../Bcmlib -j
	make -C ../qd 
	make -C ../Bslib -j
	make -C ../lapacksm 

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

run: $(RUN)
	./$(RUN)

mpirun: $(RUN)
	$(mpirun) ./$(RUN)

