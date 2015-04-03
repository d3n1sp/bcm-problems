f       = f
c       = cpp
o       = o
a       = a
r       =

.SUFFIXES: .$f .$c .$o

CC      = mpicxx -O3 -c 
CF      = mpif77 -O3 -c 
LINK    =

LIB     = qd.$a
LIBOBJS = \
	./bits.$o \
	./c_dd.$o \
	./c_qd.$o \
	./dd_const.$o \
	./dd_real.$o \
	./fpu.$o \
	./qd_const.$o \
	./qd_real.$o \
	./util.$o \

$(LIB): $(LIBOBJS)
	ar r $(LIB) $(LIBOBJS)
	ranlib $(LIB)

.$f.$o:
	$(CF) $<

.$c.$o:
	$(CC) $<

clean:
	-rm *.o
	-rm $(LIB)

