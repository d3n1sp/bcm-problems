f       = f
c       = cpp
o       = o
a       = a
r       =

.SUFFIXES: .$f .$c .$o

CC      = mpicxx -O3 -c # -D__PMPITRACE__ # -DLAM_BUILDING
CF      = mpif77 -O3 -c
LINK    = mpicxx -m32

LIB     = bslib.$a
LIBOBJS = \
	atac.$o cg.$o corr.$o Decomp.$o DecompBnd.$o DecompSplit.$o ExchangeMPI.$o fct.$o fourier.$o func.$o \
	globals.$o gmres.$o gps.$o gsmatrix.$o gsmatrixc.$o hyst.$o ich2browc.$o \
	ich2gcs.$o ich2r.$o ich2rb.$o ich2rs.$o ilu2browc.$o \
	ilu2browr.$o ilu2browr2ind.$o ilu2cs.$o \
	ilu2gcs.$o ilu2grs.$o ilu2grs2ind.$o \
	ilu2r.$o ilu2rs.$o lanczos.$o lanczosc.$o lapackd.$o lapackdm.$o \
	mvm.$o mvmc.$o newton.$o ordmatr.$o ordmatrnd.$o parfbss.$o qrd.$o rcm.$o sendrecv.$o smatrix.$o smatrixc.$o \
	smatrixcs.$o smatrixr.$o smatrixrb.$o smatrixrs.$o sofgmres.$o solvers.$o sor.$o svector.$o \
	svectorc.$o tree.$o unroll.$o

$(LIB): $(LIBOBJS)
	ar r $(LIB) $(LIBOBJS)
	ranlib $(LIB)

.$f.$o:
	$(CF) $<

.$c.$o:
	$(CC) $<

.c.$o:
	$(CC) $*.$c

clean:
	-rm *.o
	-rm $(LIB)
