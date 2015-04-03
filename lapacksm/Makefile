f       = f
c       = cpp
o       = o
a       = a
r       =

.SUFFIXES: .$f .$c .$o

CC      = gCC -O3 -c
CF      = mpif77 -O3 -c
CFNO    = mpif77 -O0 -c  #-fno-unroll-loops -ffloat-store -c
LINK    = gCC

LIB     = lapacksm.$a
LIBOBJS = \
	daxpy.$o dbdsqr.$o dcopy.$o ddot.$o dgebd2.$o dgebrd.$o \
	dgelq2.$o dgelqf.$o dgemm.$o dgemv.$o dgeqr2.$o dgeqrf.$o \
	dger.$o dgesvd.$o dlabrd.$o dlacpy.$o dladiv.$o \
	dlae2.$o dlaev2.$o dlamch.$o dlange.$o dlanst.$o \
	dlansy.$o dlapy2.$o dlapy3.$o dlarf.$o dlarfb.$o \
	dlarfg.$o dlarft.$o dlartg.$o dlas2.$o dlascl.$o \
	dlaset.$o \
	dlasq1.$o dlasq2.$o dlasq3.$o dlasq4.$o dlasq5.$o dlasq6.$o \
	dlasr.$o dlasrt.$o dlassq.$o \
	dlasv2.$o dlatrd.$o dnrm2.$o dorg2l.$o dorg2r.$o \
	dorgbr.$o dorgl2.$o dorglq.$o dorgql.$o dorgqr.$o \
	dorgtr.$o dorm2r.$o dormbr.$o dorml2.$o dormlq.$o \
	dormqr.$o drot.$o dscal.$o dsteqr.$o dsterf.$o \
	dswap.$o dsyev.$o dsymv.$o dsyr2.$o dsyr2k.$o \
	dsytd2.$o dsytrd.$o dtrmm.$o dtrmv.$o dznrm2.$o \
	ilaenv.$o lsame.$o xerbla.$o zbdsqr.$o zcopy.$o \
	zdrot.$o zdscal.$o zgebd2.$o zgebrd.$o zgelq2.$o \
	zgelqf.$o zgemm.$o zgemv.$o zgeqr2.$o zgeqrf.$o \
	zgerc.$o zgesvd.$o zlabrd.$o zlacgv.$o zlacpy.$o \
	zladiv.$o zlange.$o zlarf.$o zlarfb.$o zlarfg.$o \
	zlarft.$o zlascl.$o zlaset.$o zlasr.$o zlassq.$o \
	zscal.$o zswap.$o ztrmm.$o ztrmv.$o zung2r.$o \
	zungbr.$o zungl2.$o zunglq.$o zungqr.$o zunm2r.$o \
	zunmbr.$o zunml2.$o zunmlq.$o zunmqr.$o \

$(LIB): $(LIBOBJS)
	ar r $(LIB) $(LIBOBJS)
	ranlib $(LIB)

dlamch.$f dlamch.$o:
	$(CFNO) dlamch.$f

.$f.$o:
	$(CF) $<

.$c.$o:
	$(CC) $<

clean:
	-rm *.o
	-rm $(LIB)
