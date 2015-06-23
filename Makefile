
MPIRUN := mpirun -np 4

mpirun:
	make -C bcm-start mpirun="$(MPIRUN)" mpirun

clean:
	make -C bcm-start clean

distclean:
	git clean -d -f -x
