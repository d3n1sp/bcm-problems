
mpirun: 
	make -C bcm-start mpirun="mpirun -np 4 -H central,node12" mpirun

clean:
	make -C bcm-start clean

distclean:
	git clean -d -f -x
