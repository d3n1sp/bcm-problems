
mpirun: 
	make -C Test_var3 mpirun="mpirun -np 4 -H central,node12" mpirun

clean:
	make -C Test_var3 clean

distclean:
	git clean -d -f -x
