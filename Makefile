all:
	cd src; make
clean:
	cd src; make clean; 
expand:
	cd src/expand; make expand
distclean:
	rm gallope expand *.nc *.dat *.tmp restart/* out2d/* out3d/*; make clean
print:
	cd src; make print
