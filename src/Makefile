all: Rfuncs

Rfuncs: 
	R CMD SHLIB Rfuncs.cpp gen_lasso.cpp filtering.cpp util.cpp def.c dmatrix.c -o structreg.so

clean: 
	rm -f *.o *.so

