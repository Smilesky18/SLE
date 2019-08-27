all: main.cpp
	g++ -O2 main.cpp Source/lu_gp_sparse.cpp -fpermissive Lib/libcxsparse.a Lib/libsuitesparseconfig.a -fPIC -lsuperlu -lblas -lm -o sparse_lu
	
debug: main.cpp
	g++ -g main.cpp Source/lu_gp_sparse.cpp -fpermissive Lib/libcxsparse.a Lib/libsuitesparseconfig.a -fPIC -lsuperlu -lblas -lm -o sparse_lu.g
	
test: main.cpp
	g++ -O3 -fexceptions -fPIC -fopenmp main.cpp Source/lu_gp_sparse.cpp -fpermissive Lib/libcxsparse.a Lib/libsuitesparseconfig.a -lsuperlu -lblas -lm -lrt -W -o sparse_lu
	
clean:
	rm -rf  sparse_lu sparse_lu.g
