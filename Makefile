all: main.cpp
	g++ -O2 main.cpp Source/lu_gp_sparse.cpp -fpermissive Lib/libcxsparse.a Lib/libsuitesparseconfig.a -fPIC -lsuperlu -lblas -lm -o main
	
debug: main.cpp
	g++ -g main.cpp Source/lu_gp_sparse.cpp -fpermissive Lib/libcxsparse.a Lib/libsuitesparseconfig.a -fPIC -lsuperlu -lblas -lm -o main.g
clean:
	rm -rf  main main.g
