#all: main
#main: main.o lu_dense.o lu_sparse.o readMatrix.o
#	g++ main.o lu_dense.o lu_sparse.o readMatrix.o -lsuperlu -lblas -o main
#lu_dense.o: lib/lu_dense.cpp
#	gcc -c lib/lu_dense.cpp
#lu_sparse.o: lib/lu_sparse.cpp
#	gcc -c lib/lu_sparse.cpp
#readMatrix.o: lib/readMatrix.cpp
#	g++ -c lib/readMatrix.cpp
#main.o: main.cpp
#	gcc -c main.cpp 
#clean:
#	rm -rf *.o main



all: main.cpp
	#g++ -g lib/lu_dense.cpp lib/lu_sparse.cpp lib/readMatrix.cpp lib/super_lu.cpp main.cpp lib/microtime.cpp -fpermissive /home/robot/my_solver/lu_solver_yzy/build/SRC/libsuperlu.a -lblas -lm -o main -O3
	#g++ lib/lu_dense.cpp lib/lu_sparse.cpp lib/readMatrix.cpp lib/super_lu.cpp main.cpp lib/microtime.cpp -fpermissive /usr/local/lib/libsuperlu.a -lblas -lm -o main -O3
	#g++ -g lib/lu_dense.cpp lib/lu_sparse.cpp lib/lu_sparse_v2.cpp lib/readMatrix.cpp lib/super_lu.cpp main.cpp lib/microtime.cpp -fpermissive -lsuperlu -lblas -lm -o main -O3
debug: main_compare.cpp
	g++ -g lib/super_lu.cpp lib/lu_dense.cpp lib/lu_gp.cpp lib/lu_gp_amd.cpp lib/lu_gp_v2.cpp lib/lu_gp_v3.cpp lib/readMatrix.cpp main_compare.cpp lib/microtime.cpp lib/lu_gp_sparse.cpp -fpermissive /usr/local/lib/libcxsparse.a /usr/local/lib/libamd.a /usr/local/lib/libsuitesparseconfig.a -fPIC -lsuperlu -lblas -lm -o  main_compare.g
clean:
	rm -rf  main main_compare.g
