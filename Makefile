#all: main
#main: main.o lu_dense.o lu_sparse.o readMatrix.o
#	g++ main.o lu_dense.o lu_sparse.o readMatrix.o -lsuperlu -lblas -o main
#lu_dense.o: lib/lu_dense.cpp
#	g++ -c lib/lu_dense.cpp
#lu_sparse.o: lib/lu_sparse.cpp
#	g++ -c lib/lu_sparse.cpp
#readMatrix.o: lib/readMatrix.cpp
#	g++ -c lib/readMatrix.cpp
#main.o: main.cpp
#	g++ -c main.cpp 
#clean:
#	rm -rf *.o main



all: main.cpp
	g++ lib/lu_dense.cpp lib/lu_sparse.cpp lib/readMatrix.cpp lib/super_lu.cpp main.cpp -lsuperlu -lblas -o main
debug: main.cpp
	g++ -g lib/lu_dense.cpp lib/lu_sparse.cpp lib/readMatrix.cpp lib/super_lu.cpp main.cpp -lsuperlu -lblas -o main.g
clean:
	rm -rf *.o main
