all: main
main: main.o lu_dense.o lu_sparse.o readMatrix.o
	g++ main.o lu_dense.o lu_sparse.o readMatrix.o -o main
lu_dense.o: lib/lu_dense.cpp
	g++ -c lib/lu_dense.cpp
lu_sparse.o: lib/lu_sparse.cpp
	g++ -c lib/lu_sparse.cpp
readMatrix.o: lib/readMatrix.cpp
	g++ -c lib/readMatrix.cpp
main.o: main.cpp
	g++ -c main.cpp 
clean:
	rm -rf *.o main
