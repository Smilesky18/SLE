all: main
main: main.o lu_dense.o lu_sparse.o readMatrix.o
	g++ main.o lu_dense.o lu_sparse.o readMatrix.o -o main -L ./
lu_dense.o: lu_dense.cpp
	g++ -c lu_dense.cpp
lu_sparse.o: lu_sparse.cpp
	g++ -c lu_sparse.cpp
readMatrix.o: readMatrix.cpp
	g++ -c readMatrix.cpp
clean:
	rm -rf *.o main
