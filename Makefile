all: main
main: main.o lu_dense.o lu_sparse.o readMatrix.o
	g++ main.o lu_dense.o lu_sparse.o readMatrix.o -o main
lu_dense.o: lu_dense.cpp
	g++ -c lu_dense.cpp
lu_sparse.o: lu_sparse.cpp
	g++ -c lu_sparse.cpp
readMatrix.o: readMatrix.cpp
	g++ -c readMatrix.cpp
main.o: main.cpp
	g++ -c main.cpp 
clean:
	rm -rf *.o main
