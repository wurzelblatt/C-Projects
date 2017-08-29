calc: main.o add.o subtract.o multiply.o
	g++ main.o add.o subtract.o multiply.o -o calc
main.o: main.cpp
	g++ -c main.cpp
add.o: add.cpp
	g++ -c add.cpp
subtract.o: subtract.cpp
	g++ -c subtract.cpp
multiply.o: multiply.cpp
	g++ -c multiply.cpp
clean:
	rm *.o calc
