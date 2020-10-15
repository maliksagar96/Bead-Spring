all: simulate

simulate: header.o molcDyn.o init.o calc.o publish.o beadSpring.o 
	g++ beadSpring.cpp molcDyn.cpp init.cpp calc.cpp publish.cpp -o main.exe

header.o: molcDyn.h
	g++ -c molcDyn.h

molcDyn.o: molcDyn.cpp
	g++ -c molcDyn.cpp

init.o: init.cpp
	g++ -c init.cpp

calc.o: calc.cpp
	g++ -c calc.cpp

publish.o: publish.cpp
	g++ -c publish.cpp

beadSpring.o: beadSpring.cpp
	g++ -c beadSpring.cpp	
