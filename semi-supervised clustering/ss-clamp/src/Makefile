all:
	g++ main.cpp -Wall -O2 -lsqlite3 -std=c++0x -pthread -D CHECK_NORMALIZED

lib:
	g++ -c lib.h -o lib.h.gch -Wall -O2 -lsqlite3 -std=c++0x -pthread

clean:
	rm *~
