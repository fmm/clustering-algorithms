all:
	g++ main.cpp -Wall -O2 -lsqlite3 -std=c++0x -pthread -D CHECK_NORMALIZED

ratio:
	g++ main.cpp -Wall -O2 -lsqlite3 -std=c++0x -pthread -D CHECK_NORMALIZED -D RATIO_ALPHA

beta:
	g++ main.cpp -Wall -O2 -lsqlite3 -std=c++0x -pthread -D CHECK_NORMALIZED -D BETA=0.01

lib:
	g++ -c lib.h -o lib.h.gch -Wall -O2 -lsqlite3 -std=c++0x -pthread

clean:
	rm *~
