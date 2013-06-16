#
#if 0
#	define EXTREMELY_DETAILED
#	define DETAILED
#	define PROTOTYPE_VARIANCE
#	define MAGIC_VARIANCE
#	define RUN_ANOTHER
#endif
#

tmp:
	g++ -Wall -O2 -g -W Main.cpp -o Run.out

all:
	g++ -Wall -O2 -g -W Main.cpp -o Run.out -D RUN_ANOTHER -D DETAILED

extremely_detailed:
	g++ -Wall -O2 -g -W Main.cpp -o Run.out -D RUN_ANOTHER -D EXTREMELY_DETAILED
	
prototype:
	g++ -Wall -O2 -g -W Main.cpp -o Run.out -D PROTOTYPE_VARIANCE
	
magic:
	g++ -Wall -O2 -g -W Main.cpp -o Run.out -D MAGIC_VARIANCE

prototype_magic:
	g++ -Wall -O2 -g -W Main.cpp -o Run.out -D MAGIC_VARIANCE -D PROTOTYPE_VARIANCE
	
clean:
	rm -rf Run.out

go:
	./Run.out

run:
	./Run.out config

pdf:
	pdflatex main.tex > log

img:
	javac Visualization.java
	java Visualization

openpdf:
	evince main.pdf &
	
openimg:
	eog data.png &

