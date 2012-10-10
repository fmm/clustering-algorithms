CC=g++
CFLAGS=-c -Wall -O2 -g -W
LDFLAGS=
SOURCES=Main.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=Run.out

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf *.o
	rm -rf Run.out

run:
	./Run.out config

pdf:
	pdflatex main_glass.tex > log

img:
	javac Visualization.java
	java Visualization

openpdf:
	evince main.pdf &
	
openimg:
	eog data.png &

