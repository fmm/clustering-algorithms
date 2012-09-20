CC=g++
CFLAGS=-c -Wall -O2 -g -W
LDFLAGS=
SOURCES=Main.cpp src/Algoritmo.cpp src/Resultado.cpp src/Array.cpp src/Cluster.cpp src/Dados.cpp src/Imprime.cpp src/Includes.cpp src/Repositorio.cpp src/Tabela.cpp sodas/MonteCarlo.cpp sodas/Pattern.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=Run.out

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf *.o
	rm -rf src/*.o
	rm -rf sodas/*.o

run:
	./Run.out config

