all: build

build: mandel-tiles.cpp
	g++ -o prog mandel-tiles.cpp -pthread

run: build
	./prog $(ARGS)

clean: 
	/bin/rm -f *.o prog
