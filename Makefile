CC=g++

all: build

build: cppga.cpp
	$(CC) cppga.cpp -o cppga

test: test.cpp
	$(CC) test.cpp -o test

clean:
	rm a.out* cppga logs.txt fitfile*
