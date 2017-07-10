CC=g++
CFLAGS=-O3 -fpermissive -I fftw++-1.13/ -L /usr/lib/ -lfftw3_threads -lfftw3 -lm #-I boost_1_53_0/

all: phi4model

phi4model: main.cpp
	 $(CC) main.cpp $(CFLAGS) -o phi4model

clear: 
	rm phi4model
