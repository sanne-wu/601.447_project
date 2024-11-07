CC = g++
CONSERVATIVE_FLAGS = -std=c++14 -Wall -Wextra -pedantic
DEBUGGING_FLAGS = -g -O0
CFLAGS = $(CONSERVATIVE_FLAGS) $(DEBUGGING_FLAGS)

minhash: main.o minhash.o
	$(CC) -o minhash main.o minhash.o

main.o: main.cpp minhash.h
	$(CC) $(CFLAGS) -c main.cpp

minhash.o: minhash.cpp minhash.h
	$(CC) $(CFLAGS) -c minhash.cpp

clean:
	rm -f *.o minhash
