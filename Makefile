CC = clang++
FLAG = -std=c++2b -Wextra -Wall
LIBS = -larmadillo
IMPLIB = -L/opt/homebrew/opt/armadillo/lib
IMPINC = -I./include -I/opt/homebrew/opt/armadillo/include

TARGET = out
SOURCE = examples/test.cpp src/mfdxa.cpp

.PHONY: all

all: compile run

compile:
	$(CC) $(SOURCE) -o $(TARGET) $(FLAG) $(LIBS) $(IMPLIB) $(IMPINC)

run:
	./$(TARGET) && rm -f $(TARGET)
