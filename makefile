# makefile
CC = g++ -std=c++11
CFLAGS  = -g -O2 -Wall -Werror
TARGET = tonelli_shanks
OUTPUT = tonelli_shanks
DEPS = -lgmp

all: $(TARGET)

$(TARGET): $(TARGET).cpp
	$(CC) $(CFLAGS) $(TARGET).cpp -o $(OUTPUT) $(DEPS)

clean:
	rm $(OUTPUT)
