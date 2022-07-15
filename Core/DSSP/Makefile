CC = g++

# gcc options
	CFLAGS = -O3 -lm 

	FLAGS = -O3 
LDFLAGS= -shared
SOURCES = $(shell echo *.cpp)
HEADERS = $(shell echo *.h)
OBJECTS = $(SOURCES:.cpp=.o)

TARGET = libDSSP.so

all: $(TARGET)

clean:
	rm -f $(OBJECTS) $(TARGET)



$(TARGET) : mas.o mkdssp.o primitives-3d.o structure.o utils.o
	$(CC) -shared $(CFLAGS) $(OBJECTS) -o $@ $(LDFLAGS)
	cp $(TARGET) ..

mas.o: mas.cpp
	g++  -c -fPIC $(CFLAGS) mas.cpp 
mkdssp.o: mkdssp.cpp
	g++ -c -fPIC $(CFLAGS) mkdssp.cpp 
primitives-3d.o: primitives-3d.cpp
	g++ -c -fPIC $(CFLAGS)  primitives-3d.cpp 
structure.o: structure.cpp
	g++ -c -fPIC $(CFLAGS) structure.cpp 
utils.o: utils.cpp
	g++ -c -fPIC $(CFLAGS) utils.cpp 
	
