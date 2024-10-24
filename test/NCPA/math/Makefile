
GTESTLIBS=-lgtest_main -lgtest
LDFLAGS=$(LIBRARY_FLAGS) $(GTESTLIBS)
INCFLAGS=-I../include -I../../include $(INCLUDE_FLAGS)
CXXFLAGS=-c -fPIC -ggdb -std=c++20 $(INCFLAGS)

SOURCES=$(wildcard test_*.cpp)
OBJECTS=$(SOURCES:.cpp=.o)
TARGET=mathTest

test: build clean

build: $(OBJECTS)
	g++ -o $(TARGET) $^ $(LDFLAGS) && ./$(TARGET)
	
%.o: %.cpp
	${CXX} $^ $(CXXFLAGS) -o $@
	
clean:
	-rm $(OBJECTS) $(TARGET)
