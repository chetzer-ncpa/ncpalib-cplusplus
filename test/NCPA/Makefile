GTESTLIBS=-lgtest_main -lgtest
LDFLAGS=${LIBRARY_PATH_FLAGS} $(GTESTLIBS)
INCFLAGS=-I../include -I../../include ${INCLUDE_PATH_FLAGS}
CXXFLAGS=-c -fPIC -std=${STANDARD} $(INCFLAGS)

SOURCES=$(wildcard test_*.cpp)
OBJECTS=$(SOURCES:.cpp=.o)
TARGET=ncpalibTestSuite

test: build

build: $(OBJECTS)
	g++ -o $(TARGET) $^ $(LDFLAGS) && ./$(TARGET)
	
%.o: %.cpp
	${CXX} $^ $(CXXFLAGS) -o $@
	
clean:
	-rm $(OBJECTS) $(TARGET)
