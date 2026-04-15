include test_order.make

GTESTLIBS=-lgtest_main -lgtest
LDFLAGS=${LIBRARY_PATH_FLAGS} $(GTESTLIBS) -lgslcblas -lgsl
INCFLAGS=-I../include -I../../include ${INCLUDE_PATH_FLAGS}
CXXFLAGS=-c -fPIC -ggdb -std=${STANDARD} $(INCFLAGS)

SOURCES=test_atmosphere.cpp
# SOURCES=test_atmosphere.cpp
# SOURCES=test_ndvector.cpp
OBJECTS=$(SOURCES:.cpp=.o)
TARGET=ncpalibTestSuite

test: build

build: $(OBJECTS)
	g++ -o $(TARGET) $^ $(LDFLAGS) && ./$(TARGET)
	
%.o: %.cpp
	${CXX} $^ $(CXXFLAGS) -o $@
	
clean:
	-rm $(OBJECTS) $(TARGET)
