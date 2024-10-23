
GTESTLIBS=-lgtest_main -lgtest
LDFLAGS=$(GTESTLIBS)
INCFLAGS=-I../include -I../../include $(INCLUDE_FLAGS)
CXXFLAGS=-c -fPIC -ggdb -std=c++20 $(INCFLAGS)

SOURCES=$(wildcard test_*.cpp)
OBJECTS=$(SOURCES:.cpp=.o)
TARGET=commonTest

test: build clean

build: $(OBJECTS)
	#g++ -o $(TARGET) $^ $(patsubst test_%.o,%.o,$^) $(LDFLAGS) && ./$(TARGET)
	g++ -o $(TARGET) $^ $(LDFLAGS) && ./$(TARGET)
	
%.o: %.cpp
	${CXX} $^ $(CXXFLAGS) -o $@
	#${CXX} $(patsubst test_%.cpp,../src/%.cpp,$^) $(CXXFLAGS) -o $(patsubst test_%.cpp,%.o,$^)
	
clean:
	-rm $(OBJECTS) $(patsubst test_%.o,%.o,$(OBJECTS)) $(TARGET)
