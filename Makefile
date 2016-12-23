UNAME := $(shell uname)
ifeq ($(UNAME), Darwin) # Mac
	LDFLAGS := -framework GLUT -framework OpenGL -L/usr/local/opt/opencv3/lib -lopencv_videoio -lopencv_core -lpartio -lz -L/usr/local/lib
	CXXFLAGS := -Wall -g -Wno-sign-compare -I/usr/include/eigen3 -I/usr/local/include -I/usr/local/include/eigen3 -O3 -I/usr/local/opt/opencv3/include `pkg-config --cflags opencv` 
	CXX := c++ -std=c++11

else # Linux?
	LDFLAGS := -lglut -lGLU -lGL -fopenmp `pkg-config --libs opencv`
	CXXFLAGS := -Wall -g -Wno-sign-compare -I/usr/include/eigen3 -I/opt/local/include/eigen3 -fopenmp -O3 `pkg-config --cflags opencv`
	CXX := g++ -std=c++11
endif

OBJ := build/world.o build/jsoncpp.o

mpm: $(OBJ) build/main.o
	$(CXX) $^ -o $@ $(LDFLAGS)

mpm_graphics: $(OBJ) build/main_graphics.o
	$(CXX) $^ -o $@ $(LDFLAGS)

build/%.o: src/%.cpp
	$(CXX) -c $(CXXFLAGS) $< -o $@

# Nicked from http://www.gnu.org/software/make/manual/make.html#Automatic-Prerequisites
%.d: src/%.cpp
	@set -e; rm -f $@; \
	$(CXX) -MM $(CXXFLAGS) $< -o $@.tmp; \
	sed 's,\($*\)\.o[ :]*,\1.o: ,g' < $@.tmp > $@; \
	rm -f $@.tmp

-include $(OBJ:.o=.d)

clean:
	rm -f build/*.o build/*.d mpm
