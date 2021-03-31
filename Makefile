PROG = MD
#PROG1 = RDF

SRC = main.cpp system.cpp timecalc.cpp collision.cpp rdf.cpp writefiles.cpp readinput.cpp readfiles.cpp
#SRC1 = rdf_fin.cpp

OBJS = ${SRC:.cpp=.o}
#OBJS1 = ${SRC1:.cpp=.o}

CXX = g++ -std=c++11 -debug -O0 -Wl,-rpath /usr/local/lib/ 
CXXFLAGS= -lgsl -lgslcblas
LIBS= 

#CXXFLAGS=-O3 -funroll-loops -DNDEBUG 

all: $(PROG) $(PROG1)

$(PROG):  $(OBJS)
	 $(CXX) $(LIBS)  $(CXXFLAGS) -o $@ $^

#$(PROG1):  $(OBJS1)
#	 $(CXX) $(LIBS)  $(CXXFLAGS) -o $@ $^


%.o:  %.cpp
	$(CXX) $(LIBS) -c -o $@ $<

clean: 
	rm -rf *.o 

distclean:
	rm -f $(PROG) *.o *.debug
