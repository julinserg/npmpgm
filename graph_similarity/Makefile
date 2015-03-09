#########################################
# why isn't main.o included in .depend? #
#########################################

SRCS=main.cc graph_gen_alg.cc chunglu_gen.cc erdosrenyi_gen.cc util_fns.cc
OBJECTS=$(SRCS:.cc=.o)

CXX = g++

CXXFLAGS = -g -Wall -I ~/workspace/newton_gmres/ -L ~/build/lib -laes -Wno-sign-compare -std=c++0x -O3

all: graph_embedding

%.o: %.c
	$(CXX) -c $<  $(CXXFLAGS)

graph_embedding: $(OBJECTS)
	$(CXX) -o $@ $^  $(CXXFLAGS)

depend: .depend

.depend: $(SRCS)
	rm -f ./.depend
	$(CXX) -MM -MT $^ $(CXXFLAGS) > ./.depend

clean:
	$(RM) *.o 

include .depend
