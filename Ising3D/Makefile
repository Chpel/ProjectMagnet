
CC=gcc
CXX=g++
RM=rm -f


CPPFLAGS=-std=c++14 -O3
LDFLAGS=-std=c++14 -O3
LDLIBS=

PROGRAM=mcmc

SRCS=mcmc3d.cpp \
     main.cpp


OBJS=$(subst .cpp,.o,$(SRCS))


all: $(PROGRAM)

mcmc.o: mcmced.cpp mcmc.h observable.h

main.o: main.cpp mcmc.h


$(PROGRAM): $(OBJS)
	$(CXX) $(LDFLAGS) -o $(PROGRAM) $(OBJS) $(LDLIBS)

clean:
	rm -rf $(OBJS) $(PROGRAM)
