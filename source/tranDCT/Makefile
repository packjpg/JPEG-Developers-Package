# Project: tranDCT
# Makefile created by Matthias Stirner 28.03.2007
# Working with GCC v3.4.5

CC      = gcc
CPP     = g++
RC      = windres
CFLAGS  = -I. -O3 -Wall -pedantic -funroll-loops -ffast-math -fsched-spec-load -fomit-frame-pointer
LDFLAGS = -s
DEPS    = 
OBJ     = trandct.o
RES     = 
BIN     = tranDCT

%.o: %.cpp $(DEPS)
	$(CPP) -c -o $@ $< $(CFLAGS)

$(BIN): $(OBJ) $(RES)
	$(CPP) -o $@ $^ $(LDFLAGS)
