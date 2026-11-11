# This file is adapted from the makefile of LoopTools
# Last modification: Duc Ninh LE (Aug 2026)

SRC = .
PREFIX = x86_64-linux

EXE = 
DEF = -D
UNDERSCORE = 1

FC = ifx
FFLAGS = -O2 -g -extend-source -warn truncated_source -assume bscc -fpe0   $(DEF)QUAD -DU77EXT=1
#FFLAGS = -pg -extend-source -warn truncated_source -assume bscc -fpe0   $(DEF)QUAD -DU77EXT=1
FFLAGS_QP = -r16 -DQUADPREC -DDBLE=QEXT -DDIMAG=QIMAG -DDCONJG=QCONJG -DDCMPLX=QCMPLX

CC = gcc
CFLAGS = -O2 -g -fomit-frame-pointer -ffast-math -Wall  -m64 -DUNDERSCORE=$(UNDERSCORE)
#CFLAGS = -pg -Wall  -m64 -DUNDERSCORE=$(UNDERSCORE)

SUFFIX =
BLD = build$(SUFFIX)

LIBDIR = $(PREFIX)/lib
INCLUDEDIR = $(PREFIX)/include
BINDIR = $(PREFIX)/bin

LIB = $(BLD)/liboopints$(SUFFIX).a

INCLUDE = $(SRC)/include/loopints.h

default: $(TS) $(LIB)

all: default

install: default
	-mkdir $(PREFIX)
	-mkdir $(LIBDIR) $(BINDIR) $(INCLUDEDIR)
	cp -p $(LIB) $(LIBDIR)
	cp -p $(INCLUDE) $(INCLUDEDIR)


OBJS = $(QOBJS) \
  $(BLD)/A0.o $(BLD)/A0C.o \
  $(BLD)/A0s.o $(BLD)/A0sQP.o \
  $(BLD)/A0QP.o $(BLD)/A0CQP.o \
  $(BLD)/Bget.o $(BLD)/BgetC.o \
  $(BLD)/Bcoeffs.o $(BLD)/BcoeffsQP.o \
  $(BLD)/BgetQP.o $(BLD)/BgetCQP.o \
  $(BLD)/Cget.o $(BLD)/CgetC.o \
  $(BLD)/Ccoeffs.o $(BLD)/CcoeffsQP.o \
  $(BLD)/CgetQP.o $(BLD)/CgetCQP.o \
  $(BLD)/Dget.o $(BLD)/DgetC.o \
  $(BLD)/Dcoeffs.o $(BLD)/DcoeffsQP.o \
  $(BLD)/D0div.o $(BLD)/D0divQP.o \
  $(BLD)/DgetQP.o $(BLD)/DgetCQP.o \
  $(BLD)/cache.o $(BLD)/cacheQP.o \
  $(BLD)/func_basic.o $(BLD)/func_basicQP.o \
  $(BLD)/default.o $(BLD)/defaultQP.o \
  $(BLD)/fTn.o $(BLD)/fTnQP.o \
  $(BLD)/fS3n.o $(BLD)/fS3nQP.o \
  $(BLD)/solve-LU.o $(BLD)/solve-LUQP.o \
  $(BLD)/ini.o $(BLD)/iniQP.o \
  $(BLD)/Dump.o $(BLD)/DumpQP.o \
  $(BLD)/ffinit.o

TS = $(BLD)/timestamp

FFINC = $(SRC)/include/ff.h $(TS)

LTINC = $(SRC)/include/defs.h $(SRC)/include/lt.h $(FFINC)

# important: how to compile
XFC = $(FC) $(FFLAGS) $(FFLAGS_QUAD) -I$(SRC)/include 
CFC = $(XFC) $(DEF)COMPLEXPARA
XFCQP = $(XFC) $(FFLAGS_QP)
CFCQP = $(CFC) $(FFLAGS_QP)

$(BLD)/A0.o: $(SRC)/A/A0.F $(LTINC)
	$(XFC) -c -o $(BLD)/A0.o $(SRC)/A/A0.F
$(BLD)/A0C.o: $(SRC)/A/A0.F $(LTINC)
	$(CFC) -c -o $(BLD)/A0C.o $(SRC)/A/A0.F
$(BLD)/A0QP.o: $(SRC)/A/A0.F $(LTINC)
	$(XFCQP) -c -o $(BLD)/A0QP.o $(SRC)/A/A0.F
$(BLD)/A0CQP.o: $(SRC)/A/A0.F $(LTINC)
	$(CFCQP) -c -o $(BLD)/A0CQP.o $(SRC)/A/A0.F
$(BLD)/A0s.o: $(SRC)/A/A0s.F $(LTINC)
	$(XFC) -c -o $(BLD)/A0s.o $(SRC)/A/A0s.F
$(BLD)/A0sQP.o: $(SRC)/A/A0s.F $(LTINC)
	$(XFCQP) -c -o $(BLD)/A0sQP.o $(SRC)/A/A0s.F

$(BLD)/Bget.o: $(SRC)/B/Bget.F $(LTINC)
	$(XFC) -c -o $(BLD)/Bget.o $(SRC)/B/Bget.F
$(BLD)/BgetC.o: $(SRC)/B/Bget.F $(LTINC)
	$(CFC) -c -o $(BLD)/BgetC.o $(SRC)/B/Bget.F
$(BLD)/BgetQP.o: $(SRC)/B/Bget.F $(LTINC)
	$(XFCQP) -c -o $(BLD)/BgetQP.o $(SRC)/B/Bget.F
$(BLD)/BgetCQP.o: $(SRC)/B/Bget.F $(LTINC)
	$(CFCQP) -c -o $(BLD)/BgetCQP.o $(SRC)/B/Bget.F

$(BLD)/Bcoeffs.o: $(SRC)/B/Bcoeffs.F $(LTINC)
	$(CFC) -c -o $(BLD)/Bcoeffs.o $(SRC)/B/Bcoeffs.F
$(BLD)/BcoeffsQP.o: $(SRC)/B/Bcoeffs.F $(LTINC)
	$(CFCQP) -c -o $(BLD)/BcoeffsQP.o $(SRC)/B/Bcoeffs.F

$(BLD)/Cget.o: $(SRC)/C/Cget.F $(LTINC)
	$(XFC) -c -o $(BLD)/Cget.o $(SRC)/C/Cget.F
$(BLD)/CgetC.o: $(SRC)/C/Cget.F $(LTINC)
	$(CFC) -c -o $(BLD)/CgetC.o $(SRC)/C/Cget.F
$(BLD)/CgetQP.o: $(SRC)/C/Cget.F $(LTINC)
	$(XFCQP) -c -o $(BLD)/CgetQP.o $(SRC)/C/Cget.F
$(BLD)/CgetCQP.o: $(SRC)/C/Cget.F $(LTINC)
	$(CFCQP) -c -o $(BLD)/CgetCQP.o $(SRC)/C/Cget.F

$(BLD)/Ccoeffs.o: $(SRC)/C/Ccoeffs.F $(LTINC)
	$(CFC) -c -o $(BLD)/Ccoeffs.o $(SRC)/C/Ccoeffs.F
$(BLD)/CcoeffsQP.o: $(SRC)/C/Ccoeffs.F $(LTINC)
	$(CFCQP) -c -o $(BLD)/CcoeffsQP.o $(SRC)/C/Ccoeffs.F

$(BLD)/Dget.o: $(SRC)/D/Dget.F $(LTINC)
	$(XFC) -c -o $(BLD)/Dget.o $(SRC)/D/Dget.F
$(BLD)/DgetC.o: $(SRC)/D/Dget.F $(LTINC)
	$(CFC) -c -o $(BLD)/DgetC.o $(SRC)/D/Dget.F
$(BLD)/DgetQP.o: $(SRC)/D/Dget.F $(LTINC)
	$(XFCQP) -c -o $(BLD)/DgetQP.o $(SRC)/D/Dget.F
$(BLD)/DgetCQP.o: $(SRC)/D/Dget.F $(LTINC)
	$(CFCQP) -c -o $(BLD)/DgetCQP.o $(SRC)/D/Dget.F

$(BLD)/Dcoeffs.o: $(SRC)/D/Dcoeffs.F $(LTINC)
	$(CFC) -c -o $(BLD)/Dcoeffs.o $(SRC)/D/Dcoeffs.F
$(BLD)/DcoeffsQP.o: $(SRC)/D/Dcoeffs.F $(LTINC)
	$(CFCQP) -c -o $(BLD)/DcoeffsQP.o $(SRC)/D/Dcoeffs.F

$(BLD)/D0div.o: $(SRC)/D/D0div.F $(LTINC)
	$(CFC) -c -o $(BLD)/D0div.o $(SRC)/D/D0div.F
$(BLD)/D0divQP.o: $(SRC)/D/D0div.F $(LTINC)
	$(CFCQP) -c -o $(BLD)/D0divQP.o $(SRC)/D/D0div.F

$(BLD)/Dump.o: $(SRC)/util/Dump.F $(LTINC)
	$(CFC) -c -o $(BLD)/Dump.o $(SRC)/util/Dump.F
$(BLD)/DumpQP.o: $(SRC)/util/Dump.F $(LTINC)
	$(CFCQP) -c -o $(BLD)/DumpQP.o $(SRC)/util/Dump.F

$(BLD)/ffinit.o: $(SRC)/util/ffinit.F $(LTINC)
	$(XFC) -c -o $(BLD)/ffinit.o $(SRC)/util/ffinit.F

$(BLD)/func_basic.o: $(SRC)/util/func_basic.F $(LTINC)
	$(CFC) -c -o $(BLD)/func_basic.o $(SRC)/util/func_basic.F
$(BLD)/func_basicQP.o: $(SRC)/util/func_basic.F $(LTINC)
	$(CFCQP) -c -o $(BLD)/func_basicQP.o $(SRC)/util/func_basic.F

$(BLD)/default.o: $(SRC)/util/default.F $(LTINC)
	$(XFC) -c -o $(BLD)/default.o $(SRC)/util/default.F
$(BLD)/defaultQP.o: $(SRC)/util/default.F $(LTINC)
	$(XFCQP) -c -o $(BLD)/defaultQP.o $(SRC)/util/default.F

$(BLD)/ini.o: $(SRC)/util/ini.F $(LTINC)
	$(XFC) -c -o $(BLD)/ini.o $(SRC)/util/ini.F
$(BLD)/iniQP.o: $(SRC)/util/ini.F $(LTINC)
	$(XFCQP) -c -o $(BLD)/iniQP.o $(SRC)/util/ini.F

$(BLD)/fTn.o: $(SRC)/util/fTn.F $(LTINC)
	$(CFC) -c -o $(BLD)/fTn.o $(SRC)/util/fTn.F
$(BLD)/fTnQP.o: $(SRC)/util/fTn.F $(LTINC)
	$(CFCQP) -c -o $(BLD)/fTnQP.o $(SRC)/util/fTn.F

$(BLD)/fS3n.o: $(SRC)/util/fS3n.F $(LTINC)
	$(CFC) -c -o $(BLD)/fS3n.o $(SRC)/util/fS3n.F
$(BLD)/fS3nQP.o: $(SRC)/util/fS3n.F $(LTINC)
	$(CFCQP) -c -o $(BLD)/fS3nQP.o $(SRC)/util/fS3n.F

$(BLD)/solve-LU.o: $(SRC)/util/solve-LU.F $(LTINC)
	$(CFC) -c -o $(BLD)/solve-LU.o $(SRC)/util/solve-LU.F
$(BLD)/solve-LUQP.o: $(SRC)/util/solve-LU.F $(LTINC)
	$(CFCQP) -c -o $(BLD)/solve-LUQP.o $(SRC)/util/solve-LU.F

$(BLD)/cache.o: $(SRC)/util/cache.c $(LTINC)
	$(CC) $(CFLAGS) -c -o $(BLD)/cache.o $(SRC)/util/cache.c
$(BLD)/cacheQP.o: $(SRC)/util/cacheQP.c $(LTINC)
	$(CC) $(CFLAGS) -c -o $(BLD)/cacheQP.o $(SRC)/util/cacheQP.c

$(LIB): $(OBJS)
	ar cUru $(LIB) $?
	-ranlib $(LIB)

$(TS):
	-mkdir $(BLD)
	touch $(TS)

clean:
	rm -fr $(BLD)

