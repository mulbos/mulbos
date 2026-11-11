# This file is adapted from the makefile of FormCalc-6.0 (https://feynarts.de/formcalc/)

PROGNAME = mulbos

# dont change these paths here; change it in run.py
HOME2 = /home/ldninh/projects/MulBos/MulBos_v100
PROC = processes/ppZZ/UbarUZZ-virtEW-decayZ2
LHAPDF_DIR = /home/ldninh/codes_hep/lib/lhapdf-6.5.6

SRC = .

LT = $(HOME2)/tools/loops/LoopInts/build
LoopLib = oopints

BS = $(HOME2)/tools/integrators/basesv5.1
UT = $(HOME2)/tools/utilities
FCT = $(HOME2)/tools/FormCalc_Tools
AM1 = $(HOME2)/$(PROC)
FC = ifx
 
LHAPDF_libX = lib
LHAPDF_LDFLAGS = -L$(LHAPDF_DIR)/$(LHAPDF_libX) -Wl,-rpath,$(LHAPDF_DIR)/$(LHAPDF_libX)
LHAPDF_LIBS = -lLHAPDF

STDLIBS = -L$(AM1) -lsquaredme -lsquaredme0 -lsquaredme1 -L$(FCT) -lutil -L$(UT) -lutilities -I$(UT) -L$(BS) -lbases $(LHAPDF_LDFLAGS) $(LHAPDF_LIBS)
STDLIBS += -L$(LT) -l$(LoopLib)

ifeq ($(FC),ifx)
FFLAGS = -extend-source
else ifeq ($(FC),gfortran)
# debugging:
#1: FFLAGS = -ffixed-line-length-none -g -fcheck=all -Wall -O0 -fcheck=bounds
#2: ./myprog (run the program)
#3: Use gdb to locate the crash:
#3a: gdb ./myprog
#3b: (gdb) run
#3c: (gdb) backtrace
# You’ll see where it crashed and on which line.
#FFLAGS = -ffixed-line-length-none -g -fcheck=all -Wall -O0 -fcheck=bounds
#FFLAGS = -ffixed-line-length-none -g -fcheck=all -Wall -O0
# default:
FFLAGS = -ffixed-line-length-none
else
$(error Unknown FC '$(FC)'; use ifx or gfortran)
endif

#CC = gcc
#CFLAGS = -O2 -g  -DUNDERSCORE -m64

.PHONY: clean

#VPATH := $(SRC):$(SRC)/renconst:$(LT)/../include
VPATH := $(SRC):$(UT):$(LT)/../include

INCLUDE := $(patsubst %,-I%,$(subst :, ,$(VPATH)))

PREFIX :=

FFLAGS += $(INCLUDE)

SRCF = main.F

$(PROGNAME): $(SRCF)
	$(FC) $(FFLAGS) -o $@ $(SRCF) $(STDLIBS)
	-$(RM) $*.o
clean:
	rm -f *.o $(PROGNAME)
