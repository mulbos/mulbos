# This file is adapted from the makefile of FormCalc
# User input path to the program main dir:
HOME1 = ../../..
LT = $(HOME1)/tools/loops/LoopInts/include
UT = $(HOME1)/tools/utilities
FCUT = $(HOME1)/tools/FormCalc_Tools

SRC = .
FC = ifx

# to control how libraries are built
ARFLAGS = rUv

ifeq ($(FC),ifx)
FFLAGS = -extend-source
else ifeq ($(FC),gfortran)
# debug segmentation fault:
#FFLAGS = -ffixed-line-length-none -g -fcheck=all -Wall -O0
# debug NAN or Infinity:
#FFLAGS = -ffixed-line-length-none -g -fbacktrace -ffpe-trap=zero,overflow,invalid
# default:
FFLAGS = -ffixed-line-length-none
else
$(error Unknown FC '$(FC)'; use ifx or gfortran)
endif

CC = gcc
CFLAGS = -O2 -g  -DUNDERSCORE -m64

.PHONY: install clean

FOLDER1 := squaredme
LIBINSTALL1 := libsquaredme.a

FOLDER2 := squaredme0
LIBINSTALL2 := libsquaredme0.a

FOLDER3 := squaredme1
LIBINSTALL3 := libsquaredme1.a

MAKECMDGOALS ?= $(FOLDER1).a $(FOLDER2).a $(FOLDER3).a install

all: $(MAKECMDGOALS)

VPATH := $(UT):$(FCUT):$(LT)
INCLUDE := $(patsubst %,-I%,$(subst :, ,$(VPATH)))
FFLAGS += $(INCLUDE)

export FC
export FFLAGS
export CC
export CFLAGS

# generic def:
DIR = $(SRC)/$(FOLDER)
LIB = $(FOLDER).a

LIBS :=

FOLDER := squaredme
-include $(DIR)/makefile

FOLDER := squaredme0
-include $(DIR)/makefile

FOLDER := squaredme1
-include $(DIR)/makefile

install:
	mv $(FOLDER1).a $(LIBINSTALL1)
	mv $(FOLDER2).a $(LIBINSTALL2)
	mv $(FOLDER3).a $(LIBINSTALL3)

clean:
	$(RM) *.a
