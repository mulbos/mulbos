# This file is adapted from the makefile of FormCalc

SRC = .
FC = ifx

# to control how libraries are built
#ARFLAGS = rUv

ifeq ($(FC),ifx)
FFLAGS = -extend-source
else ifeq ($(FC),gfortran)
# debug segmentation fault:
#FFLAGS = -ffixed-line-length-none -g -fcheck=all -Wall -O0
# debug NAN or Infinity:
#FFLAGS = -ffixed-line-length-none -g -fbacktrace -ffpe-trap=zero,overflow,underflow,invalid
# default:
FFLAGS = -ffixed-line-length-none
else
$(error Unknown FC '$(FC)'; use ifx or gfortran)
endif

CC = gcc
CFLAGS = -O2 -g  -DUNDERSCORE -m64

.PHONY: install clean

FOLDER := util
LIBINSTALL := libutil.a

MAKECMDGOALS ?= $(FOLDER).a install

all: $(MAKECMDGOALS)

VPATH := .
INCLUDE := $(patsubst %,-I%,$(subst :, ,$(VPATH)))
FFLAGS += $(INCLUDE)

export FC
export FFLAGS
export CC
export CFLAGS

DIR = $(SRC)/$(FOLDER)
LIB = $(FOLDER).a

-include $(DIR)/makefile

install:
	mv $(FOLDER).a $(LIBINSTALL)

clean:
	$(RM) $(LIBINSTALL)
