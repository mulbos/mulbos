kin := $(DIR)/kin

.SUFFIXES: .F .F90 .o

.F.o:
	$(FC) $(FFLAGS) -c $< -o $@

.F90.o:
	$(FC) $(FFLAGS) -c $< -o $@

OBJS := \
  $(kin)/VecSet.o \
  $(kin)/SInvariant.o \
  $(kin)/TInvariant.o \
  $(kin)/basic_funcs_FC6.o
#  $(kin)/functions_FC9.o \

$(LIB): $(LIB)($(OBJS))

LIBS += $(LIB)

