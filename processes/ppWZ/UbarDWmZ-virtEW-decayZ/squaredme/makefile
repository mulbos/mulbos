OBJS := \
  $(DIR)/decayWmll_NLO.o \
  $(DIR)/decayZll_NLO.o \
  $(DIR)/abbr0_s.o \
  $(DIR)/abbr1_s.o \
  $(DIR)/abbr0_hel.o \
  $(DIR)/born.o \
  $(DIR)/SquaredME.o

$(LIB): $(LIB)($(OBJS))

$(LIB)($(OBJS)): $(DIR)/vars.h $(DECL_H)

LIBS += $(LIB)

