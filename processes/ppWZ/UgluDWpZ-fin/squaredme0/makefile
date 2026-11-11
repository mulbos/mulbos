OBJS := \
  $(DIR)/decayWll.o \
  $(DIR)/decayVll.o \
  $(DIR)/abbr0_s.o \
  $(DIR)/abbr0_hel.o \
  $(DIR)/born.o \
  $(DIR)/SquaredME.o

$(LIB): $(LIB)($(OBJS))

$(LIB)($(OBJS)): $(DIR)/vars.h $(DECL_H)

LIBS += $(LIB)

