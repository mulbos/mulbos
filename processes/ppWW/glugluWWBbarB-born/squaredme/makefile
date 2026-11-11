OBJS := \
  $(DIR)/decayWll.o \
  $(DIR)/decayWmll.o \
  $(DIR)/abbr0_s.o \
  $(DIR)/born.o \
  $(DIR)/abbr0_angle.o \
  $(DIR)/abbr0_hel1.o \
  $(DIR)/abbr0_hel2.o \
  $(DIR)/abbr0_hel3.o \
  $(DIR)/SquaredME.o

$(LIB): $(LIB)($(OBJS))

$(LIB)($(OBJS)): $(DIR)/vars.h $(DECL_H)

LIBS += $(LIB)

