OBJS := \
  $(DIR)/decayWmll.o \
  $(DIR)/decayVll.o \
  $(DIR)/abbr0_s.o \
  $(DIR)/abbr1_s.o \
  $(DIR)/abbr1_angle.o \
  $(DIR)/abbr0_hel.o \
  $(DIR)/abbr1_hel.o \
  $(DIR)/born.o \
  $(DIR)/self.o \
  $(DIR)/vert.o \
  $(DIR)/box.o \
  $(DIR)/SquaredME.o

$(LIB): $(LIB)($(OBJS))

$(LIB)($(OBJS)): $(DIR)/vars.h $(DECL_H)

LIBS += $(LIB)

