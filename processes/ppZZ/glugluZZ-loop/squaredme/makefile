OBJS := \
  $(DIR)/decayVll.o \
  $(DIR)/abbr1_s.o \
  $(DIR)/abbr1_angle.o \
  $(DIR)/abbr1_hel.o \
  $(DIR)/vert.o \
  $(DIR)/box_Gen5.o \
  $(DIR)/SquaredME.o

$(LIB): $(LIB)($(OBJS))

$(LIB)($(OBJS)): $(DIR)/vars.h $(DECL_H)

LIBS += $(LIB)

