OBJS := \
  $(DIR)/decayVll.o \
  $(DIR)/abbr0_s.o \
  $(DIR)/abbr1_s.o \
  $(DIR)/abbr0_angle.o \
  $(DIR)/abbr1_angle.o \
  $(DIR)/abbr0_hel.o \
  $(DIR)/abbr1_hel.o \
  $(DIR)/born.o \
  $(DIR)/self.o \
  $(DIR)/vert1.o \
  $(DIR)/vert2.o \
  $(DIR)/vert3.o \
  $(DIR)/vert3_Gen5.o \
  $(DIR)/box3.o \
  $(DIR)/box4.o \
  $(DIR)/box5.o \
  $(DIR)/SquaredME.o

$(LIB): $(LIB)($(OBJS))

$(LIB)($(OBJS)): $(DIR)/vars.h $(DECL_H)

LIBS += $(LIB)

