OBJS := \
  $(DIR)/decayWll.o \
  $(DIR)/decayWmll.o \
  $(DIR)/abbr1_s.o \
  $(DIR)/abbr1_angle.o \
  $(DIR)/abbr1_hel1.o \
  $(DIR)/abbr1_hel2.o \
  $(DIR)/vertH.o \
  $(DIR)/vertZ_Gen5.o \
  $(DIR)/box_Gen5.o \
  $(DIR)/SquaredME.o

$(LIB): $(LIB)($(OBJS))

$(LIB)($(OBJS)): $(DIR)/vars.h $(DECL_H)

LIBS += $(LIB)

