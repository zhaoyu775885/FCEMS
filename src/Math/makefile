ROOT = ../..
OBJ_DIR = $(ROOT)/obj/Math

include $(ROOT)/build.variables

$(shell test -d $(OBJ_DIR) || mkdir $(OBJ_DIR))

OBJS = $(OBJ_DIR)/GaussQuad.o\
	   $(OBJ_DIR)/Geom.o\
	   $(OBJ_DIR)/Point3d.o\
	   $(OBJ_DIR)/SingularIntegral.o\
	   $(OBJ_DIR)/Vector.o\

all: $(OBJS)

$(OBJ_DIR)/%.o: %.cpp
	$(CC) $(FLAGS) $(OPT_FLAGS) -o $@ -c $<

clean:
	-rm -rf $(OBJ_DIR) *~ .*~
