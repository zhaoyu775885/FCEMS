ROOT = ../..
OBJ_DIR = $(ROOT)/obj/Mesh

include $(ROOT)/build.variables

$(shell test -d $(OBJ_DIR) || mkdir $(OBJ_DIR))

OBJS = $(OBJ_DIR)/Mesh.o\
	   $(OBJ_DIR)/Patch.o

all: $(OBJS)

$(OBJ_DIR)/%.o: %.cpp
	$(CC) $(FLAGS) $(OPT_FLAGS) -o $@ -c $<

clean:
	-rm -rf $(OBJ_DIR) *~ .*~
