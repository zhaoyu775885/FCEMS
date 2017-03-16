ROOT=../..
OBJ_DIR=$(ROOT)/obj/Basics

include $(ROOT)/build.variables

$(shell test -d $(OBJ_DIR) || mkdir -p $(OBJ_DIR))

OBJS = $(OBJ_DIR)/routine.o

all: $(OBJS)

$(OBJ_DIR)/%.o: %.cpp 
	$(CC) $(FLAGS) $(OPT_FLAGS) -o $@ -c $<

clean:
	-rm -rf $(OBJ_DIR) *~ .*~
