ROOT = ../..
OBJ_DIR = $(ROOT)/obj/h2lib

include $(ROOT)/build.variables

$(shell test -d $(OBJ_DIR) || mkdir $(OBJ_DIR))

OBJS = $(OBJ_DIR)/amatrix.o\
	   $(OBJ_DIR)/avector.o\
	   $(OBJ_DIR)/block.o\
	   $(OBJ_DIR)/cluster.o\
	   $(OBJ_DIR)/clustergeometry.o\
	   $(OBJ_DIR)/eigensolvers.o\
	   $(OBJ_DIR)/factorizations.o\
	   $(OBJ_DIR)/h2basic.o\
	   $(OBJ_DIR)/h2settings.o\
	   $(OBJ_DIR)/harith.o\
	   $(OBJ_DIR)/hmatrix.o\
	   $(OBJ_DIR)/matrixplot.o\
	   $(OBJ_DIR)/realavector.o\
	   $(OBJ_DIR)/rkmatrix.o\
	   $(OBJ_DIR)/truncation.o

all: $(OBJS)

$(OBJ_DIR)/%.o: %.cpp
	$(CC) $(FLAGS) $(OPT_FLAGS) -o $@ -c $<

clean:
	-rm -rf $(OBJ_DIR) *~ .*~
