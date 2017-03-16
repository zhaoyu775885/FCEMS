ROOT = ..
OBJ_DIR = $(ROOT)/obj
include $(ROOT)/build.variables

PROG = fcems
DIRS = Basics\
	   h2lib\
	   Math\
	   Mesh

LIBS = $(MKL_LIBS)
LIBS += $(OPENGL_LIBS)

$(shell test -d $(OBJ_DIR) || mkdir -p $(OBJ_DIR))

OBJS = $(OBJ_DIR)/BasicBEM.o\
	   $(OBJ_DIR)/BasisFunction.o\
	   $(OBJ_DIR)/FastBEM.o\
	   $(OBJ_DIR)/Green.o\
	   $(OBJ_DIR)/Intg.o\
	   $(OBJ_DIR)/LinAlg.o\
	   $(OBJ_DIR)/main.o\
	   $(OBJ_DIR)/PostProcess.o

OBJS_BASICS = $(ROOT)/obj/Basics/*.o
OBJS_H2LIB = $(ROOT)/obj/h2lib/*.o
OBJS_MATH = $(ROOT)/obj/Math/*.o
OBJS_MESH = $(ROOT)/obj/Mesh/*.o

all: $(ROOT)/$(PROG)

$(ROOT)/$(PROG): $(OBJS)
	$(foreach i, $(DIRS), make -C $(i);)
	$(CC) $(CFLAGS) $(OPT_FLAGS) $(OBJS) $(OBJS_H2LIB) $(OBJS_BASICS) $(OBJS_MATH) $(OBJS_MESH) -o $(ROOT)/$(PROG) $(LIBS)

$(OBJ_DIR)/%.o:%.cpp
	$(CC) $(FLAGS) $(OPT_FLAGS) -o $@ -c $< 

clean:
	$(foreach i, $(DIRS), make -C $(i) clean;)
	-rm -f $(OBJS) $(ROOT)/$(PROG) *~ .*~
