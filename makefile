ROOT= .
SRC= src

all:
	(test -d $(ROOT)/lib || mkdir -p $(ROOT)/lib;)
	(test -d $(ROOT)/obj || mkdir -p $(ROOT)/obj;)
	$(foreach i, $(SRC), make -C $(i);)

clean:
	$(foreach i, $(SRC), make -C $(i) clean;)
	-rm -rf $(ROOT)/lib $(ROOT)/obj fcems *~ .*~
