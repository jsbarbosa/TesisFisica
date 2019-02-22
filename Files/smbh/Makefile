export OPENGL = 0
#SRC_DIR = ../../../../rebound/src
SRC_DIR = src
FILE_NAME = core.c

include $(SRC_DIR)/Makefile.defs

rebound.so : $(FILE_NAME) $(SRC_DIR)/librebound.so headers.h
	rm -f librebound.so
	ln -s $(SRC_DIR)/librebound.so .
	$(CC) -I $(SRC_DIR) -I . -Wl,-rpath,./ $(OPT) $(PREDEF) $(FILE_NAME) -L. -lrebound $(LIB) -fPIC -shared -o rebound.so

$(SRC_DIR)/librebound.so : $(SCR_DIR)/*
	$(MAKE) -C $(SRC_DIR)

clean:
	rm -f *.so *.run *.dat
	$(MAKE) -C $(SRC_DIR) clean
