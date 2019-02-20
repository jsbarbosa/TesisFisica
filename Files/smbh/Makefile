export OPENGL = 0
#SRC_DIR = ../../../../rebound/src
SRC_DIR = src/
FILE_NAME = core.c

include $(SRC_DIR)/Makefile.defs

rebound.so : $(FILE_NAME)
	@echo "Compiling problem file ..."
	rm -f librebound.so
	ln -s $(SRC_DIR)/librebound.so .
	$(CC) -I $(SRC_DIR) -I . -Wl,-rpath,./ $(OPT) $(PREDEF) $(FILE_NAME) -L. -lrebound $(LIB) -fPIC -shared -o rebound.so
	@echo ""
	@echo "REBOUND compiled successfully."

clean:
	@echo "Cleaning up shared library librebound.so ..."
	rm -f *.so *.run *.dat
