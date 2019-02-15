export OPENGL = 0
SRC_DIR = ../../../../rebound/src
FILE_NAME = c_simple.c

include $(SRC_DIR)/Makefile.defs

*.dat : plots.py rebound.run
	python3 plots.py

rebound.run : c_simple.c
	@echo "Compiling problem file ..."
	rm -f librebound.so
	ln -s $(SRC_DIR)/librebound.so .
	$(CC) -I $(SRC_DIR) -Wl,-rpath,./ $(OPT) $(PREDEF) $(FILE_NAME) -L. -lrebound $(LIB) -o rebound.run
	@echo ""
	@echo "REBOUND compiled successfully."

clean:
	@echo "Cleaning up shared library librebound.so ..."
	rm -f librebound.so rebound.run *.dat
