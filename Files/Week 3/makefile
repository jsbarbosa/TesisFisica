floating.png : floating.py floating.dat
	python3 floating.py

floating.dat : a.out
	./a.out > floating.dat
a.out : ieee754.c
	gcc ieee754.c -lm -Wall

clean :
	rm a.out floating.dat
