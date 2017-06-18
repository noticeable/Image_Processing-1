CFLAGS=-c -g -std=gnu99
LDFLAGS=-lm

essai1: essai1.o pgm.o fft.o
	clang $(LDFLAGS)  $^  -o $@

essai2: essai2.o pgm.o fft.o
	gcc $(LDFLAGS)  $^  -o $@

part1: part1.o pgm.o fft.o divers.o
	clang $(LDFLAGS) $^  -o $@

contours: contoursgradient.o pgm.o fft.o divers.o
	clang $(LDFLAGS) $^  -o $@

contourslaplacien: contourslaplacien.o pgm.o fft.o divers.o
	clang $(LDFLAGS) $^  -o $@

%.o:%.c
	clang $(CFLAGS) $<

.PHONY: clean

clean:
	rm -rf *.o part1 essai2 essai1 contours contourslaplacien tp2 res.pgm
