CFLAGS = -O3 -march=native -fopenmp -Wall -Wextra
LDLIBS = -lm
all: libmindist.so qtsearch
libmindist.so: mindist.c
	gcc $(CFLAGS) -shared -fPIC -o $@ $<
qtsearch: qtsearch.c mindist.c
	gcc $(CFLAGS) -o $@ qtsearch.c mindist.c $(LDLIBS)
