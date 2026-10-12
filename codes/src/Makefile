CFLAGS = -O3 -march=native -fopenmp -Wall -Wextra
LDLIBS = -lm
all: libmindist.so qtsearch qtclimb qtx qtX
libmindist.so: mindist.c qmd.c
	gcc $(CFLAGS) -shared -fPIC -o $@ mindist.c qmd.c
qtsearch: qtsearch.c mindist.c
	gcc $(CFLAGS) -o $@ qtsearch.c mindist.c $(LDLIBS)
qtclimb: qtclimb.c mindist.c
	gcc $(CFLAGS) -o $@ qtclimb.c mindist.c $(LDLIBS)
qtx: qtx.c mindist.c
	gcc $(CFLAGS) -o $@ qtx.c mindist.c $(LDLIBS)
qtX: qtX.c mindist.c
	gcc $(CFLAGS) -o $@ qtX.c mindist.c $(LDLIBS)
