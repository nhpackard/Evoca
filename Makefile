CC     = gcc
CFLAGS = -O2 -Wall -fPIC

UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
    SHARED_FLAG = -dynamiclib
    LIB         = C/libevoca.dylib
else
    SHARED_FLAG = -shared
    LIB         = C/libevoca.so
endif

.PHONY: all clean

all: $(LIB)

$(LIB): C/evoca.c C/evoca.h
	$(CC) $(CFLAGS) $(SHARED_FLAG) -o $@ C/evoca.c

clean:
	rm -f C/libevoca.so C/libevoca.dylib
