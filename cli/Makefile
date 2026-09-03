CC      ?= cc
CFLAGS  ?= -O2 -std=c99 -Wall -Wextra -pedantic
LDLIBS  ?= -lz
PREFIX  ?= /usr/local

physmerge: physmerge.c
	$(CC) $(CFLAGS) -o $@ physmerge.c $(LDLIBS)

# build without zlib (no .gz input support)
nozlib: physmerge.c
	$(CC) $(CFLAGS) -DPHYSMERGE_NO_ZLIB -o physmerge physmerge.c

install: physmerge
	install -d $(DESTDIR)$(PREFIX)/bin
	install -m 755 physmerge $(DESTDIR)$(PREFIX)/bin/physmerge

clean:
	rm -f physmerge

.PHONY: install clean nozlib
