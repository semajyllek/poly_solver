VERSION := $(shell cat VERSION)

# compiler and GMP: GMP from pkg-config when it knows it, else Homebrew's prefix, else the system paths
CC ?= cc
GMP_CFLAGS := $(shell pkg-config --cflags gmp 2>/dev/null || { p=$$(brew --prefix 2>/dev/null) && echo -I$$p/include; })
GMP_LIBS := $(shell pkg-config --libs gmp 2>/dev/null || { p=$$(brew --prefix 2>/dev/null) && echo -L$$p/lib -lgmp || echo -lgmp; })
# CFLAGS and LDFLAGS are the user's (make CFLAGS=-O3); what the build needs is added separately
CFLAGS ?= -Wall -Wextra -g -O2
ALL_CFLAGS = $(CFLAGS) $(GMP_CFLAGS) -MMD -MP
LIBS = $(GMP_LIBS) -lm

# install locations (make install PREFIX=... DESTDIR=...)
PREFIX ?= /usr/local

# library sources: every src/*.c (token_list.c and token_stack.c are unused)
LIB_NAMES = $(filter-out token_list.c token_stack.c,$(notdir $(wildcard src/*.c)))
LIB_OBJS = $(addprefix build/,$(LIB_NAMES:.c=.o))
# library, tests and CLI see the internal headers, with every global renamed to ps_/psi_ by the
# generated src/namespace.h; only the ps_ API is visible outside the shared library
INTERNAL = -Isrc -Iinclude -include src/namespace.h -fvisibility=hidden -fPIC

# test binary (tests/) and CLI (cli/), both at the top level so tools/golden.py and the README run ./polysolver
TEST_OBJS = $(LIB_OBJS) build/test_symbolic_solver.o build/test_properties.o build/test_subres_ref.o
TEST_TARGET = test_symbolic_solver
CLI_TARGET = polysolver

all: $(TEST_TARGET)

$(TEST_TARGET): $(TEST_OBJS)
	$(CC) $(ALL_CFLAGS) $(LDFLAGS) -o $@ $(TEST_OBJS) $(LIBS)

cli: $(CLI_TARGET)

$(CLI_TARGET): $(LIB_OBJS) build/main.o
	$(CC) $(ALL_CFLAGS) $(LDFLAGS) -o $@ $(LIB_OBJS) build/main.o $(LIBS)

build/%.o: src/%.c | build
	$(CC) $(ALL_CFLAGS) $(INTERNAL) -c $< -o $@
build/%.o: tests/%.c | build
	$(CC) $(ALL_CFLAGS) $(INTERNAL) -c $< -o $@
build/main.o: cli/main.c | build
	$(CC) $(ALL_CFLAGS) $(INTERNAL) -c $< -o $@
build:
	mkdir -p build

# the rename map and the flags reach every object
$(TEST_OBJS) build/main.o: src/namespace.h Makefile

# regenerate src/namespace.h after adding or removing a global (symcheck says when): the library
# compiled without it, its globals listed by nm, each mapped to tools/public_api.txt's name or psi_
namespace: | build
	mkdir -p build/raw
	for f in $(LIB_NAMES); do $(CC) $(ALL_CFLAGS) -Isrc -Iinclude -c src/$$f -o build/raw/$${f%.c}.o || exit 1; done
	nm -g build/raw/*.o > build/raw/globals.txt
	python3 tools/gen_namespace.py build/raw/globals.txt > src/namespace.h.new && mv src/namespace.h.new src/namespace.h
	rm -rf build/raw

# static and shared libraries of everything but the entry points
LIB_TARGET = build/libpolysolver.a
SOVERSION = 0
ifeq ($(shell uname -s),Darwin)
SHARED_NAME = libpolysolver.$(SOVERSION).dylib
SHARED_LINK = libpolysolver.dylib
SHARED_FLAGS = -dynamiclib -headerpad_max_install_names -install_name @rpath/$(SHARED_NAME) -compatibility_version $(SOVERSION) -current_version $(VERSION)
else
SHARED_NAME = libpolysolver.so.$(SOVERSION)
SHARED_LINK = libpolysolver.so
SHARED_FLAGS = -shared -Wl,-soname,$(SHARED_NAME)
endif
SHARED_TARGET = build/$(SHARED_NAME)

lib: $(LIB_TARGET) $(SHARED_TARGET)

$(LIB_TARGET): $(LIB_OBJS)
	rm -f $@
	ar rcs $@ $(LIB_OBJS)

$(SHARED_TARGET): $(LIB_OBJS)
	$(CC) $(SHARED_FLAGS) $(LDFLAGS) -o $@ $(LIB_OBJS) $(LIBS)
	ln -sf $(SHARED_NAME) build/$(SHARED_LINK)

# every global of the static library is ps_ or psi_, and the shared one exports only ps_
symcheck: $(LIB_TARGET) $(SHARED_TARGET)
	@bad=$$(nm -g $(LIB_TARGET) | awk '$$2 ~ /^[TDBSC]$$/ {print $$3}' | sed 's/^_//' | grep -vE '^psi?_' ); \
	 if [ -n "$$bad" ]; then echo "globals outside ps_/psi_ (run make namespace):"; echo "$$bad"; exit 1; fi
	@bad=$$(nm -g $(SHARED_TARGET) | awk '$$2 ~ /^[TDBSC]$$/ {print $$3}' | sed 's/^_//' | grep -vE '^(ps_|_?init$$|_?fini$$)'); \
	 if [ -n "$$bad" ]; then echo "shared library exports beyond ps_:"; echo "$$bad"; exit 1; fi
	@echo "symbols: $$(nm -g $(SHARED_TARGET) | awk '$$2 == "T"' | wc -l | tr -d ' ') exported, all ps_"

# a program with its own Expr, integrate, poly_add and RISCH_OK links and runs against both libraries
PUBLIC_CFLAGS = $(CFLAGS) $(GMP_CFLAGS) -Iinclude
clash: $(LIB_TARGET) $(SHARED_TARGET)
	$(CC) $(PUBLIC_CFLAGS) -o build/clash_static tests/clash.c $(LIB_TARGET) $(LDFLAGS) $(LIBS)
	$(CC) $(PUBLIC_CFLAGS) -o build/clash_shared tests/clash.c -Lbuild -lpolysolver -Wl,-rpath,$(CURDIR)/build $(LDFLAGS) $(LIBS)
	./build/clash_static && ./build/clash_shared

# the README's library examples, built against the library and run
EXAMPLES = rationals polynomials factor_solve expressions integrate definite linear big_numbers
EXAMPLE_BINS = $(addprefix examples/,$(EXAMPLES))

examples/%: examples/%.c $(LIB_TARGET) include/polysolver.h
	$(CC) $(PUBLIC_CFLAGS) $(LDFLAGS) -o $@ $< $(LIB_TARGET) $(LIBS)

examples: $(EXAMPLE_BINS)
	@for e in $(EXAMPLE_BINS); do echo "== $$e"; ./$$e || exit 1; done

# the fast lane after a change: the suite (with its CPU time) and the golden CLI corpus
check: $(TEST_TARGET) $(CLI_TARGET) symcheck clash
	@grep -q '#define PS_VERSION "$(VERSION)"' include/polysolver.h || { echo "PS_VERSION in include/polysolver.h is not VERSION ($(VERSION))"; exit 1; }
	./$(TEST_TARGET) | tail -1
	python3 tools/golden.py

# install: headers, both libraries, the CLI and a pkg-config file under $(DESTDIR)$(PREFIX)
install: $(LIB_TARGET) $(SHARED_TARGET) $(CLI_TARGET)
	sed -e 's|@PREFIX@|$(PREFIX)|' -e 's|@VERSION@|$(VERSION)|' polysolver.pc.in > build/polysolver.pc
	install -d $(DESTDIR)$(PREFIX)/include $(DESTDIR)$(PREFIX)/lib/pkgconfig $(DESTDIR)$(PREFIX)/bin
	install -m 644 include/polysolver.h $(DESTDIR)$(PREFIX)/include/
	install -m 644 $(LIB_TARGET) $(DESTDIR)$(PREFIX)/lib/
	install -m 755 $(SHARED_TARGET) $(DESTDIR)$(PREFIX)/lib/
	if [ "$$(uname -s)" = Darwin ]; then install_name_tool -id $(PREFIX)/lib/$(SHARED_NAME) $(DESTDIR)$(PREFIX)/lib/$(SHARED_NAME); fi
	ln -sf $(SHARED_NAME) $(DESTDIR)$(PREFIX)/lib/$(SHARED_LINK)
	install -m 755 $(CLI_TARGET) $(DESTDIR)$(PREFIX)/bin/
	install -m 644 build/polysolver.pc $(DESTDIR)$(PREFIX)/lib/pkgconfig/

uninstall:
	rm -f $(DESTDIR)$(PREFIX)/include/polysolver.h $(DESTDIR)$(PREFIX)/lib/libpolysolver.a \
	      $(DESTDIR)$(PREFIX)/lib/$(SHARED_NAME) $(DESTDIR)$(PREFIX)/lib/$(SHARED_LINK) \
	      $(DESTDIR)$(PREFIX)/bin/$(CLI_TARGET) $(DESTDIR)$(PREFIX)/lib/pkgconfig/polysolver.pc

# a source tarball of the tracked files: build/polysolver-$(VERSION).tar.gz
dist: | build
	git archive --format=tar.gz --prefix=polysolver-$(VERSION)/ -o build/polysolver-$(VERSION).tar.gz HEAD

# AddressSanitizer + UndefinedBehaviorSanitizer build of the test binary
sanitize: clean
	$(MAKE) CFLAGS="$(CFLAGS) -O1 -fsanitize=address,undefined -fno-omit-frame-pointer" \
	        LDFLAGS="$(LDFLAGS) -fsanitize=address,undefined"

clean:
	rm -rf build examples/*.dSYM
	rm -f $(TEST_TARGET) $(CLI_TARGET) $(EXAMPLE_BINS)

-include $(wildcard build/*.d)

.PHONY: all cli lib examples clean sanitize check namespace symcheck clash install uninstall dist
