CC = gcc
CFLAGS = -Wall -Wextra -g -I/opt/homebrew/include -MMD -MP
LDFLAGS = -L/opt/homebrew/lib -lgmp -lm

# shared library sources (everything except the two entry points)
LIB_SRCS = strbuf.c rational_num.c polynomial.c poly_zp.c poly_zx.c zx_factor.c factor.c roots.c solver.c subres.c logpart.c surd.c rootform.c definite.c constfold.c qsurd.c modgcd.c expfield.c risch.c kpoly.c residue.c cfield.c expsum.c tanfield.c deftower.c risch2.c k2poly.c \
           expression.c integrate.c matrix.c partial_frac.c token.c parser.c

# Test binary
TEST_SRCS = $(LIB_SRCS) test_symbolic_solver.c test_properties.c test_subres_ref.c
TEST_OBJS = $(TEST_SRCS:.c=.o)
TEST_TARGET = test_symbolic_solver

# CLI binary
CLI_SRCS = $(LIB_SRCS) main.c
CLI_OBJS = $(CLI_SRCS:.c=.o)
CLI_TARGET = poly_solver

all: $(TEST_TARGET)

$(TEST_TARGET): $(TEST_OBJS)
	$(CC) $(CFLAGS) -o $@ $(TEST_OBJS) $(LDFLAGS)

cli: $(CLI_TARGET)

$(CLI_TARGET): $(CLI_OBJS)
	$(CC) $(CFLAGS) -o $@ $(CLI_OBJS) $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# static library of everything but the entry points: cc -I. prog.c libpolysolver.a -lgmp -lm
LIB_OBJS = $(LIB_SRCS:.c=.o)
LIB_TARGET = libpolysolver.a

lib: $(LIB_TARGET)

$(LIB_TARGET): $(LIB_OBJS)
	ar rcs $@ $(LIB_OBJS)

# the README's library examples, built against the library and run
EXAMPLES = rationals polynomials factor_solve expressions integrate definite linear big_numbers
EXAMPLE_BINS = $(addprefix examples/,$(EXAMPLES))

examples/%: examples/%.c $(LIB_TARGET) poly_solver.h
	$(CC) $(CFLAGS) -I. -o $@ $< $(LIB_TARGET) $(LDFLAGS)

examples: $(EXAMPLE_BINS)
	@for e in $(EXAMPLE_BINS); do echo "== $$e"; ./$$e || exit 1; done

# the fast lane after a change: the suite (with its CPU time) and the golden CLI corpus
check: $(TEST_TARGET) $(CLI_TARGET)
	./$(TEST_TARGET) | tail -1
	python3 golden.py

# AddressSanitizer + UndefinedBehaviorSanitizer build of the test binary
sanitize: clean
	$(MAKE) CFLAGS="$(CFLAGS) -fsanitize=address,undefined -fno-omit-frame-pointer" \
	        LDFLAGS="$(LDFLAGS) -fsanitize=address,undefined"

clean:
	rm -f *.o *.d examples/*.d $(TEST_TARGET) $(CLI_TARGET) $(LIB_TARGET) $(EXAMPLE_BINS)

-include $(TEST_SRCS:.c=.d) main.d

.PHONY: all cli lib examples clean sanitize check
