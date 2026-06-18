CC = gcc
CFLAGS = -Wall -Wextra -g
LDFLAGS = -lm

# Test binary
TEST_SRCS = rational_num.c polynomial.c expression.c factor.c solver.c matrix.c partial_frac.c token.c parser.c test_symbolic_solver.c
TEST_OBJS = $(TEST_SRCS:.c=.o)
TEST_TARGET = test_symbolic_solver

# CLI binary
CLI_SRCS = rational_num.c polynomial.c expression.c factor.c solver.c matrix.c partial_frac.c token.c parser.c main.c
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

clean:
	rm -f $(TEST_OBJS) $(CLI_OBJS) $(TEST_TARGET) $(CLI_TARGET)
