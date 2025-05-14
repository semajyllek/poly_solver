#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

typedef struct Term {
    double coefficient; // Coefficient of the term
    int exponent;       // Exponent of the term
    struct Term* next;  // Pointer to the next term in the polynomial
} Term;

typedef struct Polynomial {
    Term* head; // Pointer to the head of the linked list of terms
} Polynomial;

// Function prototypes
Polynomial* poly_add(Polynomial* poly1, Polynomial* poly2);
Polynomial* poly_subtract(Polynomial* poly1, Polynomial* poly2);
Polynomial* poly_multiply(Polynomial* poly1, Polynomial* poly2);
Polynomial* poly_divide(Polynomial* dividend, Polynomial* divisor);
void free_polynomial(Polynomial* poly);
void print_polynomial(const Polynomial* poly);

#endif // POLYNOMIAL_H
