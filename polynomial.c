#include "polynomial.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Function to create a term
Term* create_term(double coefficient, int exponent) {
    Term* term = malloc(sizeof(Term));
    term->coefficient = coefficient;
    term->exponent = exponent;
    term->next = NULL;
    return term;
}

// Function to create a new polynomial
Polynomial* create_polynomial() {
    Polynomial* poly = malloc(sizeof(Polynomial));
    poly->head = NULL;
    return poly;
}

// Function to add a term to a polynomial
void add_term(Polynomial* poly, double coefficient, int exponent) {
    if (coefficient == 0) return; // Do not add zero terms

    Term* new_term = create_term(coefficient, exponent);
    if (poly->head == NULL) {
        poly->head = new_term; // If polynomial is empty, add the first term
        return;
    }

    Term* current = poly->head;
    Term* previous = NULL;
    
    // Traverse to find the correct position to insert the new term
    while (current != NULL && current->exponent > exponent) {
        previous = current;
        current = current->next;
    }
    
    if (current != NULL && current->exponent == exponent) {
        // Combine like terms
        current->coefficient += coefficient;
        if (current->coefficient == 0) {
            // Remove term if the coefficient becomes zero
            if (previous == NULL) {
                poly->head = current->next; // Remove head term
            } else {
                previous->next = current->next; // Remove the term
            }
            free(current);
        }
    } else {
        // Insert the new term in the correct position in the linked list
        new_term->next = current; 
        if (previous == NULL) {
            poly->head = new_term; // If inserting at the head
        } else {
            previous->next = new_term; // Link previous to new term
        }
    }
}

// Function for polynomial addition
Polynomial* poly_add(Polynomial* poly1, Polynomial* poly2) {
    Polynomial* result = create_polynomial();

    // Add terms from the first polynomial
    Term* current = poly1->head;
    while (current != NULL) {
        add_term(result, current->coefficient, current->exponent);
        current = current->next;
    }

    // Add terms from the second polynomial
    current = poly2->head;
    while (current != NULL) {
        add_term(result, current->coefficient, current->exponent);
        current = current->next;
    }

    return result; // Return the resulting polynomial
}

// Function for polynomial subtraction
Polynomial* poly_subtract(Polynomial* poly1, Polynomial* poly2) {
    // Negate the second polynomial
    Polynomial* neg_poly2 = create_polynomial();
    Term* current = poly2->head;
    while (current != NULL) {
        add_term(neg_poly2, -current->coefficient, current->exponent);
        current = current->next;
    }

    // Use poly_add to sum the first polynomial with the negated second polynomial
    Polynomial* result = poly_add(poly1, neg_poly2);

    // Free the negated polynomial if you no longer need it
    free_polynomial(neg_poly2);

    return result; // Return the resulting polynomial
}

// Function for polynomial multiplication
Polynomial* poly_multiply(Polynomial* poly1, Polynomial* poly2) {
    Polynomial* result = create_polynomial();

    // Iterate over each term in the first polynomial
    Term* term1 = poly1->head;
    while (term1 != NULL) {
        // Iterate over each term in the second polynomial
        Term* term2 = poly2->head;
        while (term2 != NULL) {
            // Multiply coefficients and add exponents
            double new_coefficient = term1->coefficient * term2->coefficient;
            int new_exponent = term1->exponent + term2->exponent;

            // Add the resulting term to the result polynomial
            add_term(result, new_coefficient, new_exponent);

            term2 = term2->next; // Move to the next term in poly2
        }
        term1 = term1->next; // Move to the next term in poly1
    }
    return result; // Return the resulting polynomial
}


// Function for polynomial division
Polynomial* poly_divide(Polynomial* dividend, Polynomial* divisor) {
    // Check if divisor is zero (no division by zero)
    if (divisor == NULL || divisor->head == NULL) {
        fprintf(stderr, "Error: Division by zero polynomial.\n");
        return NULL; // Handle error appropriately
    }

    Polynomial* quotient = create_polynomial(); // Polynomial to hold the quotient
    Polynomial* remainder = create_polynomial(); // Polynomial to hold the remainder

    // Copy the dividend to the remainder
    Term* current = dividend->head;
    while (current != NULL) {
        add_term(remainder, current->coefficient, current->exponent);
        current = current->next;
    }

    // Perform division until degree of remainder < degree of divisor
    while (remainder->head != NULL && remainder->head->exponent >= divisor->head->exponent) {
        // Leading term of the quotient
        double coeff = remainder->head->coefficient / divisor->head->coefficient; // Determine the leading coefficient
        int exp = remainder->head->exponent - divisor->head->exponent; // Calculate new exponent

        // Add leading term to the quotient
        add_term(quotient, coeff, exp);

        // Create a temporary polynomial for the product of the leading term with the divisor
        Polynomial* temp = create_polynomial();
        add_term(temp, coeff, exp); // New term for temp

        // Multiply the divisor by the temp polynomial
        Polynomial* multiplied_terms = poly_multiply(temp, divisor);

        // Subtract the multiplied terms from the remainder
        Polynomial* new_remainder = poly_subtract(remainder, multiplied_terms);

        // Free the resources
        free_polynomial(remainder);
        remainder = new_remainder;

        free_polynomial(multiplied_terms);
        free_polynomial(temp);
    }

    // Clean up the remainder when done
    free_polynomial(remainder);
    
    return quotient; // Return the resulting quotient
}



// Function to calculate the GCD of two polynomials
Polynomial* poly_gcd(Polynomial* poly1, Polynomial* poly2) {
    // Ensure both polynomials are valid
    if (!poly1 || !poly2 || (poly1->head == NULL && poly2->head == NULL)) {
        return NULL; // Return NULL if either polynomial is null or both are empty
    }

    // Use a temporary polynomial to store the dividend
    Polynomial* remainder;

    // Make sure poly1 is the polynomial with a non-zero leading term
    if (poly1->head == NULL) { // If poly1 is empty, return poly2
        return poly2; 
    }
    if (poly2->head == NULL) { // If poly2 is empty, return poly1
        return poly1; 
    }

    // While the second polynomial is not zero, continue the process
    while (poly2->head != NULL) {
        // Calculate the remainder of poly1 divided by poly2
        remainder = poly_divide(poly1, poly2); // Assuming poly_divide returns the quotient
        Polynomial* new_remainder = poly_mod(poly1, poly2); // Get the remainder of the division
        free_polynomial(poly1); // Free the previous dividend
        poly1 = poly2; // Assign poly1 as poly2 for the next iteration
        poly2 = new_remainder; // Update poly2 to the remainder
    }
    // At the end, poly1 represents the GCD
    return poly1; // Return GCD polynomial
}


// Function to calculate the LCM of two polynomials
Polynomial* poly_lcm(Polynomial* poly1, Polynomial* poly2) {
	// Calculate the product of the two polynomials
	Polynomial* product = poly_multiply(poly1, poly2);
	// Calculate the GCD of the two polynomials
	Polynomial* gcd = poly_gcd(poly1, poly2);
	// Divide the product by the GCD to get the LCM
	Polynomial* lcm = poly_divide(product, gcd);
	// Free the temporary product polynomial
	free_polynomial(product);
	return lcm; // Return the LCM polynomial
}





double evaluate_polynomial(Polynomial* poly, double x) {
    double result = 0.0;
    Term* current = poly->head;

    while (current != NULL) {
        result += current->coefficient * pow(x, current->exponent);
        current = current->next;
    }

    return result; // Return the evaluation result
}


int polynomial_degree(Polynomial* poly) {
    if (poly->head == NULL) { // Empty polynomial
        return -1; // Degree is undefined
    }
    return poly->head->exponent; // Highest exponent is the degree
}




// Function to free a polynomial
void free_polynomial(Polynomial* poly) {
	Term* current = poly->head;
	while (current != NULL) {
		Term* temp = current;
		current = current->next;
		free(temp);
	}
	free(poly);
}

// Function to print a polynomial
void print_polynomial(const Polynomial* poly) {
	Term* current = poly->head;
	while (current != NULL) {
		if (current->coefficient != 0) {
			if (current != poly->head && current->coefficient > 0) {
				printf(" + ");
			}
			if (current->coefficient < 0) {
				printf(" - ");
			}
			if (fabs(current->coefficient) != 1 || current->exponent == 0) {
				printf("%g", fabs(current->coefficient));
			}
			if (current->exponent > 0) {
				printf("x");
				if (current->exponent > 1) {
					printf("^%d", current->exponent);
				}
			}
		}
		current = current->next;
	}
	printf("\n");
}


// Function to calculate the derivative of a polynomial
Polynomial* poly_derivative(Polynomial* poly) {
	Polynomial* derivative = create_polynomial(); // Create an empty polynomial for the derivative
	Term* current = poly->head; // Start from the head of the polynomial

	while (current != NULL) {
		if (current->exponent > 0) {
			double new_coefficient = current->coefficient * current->exponent; // Differentiate the term
			int new_exponent = current->exponent - 1; // Reduce the exponent by 1
			add_term(derivative, new_coefficient, new_exponent); // Add the term to the derivative
		}
		current = current->next; // Move to the next term
	}

	return derivative; // Return the resulting derivative
}


// Function to calculate the integral of a polynomial
Polynomial* poly_integral(Polynomial* poly) {
	Polynomial* integral = create_polynomial(); // Create an empty polynomial for the integral
	Term* current = poly->head; // Start from the head of the polynomial

	while (current != NULL) {
		double new_coefficient = current->coefficient / (current->exponent + 1); // Integrate the term
		int new_exponent = current->exponent + 1; // Increase the exponent by 1
		add_term(integral, new_coefficient, new_exponent); // Add the term to the integral
		current = current->next; // Move to the next term
	}

	return integral; // Return the resulting integral
}

// Function to calculate the definite integral of a polynomial
double poly_definite_integral(Polynomial* poly, double lower, double upper) {
	Polynomial* integral = poly_integral(poly); // Calculate the indefinite integral
	double result = evaluate_polynomial(integral, upper) - evaluate_polynomial(integral, lower); // Evaluate the integral
	free_polynomial(integral); // Free the integral polynomial
	return result; // Return the definite integral
}

// Function to calculate the remainder of polynomial divi