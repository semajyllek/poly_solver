#include "token.h"
#include <stdlib.h>

typedef struct StackNode {
    Token* token;             // Token contained in the node
    struct StackNode* next;   // Pointer to the next node
} StackNode;

typedef struct Stack {
    StackNode* top; // The top node of the stack
} Stack;

// Function to initialize the stack
void init_stack(Stack* stack) {
    stack->top = NULL;
}

// Function to check if the stack is empty
bool is_empty(Stack* stack) {
    return stack->top == NULL;
}

// Function to push a token onto the stack
void push(Stack* stack, Token* token) {
    StackNode* new_node = (StackNode*)malloc(sizeof(StackNode));
    new_node->token = token;
    new_node->next = stack->top;
    stack->top = new_node;
}

// Function to pop a token from the stack
Token* pop(Stack* stack) {
    if (is_empty(stack)) return NULL; // Return NULL if stack is empty
    StackNode* top_node = stack->top;
    Token* token = top_node->token;
    stack->top = top_node->next; // Update the top pointer
    free(top_node);
    return token;
}

// Function to peek at the top token (without removing it)
Token* peek(Stack* stack) {
    if (is_empty(stack)) return NULL;
    return stack->top->token; // Return the top token without popping it
}
