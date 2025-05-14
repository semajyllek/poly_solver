#include <stdlib.h>
#include "token.h"



typedef struct ListNode {
    Token* token;              // The token in this node
    struct ListNode* next;     // Pointer to the next node
} ListNode;

typedef struct List {
    ListNode* head; // Head of the list
} List;

// Function to initialize the list
void initialize_list(List* list) {
    list->head = NULL;
}

// Function to add a token to the end of the list
void add_to_list(List* list, Token* token) {
    ListNode* new_node = (ListNode*)malloc(sizeof(ListNode));
    new_node->token = token;
    new_node->next = NULL;

    if (list->head == NULL) {
        list->head = new_node; // If the list is empty, set new node as head
    } else {
        ListNode* current = list->head;
        while (current->next != NULL) {
            current = current->next; // Traverse to the end of the list
        }
        current->next = new_node; // Add new node to the end of the list
    }
}

// Function to free the list
void free_list(List* list) {
    ListNode* current = list->head;
    while (current != NULL) {
        ListNode* next = current->next;
        free(current->token); // Free token if dynamically allocated
        free(current);
        current = next;
    }
    list->head = NULL; // Clear the list head
}
