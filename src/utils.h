#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char* append_strings(const char* str1, const char* str2) {
    // Calculate the length of both strings
    size_t len1 = strlen(str1);
    size_t len2 = strlen(str2);
    
    // Allocate memory for the new string
    // +1 for the null terminator
    char* result = malloc(len1 + len2 + 1);
    
    // Check if malloc succeeded
    if (result == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        exit(1);
    }
    
    // Copy the first string into the result
    strcpy(result, str1);
    
    // Append the second string to the result
    strcat(result, str2);
    
    return result; // Return the new concatenated string
}