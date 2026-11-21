#ifndef STRBUF_H
#define STRBUF_H

#include <stddef.h>

// growable, always NUL-terminated string builder (no fixed-size buffers anywhere)
typedef struct {
    char* buf;
    size_t len;
    size_t cap;
} StrBuf;

StrBuf sb_new(void);
void   sb_append(StrBuf* sb, const char* s);
void   sb_append_owned(StrBuf* sb, char* s);         // appends s, then frees it
void   sb_appendf(StrBuf* sb, const char* fmt, ...); // printf-style append of any length
char*  sb_finish(StrBuf* sb);                        // hands the buffer to the caller

#endif
