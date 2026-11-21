#include "strbuf.h"
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

StrBuf sb_new(void) {
    StrBuf sb = { malloc(64), 0, 64 };
    sb.buf[0] = '\0';
    return sb;
}

// grow capacity so that `extra` more bytes plus the terminator fit
static void sb_reserve(StrBuf* sb, size_t extra) {
    while (sb->len + extra + 1 > sb->cap) sb->cap *= 2;
    sb->buf = realloc(sb->buf, sb->cap);
}

void sb_append(StrBuf* sb, const char* s) {
    size_t slen = strlen(s);
    sb_reserve(sb, slen);
    memcpy(sb->buf + sb->len, s, slen + 1);
    sb->len += slen;
}

void sb_append_owned(StrBuf* sb, char* s) {
    sb_append(sb, s);
    free(s);
}

void sb_appendf(StrBuf* sb, const char* fmt, ...) {
    va_list ap, ap2;
    va_start(ap, fmt);
    va_copy(ap2, ap);
    int n = vsnprintf(NULL, 0, fmt, ap); // measure first, then write exactly
    va_end(ap);
    if (n > 0) {
        sb_reserve(sb, (size_t)n);
        vsnprintf(sb->buf + sb->len, (size_t)n + 1, fmt, ap2);
        sb->len += (size_t)n;
    }
    va_end(ap2);
}

char* sb_finish(StrBuf* sb) {
    char* out = sb->buf;
    sb->buf = NULL;
    sb->len = sb->cap = 0;
    return out;
}
