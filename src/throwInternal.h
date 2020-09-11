#pragma once
#define THROWINTERNAL() throwInternal(__FILE__, __LINE__)
void throwInternal(const char* file, int line);
