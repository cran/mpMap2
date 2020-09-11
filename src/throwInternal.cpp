#include "throwInternal.h"
#include <sstream>
#include <Rcpp.h>
void throwInternal(const char* file, int line)
{
	std::stringstream ss;
	ss << file << ":" << line << " Internal error";
	throw std::runtime_error(ss.str().c_str());
}
