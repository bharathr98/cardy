#include "utils.hpp"

unsigned long long int read_urandom()
{
	union {
		unsigned long long int value;
		char cs[sizeof(unsigned long long int)];
	} u;

	std::ifstream rfin("/dev/urandom");
	rfin.read(u.cs, sizeof(u.cs));
	rfin.close();

	return u.value;
}