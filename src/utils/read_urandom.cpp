#include "utils.hpp"
#include<gmpxx.h>

uint64_t read_urandom()
{
	log("Read random called");
	union {
		uint64_t value;
		char cs[sizeof(uint64_t)];
	} u;

	std::ifstream rfin("/dev/urandom");
	rfin.read(u.cs, sizeof(u.cs));
	rfin.close();

	return u.value;
}