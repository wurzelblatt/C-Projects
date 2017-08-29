#include <cstdio>
#include <iostream>


union Udoub{
	double d;
	unsigned char c[8];
};


void floatvalue(double input) {
	Udoub u;
	u.d = input;
	//u.d = 6.5;
	for ( int i=7; i>=0; i--) printf("%02x", u.c[i]);
	printf("\n");
}
