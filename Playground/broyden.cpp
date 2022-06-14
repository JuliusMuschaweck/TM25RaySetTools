#include<nr3.h>
#include<ludcmp.h>
#include<qrdcmp.h>
#include<roots_multidim.h>
#include<iostream>

VecDoub func(const VecDoub& x)
{
	VecDoub rv(2);
	rv[0] = exp(x[0]) + x[1];
	rv[1] = 2 + x[1] - x[0] * x[0];
	return rv;
}

void TestBroyden()
{
	VecDoub x0(2);
	x0[0] = 2;
	x0[1] = -3;
	Bool check;
	broydn(x0, check, func);
}