#include "LinAlg3.h"
#include<limits>

void TestLinAlg3()
	{
	using R = float;
	using V3 = TVec3<R>;
	using M3 = TMat3<R>;
	R eps = std::numeric_limits<R>::epsilon();

	V3 v{ 1,2,3 };
	V3 w{ 4,5, R(6.1) };
	V3 v0{ 0,0,0 };

	//template< typename R>
	//TVec3<R> operator+(const TVec3<R>& lhs, const TVec3<R>& rhs);
	assert((v + w) == (V3{ 5,7,R(9.1) }));

	//template< typename R>
	//TVec3<R> operator-(const TVec3<R>& lhs, const TVec3<R>& rhs);
	assert((v - w) == (V3{ -3,-3,R(-3.1) }));
	//
	//// scalar product
	//template< typename R>
	//R operator*(const TVec3<R>& lhs, const TVec3<R>& rhs);
	assert((v * w) == 4 + 10 + 3 * R(6.1));
	//
	//// multiplication with scalar from left and right
	//template< typename R>
	//TVec3<R> operator*(R lhs, const TVec3<R>& rhs);
	assert((2 * w) == (V3{ 8, 10, R(12.2) }));
	//template< typename R>
	//TVec3<R> operator*(const TVec3<R>& lhs, R rhs);
	assert((2 * w) == w * 2);
	//
	//// divide by scalar
	//template< typename R>
	//TVec3<R> operator/(const TVec3<R>& lhs, R rhs);
	assert((w / 2) == w * 0.5);
	//
	//// add subtract vector assignment
	//template< typename R>
	//const TVec3<R>& operator +=(TVec3<R>& lhs, const TVec3<R>& rhs);
	{
	V3 tmp = v;
	tmp += w;
	assert(tmp == v + w);
	}
	//template< typename R>
	//const TVec3<R>& operator -=(TVec3<R>& lhs, const TVec3<R>& rhs);
	{
	V3 tmp = v;
	tmp -= w;
	assert(tmp == v - w);
	}
	//
	//// multiply divide scalar assignment
	//template< typename R>
	//const TVec3<R>& operator *=(TVec3<R>& lhs, R rhs);
	{
	V3 tmp = v;
	tmp *= 2;
	assert(tmp == (v*2.0));
	}
	//template< typename R>
	//const TVec3<R>& operator /=(TVec3<R>& lhs, R rhs);
	{
	V3 tmp = v;
	tmp /= 2.0;
	assert(tmp == (v / 2));
	}
	//
	//// Norm
	//template< typename R>
	//R Norm(const TVec3<R>& rhs);
	assert(Norm(V3({ 1,2,3 })) == sqrtf(1.0f + 2 * 2 + 3 * 3));
	//
	//// ***************  matrix operations *****************************
	//
	M3 M{ V3({1,2,3}), {4,5,2}, {-4,0,4} };
	M3 Mx2{ V3({2,4,6}), {8,10,4}, {-8,0,8} };
	M3 Mx3{ V3({3,6,9}), {12,15,6}, {-12,0,12} };
	assert(Det(M) != 0);
	//// linear element access for matrix
	//template< typename R, size_t i>
	//R Elem(const TMat3<R>& lhs);
	assert((Elem<float, 6>(M)) == R(-4));
	//template< typename R>
	//R Elem(const TMat3<R>& lhs, size_t i); // i < 9, but this is not checked;
	assert((Elem(M, 6)) == -4);
	//template< typename R, size_t i>
	//TVec3<R> Column(const TMat3<R>& rhs);
	assert((Column<float, 1>(M)) == (V3{ 2,5,0 }));
	//
	//
	//// matrix-scalar operations
	//// matrix-scalar multiplication from left and right
	//template< typename R>
	//TMat3<R> operator*(R lhs, const TMat3<R>& rhs);
	assert((2 * M) == Mx2);
	//template< typename R>
	//TMat3<R> operator*(const TMat3<R>& lhs, R rhs);
	assert((M * 2) == Mx2);
	//
	//// matrix-scalar division
	//template< typename R>
	//TMat3<R> operator/(const TMat3<R>& lhs, R rhs);
	assert((Mx2 / 2) == M);
	//
	//// matrix-scalar multiplication and division assigment
	//template< typename R>
	//const TMat3<R>& operator *= (TMat3<R>& lhs, R rhs);
	{
	M3 tmp(M);
	tmp *= 2;
	assert(tmp == Mx2);
	}
	//template< typename R>
	//const TMat3<R>& operator /= (TMat3<R>& lhs, R rhs);
	{
	M3 tmp(Mx2);
	tmp /= 2;
	assert(tmp == M);
	}
	//
	//
	//// matrix-vector operations
	//// matrix-vector multiplication from left and right
	//template< typename R>
	//TVec3<R> operator* (const TVec3<R>& lhs, const TMat3<R>& rhs);
		// 	V3 v{ 1,2,3 };
	assert((v * M) == (V3{ 1 + 8 - 12,2 + 10,3 + 4 + 12 }));
	//template< typename R>
	//TVec3<R> operator* (const TMat3<R>& lhs, const TVec3<R>& rhs);
	assert((M * v) == (V3{ 1 + 4 + 9, 4 + 10 + 6, -4 + 12 }));
	//
	//
	//// matrix-matrix operations
	//// matrix addition and subtraction
	//template< typename R>
	//TMat3<R> operator+ (const TMat3<R>& lhs, const TMat3<R>& rhs);
	assert((M + Mx2) == (Mx3));
	assert((M + Mx2) == (Mx2 + M));
	//template< typename R>
	//TMat3<R> operator- (const TMat3<R>& lhs, const TMat3<R>& rhs);
	assert((Mx3 - Mx2) == M);
	//
	//// matrix-matrix multiplication
	//template< typename R>
	//TMat3<R> operator*(const TMat3<R>& lhs, const TMat3<R>& rhs);
	{
	M3 tmp{ V3{1,2,3}, {4,5,6}, {7,8,9 } };
	M3 tmp2{ V3{3,1,4}, {2,3,1}, { 3,4,9 } };
	M3 test{ V3{16,19,33}, {40,43,75}, {64,67,117  } };
	assert(test == (tmp * tmp2));
	}
	//
	//// add subtract assigment
	//template< typename R>
	//const TMat3<R>& operator += (TMat3<R>& lhs, const TMat3<R>& rhs);
	{
	M3 tmp = Mx2;
	tmp += M;
	assert(tmp == Mx3);
	}
	//template< typename R>
	//const TMat3<R>& operator -= (TMat3<R>& lhs, const TMat3<R>& rhs);
	{
	M3 tmp = Mx2;
	tmp -= M;
	assert(tmp == M);
	}
	//
	//// other operations: outer product of vectors
	//template< typename R>
	//TMat3<R> Outer(const TVec3<R>& lhs, const TVec3<R>& rhs);
	{
	M3 tmp = Outer(V3{ 2, 3, 5 }, V3{ 7,11,13 });
	assert(tmp == (M3{ V3{14,22,26}, {21,33,39}, {35,55,65} }));
	}
	//
	//// matrix inverse, transpose, determinant, trace, Frobenius norm, diagonal
	//template< typename R>
	//TMat3<R> Inverse(const TMat3<R>& rhs);
	{
	M3 a{
		V3{16,19,33},
		{40,43,75},
		{40,35,45}
		};
	M3 ia{
		V3{R(-0.575), 0.25, R(0.005)},
		{1, -0.5, R(0.1)},
		{R(-0.2666666666666666666666666), R(0.16666666666666666666666), R(-0.06)}
		};
	M3 tmp = a * ia;
	R rzero = FrobeniusNorm(Eye3<R>() - tmp);
	tmp = Inverse(a);
	R d = Det(a);
	M3 mzero = Inverse(a) - ia;
	rzero = abs(eps - FrobeniusNorm(mzero));

	M3 tmp2 = RotX3<R>(0.1);
	M3 itmp2 = Inverse(tmp2);
	M3 tmp3 = tmp2 - Transpose(itmp2);

	M3 tt = Inverse(Eye3<R>());
	assert(rzero < 5 * eps);
	}
	//template< typename R>
	//TMat3<R> Transpose(const TMat3<R>& rhs);
	assert(
		(
		Transpose(M3{
		V3{ 11, 12, 13 },
		{21,22,23 },
		{31,32,33 }
		})) ==
		(M3{
		V3{ 11, 21, 31 },
		{12,22,32 },
		{13,23,33 }
		}));
	//template< typename R>
	//void DoTranspose(TMat3<R>& rhs);
	{
	M3 tmp{
		V3{ 11, 12, 13 },
		{21,22,23 },
		{31,32,33 } };
	M3 tmp2 = tmp;
	DoTranspose(tmp2);
	assert(tmp2 == (Transpose(tmp)));
	}
	//template< typename R>
	//R Det(const TMat3<R> rhs);
	{
	M3 a{
	V3{16,19,33},
	{40,43,75},
	{40,35,45} };
	assert(Det(a) == 1200);
	}
	//template< typename R>
	//R Trace(const TMat3<R> rhs);
	{
	M3 a{
	V3{16,19,33},
	{40,43,75},
	{40,35,45} };
	assert((Trace(a)) == (16 + 43 + 45));
	}
	//template< typename R>
	//TVec3<R> Diag(const TMat3<R> rhs);
	{
	M3 a{
	V3{16,19,33},
	{40,43,75},
	{40,35,45} };
	assert((Diag(a)) == (V3{ 16 , 43 , 45 }));
	}
	//template< typename R>
	//R FrobeniusNorm(const TMat3<R> rhs);
	{
	M3 a{
	V3{16,19,33},
	{40,43,75},
	{40,35,45} };
	assert((FrobeniusNorm(a)) == (sqrtf(R(15630))));
	}
	//
	//// generators
	//template< typename R>
	//TMat3<R> Zero3(); // all zeros
	assert((Zero3<R>()) == (M3{ V3{0,0,0},{0,0,0},{0,0,0} }));
	//template< typename R>
	//TMat3<R> Ones3(); // all ones
	assert((Ones3<R>()) == (M3{ V3{1,1,1},{1,1,1},{1,1,1} }));
	//template< typename R>
	//TMat3<R> Eye3(); // identity matrix
	assert((Eye3<R>()) == (M3{ V3{1,0,0},{0,1,0},{0,0,1} }));
	//
	//// rotations, see https://en.wikipedia.org/wiki/Rotation_matrix
	//template< typename R>
	//TMat3<R> RotX3(R theta_rad);
	//template< typename R>
	//TMat3<R> RotY3(R theta_rad);
	//template< typename R>
	//TMat3<R> RotZ3(R theta_rad);
	//template< typename R>
	R theta = asin(0.5); // 30°
	M3 Rx = RotX3<R>(theta);
	assert((Rx[0][0] == 1) && (FrobeniusNorm(Transpose(Rx) - Inverse(Rx)) < 5 * eps) && (abs(Det(Rx) - 1) < 5 * eps));
	M3 Ry = RotY3<R>(theta);
	assert((Ry[1][1] == 1) && (FrobeniusNorm(Transpose(Ry) - Inverse(Ry)) < 5 * eps) && (abs(Det(Ry) - 1) < 5 * eps));
	M3 Rz = RotZ3<R>(theta);
	assert((Rz[2][2] == 1) && (FrobeniusNorm(Transpose(Rx) - Inverse(Rx)) < 5 * eps) && (abs(Det(Rx) - 1) < 5 * eps));
	//TMat3<R> Rot3(const TVec3<R>& u, R theta_rad);
	M3 Rxyz = Rot3<R>(V3{ 1,2,3 }, theta);
	assert((FrobeniusNorm(Transpose(Rxyz) - Inverse(Rxyz)) < 5 * eps) && (abs(Det(Rxyz) - 1) < 5 * eps));
	V3 test = Rxyz * V3{ 1,2,3 };
	assert(Norm(test - V3{ 1,2,3 }) < 5 * eps);
	{
	M3 a{
		V3{16,19,33},
		{40,43,75},
		{40,35,45} };
	V3 b{ 1,2,5 };
	V3 x = Solve(a, b);
	V3 diff = a * x - b;
	R nd = Norm(diff) / FrobeniusNorm(a);
	assert(nd < 5 * eps);
	}

}