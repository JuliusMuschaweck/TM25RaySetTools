// linear algebra in 3D
// J. Muschaweck, JMO GmbH 2019
#ifndef LINALG3_H
#define LINALG3_H

#include <array>
#include <assert.h>
#include <math.h>
// #include <algorithm>
#include <utility> // for swap

// type definitions for 3D vectors and matrices
template< typename R>
using TVec3 = std::array<R, 3>;
using TVec3f = TVec3<float>;
using TVec3d = TVec3<double>;

template< typename R>
using TMat3 = std::array<std::array<R, 3>, 3>;
using TMat3f = TMat3<float>;
using TMat3d = TMat3<double>;

// ************  vector operations  ********************** 

//vector addition and subtraction
template< typename R>
TVec3<R> operator+(const TVec3<R>& lhs, const TVec3<R>& rhs);
template< typename R>
TVec3<R> operator-(const TVec3<R>& lhs, const TVec3<R>& rhs);

// scalar product
template< typename R>
R operator*(const TVec3<R>& lhs, const TVec3<R>& rhs);
template< typename R>
R Sqr(const TVec3<R>& rhs); // rhs*rhs

// cross product
template< typename R>
TVec3<R> Cross (const TVec3<R>& lhs, const TVec3<R>& rhs);

// multiplication with scalar from left and right
template< typename R, typename U>
TVec3<R> operator*(U lhs, const TVec3<R>& rhs);
template< typename R, typename U>
TVec3<R> operator*(const TVec3<R>& lhs, U rhs);

// divide by scalar
template< typename R, typename U>
TVec3<R> operator/(const TVec3<R>& lhs, U rhs);

// add subtract vector assignment
template< typename R>
const TVec3<R>& operator +=(TVec3<R>& lhs, const TVec3<R>& rhs);
template< typename R>
const TVec3<R>& operator -=(TVec3<R>& lhs, const TVec3<R>& rhs);

// multiply divide scalar assignment
template< typename R, typename U>
const TVec3<R>& operator *=(TVec3<R>& lhs, U rhs);
template< typename R, typename U>
const TVec3<R>& operator /=(TVec3<R>& lhs, U rhs);

// Norm
template< typename R>
R Norm(const TVec3<R>& rhs);

// Unit vector
template<typename R>
TVec3<R> UnitVector(const TVec3<R>& rhs);

// ***************  matrix operations *****************************

// linear element access for matrix
template< typename R, size_t i>
R Elem(const TMat3<R>& lhs);
template< typename R>
R Elem(const TMat3<R>& lhs, size_t i); // i < 9, but this is not checked;
template< typename R, size_t i>
TVec3<R> Column(const TMat3<R>& rhs);


// matrix-scalar operations
// matrix-scalar multiplication from left and right
template< typename R, typename U>
TMat3<R> operator*(U lhs, const TMat3<R>& rhs);
template< typename R, typename U>
TMat3<R> operator*(const TMat3<R>& lhs, U rhs);

// matrix-scalar division
template< typename R, typename U>
TMat3<R> operator/(const TMat3<R>& lhs, U rhs);

// matrix-scalar multiplication and division assigment
template< typename R, typename U>
const TMat3<R>& operator *= (TMat3<R>& lhs, U rhs);
template< typename R, typename U>
const TMat3<R>& operator /= (TMat3<R>& lhs, U rhs);


// matrix-vector operations
// matrix-vector multiplication from left and right
template< typename R>
TVec3<R> operator* (const TVec3<R>& lhs, const TMat3<R>& rhs);
template< typename R>
TVec3<R> operator* (const TMat3<R>& lhs, const TVec3<R>& rhs);


// matrix-matrix operations
// matrix addition and subtraction
template< typename R>
TMat3<R> operator+ (const TMat3<R>& lhs, const TMat3<R>& rhs);
template< typename R>
TMat3<R> operator- (const TMat3<R>& lhs, const TMat3<R>& rhs);

// matrix-matrix multiplication
template< typename R>
TMat3<R> operator*(const TMat3<R>& lhs, const TMat3<R>& rhs);

// add subtract assigment
template< typename R>
const TMat3<R>& operator += (TMat3<R>& lhs, const TMat3<R>& rhs);
template< typename R>
const TMat3<R>& operator -= (TMat3<R>& lhs, const TMat3<R>& rhs);

// other operations: outer product of vectors
template< typename R>
TMat3<R> Outer(const TVec3<R>& lhs, const TVec3<R>& rhs);

// matrix inverse, transpose, determinant, trace, Frobenius norm, diagonal
template< typename R>
TMat3<R> Inverse(const TMat3<R>& rhs);
template< typename R>
TMat3<R> Transpose(const TMat3<R>& rhs);
template< typename R>
void DoTranspose(TMat3<R>& rhs);
template< typename R>
R Det(const TMat3<R> rhs);
template< typename R>
R Trace(const TMat3<R> rhs);
template< typename R>
TVec3<R> Diag(const TMat3<R> rhs);
template< typename R>
R FrobeniusNorm(const TMat3<R> rhs);


// generators
template< typename R>
TMat3<R> Zero3(); // all zeros
template< typename R>
TMat3<R> Ones3(); // all ones
template< typename R>
TMat3<R> Eye3(); // identity matrix

// rotations, see https://en.wikipedia.org/wiki/Rotation_matrix
template< typename R, typename U>
TMat3<R> RotX3(U theta_rad); 
template< typename R, typename U>
TMat3<R> RotY3(U theta_rad);
template< typename R, typename U>
TMat3<R> RotZ3(U theta_rad);
template< typename R, typename U>
TMat3<R> Rot3(const TVec3<R>& u, U theta_rad);

// solve linear equation 
template< typename R>  // solve Ax=b
TVec3<R> Solve(const TMat3<R>& A, const TVec3<R>& b);


// testing
void TestLinAlg3();

// *********************************************************************************************
// template definitions
// *********************************************************************************************

// ************  vector operations  ********************** 

//vector addition and subtraction
template< typename R>
TVec3<R> operator+(const TVec3<R>& lhs, const TVec3<R>& rhs)
	{
	return TVec3<R>{lhs[0] + rhs[0], lhs[1] + rhs[1], lhs[2] + rhs[2]};
	}

template< typename R>
TVec3<R> operator-(const TVec3<R>& lhs, const TVec3<R>& rhs)
	{
	return TVec3<R>{lhs[0] - rhs[0], lhs[1] - rhs[1], lhs[2] - rhs[2]};
	};

// scalar product
template< typename R>
R operator*(const TVec3<R>& lhs, const TVec3<R>& rhs)
	{
	return lhs[0] * rhs[0] + lhs[1] * rhs[1] + lhs[2] * rhs[2];
	}
template<typename R>
inline R Sqr(const TVec3<R>& rhs)
	{
	return rhs[0] * rhs[0] + rhs[1] * rhs[1] + rhs[2] * rhs[2];
	}
template<typename R>
inline TVec3<R> Cross(const TVec3<R>& lhs, const TVec3<R>& rhs)
	{
	return { lhs[1] * rhs[2] - lhs[2] * rhs[1], lhs[2] * rhs[0] - lhs[0] * rhs[2], lhs[0] * rhs[1] - lhs[1] * rhs[0] };
	}
;

// multiplication with scalar from left and right
template< typename R, typename U>
TVec3<R> operator*(U lhs, const TVec3<R>& rhs)
	{
	R l = static_cast<R>(lhs);
	return TVec3<R>{l*rhs[0], l*rhs[1], l * rhs[2]};
	};

template< typename R, typename U>
TVec3<R> operator*(const TVec3<R>& lhs, U rhs)
	{
	return rhs * lhs;
	};

// divide by scalar
template< typename R, typename U>
TVec3<R> operator/(const TVec3<R>& lhs, U rhs)
	{
	R r = static_cast<R>(rhs);
	assert(r != 0);
	return TVec3<R>{lhs[0] / r, lhs[1] / r, lhs[2] / r};
	}


// add subtract vector assignment
template< typename R>
const TVec3<R>& operator +=(TVec3<R>& lhs, const TVec3<R>& rhs)
	{
	lhs[0] += rhs[0]; lhs[1] += rhs[1]; lhs[2] += rhs[2];
	return lhs;
	}

template< typename R>
const TVec3<R>& operator -=(TVec3<R>& lhs, const TVec3<R>& rhs)
	{
	lhs[0] -= rhs[0];
	lhs[1] -= rhs[1];
	lhs[2] -= rhs[2];
	return lhs;
	}

// multiply divide scalar assignment
template< typename R, typename U>
const TVec3<R>& operator *=(TVec3<R>& lhs, U rhs)
	{
	R r = static_cast<R>(rhs);
	lhs[0] *= r;
	lhs[1] *= r;
	lhs[2] *= r;
	return lhs;
	}

template< typename R, typename U>
const TVec3<R>& operator /=(TVec3<R>& lhs, U rhs)
	{
	R r = static_cast<R>(rhs);
	assert(r != 0);
	lhs[0] /= r;
	lhs[1] /= r;
	lhs[2] /= r;
	return lhs;
	}


// Norm
template< typename R>
R Norm(const TVec3<R>& rhs)
	{
	return sqrt(rhs*rhs);
	}

// Unit vector
template<typename R>
inline TVec3<R> UnitVector(const TVec3<R>& rhs)
	{
	return rhs / Norm(rhs);
	}
;



// ***************  matrix operations *****************************

// linear element access for matrix
template< typename R, size_t i>
R Elem(const TMat3<R>& lhs)
	{
	static_assert(i < 9, "Elem: i >= 9");
	return lhs[i / 3][i % 3];
	}

template< typename R>
R Elem(const TMat3<R>& lhs, size_t i) // i < 9, but this is not checked;
	{
	assert(i < 9);
	return lhs[i / 3][i % 3];
	}
template< typename R, size_t i>
TVec3<R> Column(const TMat3<R>& rhs)
	{
	return TVec3<R>{rhs[0][i], rhs[1][i], rhs[2][i]};
	}


// matrix-scalar operations
// matrix-scalar multiplication from left and right
template< typename R, typename U>
TMat3<R> operator*(U lhs, const TMat3<R>& rhs)
	{
	R l = static_cast<R>(lhs);
	return TMat3<R>{
			TVec3<R>{rhs[0][0] * l, rhs[0][1] * l, rhs[0][2] * l},
			{rhs[1][0] * l, rhs[1][1] * l, rhs[1][2] * l},
			{rhs[2][0] * l, rhs[2][1] * l, rhs[2][2] * l } 
		};
	}
template<typename R, typename U>
TMat3<R> operator*(const TMat3<R>& lhs, U rhs)
	{
	return rhs * lhs;
	}

// matrix-scalar division
template< typename R, typename U>
TMat3<R> operator/(const TMat3<R>& lhs, U rhs)
	{
	assert(rhs != 0);
	return lhs * (1 / static_cast<R>(rhs));
	}

// matrix-scalar multiplication and division assigment
template< typename R, typename U>
const TMat3<R>& operator *= (TMat3<R>& lhs, U rhs)
	{
	R r = static_cast<R>(rhs);
	lhs[0][0] *= r; lhs[0][1] *= r; lhs[0][2] *= r;
	lhs[1][0] *= r; lhs[1][1] *= r; lhs[1][2] *= r;
	lhs[2][0] *= r; lhs[2][1] *= r; lhs[2][2] *= r;
	return lhs;
	}

template< typename R, typename U>
const TMat3<R>& operator /= (TMat3<R>& lhs, U rhs)
	{
	assert(rhs != 0);
	return lhs *= (1 / static_cast<R>(rhs));
	}


// matrix-vector operations
// matrix-vector multiplication from left and right
template< typename R>
TVec3<R> operator* (const TVec3<R>& lhs, const TMat3<R>& rhs)
	{
	return TVec3<R>{
		lhs[0] * rhs[0][0] + lhs[1] * rhs[1][0] + lhs[2] * rhs[2][0],
		lhs[0] * rhs[0][1] + lhs[1] * rhs[1][1] + lhs[2] * rhs[2][1],
		lhs[0] * rhs[0][2] + lhs[1] * rhs[1][2] + lhs[2] * rhs[2][2]
		};
	};

template< typename R>
TVec3<R> operator* (const TMat3<R>& lhs, const TVec3<R>& rhs)
	{
	return TVec3<R>{ lhs[0] * rhs, lhs[1] * rhs, lhs[2] * rhs};
	};


// matrix-matrix operations
// matrix addition and subtraction
template< typename R>
TMat3<R> operator+ (const TMat3<R>& lhs, const TMat3<R>& rhs)
	{
	return TMat3<R>{ lhs[0] + rhs[0], lhs[1] + rhs[1], lhs[2] + rhs[2] };
	}

template< typename R>
TMat3<R> operator- (const TMat3<R>& lhs, const TMat3<R>& rhs)
	{
	return TMat3<R>{ lhs[0] - rhs[0], lhs[1] - rhs[1], lhs[2] - rhs[2] };
	}

// matrix-matrix multiplication
template< typename R>
TMat3<R> operator*(const TMat3<R>& lhs, const TMat3<R>& rhs)
	{
	return TMat3<R>{
		TVec3<R>{ lhs[0] * Column<R, 0>(rhs), lhs[0] * Column<R, 1>(rhs), lhs[0] * Column<R, 2>(rhs)},
			{ lhs[1] * Column<R, 0>(rhs), lhs[1] * Column<R, 1>(rhs), lhs[1] * Column<R, 2>(rhs) },
			{ lhs[2] * Column<R, 0>(rhs), lhs[2] * Column<R, 1>(rhs), lhs[2] * Column<R, 2>(rhs) }
		};
	}

// add subtract assigment
template< typename R>
const TMat3<R>& operator += (TMat3<R>& lhs, const TMat3<R>& rhs)
	{
	lhs[0] += rhs[0]; lhs[1] += rhs[1]; lhs[2] += rhs[2];
	return lhs;
	}
template< typename R>
const TMat3<R>& operator -= (TMat3<R>& lhs, const TMat3<R>& rhs)
	{
	lhs[0] -= rhs[0]; lhs[1] -= rhs[1]; lhs[2] -= rhs[2];
	return lhs;
	}

// other operations: outer product of vectors
template< typename R>
TMat3<R> Outer(const TVec3<R>& lhs, const TVec3<R>& rhs)
	{
	return TMat3<R>{
		TVec3<R>{ lhs[0] * rhs[0], lhs[0] * rhs[1], lhs[0] * rhs[2]},
			{ lhs[1] * rhs[0], lhs[1] * rhs[1], lhs[1] * rhs[2] },
			{ lhs[2] * rhs[0], lhs[2] * rhs[1], lhs[2] * rhs[2] }
		};
	}

// matrix inverse, transpose, determinant, trace, Frobenius norm, diagonal
template< typename R>
TMat3<R> Inverse(const TMat3<R>& rhs)
	{
	R det = Det(rhs);
	assert(det != 0);
	R a = rhs[0][0]; R b = rhs[0][1]; R c = rhs[0][2];
	R d = rhs[1][0]; R e = rhs[1][1]; R f = rhs[1][2];
	R g = rhs[2][0]; R h = rhs[2][1]; R i = rhs[2][2];
	return TMat3<R> {
				TVec3<R>{e*i-f*h, c*h-b*i, b*f-c*e }, 
				TVec3<R>{f*g-d*i, a*i-c*g, c*d-a*f },
				TVec3<R>{d*h-e*g, b*g-a*h, a*e-b*d }
		} / det;
	}

template< typename R>
TMat3<R> Transpose(const TMat3<R>& rhs)
	{
	return TMat3<R> {
			TVec3<R>{rhs[0][0], rhs[1][0], rhs[2][0] },
			{rhs[0][1], rhs[1][1], rhs[2][1] },
			{rhs[0][2], rhs[1][2], rhs[2][2] }
		};
	}

template< typename R>
void DoTranspose(TMat3<R>& rhs)
	{
	std::swap(rhs[1][0], rhs[0][1]);
	std::swap(rhs[2][0], rhs[0][2]);
	std::swap(rhs[2][1], rhs[1][2]);
	}

template< typename R>
R Det(const TMat3<R> rhs)
	{
	auto SubDet = [rhs](size_t i, size_t j)
		{
		return rhs[(i + 1) % 3][(j + 1) % 3] * rhs[(i + 2) % 3][(j + 2) % 3]
			- rhs[(i + 2) % 3][(j + 1) % 3] * rhs[(i + 1) % 3][(j + 2) % 3];
		};
	return rhs[0][0] * SubDet(0, 0) + rhs[1][0] * SubDet(1, 0) + rhs[2][0] * SubDet(2, 0);
	}

template< typename R>
R Trace(const TMat3<R> rhs)
	{
	return rhs[0][0] + rhs[1][1] + rhs[2][2];
	}

template< typename R>
TVec3<R> Diag(const TMat3<R> rhs)
	{
	return TVec3<R>{rhs[0][0], rhs[1][1], rhs[2][2]};
	}

template< typename R>
R FrobeniusNorm(const TMat3<R> rhs)
	{
	auto Sq = [](R r) {return r * r; };
	return sqrt(
		Sq(rhs[0][0]) + Sq(rhs[0][1]) + Sq(rhs[0][2]) +
		Sq(rhs[1][0]) + Sq(rhs[1][1]) + Sq(rhs[1][2]) +
		Sq(rhs[2][0]) + Sq(rhs[2][1]) + Sq(rhs[2][2])
	);
	}

// generators
template< typename R>
TMat3<R> Zero3() // all zeros
	{
	return TMat3<R> { TVec3<R>{0, 0, 0}, { 0, 0, 0 }, { 0, 0, 0 }};
	}
template< typename R>
TMat3<R> Ones3() // all ones
	{
	return TMat3<R> { TVec3<R>{1, 1, 1}, { 1, 1, 1 }, { 1, 1, 1 }};
	}
template< typename R>
TMat3<R> Eye3() // identity matrix
	{
	return TMat3<R> { TVec3<R>{1, 0, 0}, { 0, 1, 0 }, { 0, 0, 1 }};
	}

// rotations, see https://en.wikipedia.org/wiki/Rotation_matrix
template< typename R, typename U>
TMat3<R> RotX3(U theta_rad)
	{
	R c = cos(theta_rad);
	R s = sin(theta_rad);
	return TMat3<R>{
		TVec3<R>{1, 0, 0},
			{0, c, -s},
			{0, s, c}
		};
	}

template< typename R, typename U>
TMat3<R> RotY3(U theta_rad)
	{
	R c = cos(theta_rad);
	R s = sin(theta_rad);
	return TMat3<R>{
		TVec3<R>{ c, 0, s },
			{ 0, 1, 0 },
			{-s, 0, c }
		};
	}

template< typename R, typename U>
TMat3<R> RotZ3(U theta_rad)
	{
	R c = cos(theta_rad);
	R s = sin(theta_rad);
	return TMat3<R>{
		TVec3<R>{ c, -s, 0 },
			{ s, c, 0 },
			{ 0, 0, 1 }
		};
	}

template< typename R, typename U>
TMat3<R> Rot3(const TVec3<R>& u, U theta_rad)
	{
	R c = cos(theta_rad);
	R s = sin(theta_rad);
	R mc = 1 - c;
	TVec3<R> uu = u / Norm(u);
	R ux = uu[0];
	R uy = uu[1];
	R uz = uu[2];
	return TMat3<R>{
		TVec3<R>{ c + ux*ux*mc,    ux*uy*mc - uz*s, ux*uz*mc + uy*s },
			{ uy*ux*mc + uz*s, c + uy*uy*mc,    uy*uz*mc - ux*s },
			{ uz*ux*mc - uy*s, uz*uy*mc + ux*s, c + uz*uz*mc }
		};
	}

template< typename R>  // solve Ax=b
TVec3<R> Solve(const TMat3<R>& A, const TVec3<R>& b)
	{
	return Inverse(A) * b;
	}


#endif // include guard