#ifndef __PhaseSpace_h
#define __PhaseSpace_h

#include "..\TM25ReadWrite\LinAlg3.h"
#include <tuple>
#include <limits>

template<typename R>
struct TLocR
	{
	using Real = R;
	Real x0_;
	Real x1_;
	};

template<typename R>
struct TDirR
	{
	using Real = R;
	Real k0_;
	Real k1_;
	};



// returns the point where the ray r0,k leaves the sphere at center with radius
// precondition: k is unit vector
template<typename R>
std::pair<bool, TVec3<R>> SphereIntersect(const TVec3<R>& center, R radius, const TVec3<R>& r0, const TVec3<R>& k);

// Phase space on a sphere given by center and radius.
// For location: Coordinates are given by stereographic projection 
// see https://en.wikipedia.org/wiki/Stereographic_projection
// Projection point = negative z axis pole of the sphere, i.e. origin - {0,0,1}*radius
// Projection plane = x-y plane through origin (not through the other pole)
// For direction: Stereographic projection preserves angles
// => back projecting the x-y directions onto the sphere creates orthogonal k0, k1.
// AreaScale returns the factor by which the area element on the sphere is magnified by back projection

template<typename R>
class TZAxisStereographicSphericalPhaseSpace
	{
	public:
		using Real = typename R;
		using TLoc = TLocR<Real>;
		using TDir = TDirR<Real>;
		using TV3 = TVec3<Real>;

		TZAxisStereographicSphericalPhaseSpace();
		TZAxisStereographicSphericalPhaseSpace(TV3 center, Real radius);

		TLoc Loc(const TV3& r, const TV3& k) const;
		TDir Dir(const TLoc& loc, const TV3& k) const;
		std::pair<bool, TLoc> LocIf(const TV3& r, const TV3& k) const;
		std::tuple<TLoc, TDir> LocDir(const TV3& r, const TV3& k) const;

		std::array<Real,4> PhaseSpacePoint(const TV3& r, const TV3& k) const;
		std::pair<bool, std::array<Real,4>> PhaseSpacePointIf(const TV3& r, const TV3& k) const;

		TV3 Loc_to_V3(const TLoc& loc) const;
		TV3 Dir_to_V3(const TLoc& loc, const TDir& dir) const;

		Real AreaScale(const TLoc& loc) const;

	private:
		TV3 center_;
		Real radius_;
	};


// Phase space on infinite plane, given by origin and the two e0, e1 axes
// Construction ensures that e0 and e1 are orthogonal unit vectors (e1 is rotated if needed)
// the two location coordinates are the e0 e1 components
// for direction, k0 and k1 are aligned with e0 and e1
// AreaScale always returns 1
template<typename R>
class TPlanarPhaseSpace
	{
	public:
		using Real = R;
		using TLoc = TLocR<Real>;
		using TDir = TDirR<Real>;
		using TV3 = TVec3<Real>;
		TPlanarPhaseSpace();
		TPlanarPhaseSpace(TV3 origin, TV3 e0, TV3 e1);
		TLoc Loc(const TV3& r, const TV3& k) const;
		TDir Dir(const TLoc& loc, const TV3& k) const;
		std::pair<bool, TLoc> LocIf(const TV3& r, const TV3& k) const;
		std::tuple<TLoc, TDir> LocDir(const TV3& r, const TV3& k) const;

		TV3 Loc_to_V3(const TLoc& loc) const;
		TV3 Dir_to_V3(const TLoc& loc, const TDir& dir) const;

		Real AreaScale(const TLoc& loc) const;

	private:
		TV3 origin_;
		TV3 e0_;
		TV3 e1_;
	};



void TestPhaseSpace();

// template definitions


//template<typename R>
//TPlanarPhaseSpace<R>::TPlanarPhaseSpace()
//	: origin_({ 0,0,0 }), e0_({ 1,0,0 }), e1_({ 0,1,0 })
//	{
//	};
//
//template<typename R>
//TPlanarPhaseSpace<R>::TPlanarPhaseSpace(TV3 origin, TV3 e0, TV3 e1)
//	: origin_(origin), e0_(UnitVector(e0)), e1_(UnitVector(e1)
//	{
//	e1_ = UnitVector(e1_ - (e0_*e1_) * e0_);
//	}


template<typename R>
std::pair<bool, TVec3<R>> SphereIntersect(const TVec3<R>& center, R radius, const TVec3<R>& r0, const TVec3<R>& k)
	{
	using Real = R;
	using TV3 = TVec3<Real>;
	TV3 r0c = r0 - center;
	Real b = k * r0c;
	Real c = Sqr(r0c) - radius * radius;
	Real D = b * b - c;
	Real NaN = std::numeric_limits<Real>::quiet_NaN();
	if (D < 0) return std::make_pair(false, TV3{ NaN,NaN,NaN });
	Real lambda = -b + sqrt(D);
	return std::make_pair(true, r0 + lambda * k);
	}



template<typename R>
TZAxisStereographicSphericalPhaseSpace<R>::TZAxisStereographicSphericalPhaseSpace()
	: center_{0,0,0}, radius_(1)
	{
	};

template<typename R>
TZAxisStereographicSphericalPhaseSpace<R>::TZAxisStereographicSphericalPhaseSpace(TV3 center, Real radius)
	: center_(center), radius_(abs(radius))
	{
	}

template<typename R>
typename TZAxisStereographicSphericalPhaseSpace<R>::TLoc TZAxisStereographicSphericalPhaseSpace<R>::Loc(
	const TV3 & r, const TV3 & k) const
	{
	std::pair<bool, TV3> p = SphereIntersect<Real>(center_, radius_, r, k);
	if (p.first == false)
		throw std::runtime_error("TZAxisStereographicSphericalPhaseSpace<R>::Loc: No intersection of ray with sphere");
	Real zpole = center_[2] - radius_;
	Real zp = p.second[2] - zpole;
	if (zp <= 0)
		throw std::runtime_error("TZAxisStereographicSphericalPhaseSpace<R>::Loc: Intersection is at - z pole");
	Real projFac = radius_ / zp;
	return TLoc{ projFac * (p.second[0] - center_[0]), projFac * (p.second[1] - center_[1]) };
	}

template<typename R>
typename TZAxisStereographicSphericalPhaseSpace<R>::TDir TZAxisStereographicSphericalPhaseSpace<R>::Dir(const TLoc & loc, const TV3 & k) const
	{
	// compute the two rotated unit vectors
	Real loc2 = loc.x0_*loc.x0_ + loc.x1_*loc.x1_;
	Real r2 = radius_ * radius_;
	if (loc2 < (10 * r2 * std::numeric_limits<Real>::epsilon())) // on or very near axis -> k is correctly oriented
		return { k[0],k[1] };
	Real c = (r2 - loc2) / (r2 + loc2); // cosine of rotation angle:
	Real mc = 1 - c;
	Real s = sqrt(loc2) * (1+c) / radius_;
	// 	s = sqrt(1 - c * c); // sine of rotation angle
	Real ux = -loc.x1_;// the rotation axis u
	Real uy = loc.x0_;
	Real unorm = sqrt(ux*ux + uy * uy);
	ux /= unorm;
	uy /= unorm; 
	// now apply rotation matrix, using uz = 0 on {1,0,0} and {0,1,0};
	Real uxuymc = ux * uy*mc;
	TV3 exRot{ c + ux * ux*mc, uxuymc,       -uy * s };
	TV3 eyRot{ uxuymc,         c + uy * uy*mc, ux*s };
	return { exRot * k, eyRot * k };
	}

template<typename R>
std::pair<bool, typename TZAxisStereographicSphericalPhaseSpace<R>::TLoc> TZAxisStereographicSphericalPhaseSpace<R>::LocIf(const TV3 & r, const TV3 & k) const
	{
	std::pair<bool, TV3> p = SphereIntersect<Real>(center_, radius_, r, k);
	if (p.first == false)
		return std::make_pair(false, TLoc{ std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::quiet_NaN() });
	Real zpole = center_[2] - radius_;
	Real zp = p.second[2] - zpole;
	if (zp <= 0)
		return std::make_pair(false, TLoc{ std::numeric_limits<Real>::quiet_NaN(), std::numeric_limits<Real>::quiet_NaN() });
	Real projFac = zp / radius_;
	return std::make_pair(true, TLoc{ projFac * (p.second[0] - center_[0]), projFac * (p.second[1] - center_[1]) });
	}

template<typename R>
std::tuple<
	typename TZAxisStereographicSphericalPhaseSpace<R>::TLoc, 
	typename TZAxisStereographicSphericalPhaseSpace<R>::TDir> 
	TZAxisStereographicSphericalPhaseSpace<R>::LocDir(const TV3 & r, const TV3 & k) const
	{
	TLoc loc = Loc(r, k);
	return std::make_tuple(loc, Dir(loc, k));
	}

template<typename R>
std::array<R, 4> TZAxisStereographicSphericalPhaseSpace<R>::PhaseSpacePoint(const TV3 & r, const TV3 & k) const
	{
	TLoc loc = Loc(r, k);
	TDir dir = Dir(loc, k);
	return { loc.x0_, loc.x1_, dir.k0_, dir.k1_ };
	}

template<typename R>
std::pair<bool, std::array<R, 4>> TZAxisStereographicSphericalPhaseSpace<R>::PhaseSpacePointIf(const TV3 & r, const TV3 & k) const
	{
	auto locIf = LocIf(r, k);
	if (locIf.first == false)
		{
		Real mynan = std::numeric_limits<Real>::quiet_nan();
		return std::make_pair(false, std::array<Real, 4>{mynan, mynan, mynan, mynan});
		}
	TDir dir = Dir(locIf.second, k);
	return std::make_pair(true, std::array<Real, 4>{locIf.second.x0_, locIf.second.x1_, dir.k0_, dir.k1_});
	}

template<typename R>
typename TZAxisStereographicSphericalPhaseSpace<R>::TV3 TZAxisStereographicSphericalPhaseSpace<R>::Loc_to_V3(const TLoc & loc) const
	{
	// compute the two rotated unit vectors
	Real loc2 = loc.x0_*loc.x0_ + loc.x1_*loc.x1_;
	Real r2 = radius_ * radius_;
	Real c = (r2 - loc2) / (r2 + loc2); // cosine of rotation angle:
	Real fac = 1 + c;
	return center_ + TV3{ fac * loc.x0_, fac * loc.x1_, c * radius_};
	}

template<typename R>
typename TZAxisStereographicSphericalPhaseSpace<R>::TV3 TZAxisStereographicSphericalPhaseSpace<R>::Dir_to_V3(const TLoc & loc, const TDir & dir) const
	{
	// compute sine and cosine of back rotation angle
	Real loc2 = loc.x0_*loc.x0_ + loc.x1_*loc.x1_;
	Real r2 = radius_ * radius_;
	// make sure dir is inside unit circle * 0.99
	constexpr Real _099 = static_cast<Real>(0.99);
	Real k0 = dir.k0_;
	Real k1 = dir.k1_;
	Real test = k0*k0 + k1*k1;
	if (test > _099)
		{
		Real fac = _099 / sqrt(test);
		k0 *= fac;
		k1 *= fac;
		}
	Real k2 = sqrt(1 - (k0*k0 + k1*k1));
	TV3 k{k0, k1, k2 };
	if (loc2 < (10 * r2 * std::numeric_limits<Real>::epsilon())) // on or very near axis -> k is correctly oriented
		return k;
	Real c = (r2 - loc2) / (r2 + loc2); // cosine of rotation angle:
	Real mc = 1 - c;
	Real s = sqrt(loc2) * (1 + c) / radius_;
	// 	s = sqrt(1 - c * c); // sine of rotation angle
	Real ux = -loc.x1_;// the rotation axis u
	Real uy = loc.x0_;
	Real unorm = sqrt(ux*ux + uy * uy);
	ux /= unorm;
	uy /= unorm;
	Real uxuymc = ux * uy*mc;
	// the rotation matrix from 3D to tangent disk
	TMat3<Real> A{	TV3{c + ux*ux*mc, uxuymc, uy*s},
					TV3{uxuymc, c + uy*uy*mc, -ux * s},
					TV3{-uy*s, ux*s, c} };
	//Real d = Det(A);
	//auto testEye = A * Transpose(A);
	// the columns of A are the unit vectors of k on tangent disk
	TV3 rv = A * k;
	//auto test1 = Norm(rv);
	//TV3 pointOnSphere{ loc.x0_ * (1 + c), loc.x1_*(1 + c), c * radius_ };
	//TV3 testSouthPole = Transpose(A) * pointOnSphere;
	return rv;
	}

template<typename R>
R TZAxisStereographicSphericalPhaseSpace<R>::AreaScale(const TLoc & loc) const
	{
	auto sqr = [](R x) {return x * x; };
	Real loc2 = sqr(loc.x0_) + sqr(loc.x1_);
	Real r2 = sqr(radius_);
	return 4 / sqr(1 + loc2 / r2);
	}

#endif
// include guard
