#ifndef __PhaseSpace_h
#define __PhaseSpace_h

#include <LinAlg3.h>
#include <tuple>
#include <limits>
#include <stdexcept>
#include <vector>

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




// Phase space on infinite plane perpendicular to z axis
// for direction, k0 and k1 are aligned with x and y
// AreaScale always returns 1
// only rays with kz > 0 or (kz == 0 and z == z0) are allowed
template<typename R>
class TPlanarZPhaseSpace
	{
	public:
		using Real = R;
		using TLoc = TLocR<Real>;
		using TDir = TDirR<Real>;
		using TV3 = TVec3<Real>;
		TPlanarZPhaseSpace(Real z0 = 0) : z0_(z0) {};
		TLoc Loc(const TV3& r, const TV3& k) const;  // precondition: k is unit vector
		TDir Dir(const TLoc& loc, const TV3& k) const; // precondition: k is unit vector
		std::pair<bool, TLoc> LocIf(const TV3& r, const TV3& k) const;
		std::tuple<TLoc, TDir> LocDir(const TV3& r, const TV3& k) const;

		std::array<Real, 4> PhaseSpacePoint(const TV3& r, const TV3& k) const;
		std::pair<bool, std::array<Real, 4>> PhaseSpacePointIf(const TV3& r, const TV3& k) const;

		TV3 Loc_to_V3(const TLoc& loc) const;
		TV3 Dir_to_V3(const TLoc& loc, const TDir& dir) const;

		Real AreaScale(const TLoc& loc) const;

	private:
		Real z0_;
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



// Phase space on infinite cylinder around the z axis
// x0 is aligned with z, x1 in [- r pi, r pi], x1=0 at x axis 
// for direction, k0 is aligned with z
// AreaScale always returns 1
template<typename R>
class TCylinderZPhaseSpace
	{
	public:
		using Real = R;
		using TLoc = TLocR<Real>;
		using TDir = TDirR<Real>;
		using TV3 = TVec3<Real>;
		TCylinderZPhaseSpace(Real radius = 1) : radius_(radius) {};
		TLoc Loc(const TV3& r, const TV3& k) const;
		TDir Dir(const TLoc& loc, const TV3& k) const;
		std::pair<bool, TLoc> LocIf(const TV3& r, const TV3& k) const;
		std::tuple<TLoc, TDir> LocDir(const TV3& r, const TV3& k) const;

		std::array<Real, 4> PhaseSpacePoint(const TV3& r, const TV3& k) const;
		std::pair<bool, std::array<Real, 4>> PhaseSpacePointIf(const TV3& r, const TV3& k) const;

		TV3 Loc_to_V3(const TLoc& loc) const;
		TV3 Dir_to_V3(const TLoc& loc, const TDir& dir) const;

		Real AreaScale(const TLoc& loc) const;

	private:
		Real radius_;
	};





// Phase space on infinite plane, given by origin and the two e0, e1 axes
// Construction ensures that e0 and e1 are orthogonal unit vectors (e1 is rotated if needed)
// the two location coordinates are the e0 e1 components
// for direction, k0 and k1 are aligned with e0 and e1
// AreaScale always returns 1
// NOT YET IMPLEMENTED
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




#include "PhaseSpaceImpl.h"

#endif