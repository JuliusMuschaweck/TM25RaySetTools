#ifndef __PHASESPACEGRID_H
#define __PHASESPACEGRID_H

#include<array>
#include<vector>
#include "PhaseSpace.h"

// Regular 4D phase space grid
template<typename R>
class TPhaseSpaceGrid
	{
	public:
		using Real = R;
		using TLoc = TLocR<Real>;
		using TDir = TDirR<Real>;
		using T4DPoint = std::array<Real, 4>;
		using T4DIndex = std::array<size_t, 4>;
		using TBox = std::pair<T4DPoint, T4DPoint>;

		TPhaseSpaceGrid(const T4DPoint& bbmin, const T4DPoint& bbmax, const T4DIndex& nn);

		Real  Value(const T4DIndex& nn) const;
		Real& Value(const T4DIndex& nn);

//	private:
		// the bounding box values
		T4DPoint bbmin_;
		T4DPoint bbmax_;
		// the number of points in each dimension
		size_t nx0_;
		size_t nx1_;
		size_t nk0_;
		size_t nk1_;
		// the values on the grid, in row major x0, x1, k0, k1 index
		std::vector<Real> v_;
	};

struct TInterpolateRaySetData;

TPhaseSpaceGrid<float> FillDefaultPhaseSpaceGrid(const TInterpolateRaySetData& raySetData, const size_t nBins);

TPhaseSpaceGrid<float> FillPhaseSpaceGrid(const TPhaseSpaceGrid<float>::TBox& boundingBox,
		const TInterpolateRaySetData& raySetData, const TPhaseSpaceGrid<float>::T4DIndex& nBins);


// **********************************************************************
// template definitions

template<typename R>
TPhaseSpaceGrid<R>::TPhaseSpaceGrid(const T4DPoint& bbmin, const T4DPoint& bbmax, const T4DIndex& nn)
	: bbmin_(bbmin), bbmax_(bbmax), nx0_(nn[0]), nx1_(nn[1]), nk0_(nn[2]), nk1_(nn[3])
	{
	v_.resize(nx0_ * nx1_ * nk0_ * nk1_);
	std::fill(v_.begin(), v_.end(), static_cast<R>(0.0));
	}

template<typename R>
R TPhaseSpaceGrid<R>::Value(const T4DIndex& nn) const
	{
	#ifndef NDEBUG
	if (nn[0] >= nx0_ || nn[1] >= nx1_ || nn[2] >= nk0_ || nn[3] >= nk1_)
		throw std::runtime_error("TPhaseSpaceGrid<R>::Real: index out of bounds");
	#endif // !NDEBUG
	size_t i = nn[0] + nx0_ * (nn[1] + nx1_ * (nn[2] + nk0_ * nn[3]));
	return v_[i];
	}

template<typename R>
R&  TPhaseSpaceGrid<R>::Value(const T4DIndex& nn)
	{
	#ifndef NDEBUG
	if (nn[0] >= nx0_ || nn[1] >= nx1_ || nn[2] >= nk0_ || nn[3] >= nk1_)
		throw std::runtime_error("TPhaseSpaceGrid<R>::Real: index out of bounds");
	#endif // !NDEBUG
	size_t i = nn[0] + nx0_ * (nn[1] + nx1_ * (nn[2] + nk0_ * nn[3]));
	return v_[i];
	}

#endif
