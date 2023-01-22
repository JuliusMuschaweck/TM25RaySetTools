#include "InterpolateRaySetData.h"
#include <sstream>



TInterpolateRaySetData::TThreeVecs TInterpolateRaySetData::ComputeLumRange(size_t ibegin, size_t iend, KDTree::Def::TIdx nNeighbors) const
	{
	TVec luminances = TVec(iend - ibegin);
	TVec volumes = TVec(iend - ibegin);
	TVec cellFluxes = TVec(iend - ibegin);
	std::stringstream s;
	if (ibegin == 0)
		{
		s << "ComputeLumRange(" << ibegin << ", " << iend << ")\n";
		safeout_.Write(s.str());
		};
	for (size_t i = ibegin; i < iend; ++i)
		{
		if (ibegin == 0)
			{
			size_t nraysPercent = (iend > 100) ? (iend / 100) : 1;
			if ((iend > 100) && (i % nraysPercent == 0))
				{
				s.str("");
				s << round(double(i) / iend * 100.0) << "% ";
				safeout_.Write(s.str());
				}
			}

		const KDTree::TKDTree::TPointIdx pi{ static_cast<KDTree::Def::TIdx>(i) };
		const KDTree::TKDTree::TNearestNeighbors inbs = kdtree_->NearestNeighborsOfPoint(pi, nNeighbors);
		// compute total volume and avg luminance of neighbors -- no write except local
		const float vol = kdtree_->TotalVolume(inbs.i_nodes_);
		float flux = 0;
		for (auto& ii : inbs.i_points_)
			flux += rayFluxes_[ii.pi_];
		const float avgLuminance = flux / vol;
		// write out luminances, volumes and cellFluxes -- RACE DANGER
		luminances[i - ibegin] = avgLuminance;
		const KDTree::TKDTree::TNodeIdx ni = kdtree_->NodeIndex(pi);
		const float iVolume = (kdtree_->Node(ni)).Volume();
		volumes[i - ibegin] = iVolume;
		cellFluxes[i - ibegin] = (iVolume * avgLuminance);
		}
	s.str("");
	s << "ComputeLumRange(" << ibegin << ", " << iend << ") done\n";
	safeout_.Write(s.str());
	return TThreeVecs(std::move(luminances), std::move(volumes), std::move(cellFluxes));
	}

void TInterpolateRaySetData::SetTotalFlux(float newTotalFlux)
	{
	double fac = static_cast<double>(newTotalFlux) / static_cast<double>(totalFlux_);
	for (auto& f : rayFluxes_)
		f = static_cast<float>(static_cast<double>(f) * fac);
	for (auto& l : luminances_)
		l = static_cast<float>(static_cast<double>(l) * fac);
	for (auto& c : cellFluxes_)
		c = static_cast<float>(static_cast<double>(c) * fac);
	totalFlux_ = newTotalFlux;
	}

// sort v, but return index of permutation instead of modifying v
// postcondition: i < j => v[rv[i]] <= v[rv[j]]
template<typename Compare = decltype(std::less())>
std::vector<size_t> IndexSort(const std::vector<float>& v, Compare comp = std::less())
	{
	std::vector<size_t> rv(v.size());
	std::iota(rv.begin(), rv.end(), 0);
	auto pred = [&v, comp](size_t lhs, size_t rhs) -> bool {return comp(v[lhs], v[rhs]); };
	std::sort(rv.begin(), rv.end(), pred);
	return rv;
	}

#ifdef NDEBUG
#pragma optimize("",off)
#endif
TInterpolateRaySetData::TCharacteristicCurve TInterpolateRaySetData::CharacteristicCurve() const
	{
	// sort descending
	std::vector<size_t> idx = IndexSort(luminances_, std::greater());
	TCharacteristicCurve rv;
	size_t nCells = idx.size();
	rv.etendue_.reserve(nCells);
	rv.luminance_.reserve(nCells);
	for (size_t i : idx)
		{
		rv.etendue_.push_back(volumes_[i]);
		rv.luminance_.push_back(luminances_[i]);
		}
	return rv;
	}

void TInterpolateRaySetData::WriteCharacteristicCurve(const std::string& fn, const TCharacteristicCurve& cc) const
	{
	TM25::TWriteFile f(fn);
	f.Write<uint64_t>(static_cast<uint64_t>(cc.etendue_.size()));
	f.WriteVector(cc.etendue_);
	f.WriteVector(cc.luminance_);
	}

// helper function for skewness distribution. 
// precondition: v ascending (not checked!), v.size() >= nBins + 1, nBins > 0
// returns nBins + 1 indices into v such that the half open intervals have approx. same width in v
// rv starts with 0, ends with v.size()-1, strictly ascending (no empty bins)
// more precisely: Given bin size: d = (v.back()-v.front())/nBins.
// for any i such that 0 < i <= nBins, v[rv[i]] is the smallest/first element in v such that v[rv[i]] >= i * d
// except if we run out of remaining values in v, then they come one by one
std::vector<size_t> BinIndices(const std::vector<float>& v, size_t nBins)
	{
	if (nBins == 0)
		throw std::runtime_error("BinIndices: nBins == 0");
	size_t nv = v.size();
	if (nv < nBins + 1)
		throw std::runtime_error("BinIndices: too many bins / too few values");
	if (nBins == 1)
		return std::vector<size_t> {0, nv};
	std::vector<size_t> rv;
	rv.push_back(0);
	double d = (v.back() - v.front()) / nBins;
	size_t vpos = 0;
	for (size_t i = 1; i < nBins; ++i)
		{
		double nextBoundary = i * d;
		++vpos; // advance at least one
		while (v[vpos] < nextBoundary)
			{
			size_t remain_v = nv - vpos;
			size_t remain_rv = nBins + 1 - i;
			if (remain_v <= remain_rv)
				break;
			++vpos;
			}
		rv.push_back(vpos);
		}
	rv.push_back(nv - 1);
	return rv;
	}

TInterpolateRaySetData::TSkewnessDistribution::BinType TInterpolateRaySetData::BinTypeFromString(std::string s) const
	{
	if (s == "sameSkewness")
		return TSkewnessDistribution::BinType::sameSkewness;
	if (s == "sameEtendue")
		return TSkewnessDistribution::BinType::sameEtendue;
	if (s == "sameFlux")
		return TSkewnessDistribution::BinType::sameFlux;
	throw std::runtime_error("TInterpolateRaySetData::BinTypeFromString: unknown bin type: " + s);
	}

TInterpolateRaySetData::TSkewnessDistribution TInterpolateRaySetData::SkewnessDistribution_z_axis(size_t nBins, TSkewnessDistribution::BinType binType, const TM25::TTM25RaySet& rs) const
	{// ray r + lambda * u, skewness s = n r dot (u cross e_z), here n = 1
	auto s = [](const TVec3f& r, const TVec3f& u) {return r[0] * u[1] - r[1] * u[0]; };
	size_t nRays = rs.NRays();
	std::vector<float> raw_skewness;
	raw_skewness.reserve(nRays);
	for (size_t i = 0; i < nRays; ++i)
		{
		std::tuple<TVec3f, TVec3f, float> iray = rs.RayLocDirFlux(i);
		raw_skewness.push_back(s(std::get<0>(iray), std::get<1>(iray)));
		}
	std::vector<size_t> s_idx = IndexSort(raw_skewness);
	std::vector<float> sorted_skewness;
	sorted_skewness.reserve(nRays);
	for (size_t i : s_idx)
		{
		sorted_skewness.push_back(raw_skewness[i]);
		}
	{ // free memory: clear won't do it! swap with empty .. 
	std::vector<float> empty;
	raw_skewness.swap(empty);
	}  // empty goes out of scope and is destroyed
	double acc_etendue = 0;
	double acc_flux = 0;
	std::vector<float> accumulated_sorted_etendue;
	accumulated_sorted_etendue.reserve(nRays);
	std::vector<float> accumulated_sorted_flux;
	accumulated_sorted_flux.reserve(nRays);
	for (size_t i : s_idx)
		{
		accumulated_sorted_etendue.push_back(static_cast<float>(acc_etendue += volumes_[i]));
		accumulated_sorted_flux.push_back(static_cast<float>(acc_flux += rayFluxes_[i]));
		}
	// now we have three arrays: sorted_skewness, accumulated_sorted_etendue, accumulated_sorted_flux.
	std::vector<size_t> bin_idx;
	switch (binType)
		{
		case TSkewnessDistribution::BinType::sameSkewness:
			bin_idx = BinIndices(sorted_skewness, nBins); break;
		case TSkewnessDistribution::BinType::sameEtendue:
			bin_idx = BinIndices(accumulated_sorted_etendue, nBins); break;
		case TSkewnessDistribution::BinType::sameFlux:
			bin_idx = BinIndices(accumulated_sorted_flux, nBins); break;
		default:
			throw std::runtime_error("TInterpolateRaySetData::SkewnessDistribution_z_axis: unknown bin type");
		}
	TInterpolateRaySetData::TSkewnessDistribution rv;
	rv.axis_direction_ = TVec3f{ 0,0,1 };
	rv.axis_point_ = TVec3f{ 0,0,0 };
	for (size_t i = 0; i <= nBins; ++i)
		rv.skewness_.push_back(sorted_skewness[bin_idx[i]]);
	for (size_t i = 0; i < nBins; ++i)
		{
		float dU = accumulated_sorted_etendue[bin_idx[i + 1]] - accumulated_sorted_etendue[bin_idx[i]];
		float dPhi = accumulated_sorted_flux[bin_idx[i + 1]] - accumulated_sorted_flux[bin_idx[i]];
		float ds = sorted_skewness[bin_idx[i + 1]] - sorted_skewness[bin_idx[i]];
		rv.dU_ds_.push_back(dU / ds);
		rv.dPhi_ds_.push_back(dPhi / ds);
		}
	return rv;
	}

//format: 4 byte floats, 3 floats for point, 3 floats for direction,  4 byte int for nBins
		// nBins+1 (!!) floats for skewness, nBins floats for dU_ds, nBins floats for dPhi_ds
void TInterpolateRaySetData::WriteSkewnessDistribution_z_axis(const std::string& fn, const TSkewnessDistribution& cc) const
	{
	TM25::TWriteFile f(fn);
	f.WriteArray(cc.axis_point_);
	f.WriteArray(cc.axis_direction_);
	uint32_t nBins = static_cast<uint32_t>(cc.skewness_.size() - 1);
	f.Write(nBins);
	f.WriteVector(cc.skewness_);
	f.WriteVector(cc.dU_ds_);
	f.WriteVector(cc.dPhi_ds_);
	}

#ifdef NDEBUG
#pragma optimize("",on)
#endif


