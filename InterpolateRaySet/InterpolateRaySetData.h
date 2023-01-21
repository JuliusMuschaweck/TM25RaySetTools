#ifndef __INTERPOLATERAYSETDATA_H
#define __INTERPOLATERAYSETDATA_H

#include <memory>
#include <vector>
#include <LinAlg3.h>
#include "InterpolateRaySet_IO.h"
#include "KDTree.h"
#include <numeric>

#define MULTITHREAD
#ifdef MULTITHREAD
#include<thread>
#include<future>
#include<mutex>
#endif


struct TInterpolateRaySetData
	{
	template<typename PhaseSpace>
	void Init(const PhaseSpace& ps, const TM25::TTM25RaySet& rs, KDTree::Def::TIdx nNeighbors, TLogPlusCout& info);

	void SetTotalFlux(float newTotalFlux);

	std::unique_ptr<KDTree::TKDTree> kdtree_;
	std::vector<float> rayFluxes_;
	std::vector<float> volumes_;
	std::vector<float> luminances_;
	std::vector<float> cellFluxes_;
	float totalVolume_ = 0.0f;
	float totalFlux_ = 0.0f;
	float avgLuminance_ = std::numeric_limits<float>::quiet_NaN();
	mutable TThreadSafe_stdcout safeout_;

	// sort cells according to luminance, result is the characteristic curve
	struct TCharacteristicCurve // see doi:10.1117/12.615874, "Characterization of the thermodynamic quality of light sources"
		{
		std::vector<float> etendue_;
		std::vector<float> luminance_;
		};
	TCharacteristicCurve CharacteristicCurve() const;
	//format: 8 byte int for nCells, 2 by n matrix of 4 byte floats, row major, first row is etendue, second row is luminance.
	void WriteCharacteristicCurve(const std::string& fn, const TCharacteristicCurve& cc) const;

	// compute the skewness distribution of etendue and flux, w.r.t. a certain axis of rotational symmetry
	// put cells into bins of equal size, choose quantity of which equal size is desired
	struct TSkewnessDistribution // see doi:10.1364/JOSAA.14.002855 "Performance limitations of rotationally symmetric nonimaging devices"
		{
		TVec3f axis_point_{ 0.0f, 0.0f, 0.0f };
		TVec3f axis_direction_{ 0.0f, 0.0f, 1.0f };
		enum class BinType { sameSkewness, sameEtendue, sameFlux };
		std::vector<float> skewness_; // nBins + 1: limits of intervals for bar graph. Compute mean of adjacent value pairs to obtain nBins values for easy plotting
		std::vector<float> dU_ds_; // nBins
		std::vector<float> dPhi_ds_; // nBins
		};
	TSkewnessDistribution::BinType BinTypeFromString(std::string s) const;
	TSkewnessDistribution SkewnessDistribution_z_axis(size_t nBins, TSkewnessDistribution::BinType binType, const TM25::TTM25RaySet& rs) const;
	TSkewnessDistribution SkewnessDistribution(TVec3f axis_point, TVec3f axis_direction, size_t nBins, TSkewnessDistribution::BinType binType, const TM25::TTM25RaySet& rs) const;
	//format: 4 byte floats, 3 floats for point, 3 floats for direction,  4 byte int for nBins
	// nBins+1 (!!) floats for skewness, nBins floats for dU_ds, nBins floats for dPhi_ds
	void WriteSkewnessDistribution_z_axis(const std::string& fn, const TSkewnessDistribution& cc) const;

	private:
		template<typename PhaseSpace>
		void CreateTree(const PhaseSpace& ps, const TM25::TTM25RaySet& rs, std::ostream& info);

		using TVec = std::vector<float>;
		using TThreeVecs = std::tuple<TVec, TVec, TVec>;
		TThreeVecs ComputeLumRange(size_t ibegin, size_t iend, KDTree::Def::TIdx nNeighbors) const;
	};




//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

template<typename PhaseSpace>
void TInterpolateRaySetData::Init(const PhaseSpace& ps, const TM25::TTM25RaySet& rs, KDTree::Def::TIdx nNeighbors, TLogPlusCout& info)
	{
	CreateTree(ps, rs, info);
	info << "averaging luminance over " << nNeighbors << " nearest neighbors\n ";
	size_t nRays = rs.NRays();
	bool multithreaded = true; 
	if (multithreaded)
		{
		size_t min_chunkSize = 1000;
		size_t max_nChunks = nRays / min_chunkSize + 1;
		size_t nCores = std::thread::hardware_concurrency();
		size_t nThreads = std::min(nCores, max_nChunks);
		size_t chunkSize = nRays / nThreads + 1;
		info << "nRays = " << nRays << ", max_nChunks = " << max_nChunks << ", nCores = " << nCores
			<< ", nThreads = " << nThreads << ", chunkSize = " << chunkSize << "\n";
		if (nThreads == 1)
			{
			TThreeVecs tmp = ComputeLumRange(0, nRays, nNeighbors);
			luminances_ = std::move(std::get<0>(tmp));
			volumes_ = std::move(std::get<1>(tmp));
			cellFluxes_ = std::move(std::get<2>(tmp));
			}
		else
			{
			std::vector<std::future<TThreeVecs> > futures(nThreads - 1);
			size_t ibegin = 0;
			size_t iend = ibegin + chunkSize;
			info << "starting " << nThreads << " luminance threads\n";
			for (size_t i = 0; i < nThreads - 1; ++i)
				{
				iend = ibegin + chunkSize;
				//				futures[i] = std::async(std::launch::async, [&]() -> TThreeVecs {return this->ComputeLumRange(ibegin, iend, nNeighbors); });
				futures[i] = std::async(std::launch::async, &TInterpolateRaySetData::ComputeLumRange, this, ibegin, iend, nNeighbors);
				ibegin = iend;
				}
			TThreeVecs lastRes = ComputeLumRange(ibegin, nRays, nNeighbors);
			ibegin = 0;
			luminances_.resize(nRays);
			volumes_.resize(nRays);
			cellFluxes_.resize(nRays);
			info << "collecting thread results\n";
			for (size_t i = 0; i < nThreads - 1; ++i)
				{
				TThreeVecs resi = futures[i].get();
				const TVec& ilum = std::get<0>(resi);
				const TVec& ivol = std::get<1>(resi);
				const TVec& iflux = std::get<2>(resi);
				size_t iend = ibegin + chunkSize;
				info << " from thread " << i << ", ibegin = " << ibegin << "\n";
				std::copy(ilum.begin(), ilum.end(), luminances_.begin() + ibegin);
				std::copy(ivol.begin(), ivol.end(), volumes_.begin() + ibegin);
				std::copy(iflux.begin(), iflux.end(), cellFluxes_.begin() + ibegin);
				ibegin = iend;
				}
			const TVec& ilum = std::get<0>(lastRes);
			const TVec& ivol = std::get<1>(lastRes);
			const TVec& iflux = std::get<2>(lastRes);
			std::copy(ilum.begin(), ilum.end(), luminances_.begin() + ibegin);
			std::copy(ivol.begin(), ivol.end(), volumes_.begin() + ibegin);
			std::copy(iflux.begin(), iflux.end(), cellFluxes_.begin() + ibegin);
			}
		}
	else // single threaded
		{
		for (size_t i = 0; i < nRays; ++i)
			{
			// show progress
			size_t nraysPercent = (nRays > 100) ? (nRays / 100) : 1;
			if ((nRays > 100) && (i % nraysPercent == 0))
				info << round(double(i) / nRays * 100.0) << "% ";
			// compute nearest neighbors of current point -- no write except local
			KDTree::TKDTree::TPointIdx pi{ static_cast<KDTree::Def::TIdx>(i) };
			KDTree::TKDTree::TNearestNeighbors inbs = kdtree_->NearestNeighborsOfPoint(pi, nNeighbors);
			// compute total volume and avg luminance of neighbors -- no write except local
			float vol = kdtree_->TotalVolume(inbs.i_nodes_);
			float flux = 0;
			for (auto i : inbs.i_points_)
				flux += rayFluxes_[i.pi_];
			float avgLuminance = flux / vol;
			// write out luminances, volumes and cellFluxes -- RACE DANGER
			luminances_.push_back(avgLuminance);
			KDTree::TKDTree::TNodeIdx ni = kdtree_->NodeIndex(pi);
			float iVolume = (kdtree_->Node(ni)).Volume();
			volumes_.push_back(iVolume);
			cellFluxes_.push_back(iVolume * avgLuminance);
			}
		}
	info << "\n";
	totalVolume_ = std::accumulate(volumes_.begin(), volumes_.end(), 0.0f);
	totalFlux_ = std::accumulate(rayFluxes_.begin(), rayFluxes_.end(), 0.0f);
	avgLuminance_ = totalFlux_ / totalVolume_;
	}


template<typename PhaseSpace>
void TInterpolateRaySetData::CreateTree(const PhaseSpace& ps, const TM25::TTM25RaySet& rs, std::ostream& info)
	{
	info << "computing phase space points\n";
	KDTree::Def::TKDPoints pspoints;
	size_t nRays = rs.NRays();
	for (size_t i = 0; i < nRays; ++i)
		{
		using std::get;
		std::tuple<TVec3f, TVec3f, float> iray = rs.RayLocDirFlux(i);
		KDTree::Def::TKDPoint phasespacepoint = ps.PhaseSpacePoint(get<0>(iray), get<1>(iray));
		pspoints.push_back(phasespacepoint);
		rayFluxes_.push_back(get<2>(iray));
		}
	info << "creating KD tree\n";
	kdtree_.reset(new KDTree::TKDTree(std::move(pspoints)));
	kdtree_->CreateTree();
	info << "shrinking edge nodes\n";
	kdtree_->ShrinkEdgeNodes();
	info << "checking tree consistency\n";
	kdtree_->CheckConsistency();
	}



#endif
