// InterpolateRaySet.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <TM25.h>
#include <ZemaxBinary.h>
#include <TranslateZemax.h>
#include "KDTree.h"
#include <iostream>
#include "PhaseSpace.h"
#include <algorithm>
#include <numeric>
#include <memory>
#include <CfgFile.h>
#include <WriteFile.h>

#define MULTITHREAD
#ifdef MULTITHREAD
#include<thread>
#include<future>
#include<mutex>
#endif

int Test()
	{
	TestPhaseSpace();
	//if (KDTree::Def::dim == 2)
	//	KDTree::TestKDTree2D("TestKDTree.m");
	//if (KDTree::Def::dim == 4)
	//	KDTree::TestKDTree4D();
	TestConfiguration();
	return 0;
	}

using std::cout;
using std::endl;


void DisplayRaySetInfo(std::ostream& info)
	{
	}

class TRaySetControlSection : public TSection
	{
	public:
		TRaySetControlSection() : TSection("RaySetControl") {};
		virtual void AddAllowedValues()
			{
			values_.insert({ "inputRayFileName",	MakeDefaultValueTokenSequence<Token::string>(std::string("inputRayFileName missing")) });
			values_.insert({ "inputRayFileFormat",	MakeDefaultValueTokenSequence<Token::string>(std::string("TM25")) });
			values_.insert({ "logFileName",			MakeDefaultValueTokenSequence<Token::string>(std::string("InterpolateRaySet.log")) });
			values_.insert({ "consoleOutput",		MakeDefaultValueTokenSequence<Token::boolean>(true) });
			values_.insert({ "scrambleInput",		MakeDefaultValueTokenSequence<Token::boolean>(true) });
			values_.insert({ "selectByMaxNumber",	MakeEmptyTokenSequence() });
			values_.insert({ "doDiagnostics",		MakeDefaultValueTokenSequence<Token::boolean>(true) });
			values_.insert({ "restrictToRelVirtualFocusDistance",MakeEmptyTokenSequence() });
			values_.insert({ "restrictToKz",		MakeEmptyTokenSequence() });
			values_.insert({ "nOutputRays",			MakeDefaultValueTokenSequence<Token::integer>(1) });
			values_.insert({ "nNeighbors",			MakeDefaultValueTokenSequence<Token::integer>(10) });
			values_.insert({ "outputRayFileName",	MakeDefaultValueTokenSequence<Token::string>(std::string("tmp.TM25RAY")) });
			values_.insert({ "outputRayFileFormat",	MakeDefaultValueTokenSequence<Token::string>(std::string("TM25")) });
			values_.insert({ "phaseSpaceType",		MakeDefaultValueTokenSequence<Token::identifier>(std::string("VirtualFocusZSphere")) });
			values_.insert({ "ZPlane_z",			MakeDefaultValueTokenSequence<Token::real>(0.0) });
			values_.insert({ "ZCylinderRadius",		MakeDefaultValueTokenSequence<Token::real>(1.0) });
			values_.insert({ "characteristicCurveFileName",	MakeEmptyTokenSequence() });
			values_.insert({ "skewness_z_FileName",	MakeEmptyTokenSequence() });
			values_.insert({ "skewness_nBins",	MakeDefaultValueTokenSequence<Token::integer>(100) });
			values_.insert({ "skewness_binType",	MakeDefaultValueTokenSequence<Token::identifier>(std::string("sameSkewness")) });
			values_.insert({ "setTotalFlux",		MakeEmptyTokenSequence() });
			}
		virtual void AddAllowedKeywords()
			{
			keywords_.insert("VirtualFocusZSphere");
			keywords_.insert("ZPlane");
			keywords_.insert("ZCylinder");

			keywords_.insert("sameSkewness");
			keywords_.insert("sameEtendue");
			keywords_.insert("sameFlux");
			}; // default: empty

	};

class TInterpolateRaySetCfg : public TConfiguration
	{
	public:
		TInterpolateRaySetCfg()
			: TConfiguration{}
			{
			TRaySetControlSection rcs;
			rcs.InitAllowed();
			AddSection(std::move(rcs));
			}
	};

TM25::TTM25RaySet ReadRaySet(const TInterpolateRaySetCfg& cfg, std::ostream& info)//const std::string& fn, std::ostream& info)
	{
	const TSection& rsc = cfg.Section("RaySetControl");
	TM25::TTM25RaySet rs;
	const std::string format = rsc.String("inputRayFileFormat");
	const std::string fn = rsc.String("inputRayFileName");
	if (format.compare("TM25") == 0)
		{
		info << "Reading TM25 ray file " << fn << endl;
		rs.Read(fn);
		// rs.Make_k_unit();
		if (rs.Warnings().empty())
			info << "No warnings" << endl;
		for (auto w : rs.Warnings())
			{
			info << "Warning: " << w << "\n";
			}
		}
	else if (format.compare("ZemaxBinary") == 0)
		{
		info << "Reading Zemax binary ray file " << fn << endl;
		TM25::TZemaxRaySet zemaxRaySet(fn);
		rs = TM25::ZemaxBinaryToTM25(zemaxRaySet);
		}
	return rs;
	}

void WriteRaySet(TM25::TTM25RaySet& rs, const TInterpolateRaySetCfg& cfg, std::ostream& info)
	{
	const TSection& rsc = cfg.Section("RaySetControl");
	const std::string format = rsc.String("outputRayFileFormat");
	const std::string fn = rsc.String("outputRayFileName");
	if (format.compare("TM25") == 0)
		{
		info << "Writing TM25 ray file " << fn << endl;
		rs.Write(fn);
		if (rs.Warnings().empty())
			info << "No warnings" << endl;
		for (auto w : rs.Warnings())
			{
			info << "Warning: " << w << "\n";
			}
		}
	else if (format.compare("ZemaxBinary") == 0)
		{
		info << "Writing  Zemax binary ray file " << fn << endl;
		TM25::TZemaxRaySet zemaxRaySet = TM25ToZemaxBinary(rs);
		zemaxRaySet.Write(fn);
		}
	}

// create an ostream that sends output to many other ostreams
class ComposeStream : public std::ostream
	{
	struct ComposeBuffer : public std::streambuf
		{
		void addBuffer(std::streambuf* buf)
			{
			bufs.push_back(buf);
			}
		virtual int overflow(int c)
			{
			auto do_sputc = [](std::streambuf* buf, int c)
				{
				buf->sputc(c);
				};
			for (auto buf : bufs)
				{
				do_sputc(buf, c);
				}
			return c;
			}
		private:
			std::vector<std::streambuf*>    bufs;

		};
	ComposeBuffer myBuffer;
	public:
		ComposeStream()
			:std::ostream(nullptr)
			{
			std::ostream::rdbuf(&myBuffer);
			}
		void LinkStream(std::ostream& out)
			{
			out.flush();
			myBuffer.addBuffer(out.rdbuf());
			}
	};

class TThreadSafe_stdcout
	{
	public:
		template<typename T>
		void Write(const T& t)
			{
			std::lock_guard<std::mutex> lock(m_);
			std::cout << t;
			}
	private:
		std::mutex m_;
	};

class TLogPlusCout : public ComposeStream
	{
	public:
		TLogPlusCout(bool doCout, const std::string& logfn)
			{
			if (doCout)
				LinkStream(std::cout);
			if (!(logfn.empty()))
				{
				log_.reset(new std::ofstream(logfn));
				if (!(log_->good()))
					{
					cout << "warning: cannot open log file " << logfn << " for writing" << endl;
					log_.reset(nullptr);
					}
				}
			if (log_)
				LinkStream(*log_);
			}
	private:
		std::unique_ptr<std::ofstream> log_;
	};

TZAxisStereographicSphericalPhaseSpace<float> CreateZAxisSphericalPhaseSpace(const TM25::TTM25RaySet& rs)
	{
	auto vf = rs.VirtualFocus();
	double dist = rs.MaxDistance(vf);
	// now there is a sphere around vf with radius dist which all selected rays intersect.
	// we make another sphere with radius 2.2*dist at vf - (0,0,dist) to erect the phase space
	// hoping rays will stay away from the -z hemisphere
//		TPhaseSpace ps(vf + TVec3f{ 0,0,-dist }, static_cast<float>(dist * 2.2));
	TZAxisStereographicSphericalPhaseSpace<float> ps(vf + TVec3f{ 0,0,0 }, static_cast<float>(dist * 1.1));
	return ps;
	}


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
	void WriteCharacteristicCurve(const std::string& fn, const TCharacteristicCurve& cc);

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
	void WriteSkewnessDistribution_z_axis(const std::string& fn, const TSkewnessDistribution& cc);

	private:
		template<typename PhaseSpace>
		void CreateTree(const PhaseSpace& ps, const TM25::TTM25RaySet& rs, std::ostream& info);
		using TVec = std::vector<float>;
		using TThreeVecs = std::tuple<TVec, TVec, TVec>;
		TThreeVecs ComputeLumRange(size_t ibegin, size_t iend, KDTree::Def::TIdx nNeighbors) const;
	};

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
		for (auto &ii : inbs.i_points_)
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


template<typename PhaseSpace>
void TInterpolateRaySetData::Init(const PhaseSpace& ps, const TM25::TTM25RaySet& rs, KDTree::Def::TIdx nNeighbors, TLogPlusCout& info)
	{ 
	CreateTree(ps, rs, info);
	info << "averaging luminance over " << nNeighbors << " nearest neighbors\n ";
	size_t nRays = rs.NRays();
	bool multithreaded = true; // DO NOT Change, multithreaded is buggy
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
			std::vector<std::future<TThreeVecs>> futures(nThreads - 1);
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

void TInterpolateRaySetData::SetTotalFlux(float newTotalFlux)
	{
	double fac = static_cast<double>(newTotalFlux) / static_cast<double>(totalFlux_);
	for (auto& f : rayFluxes_)
		f = static_cast<float>(static_cast<double>(f) * fac);
	for (auto& l : luminances_)
		l = static_cast<float>(static_cast<double>(l) * fac);
	for (auto& c : cellFluxes_)
		c = static_cast<float>(static_cast<double>(c) * fac);
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

void TInterpolateRaySetData::WriteCharacteristicCurve(const std::string& fn, const TCharacteristicCurve& cc)
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
void TInterpolateRaySetData::WriteSkewnessDistribution_z_axis(const std::string& fn, const TSkewnessDistribution& cc)
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



template<typename PhaseSpace>
TM25::TDefaultRayArray ComputeRayArray(PhaseSpace ps, const KDTree::Def::TKDPoints& interpolatedPsPoints, const std::vector<float>& interpolatedFluxes)
	{
	std::vector<float> ray(7);
	const size_t nInterpolatedRays = interpolatedPsPoints.size();
	TM25::TDefaultRayArray rayArray(nInterpolatedRays, 7); // 7 items
	for (KDTree::Def::TIdx i = 0; i < nInterpolatedRays; ++i)
		{
		KDTree::Def::TKDPoint thisPsPoint = interpolatedPsPoints[i];
		float thisFlux = interpolatedFluxes[i];
		typename PhaseSpace::TV3 loc3 = ps.Loc_to_V3({ thisPsPoint[0],thisPsPoint[1] });
		typename PhaseSpace::TV3 dir3 = ps.Dir_to_V3({ thisPsPoint[0],thisPsPoint[1] }, { thisPsPoint[2],thisPsPoint[3] });
		rayArray.SetRay<7>(i, { loc3[0], loc3[1], loc3[2], dir3[0], dir3[1], dir3[2], thisFlux });
		}
	return rayArray;
	};

int main(int argc, char* argv[])
	{
	cout << "InterpolateRaySet, by Julius Muschaweck, 2019" << endl;
	if (argc < 2)
		{
		cout << "Call me 'InterpolateRaySet cfgFn', where cfgFn is the name of the configuration file" << endl;
		return 0;
		}
	try
		{
		//return Test();
		// Test();
		std::string cfgFn{ argv[1] };
		cout << "parsing configuration file " << cfgFn << endl;
		TInterpolateRaySetCfg cfg;
		cfg.ParseCfgFile(cfgFn);
		const TSection& rsc = cfg.Section("RaySetControl");
		TLogPlusCout logS(rsc.Bool("consoleOutput"), rsc.String("logFileName"));
		logS << "config file " << cfgFn << " successfully read, file content:\n";
		logS << cfg.Content();
		logS << "%% end of configuration file\n\n";
		TM25::TTM25RaySet rs = ReadRaySet(cfg, logS);
		if (rsc.Bool("scrambleInput"))
			{
			rs.Shuffle();
			logS << "input scrambled\n";
			}
		if (!(rsc.IsEmpty("selectByMaxNumber")))
			{
			const int nmax = rsc.Int("selectByMaxNumber");
			int n = 0;
			auto pred = [&n, nmax](const float*, size_t) {return ++n <= nmax; };
			rs.SelectSubset(pred);
			}
		if (rsc.Bool("doDiagnostics"))
			{
			logS << rs.Diagnostics();
			}
		if (!(rsc.IsEmpty("restrictToRelVirtualFocusDistance")))
			{
			auto vf = rs.VirtualFocus();
			double fac = rsc.Real("restrictToRelVirtualFocusDistance");
			const float dist = static_cast<float>(rs.MaxDistance(vf) * fac);
			rs.SelectMaxDistance(dist, vf);
			logS << "# of rays after selecting relative virtual focus distance < " << dist << ": " << rs.NRays() << "\n";
			}
		if (!(rsc.IsEmpty("restrictToKz")))
			{
			const double kzmin = rsc.Real("restrictToKz");
			rs.SelectSubset([kzmin](const float* r, size_t) {float kz = *(r + 5); return kz >= kzmin; });
			logS << "# of rays after selecting kz >= " << kzmin << ": " << rs.NRays() << "\n";
			}

		// compute kd tree and information about rays, cells, etc-
		TInterpolateRaySetData interpRaySetData;

		std::string pstype = rsc.Keyword("phaseSpaceType");
		if (pstype.compare("VirtualFocusZSphere") == 0)
			{
			TZAxisStereographicSphericalPhaseSpace<float>	ps = CreateZAxisSphericalPhaseSpace(rs);
			interpRaySetData.Init(ps, rs, rsc.Int("nNeighbors"), logS);
			}
		else if (pstype.compare("ZPlane") == 0)
			{
			double z0 = rsc.Real("ZPlane_z");
			TPlanarZPhaseSpace<float> ps(static_cast<float>(z0));
			interpRaySetData.Init(ps, rs, rsc.Int("nNeighbors"), logS);
			}
		else if (pstype.compare("ZCylinder") == 0)
			{
			double radius = rsc.Real("ZCylinderRadius");
			TCylinderZPhaseSpace<float> ps(static_cast<float>(radius));
			interpRaySetData.Init(ps, rs, rsc.Int("nNeighbors"), logS);
			}
		else
			throw std::runtime_error("InterpolateRaySet: no valid phase space type");
		if (!rsc.IsEmpty("setTotalFlux"))
			interpRaySetData.SetTotalFlux(static_cast<float>(rsc.Real("setTotalFlux")));

		size_t nRays = rs.NRays();

		if (!rsc.IsEmpty("characteristicCurveFileName"))
			interpRaySetData.WriteCharacteristicCurve(rsc.String("characteristicCurveFileName"),
				interpRaySetData.CharacteristicCurve());
		if (!rsc.IsEmpty("skewness_z_FileName"))
			{
			size_t nBins = rsc.Int("skewness_nBins");
			std::string fn = rsc.String("skewness_z_FileName");
			std::string bts = rsc.Keyword("skewness_binType");
			TInterpolateRaySetData::TSkewnessDistribution::BinType binType = interpRaySetData.BinTypeFromString(bts);
			TInterpolateRaySetData::TSkewnessDistribution sd = interpRaySetData.SkewnessDistribution_z_axis(nBins, binType, rs);
			interpRaySetData.WriteSkewnessDistribution_z_axis(fn, sd);
			}


		// now we have an array of rays, a corresponding array of points in phase space, and corresponding arrays of volumes, fluxes and luminances. 
		// we also have the KD tree structure which gives us a phase space bounding box for each ray.
		// we are ready to create additional rays.
		// Two strategies: 
		// a) constant flux rays, that is many rays for high luminance regions, few rays for low luminance regions
		// b) constant etendue rays, that is high power rays for high luminance regions, low power rays for low luminance regions
		// or any mixture of both. 
		// We use strategy a) here. That is, # of rays per cell is proportional to cell flux
		// Compute the # of rays per cell, including the one that's already in there
		size_t nTotalRays = rsc.Int("nOutputRays");
		logS << "computing # of rays per cell (total target = " << nTotalRays << ")\n";
		std::vector<KDTree::Def::TIdx> raysPerCell(nRays);
		for (size_t i = 0; i < nRays; ++i)
			{
			raysPerCell[i] = static_cast<KDTree::Def::TIdx>(round(interpRaySetData.cellFluxes_[i] / interpRaySetData.totalFlux_ * nTotalRays));
			}
		size_t totalRays = std::accumulate(raysPerCell.begin(), raysPerCell.end(), 0);
		logS << "actual total = " << totalRays << "\n";
		// Again, several strategies:
		// a) Create only new rays, discarding the original ray set
		// b) Keep the original ray set, and add new rays
		//		i) Keep the original ray power, which will be sometimes too large if the cell is small
		//		ii) Keep only the original ray's location and direction, and adjust its power to match the addtl rays
		// We use b) ii) because the original ray set has some value, while the noise is reduced
		// To distribute the rays in the cell, we subdivide the cell and put each new ray at the center of its subcell
		// For e.g. 7 rays, we split the cell into two along dimension 0, with 4/7 and 3/7 volume each and continue down.
		// The original ray will reside in exactly one subcell, the others are added

		// Generate the additional rays per phase space
		KDTree::Def::TKDPoints interpolatedPsPoints;
		std::vector<float> interpolatedFluxes;
		std::default_random_engine dre;
		std::uniform_real_distribution<float> dr; // default [0;1] just fine
		auto ranGen = [&dre, &dr]() -> float {return dr(dre); };
		for (KDTree::Def::TIdx i = 0; i < nRays; ++i)
			{
			const KDTree::Def::TKDPoint& thisPsPoint = interpRaySetData.kdtree_->Point({ i });
			float thisVolume = interpRaySetData.volumes_[i];
			float thisLuminance = interpRaySetData.luminances_[i];
			float thisFlux = interpRaySetData.cellFluxes_[i];
			KDTree::Def::TIdx thisNoOfRays = raysPerCell[i];
			KDTree::TKDTree::TPointIdx pi{ i };
			const KDTree::TNode& thisNode = interpRaySetData.kdtree_->Node(interpRaySetData.kdtree_->NodeIndex(pi));
			//std::vector<KDTree::Def::TKDPoint> newPsPoints = thisNode.Partition(thisNoOfRays, thisPsPoint);
			std::vector<KDTree::Def::TKDPoint> newPsPoints = thisNode.RandomPartition(thisNoOfRays, thisPsPoint, ranGen);
			auto sqr = [](double x) {return x * x; };
			const double kzmin2 = 0.005; // sqr(rsc.Real("restrictToKz"));
			for (auto& npsp : newPsPoints)
				{
				// disregard rays with k0^2+k1^2>1
				if ((1 - sqr(npsp[2]) - sqr(npsp[3])) < kzmin2)
					continue;
				interpolatedPsPoints.push_back(npsp);
				interpolatedFluxes.push_back(thisFlux / thisNoOfRays);
				}
			}
		// now interpolatedPsPoints and interpolatedFluxes contain all the interpolated ray information on phase space level
		// move rays back to 3D space

		TM25::TDefaultRayArray rayArray;
		if (pstype.compare("VirtualFocusZSphere") == 0)
			{
			TZAxisStereographicSphericalPhaseSpace<float>	ps = CreateZAxisSphericalPhaseSpace(rs);
			rayArray = ComputeRayArray(ps, interpolatedPsPoints, interpolatedFluxes);
			}
		else if (pstype.compare("ZPlane") == 0)
			{
			double z0 = rsc.Real("ZPlane_z");
			TPlanarZPhaseSpace<float> ps(static_cast<float>(z0));
			rayArray = ComputeRayArray(ps, interpolatedPsPoints, interpolatedFluxes);
			}
		else if (pstype.compare("ZCylinder") == 0)
			{
			double radius = rsc.Real("ZCylinderRadius");
			TCylinderZPhaseSpace<float> ps(static_cast<float>(radius));
			rayArray = ComputeRayArray(ps, interpolatedPsPoints, interpolatedFluxes);
			}
		else
			throw std::runtime_error("InterpolateRaySet: no valid pahse space type");

		TM25::TTM25Header header = rs.Header();
		header.n_rays_4_7_1_6 = rayArray.NRays();
		TM25::TTM25RaySet interpolatedRaySet(header, std::move(rayArray));

		logS << "writing " << header.n_rays_4_7_1_6 << " rays to TM25 ray file " << rsc.String("outputRayFileName");
		

		interpolatedRaySet.Write(rsc.String("outputRayFileName"));

		TM25::TTM25RaySet test_rs;
		test_rs.Read(rsc.String("outputRayFileName"));
		if (test_rs.Warnings().empty())
			logS << "No warnings" << endl;
		for (auto w : test_rs.Warnings())
			{
			logS << "Warning: " << w << "\n";
			}
		auto bb = test_rs.RayArray().BoundingBox();
		logS << "Bounding Box: x in [" << bb.first[0] << ',' << bb.second[0] << "], y in ["
			<< bb.first[1] << ',' << bb.second[1] << "], z in [" << bb.first[2] << ',' << bb.second[2] << "]" << '\n';
		auto vf = test_rs.VirtualFocus();
		logS << "Virtual Focus: F = [" << vf[0] << ',' << vf[1] << ',' << vf[2] << ']' << '\n';
		logS << "Maximum distance: d = " << test_rs.MaxDistance(vf) << '\n';

		}
	catch (TM25::TM25Error e)
		{
		std::cout << e.what() << std::endl;
		}
	catch (std::runtime_error e)
		{
		std::cout << "std::runtime_error: " << e.what() << std::endl;
		}
	//catch (...)
	//	{
	//	std::cout << "unknown error" << std::endl;
	//	}
	return 0;
	}
