// InterpolateRaySet.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "../TM25ReadWrite/TM25.h"
#include "KDTree.h"
#include <iostream>
#include "PhaseSpace.h"
#include <algorithm>
#include <numeric>
#include <memory>
#include "CfgFile.h"

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


TM25::TTM25RaySet ReadRaySet(const std::string& fn, std::ostream& info)
	{
	info << "Reading ray file " << fn << endl;
	TM25::TTM25RaySet rs;
	rs.Read(fn);
	// rs.Make_k_unit();
	if (rs.Warnings().empty())
		info << "No warnings" << endl;
	for (auto w : rs.Warnings())
		{
		info << "Warning: " << w << "\n";
		}
	return rs;
	}

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
			values_.insert({ "logFileName",			MakeDefaultValueTokenSequence<Token::string>(std::string("InterpolateRaySet.log")) });
			values_.insert({ "consoleOutput",		MakeDefaultValueTokenSequence<Token::boolean>(true) });
			values_.insert({ "scrambleInput",		MakeDefaultValueTokenSequence<Token::boolean>(true) });
			values_.insert({ "selectByMaxNumber",	MakeEmptyTokenSequence()});
			values_.insert({ "doDiagnostics",		MakeDefaultValueTokenSequence<Token::boolean>(true) });
			values_.insert({ "restrictToRelVirtualFocusDistance",MakeEmptyTokenSequence() });
			values_.insert({ "restrictToKz",		MakeDefaultValueTokenSequence<Token::real>(0.05) });
			values_.insert({ "nOutputRays",			MakeDefaultValueTokenSequence<Token::integer>(1) });
			values_.insert({ "nNeighbors",			MakeDefaultValueTokenSequence<Token::integer>(10) });
			values_.insert({ "outputRayFileName",	MakeDefaultValueTokenSequence<Token::string>(std::string("tmp.TM25RAY")) });
			values_.insert({ "phaseSpaceType",		MakeDefaultValueTokenSequence<Token::identifier>(std::string("VirtualFocusZSphere")) });
			values_.insert({ "ZPlane_z",			MakeDefaultValueTokenSequence<Token::real>(0.0) });			
			
			}
		virtual void AddAllowedKeywords() 
			{
			keywords_.insert("VirtualFocusZSphere");
			keywords_.insert("ZPlane");

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
		std::unique_ptr<KDTree::TKDTree> kdtree_;
		std::vector<float> rayFluxes_;
		std::vector<float> volumes_;
		std::vector<float> luminances_;
		std::vector<float> cellFluxes_;
		float totalVolume_;
		float totalFlux_;
		float avgLuminance_;
	private:
		template<typename PhaseSpace>
		void CreateTree(const PhaseSpace& ps, const TM25::TTM25RaySet& rs, std::ostream& info);
	};

template<typename PhaseSpace>
void TInterpolateRaySetData::Init(const PhaseSpace& ps, const TM25::TTM25RaySet& rs, KDTree::Def::TIdx nNeighbors, TLogPlusCout& info)
	{
	CreateTree(ps, rs, info);
	info << "averaging luminance over " << nNeighbors << " nearest neighbors\n ";
	size_t nRays = rs.NRays();
	for (size_t i = 0; i < nRays; ++i)
		{
		size_t nraysPercent = (nRays>100) ? (nRays / 100) : 1;
		if ((nRays > 100) && (i % nraysPercent == 0))
			info << round(double(i) / nRays * 100.0) << "% ";
		KDTree::TKDTree::TPointIdx pi{ static_cast<KDTree::Def::TIdx>(i) };
		KDTree::TKDTree::TNearestNeighbors inbs = kdtree_->NearestNeighborsOfPoint(pi, nNeighbors);
		float vol = kdtree_->TotalVolume(inbs.i_nodes_);
		float flux = 0;
		for (auto i : inbs.i_points_)
			flux += rayFluxes_[i.pi_];
		float avgLuminance = flux / vol;
		luminances_.push_back(avgLuminance);
		KDTree::TKDTree::TNodeIdx ni = kdtree_->NodeIndex(pi);
		float iVolume = (kdtree_->Node(ni)).Volume();
		volumes_.push_back(iVolume);
		cellFluxes_.push_back(iVolume * avgLuminance);
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

int main(int argc, char *argv[])
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
		Test();
		std::string cfgFn{ argv[1] };
		cout << "parsing configuration file " << cfgFn << endl;
		TInterpolateRaySetCfg cfg;
		cfg.ParseCfgFile(cfgFn);
		const TSection& rsc = cfg.Section("RaySetControl");
		TLogPlusCout logS(rsc.Bool("consoleOutput"), rsc.String("logFileName"));
		logS << "config file "<<cfgFn << " successfully read, file content:\n";
		logS << cfg.Content();
		logS << "%% end of configuration file\n\n";
		TM25::TTM25RaySet rs = ReadRaySet(rsc.String("inputRayFileName"), logS);
		if (rsc.Bool("scrambleInput"))
			{
			rs.Shuffle();
			logS << "input scrambled\n";
			}
		if (!(rsc.IsEmpty("selectByMaxNumber")))
			{
			const int nmax = rsc.Int("selectByMaxNumber");
			int n = 0;
			auto pred = [&n,nmax](const float*, size_t) {return ++n <= nmax; };
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
		size_t nRays = rs.NRays();

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
			const double kzmin2 = sqr(rsc.Real("restrictToKz"));
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

		TM25::TTM25Header header = rs.Header();
		header.n_rays_4_7_1_6 = rayArray.NRays();
		TM25::TTM25RaySet interpolatedRaySet(header, std::move(rayArray));

		logS << "writing " << rayArray.NRays() << " rays to TM25 ray file " << rsc.String("outputRayFileName");
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
