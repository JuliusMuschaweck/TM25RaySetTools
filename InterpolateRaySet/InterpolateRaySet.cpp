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
#include <WriteFile.h>
#include "InterpolateRaySetConfig.h"
#include <ASCIIRayFile.h>
#include "InterpolateRaySet_IO.h"
#include "InterpolateRaySetData.h"
#include "PhaseSpaceGrid.h"

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

void PrepareRaySet(const TInterpolateRaySetCfg& cfg, const TSection& rsc, TLogPlusCout& logS, TM25::TTM25RaySet& rs)
	{
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
	if (!(rsc.IsEmpty("restrictToXYBox")))
		{
		std::vector<double> rawbox = rsc.RealVector("restrictToXYBox"); // xmin ymin xmax ymax
		double xmin = rawbox[0];
		double ymin = rawbox[1];
		double xmax = rawbox[2];
		double ymax = rawbox[3];

		auto selector = [xmin, ymin, xmax, ymax](const float* r, size_t)
			{
			float x = *r;
			float y = *(r + 1);
			bool rv = x >= xmin && x <= xmax && y >= ymin && y <= ymax;
			return rv;
			};
		rs.SelectSubset(selector);
		logS << "# of rays after selecting x/y box >= (" << "xmin <<','<< ymin << ',' << xmax << ',' << ymax" << "): " << rs.NRays() << "\n";
		}
	if (!(rsc.IsEmpty("restrictToFirstNRays")))
		{
		size_t i = 0;
		const size_t n = static_cast<size_t>(rsc.Int("restrictToFirstNRays"));
		auto selector = [&i, n](const float* r, size_t)
			{
			return i++ < n;
			};
		rs.SelectSubset(selector);
		logS << "# of rays after selecting first "<< n<<" rays: " << rs.NRays() << "\n";
		}
	}

void PrepareRaySetData(TInterpolateRaySetData& interpRaySetData, const TSection& rsc, TLogPlusCout& logS, const TM25::TTM25RaySet& rs)
	{
	std::string pstype = rsc.Keyword("phaseSpaceType");
	logS << "preparing phase space of type " << pstype << "\n";
	int nClip = rsc.Int("nClip");
	int nNeighbors = rsc.Int("nNeighbors");
	if (pstype.compare("VirtualFocusZSphere") == 0)
		{
		TZAxisStereographicSphericalPhaseSpace<float> ps = CreateZAxisSphericalPhaseSpace(rs);
		interpRaySetData.Init(ps, rs, nNeighbors, nClip, logS);
		}
	else if (pstype.compare("ZPlane") == 0)
		{
		double z0 = rsc.Real("ZPlane_z");
		TPlanarZPhaseSpace<float> ps(static_cast<float>(z0));
		interpRaySetData.Init(ps, rs, nNeighbors, nClip, logS);
		}
	else if (pstype.compare("ZCylinder") == 0)
		{
		double radius = rsc.Real("ZCylinderRadius");
		TCylinderZPhaseSpace<float> ps(static_cast<float>(radius));
		interpRaySetData.Init(ps, rs, nNeighbors, nClip, logS);
		}
	else
		throw std::runtime_error("InterpolateRaySet: no valid phase space type");
	if (!rsc.IsEmpty("setTotalFlux"))
		interpRaySetData.SetTotalFlux(static_cast<float>(rsc.Real("setTotalFlux")));
	}

void WriteCharacteristicCurve(const TInterpolateRaySetData& interpRaySetData, const TSection& rsc, TLogPlusCout& logS)
	{
	if (rsc.IsEmpty("characteristicCurveFileName"))
		{
		logS << "no characteristicCurveFileName given\n";
		return;
		}
	std::string fn = rsc.String("characteristicCurveFileName");
	logS << "Writing characteristic curve to " << fn << "\n";
	interpRaySetData.WriteCharacteristicCurve(fn,
			interpRaySetData.CharacteristicCurve());
	}

void WriteSkewnessDistribution(const TInterpolateRaySetData& interpRaySetData, const TSection& rsc, TLogPlusCout& logS, const TM25::TTM25RaySet& rs)
	{
	if (rsc.IsEmpty("skewness_z_FileName"))
		{
		logS << "no skewness_z_FileName given\n";
		return;
		}
	std::string fn = rsc.String("skewness_z_FileName");
	size_t nBins = rsc.Int("skewness_nBins");
	std::string bts = rsc.Keyword("skewness_binType");
	logS << "Writing skewness data of type "<<bts<<" to " << fn << "\n";
	TInterpolateRaySetData::TSkewnessDistribution::BinType binType = interpRaySetData.BinTypeFromString(bts);
	TInterpolateRaySetData::TSkewnessDistribution sd = interpRaySetData.SkewnessDistribution_z_axis(nBins, binType, rs);
	interpRaySetData.WriteSkewnessDistribution_z_axis(fn, sd);
	}

TM25::TDefaultRayArray ComputeInterpolatedRays(const TInterpolateRaySetData& interpRaySetData, const TSection& rsc, TLogPlusCout& logS, const TM25::TTM25RaySet& rs)
	{
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
	size_t nRays = rs.NRays();
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
	std::string pstype = rsc.Keyword("phaseSpaceType");
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
		throw std::runtime_error("InterpolateRaySet: no valid phase space type");
	return rayArray;
	}


void WriteInterpolatedRayFile(const TSection& rsc, TLogPlusCout& logS, const TM25::TTM25RaySet& rs, 
		const TInterpolateRaySetData& interpRaySetData)
	{
	if (rsc.IsEmpty("outputRayFileName"))
		{
		logS << "no outputRayFileName given\n";
		return;
		}

	std::string fn = rsc.String("outputRayFileName");
	TM25::TDefaultRayArray rayArray = ComputeInterpolatedRays(interpRaySetData, rsc, logS, rs);
	TM25::TTM25Header header = rs.Header();
	header.n_rays_4_7_1_6 = rayArray.NRays();
	TM25::TTM25RaySet interpolatedRaySet(header, std::move(rayArray));

	logS << "writing " << header.n_rays_4_7_1_6 << " rays to TM25 ray file " << fn;


	interpolatedRaySet.Write(rsc.String("outputRayFileName"));

	TM25::TTM25RaySet test_rs;
	test_rs.Read(rsc.String("outputRayFileName"));
	if (test_rs.Warnings().empty())
		logS << "No warnings" << std::endl;
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

void WriteLuminanceLookupTable(TInterpolateRaySetData& interpRaySetData, const TSection& rsc, TLogPlusCout& logS)
	{
	if (rsc.IsEmpty("LuminanceLookupTable_FileName"))
		{
		logS << "no LuminanceLookupTable_FileName given\n";
		return;
		}
	std::string fn = rsc.String("LuminanceLookupTable_FileName");
	double xmin = rsc.Real("LuminanceLookupTable_xmin");
	double xmax = rsc.Real("LuminanceLookupTable_xmax");
	double ymin = rsc.Real("LuminanceLookupTable_ymin");
	double ymax = rsc.Real("LuminanceLookupTable_ymax");
	double kxmin = rsc.Real("LuminanceLookupTable_kxmin");
	double kxmax = rsc.Real("LuminanceLookupTable_kxmax");
	double kymin = rsc.Real("LuminanceLookupTable_kymin");
	double kymax = rsc.Real("LuminanceLookupTable_kymax");

	auto flt = [](double rhs) {return static_cast<float>(rhs); };
	auto szt = [](int rhs) {return static_cast<size_t>(rhs); };

	int xPoints = rsc.Int("LuminanceLookupTable_xPoints");
	int yPoints = rsc.Int("LuminanceLookupTable_yPoints");
	int kxPoints = rsc.Int("LuminanceLookupTable_kxPoints");
	int kyPoints = rsc.Int("LuminanceLookupTable_kyPoints");
	TPhaseSpaceGrid<float>::TBox boundingBox{ {flt(xmin), flt(ymin), flt(kxmin), flt(kymin)},{flt(xmax), flt(ymax), flt(kxmax), flt(kymax)} };
	TPhaseSpaceGrid<float>::T4DIndex bins{ szt(xPoints), szt(yPoints), szt(kxPoints), szt(kyPoints) };
	

	logS << "computing luminance interpolation table with " << xPoints * yPoints * kxPoints * kyPoints << " entries\n";
	TPhaseSpaceGrid<float> grid = FillPhaseSpaceGrid(boundingBox, interpRaySetData, bins);

	logS << "writing luminance interpolation table to file " << fn << "\n";
	TM25::TWriteFile f(fn);
	f.Write<TPhaseSpaceGrid<float>::T4DPoint>(grid.bbmin_);
	f.Write<TPhaseSpaceGrid<float>::T4DPoint>(grid.bbmax_);
	f.Write<size_t>(grid.nx0_);
	f.Write<size_t>(grid.nx1_);
	f.Write<size_t>(grid.nk0_);
	f.Write<size_t>(grid.nk1_);
	for (float d : grid.v_)
		f.Write<float>(d);
	}

int main(int argc, char* argv[])
	{
	using std::cout;
	using std::endl;
	cout << "InterpolateRaySet, by Julius Muschaweck, 2019" << endl;
	if (argc < 2)
		{
		cout << "Call me 'InterpolateRaySet cfgFn', where cfgFn is the name of the configuration file" << endl;
		return 0;
		}
	try
		{
		// read and parse configuration file
		std::string cfgFn{ argv[1] };
		cout << "parsing configuration file " << cfgFn << endl;
		TInterpolateRaySetCfg cfg;
		cfg.ParseCfgFile(cfgFn);
		const TSection& rsc = cfg.Section("RaySetControl");
		TLogPlusCout logS(rsc.Bool("consoleOutput"), rsc.String("logFileName"));

		logS << "config file " << cfgFn << " successfully read, file content:\n";
		logS << cfg.Content();
		logS << "%% end of configuration file\n\n";

		// read the raw input ray file
		TM25::TTM25RaySet rs = ReadRaySet(cfg, logS);
		// apply scramble, ray selection etc. 
		// to "clean up" the raw ray file by selecting a proper subset of the incoming rays
		PrepareRaySet(cfg, rsc, logS, rs);

		// compute kd tree and information about rays, cells, etc-
		// also finds nearest neighbors for each ray and computes "moving average" luminance for each ray
		TInterpolateRaySetData interpRaySetData;
		PrepareRaySetData(interpRaySetData, rsc, logS, rs);

		// the following routines do nothing except a message if the corresponding file name is not given in cfg file.

		// interesting but not needed to generate either large ray sets or luminance lookup tables
		WriteCharacteristicCurve(interpRaySetData, rsc, logS);
		WriteSkewnessDistribution(interpRaySetData, rsc, logS, rs);

		// generate as many "interpolated" rays as you like -- not needed to generate luminance lookup tables
		// currently only in TM25, but writing to Zemax Binary is easy to hook up
		WriteInterpolatedRayFile(rsc, logS, rs, interpRaySetData);
		
		// compute average luminance values on a regular grid in phase space
		// and write to binary output file
		WriteLuminanceLookupTable(interpRaySetData, rsc, logS);
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
