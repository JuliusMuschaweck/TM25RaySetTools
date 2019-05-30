// InterpolateRaySet.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "../TM25ReadWrite/TM25.h"
#include "KDTree.h"
#include <iostream>
#include "PhaseSpace.h"
#include <algorithm>
#include <numeric>
int Test()
	{
	TestPhaseSpace();
	if (KDTree::Def::dim == 2)
		KDTree::TestKDTree2D("TestKDTree.m");
	if (KDTree::Def::dim == 4)
		KDTree::TestKDTree4D();
	return 0;
	}

int main()
	{
	using std::cout;
	using std::endl;
	try
		{
		//return Test();
		Test();
		TM25::TTM25RaySet rs;
		// std::string fn = "../rayfile_LERTDUW_S2WP_blue_100k_20161013_IES_TM25.TM25RAY";
		// this small ray file is on GitHub
		std::string fn = "../rayfile_LERTDUW_S2WP_green_20M_20161013_IES_TM25.TM25RAY";
		// this 500 MB ray file with 20 million rays can be downloaded from www.osram-os.com
		cout << "Reading ray file " << fn << endl;
		rs.Read(fn);
		// rs.Make_k_unit();
		if (rs.Warnings().empty())
			cout << "No warnings" << endl;
		for (auto w : rs.Warnings())
			{
			std::cout << "Warning: "<< w << "\n";
			}
		int n = 0;
//		auto pred = [&n](const float*, size_t) {return ++n < 10000; };
//		rs.SelectSubset(pred);
		auto bb = rs.RayArray().BoundingBox();
		cout << "Bounding Box: x in [" << bb.first[0] << ',' << bb.second[0] << "], y in ["
			<< bb.first[1] << ',' << bb.second[1] << "], z in [" << bb.first[2] << ',' << bb.second[2] << "]" << endl;
		auto vf = rs.VirtualFocus();
		cout << "Virtual Focus: F = [" << vf[0] << ',' << vf[1] << ',' << vf[2] << ']' << endl;
		cout << "Maximum distance: d = " << rs.MaxDistance(vf) << std::endl;
		auto dh = rs.DistanceHistogram(vf, 10);
		cout << "Distance histogram:\n";
		size_t i = 0;
		size_t nRaysHisto = 0;
		double fluxHisto = 0;
		for (auto idh : dh)
			{
			nRaysHisto += idh.nRays_;
			fluxHisto += idh.flux_;
			std::cout << "bin " << i++ << ": distance <= " << idh.dist_ << ", # of rays = " << idh.nRays_<<'/'<<nRaysHisto 
				<< ", flux = " << idh.flux_<<'/'<<fluxHisto << endl;
			}
		cout << "total # of rays = " << nRaysHisto << ", total flux = " << fluxHisto << endl;

		auto fh = rs.FluxHistogram(10);
		cout << "Flux histogram:\n";
		i = 0;
		nRaysHisto = 0;
		fluxHisto = 0;
		for (auto ifh : fh)
			{
			nRaysHisto += ifh.nRays_;
			fluxHisto += ifh.fluxInBin_;
			std::cout << "bin " << i++ << ": flux <= " << ifh.fluxLimit_<< ", # of rays = " << ifh.nRays_ << '/' << nRaysHisto
				<< ", flux = " << ifh.fluxInBin_ << '/' << fluxHisto << endl;
			}
		cout << "total # of rays = " << nRaysHisto << ", total flux = " << fluxHisto << endl;

		// 99% of rays are within half the max distance
		float dist = static_cast<float>(rs.MaxDistance(vf) / 2);
		rs.SelectMaxDistance(dist, vf);
		cout << "# of rays after selecting half distance: " << rs.NRays() << " dist "<< dist <<" >= " << rs.MaxDistance(vf) << endl;
		// now there is a sphere around vf with radius dist which all selected rays intersect.
		// we make another sphere with radius 2.2*dist at vf - (0,0,dist) to erect the phase space
		// hoping rays will stay away from the -z hemisphere
		using TPhaseSpace = TZAxisStereographicSphericalPhaseSpace<float>;
//		TPhaseSpace ps(vf + TVec3f{ 0,0,-dist }, static_cast<float>(dist * 2.2));
		TPhaseSpace ps(vf + TVec3f{ 0,0,0 }, static_cast<float>(dist * 1.1));
		size_t nRays = rs.NRays();
		
		// create the array of points for the KD tree
		// as well as a same size and order array of fluxes
		KDTree::Def::TKDPoints pspoints;
		std::vector<float> fluxes;
		for (size_t i = 0; i < nRays; ++i)
			{
			using std::get;
			std::tuple<TVec3f, TVec3f, float> iray = rs.RayLocDirFlux(i);
			KDTree::Def::TKDPoint phasespacepoint = ps.PhaseSpacePoint(get<0>(iray), get<1>(iray));
			pspoints.push_back(phasespacepoint);
			fluxes.push_back(get<2>(iray));
			}
		// create tree using the move constructor
		KDTree::TKDTree kdtree(std::move(pspoints));
		kdtree.CreateTree();
		kdtree.ShrinkEdgeNodes();
		kdtree.CheckConsistency();
		// for each cell, find the 10 nearest neighbors
		// compute their total volume, total power and average luminance
		// set cell luminance to avg luminance over neighbors
		std::vector<float> volumes;
		std::vector<float> luminances;
		KDTree::Def::TIdx nNeighbors = 10;
		cout << "finding nearest neighbors" << endl;
		for (size_t i = 0; i < nRays; ++i)
			{
			size_t nraysPercent = nRays / 100;
			if (i % nraysPercent == 0)
				cout << round(double(i) / nRays * 100.0) << "% ";
			cout.flush();
			KDTree::TKDTree::TPointIdx pi{ static_cast<KDTree::Def::TIdx>(i) };
			KDTree::TKDTree::TNearestNeighbors inbs = kdtree.NearestNeighborsOfPoint(pi, nNeighbors);
			float vol = kdtree.TotalVolume(inbs.i_nodes_);
			float flux = 0;
			for (auto i : inbs.i_points_)
				flux += fluxes[i.pi_];
			float avgLuminance = flux / vol;
			luminances.push_back(avgLuminance);
			volumes.push_back(kdtree.Nodes()[kdtree.NodeIndex()[pi.pi_].ni_].Volume());
			}
		cout << endl;
		float totalVolume = std::accumulate(volumes.begin(), volumes.end(), 0.0f);
		float totalFlux = std::accumulate(fluxes.begin(), fluxes.end(), 0.0f);
		float avgLuminance = totalFlux / totalVolume;
		std::vector<float> cellFluxes(nRays);
		for (size_t i = 0; i < nRays; ++i)
			cellFluxes[i] = volumes[i] * luminances[i];
		// now we have an array of rays, a corresponding array of points in phase space, and corresponding arrays of volumes, fluxes and luminances. 
		// we also have the KD tree structure which gives us a phase space bounding box for each ray.
		// we are ready to create additional rays.
		// Two strategies: 
		// a) constant flux rays, that is many rays for high luminance regions, few rays for low luminance regions
		// b) constant etendue rays, that is high power rays for high luminance regions, low power rays for low luminance regions
		// or any mixture of both. 
		// We use strategy a) here. That is, # of rays per cell is proportional to cell flux
		// Compute the # of rays per cell, including the one that's already in there
		size_t nTotalRays = 10 * nRays;
		std::vector<KDTree::Def::TIdx> raysPerCell(nRays);
		for (size_t i = 0; i < nRays; ++i)
			{
			raysPerCell[i] = static_cast<KDTree::Def::TIdx>(round(fluxes[i] / totalFlux * nTotalRays));
			}
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
		for (KDTree::Def::TIdx i = 0; i < nRays; ++i)
			{
			const KDTree::Def::TKDPoint& thisPsPoint = kdtree.Point({ i });
			float thisVolume = volumes[i];
			float thisLuminance = luminances[i];
			float thisFlux = cellFluxes[i];
			KDTree::Def::TIdx thisNoOfRays = raysPerCell[i];
			KDTree::TKDTree::TPointIdx pi{ i };
			const KDTree::TNode& thisNode = kdtree.Node(kdtree.NodeIndex(pi));
			std::vector<KDTree::Def::TKDPoint> newPsPoints = thisNode.Partition(thisNoOfRays, thisPsPoint);
			for (auto& npsp : newPsPoints)
				{
				interpolatedPsPoints.push_back(npsp);
				interpolatedFluxes.push_back(thisFlux / thisNoOfRays);
				}
			}
		// now interpolatedPsPoints and interpolatedFluxes contain all the interpolated ray information on phase space level
		// move rays back to 3D space
		size_t nInterpolatedRays = interpolatedPsPoints.size();
		TM25::TDefaultRayArray rayArray(nInterpolatedRays, 7); // 7 items
		std::vector<float> ray(7);
		for (KDTree::Def::TIdx i = 0; i < nInterpolatedRays; ++i)
			{
			KDTree::Def::TKDPoint thisPsPoint = interpolatedPsPoints[i];
			float thisFlux = interpolatedFluxes[i];
			TPhaseSpace::TV3 loc3 = ps.Loc_to_V3({ thisPsPoint[0],thisPsPoint[1] });
			TPhaseSpace::TV3 dir3 = ps.Dir_to_V3({ thisPsPoint[0],thisPsPoint[1] }, { thisPsPoint[2],thisPsPoint[3] });
			rayArray.SetRay<7>(i, { loc3[0], loc3[1], loc3[2], dir3[0], dir3[1], dir3[2], thisFlux });
			}
		TM25::TTM25Header header = rs.Header();
		header.n_rays_4_7_1_6 = rayArray.NRays();
		TM25::TTM25RaySet interpolatedRaySet(header, std::move(rayArray));
		interpolatedRaySet.Write("InterpolatedRaySet.TM25RAY");

		TM25::TTM25RaySet test_rs;
		test_rs.Read("InterpolatedRaySet.TM25RAY");
		if (test_rs.Warnings().empty())
			cout << "No warnings" << endl;
		for (auto w : test_rs.Warnings())
			{
			std::cout << "Warning: " << w << "\n";
			}
		bb = test_rs.RayArray().BoundingBox();
		cout << "Bounding Box: x in [" << bb.first[0] << ',' << bb.second[0] << "], y in ["
			<< bb.first[1] << ',' << bb.second[1] << "], z in [" << bb.first[2] << ',' << bb.second[2] << "]" << endl;
		vf = test_rs.VirtualFocus();
		cout << "Virtual Focus: F = [" << vf[0] << ',' << vf[1] << ',' << vf[2] << ']' << endl;
		cout << "Maximum distance: d = " << test_rs.MaxDistance(vf) << std::endl;

		}
	catch (TM25::TM25Error e)
		{
		std::cout << e.what() << std::endl;
		}
	catch (std::runtime_error e)
		{
		std::cout << "std::runtime_error: " << e.what() << std::endl;
		}
	catch (...)
		{
		std::cout << "unknown error" << std::endl;
		}
	return 0;
	}
