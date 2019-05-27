// InterpolateRaySet.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "../TM25ReadWrite/TM25.h"
#include "KDTree.h"
#include <iostream>
#include "PhaseSpace.h"
#include <algorithm>

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
//		return Test();
		TM25::TTM25RaySet rs;
		std::string fn = "../rayfile_LERTDUW_S2WP_blue_100k_20161013_IES_TM25.TM25RAY";
		// this small ray file is on GitHub
		// std::string fn = "../rayfile_LERTDUW_S2WP_green_20M_20161013_IES_TM25.TM25RAY";
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
		auto pred = [&n](const float*, size_t) {return ++n < 10000; };
		rs.SelectSubset(pred);
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
		TZAxisStereographicSphericalPhaseSpace<float> ps(vf + TVec3f{ 0,0,-dist }, static_cast<float>(dist * 2.2));
		size_t nr = rs.NRays();
		
		// create the array of points for the KD tree
		// as well as a same size and order array of fluxes
		KDTree::Def::TKDPoints pspoints;
		std::vector<float> fluxes;
		for (size_t i = 0; i < nr; ++i)
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
		for (size_t i = 0; i < nr; ++i)
			{
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
		// now we have an array of rays, a corresponding array of points in phase space, and corresponding arrays of volumes, fluxes and luminances. 
		// we also have the KD tree structure which gives us a phase space bounding box for each ray.
		// we are ready to create additional rays.

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
