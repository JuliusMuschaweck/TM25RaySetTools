// InterpolateRaySet.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "../TM25ReadWrite/TM25.h"
#include "KDTree.h"
#include <iostream>
#include "PhaseSpace.h"

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
		KDTree::TKDTree kdtree(std::move(pspoints));
		kdtree.CreateTree();
		std::vector<float> volumes;
		std::vector<float> luminances;
		size_t nNeighbors = 10;
		for (size_t i = 0; i < nr; ++i)
			{
			KDTree::TKDTree::TNearestNeighbors inbs = kdtree.NearestNeighborsOfPoint(i, nNeighbors);
			float vol = kdtree.
			}

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
