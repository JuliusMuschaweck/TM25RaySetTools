// CreatePointSourceRayFile.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "CreatePointSourceRayFileConfig.h"
#include <LogPlusCout.h>
#include <TM25.h>
#include <ReadASCIIMatrix.h>
#include <Sobol.h>
#include <numbers>
#include <iomanip>
#include <time.h>

#include <TM25.h>

TM25::TTM25RaySet CreateRayPointSourceRaySet(const TSection& rsc, TLogPlusCout& logS, const TM25::TReadASCIIMatrix& I_kx_ky)
	{
	TM25::TDefaultRayArray rays;
	const size_t nRays = rsc.Int("nOutputRays");
	constexpr size_t nItems = 7; // x, y, z, kx, ky kz power
	rays.Resize(nRays, nItems);
	TSobol<3> sobol; // first two: to be transformed to kx, ky. Third: intensity fraction
	
	double zmax = I_kx_ky.zMax();
	double zmin = I_kx_ky.zMin();
	double rel_threshold = rsc.Real("rel_threshold");

	size_t i = 0;
	while (i < nRays)
		{
		std::array<double, 3> sob = sobol();
		double r = sqrt(sob[0]);
		double theta = 2 * std::numbers::pi * sob[1];
		double kx = r * sin(theta);
		double ky = r * cos(theta);
		double z = I_kx_ky.Interpolate(kx, ky, NAN);
		if (std::isnan(z))
			continue;
		double zrel = (z - zmin) / (zmax - zmin);
		float fkx = static_cast<float>(kx);
		float fky = static_cast<float>(ky);
		float fkz = static_cast<float>(sqrt(1 - kx * kx - ky * ky));
		float fzrel = static_cast<float>(zrel);
		if (zrel > rel_threshold) // just use it
			{
			std::array<float, 7> ray{ 0,0,0,fkx,fky,fkz,fzrel };
			rays.SetRay<7>(i, ray);
			++i;
			}
		else if (zrel / rel_threshold > sob[2])
				{
				std::array<float, 7> ray{ 0,0,0,fkx,fky,fkz, static_cast<float>(rel_threshold)};
				rays.SetRay<7>(i, ray);
				++i;
				}
		}
	
	TM25::TTM25Header header;
	header.phi_4_7_1_5 = static_cast<float>(rsc.Real("setTotalFlux"));
	header.n_rays_4_7_1_6 = nRays;
	const std::time_t t = std::time(nullptr);
	std::tm tm;
	gmtime_s(&tm, &t);
	char timestr[100];
	std::strftime(timestr, 100, "%Y-%m-%dT%H:%M:%S%z", &tm);
	header.file_date_time_str_4_7_1_7 = timestr;
		// "2013-09-04T08:30:29+01:00"; // example from TM25
	TM25::TTM25RaySet rv(header, rays);
	return rv;
	}



int main(int argc, char* argv[])
{
	using std::cout;
	using std::endl;
	cout << "CreatePointSourceRayFile, by Julius Muschaweck, 2025" << endl;
	if (argc < 2)
		{
		cout << "Call me 'CreatePointSourceRayFile cfgFn', where cfgFn is the name of the configuration file" << endl;
		return 0;
		}
	try
		{
		// read and parse configuration file
		std::string cfgFn{ argv[1] };
		cout << "parsing configuration file " << cfgFn << endl;
		TPointSourceRayFileCfg cfg;
		cfg.ParseCfgFile(cfgFn);
		const TSection& rsc = cfg.Section("PointSourceRayFileControl");
		TLogPlusCout logS(rsc.Bool("consoleOutput"), rsc.String("logFileName"));

		logS << "config file " << cfgFn << " successfully read, file content:\n";
		logS << cfg.Content();
		logS << "%% end of configuration file\n\n";

		// read intensity file
		TM25::TReadASCIIMatrix I_kx_ky;
		std::string delims = rsc.String("delimiters");
		double eps = rsc.Real("equidistance_rel_epsilon");
		I_kx_ky.Read("farfield.txt", delims, eps);
		TM25::TTM25RaySet rayset = CreateRayPointSourceRaySet(rsc, logS, I_kx_ky);

		rayset.Write(rsc.String("outputRayFileName"));

		}
	catch (TM25::TM25Error e)
		{
		std::cout << e.what() << std::endl;
		}
	catch (std::runtime_error e)
		{
		std::cout << "std::runtime_error: " << e.what() << std::endl;
		}
	catch (std::invalid_argument e)
		{
		std::cout << "std::invalid_argument: " << e.what() << std::endl;
		}


}

