#ifndef __ASCIIRAYFILE_H
#define __ASCIIRAYFILE_H

#include "TM25.h"

namespace TM25
	{

	struct TASCIIRaySetOptions
		{
		int wavelengthColumn_ = 7; // 8th column -- zero based index
		enum class WLU {nanometer, micrometer};
		WLU wavelengthUnit_ = WLU::nanometer;
		size_t minHeaderLines = 0; // no of starting lines that will be considered header no matter what's in them.
		};

	TTM25RaySet ReadGenericASCIIRaySet(const std::string& fn, TASCIIRaySetOptions opts);
	// sequence of header lines, followed by sequence of ray lines, followed by sequence of footer lines
	// header lines and footer lines can be anything that cannot be converted to seven or eight floating point numbers
	// ray lines must contain seven or eight floating point numbers separated by blank or \t
	// all ray lines must contain same number of numbers (seven or eight)
	// ray lines contain x, y, z, kx, ky, kz, (power | (power, wavelength) | (wavelength, power))
	// returns TTM25RaySet. 
	// The ASCII header and footer content will be in additional_info_4_7_3_8
	// Spectral data identifier will be 0 (no wavelength) or 2 (wavelength per ray)



	}

#endif