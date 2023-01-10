#include "ReadASAPRaySet.h"
#include <ReadFile.h>

TM25::TTM25RaySet ReadASAPRaySet(const std::string& fn)
	{
	TM25::TReadFile f(fn);
	// record 1
	std::array<int32_t, 7> rec1 = f.Read<int32_t, 7>();
	int LENR = (rec1[0] - 247) / 256 / 4;
	if (LENR != 7)
		throw std::runtime_error("ReadASAPRaySet: LENR != 7 -> " + std::to_string(LENR));
	// record 2
	std::array<char, 24> TITLE = f.Read<char, 24>();
	std::string title(TITLE.begin(), TITLE.end());
	int32_t NUMF = f.Read<int32_t>();
	// record 3
	std::array<char, 8> ZLABEL = f.Read<char, 8>();
	std::string zlabel(ZLABEL.begin(), ZLABEL.end());
	union i32_float
		{
		int32_t i;
		float f;
		};
	i32_float ZVAL;
	ZVAL.i = f.Read<int32_t>();
	std::array<char, 16> FLABEL = f.Read<char, 16>();
	std::string flabel(FLABEL.begin(), FLABEL.end());
	// record 4
	std::array<char, 16> YLABEL = f.Read<char, 16>();
	std::string ylabel(YLABEL.begin(), YLABEL.end());
	i32_float YMIN; YMIN.i = f.Read<int32_t>();
	i32_float YMAX; YMAX.i = f.Read<int32_t>();
	i32_float NUMY; NUMY.i = f.Read<int32_t>();
	// record 5
	std::array<char, 16> XLABEL = f.Read<char, 16>();
	std::string xlabel(XLABEL.begin(), XLABEL.end());
	i32_float XMIN; XMIN.i = f.Read<int32_t>();
	i32_float XMAX; XMAX.i = f.Read<int32_t>();
	i32_float NUMX; NUMX.i = f.Read<int32_t>();

	TM25::TTM25Header h;
	h.n_rays_4_7_1_6 = - NUMY.i;
	h.file_date_time_str_4_7_1_7 = "no date string implemented";
	h.name_4_7_3_1 = std::u32string(fn.begin(), fn.end());

	size_t nItems = 7;
	size_t nRays = h.n_rays_4_7_1_6;
	TM25::TDefaultRayArray r(nRays, nItems);

	for (size_t i = 0; i < nRays; ++i)
		{
		std::array<float, 7> ray = f.Read<float, 7>();
		r.SetRay<7>(i, ray);
		}

	TM25::TTM25RaySet rv(h, std::move(r));
	return rv;
	}
