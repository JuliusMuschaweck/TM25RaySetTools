#include "ASCIIRayFile.h"
#include "ReadFile.h"
#include <charconv>
#include <limits>
#include <Timer.h>

namespace TM25
	{

	bool IsRayLine(const std::string& s, std::array<float, 8>& ray, size_t& nItems)
		{
		auto items = SplitString(s);
		nItems= items.size();
		if (nItems < 7 || nItems > 8)
			return false;
		const char* b = s.c_str();
		auto iray = ray.begin();
		for (const auto& item: items)
			{
			std::from_chars_result res = std::from_chars(b + item.first, b + item.second, *iray);
			if (res.ec != std::errc())
				return false;
			++iray;
			}
		if (nItems == 7)
			ray[7] = std::numeric_limits<float>::quiet_NaN();
		return true;
		}

	inline void AppendRawRay(std::vector<float>& data, std::array<float, 8> ray, size_t nItems)
		{
		size_t sizeBefore = data.size();
		data.resize(sizeBefore + nItems);
		float* raystart = data.data() + sizeBefore;
		std::memcpy(raystart, ray.data(), nItems * sizeof(float));
		}

	struct TRawASCIIRaySetContent
		{
		std::vector<std::string> header_;
		std::vector<std::string> footer_;
		std::vector<float> data_;
		size_t nItems_ = 0;
		size_t nRays_ = 0;

		std::pair<float, float> MinMaxColValue(size_t column)
			{
			std::pair<float, float> rv
				{
				std::numeric_limits<float>::max(),
				std::numeric_limits<float>::min() };
			for (size_t i = 0; i < nRays_; ++i)
				{
				float f = data_[i * nItems_ + column];
				if (f < rv.first) rv.first = f;
				if (f > rv.second) rv.second = f;
				}
			}
		void SwapColumns(size_t lhs, size_t rhs)
			{
			for (size_t i = 0; i < nRays_; ++i)
				{
				std::swap(data_[i * nItems_ + lhs], data_[i * nItems_ + rhs]);
				}
			}
		};


	TRawASCIIRaySetContent ReadRawASCIIRaySet(const std::string& fn)
		{
		TRawASCIIRaySetContent rv;
		TReadFile f(fn);
		float nan = std::numeric_limits<float>::signaling_NaN();
		std::array<float, 8> ray{ nan, nan, nan, nan, nan, nan, nan, nan };
		std::string buf;
		size_t nItems = 0;
		// read header lines, if any, detect first ray and its nItems
		while (!f.AtEof())
			{
			if (IsRayLine(f.ReadLine(buf), ray, nItems))
				{
				AppendRawRay(rv.data_, ray, nItems);
				rv.nRays_ = 1;
				break;
				}
			else
				{
				if (!buf.empty())
					rv.header_.push_back(buf);
				}
			}
		// now header is complete and first ray is read
		const size_t nExpectedItems = nItems;
		rv.nItems_ = nItems;
		while (!f.AtEof())
			{
			if (IsRayLine(f.ReadLine(buf), ray, nItems))
				{
				if (nItems != nExpectedItems)
					throw std::runtime_error("ReadRawASCIIRaySet:: wrong number of items in ray line");
				AppendRawRay(rv.data_, ray, nItems);
				++rv.nRays_;
				}
			else
				{
				if (!buf.empty())
					rv.footer_.push_back(buf);
				break;
				}
			}
		// now rays are complete, and first footer line is read, if any
		while (!f.AtEof())
			{
			f.ReadLine(buf);
			if (!buf.empty())
				rv.footer_.push_back(buf);
			}
		return rv;
		}

	TTM25RaySet ReadGenericASCIIRaySet(const std::string& fn, TASCIIRaySetOptions opts)
		// sequence of header lines, followed by sequence of ray lines, followed by sequence of footer lines
		// header lines and footer lines can be anything that cannot be converted to seven or eight floating point numbers
		// ray lines must contain seven or eight floating point numbers separated by blank or \t
		// all ray lines must contain same number of numbers (seven or eight)
		// ray lines contain x, y, z, kx, ky, kz, (power | (power, wavelength) | (wavelength, power))
		// returns TTM25RaySet. 
		// The ASCII header and footer content will be in additional_info_4_7_3_8
		// Spectral data identifier will be 0 (no wavelength) or 2 (wavelength per ray)
		{
		TRawASCIIRaySetContent rawdata = ReadRawASCIIRaySet(fn);
		TTM25Header header;
		header.n_rays_4_7_1_6 = rawdata.nRays_;
		header.file_date_time_str_4_7_1_7 = currentISO8601TimeUTC();
		if (rawdata.nItems_ == 7) // only power
			{
			header.spectrum_type_4_7_1_9 = 0; // no info
			header.lambda_4_7_1_10 = std::numeric_limits<float>::signaling_NaN();
			}
		else // nItems_ == 8 : individual wavelength
			{
			header.spectrum_type_4_7_1_9 = 2; // wavelength per ray
			header.lambda_4_7_1_10 = std::numeric_limits<float>::signaling_NaN();
			if (opts.wavelengthColumn_ != 7)
				rawdata.SwapColumns(opts.wavelengthColumn_, 7);
			auto minmax = rawdata.MinMaxColValue(7);
			header.lambda_min_4_7_1_11 = minmax.first;
			header.lambda_max_4_7_1_12 = minmax.second;
			std::basic_ostringstream<std::u32string::value_type> s32;
			s32 << U"Raw ASCII ray file read from " << ToU32String(fn) << std::endl;
			s32 << U"Header:" << std::endl;
			for (const auto& hline : rawdata.header_)
				s32 << ToU32String(hline) << std::endl;
			s32 << std::endl << U"Footer:" << std::endl;
			for (const auto& fline : rawdata.footer_)
				s32 << ToU32String(fline) << std::endl;
			}
		
		TDefaultRayArray rays(rawdata.nRays_, rawdata.nItems_, std::move(rawdata.data_));
		TTM25RaySet rv(header, std::move(rays));
		return rv;
		}

	}