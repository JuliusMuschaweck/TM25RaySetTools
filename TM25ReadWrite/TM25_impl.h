// template definitions
// to be included by TM25.h

#include "TM25Util.h"
#include <sstream>
namespace TM25
	{

	template<typename TRayArray>
	TBasicTM25RaySet<TRayArray>::TBasicTM25RaySet()
		{
		ray_array_.Resize(1, 7);
		header_.n_rays_4_7_1_6 = 1;
		ray_array_.SetRay(0, { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 1.0f });
		}

	template<typename TRayArray>
	TBasicTM25RaySet<TRayArray>::TBasicTM25RaySet(const TTM25Header& h, const TRayArray& r)
		: header_(h), items_(h), ray_array_(r) {};

	template<typename TRayArray>
	TBasicTM25RaySet<TRayArray>::TBasicTM25RaySet(const TTM25Header& h, TRayArray&& r)
		: header_(h), items_(h), ray_array_(std::move(r)) {};


	template<typename TRayArray>
	void TBasicTM25RaySet<TRayArray>::Read(std::string filename)
		{
		warnings_.clear();
		TReadFile f(filename);
		ReadHeader(f);
		items_ = TRaySetItems(header_);
		ReadRayData(f);
		}

	template<typename TRayArray>
	const TRaySetItems& TBasicTM25RaySet<TRayArray>::Items() const
		{
		return items_;
		}

	template<typename TRayArray>
	const TRayArray& TBasicTM25RaySet<TRayArray>::RayArray() const
		{
		return ray_array_;
		}

	template<typename TRayArray>
	const TTM25Header& TBasicTM25RaySet<TRayArray>::Header() const
		{
		return header_;
		}

	template<typename TRayArray>
	const std::vector<std::string>& TBasicTM25RaySet<TRayArray>::Warnings() const
		{
		return warnings_;
		}

	template<typename TRayArray>
	size_t TBasicTM25RaySet<TRayArray>::NRays() const
		{
		if (header_.n_rays_4_7_1_6 != ray_array_.NRays())
			{
			std::stringstream s;
			s << "TBasicTM25RaySet<TRayArray>::NRays(): Inconsistent ray data size, "
				<< "header_.n_rays_4_7_1_6 says " << header_.n_rays_4_7_1_6
				<< ", ray_array_.NRays() says " << ray_array_.NRays();
			throw TM25Error(s.str());
			}
		return header_.n_rays_4_7_1_6;
		}

	template<typename TRayArray>
	void TBasicTM25RaySet<TRayArray>::ReadHeader(TReadFile& f)
		{
		header_ = TTM25Header();
		std::string section;
		auto Is0or1 = [](int i) {return (i == 0) || (i == 1); };
		auto IsNonNegOrsNaN = [](float i) {return (i == 0.0f) || (std::isnormal(i) && i > 0.0f) ||
			(i == std::numeric_limits<float>::signaling_NaN()); };
		auto IsPosOrsNaN = [](float i) {return (i > 0.0f) ||
			(i == std::numeric_limits<float>::signaling_NaN()); };
		auto Is0to4 = [](int i) {return (i >= 0) && (i <= 4); };
		auto TestFlag01 = [Is0or1](int i, const std::string& s) {if (!Is0or1(i))
			throw TM25Error(s + " (" + std::to_string(i) + ") must be 0 or 1"); };
		try
			{
			// file type, must be "TM25"
			section = "4.7.1.1";
			std::array<char, 4> ftype = f.Read<char, 4>();
			std::array<char, 4> ftype_std{ 'T','M','2','5' };
			if (ftype != ftype_std) throw TM25Error("file type should be TM25, is "
				+ std::string(ftype.begin(), ftype.end()));
			// file version, must be 2013
			section = "4.7.1.2";
			header_.version_4_7_1_2 = f.Read<int>();
			if (header_.version_4_7_1_2 != 2013) throw TM25Error("file version should be 2013, is "
				+ std::to_string(header_.version_4_7_1_2));
			// creation method, should be 0 or 1
			section = "4.7.1.3";
			header_.creation_method_4_7_1_3 = f.Read<int>();
			if (!Is0or1(header_.creation_method_4_7_1_3))
				warnings_.push_back(section + ": creation method should be 0 or 1, is "
					+ std::to_string(header_.creation_method_4_7_1_3));
			// total lum flux
			section = "4.7.1.4";
			header_.phi_v_4_7_1_4 = f.Read<float>();
			if (!IsNonNegOrsNaN(header_.phi_v_4_7_1_4))
				throw TM25Error("luminous flux is not nonnegative or NaN: " + std::to_string(header_.phi_v_4_7_1_4));
			// total radiant flux
			section = "4.7.1.5";
			header_.phi_4_7_1_5 = f.Read<float>();
			if (!IsNonNegOrsNaN(header_.phi_4_7_1_5))
				throw TM25Error("radiant flux is not nonnegative or NaN: " + std::to_string(header_.phi_4_7_1_5));
			// number of rays
			section = "4.7.1.6";
			header_.n_rays_4_7_1_6 = f.Read<uint64_t>();
			if (header_.n_rays_4_7_1_6 == 0)
				warnings_.push_back(section + ": number of rays is zero");
			// file creation date and time
			section = "4.7.1.7";
			std::array<char, 28> tmp = f.Read<std::array<char, 28>>();
			AssignCharArrayToString<28>(tmp, header_.file_date_time_str_4_7_1_7); // no error condition
			// ray start position
			section = "4.7.1.8";
			header_.start_position_4_7_1_8 = f.Read<int>();
			if (header_.start_position_4_7_1_8 < 0 || header_.start_position_4_7_1_8 > 7)
				warnings_.push_back(section + ": header_.start_position not well defined (should be 0..7): "
					+ std::to_string(header_.start_position_4_7_1_8));
			// spectral data identifier
			section = "4.7.1.9";
			header_.spectrum_type_4_7_1_9 = f.Read<int>();
			if (!Is0to4(header_.spectrum_type_4_7_1_9))
				throw TM25Error("Spectral data identifier must be in [0;4]: " + std::to_string(header_.spectrum_type_4_7_1_9));
			// single wavelength
			section = "4.7.1.10";
			header_.lambda_4_7_1_10 = f.Read<float>();
			if (header_.spectrum_type_4_7_1_9 == 1 && !(header_.lambda_4_7_1_10 > 0))
				throw TM25Error("Single wavelength must be positive number if spectral data identifier is 1: " + std::to_string(header_.lambda_4_7_1_10));
			// minimum wavelength
			section = "4.7.1.11";
			header_.lambda_min_4_7_1_11 = f.Read<float>();
			if (header_.spectrum_type_4_7_1_9 >= 2 && !(header_.lambda_min_4_7_1_11 > 0))
				warnings_.push_back(section + ": minimum wavelength should be > 0 when spectral data is present: " + std::to_string(header_.lambda_min_4_7_1_11));
			// maximum wavelength
			section = "4.7.1.12";
			header_.lambda_max_4_7_1_12 = f.Read<float>();
			if (header_.spectrum_type_4_7_1_9 >= 2 && !(header_.lambda_max_4_7_1_12 > 0))
				warnings_.push_back(section + ": maximum wavelength should be > 0 when spectral data is present: " + std::to_string(header_.lambda_max_4_7_1_12));
			// no of spectral tables
			section = "4.7.13";
			header_.n_spectra_4_7_1_13 = f.Read<int>();
			if ((header_.spectrum_type_4_7_1_9 == 3) || (header_.spectrum_type_4_7_1_9 == 4))
				if (header_.n_spectra_4_7_1_13 <= 0)
					throw TM25Error("Spectrum identifier (" + std::to_string(header_.spectrum_type_4_7_1_9) + ") requires spectral table, but there is none, # of tables is: "
						+ std::to_string(header_.n_spectra_4_7_1_13));
			// no of addtl ray data items
			section = "4.7.14";
			header_.n_addtl_items_4_7_1_14 = f.Read<int>();
			if (header_.n_addtl_items_4_7_1_14 < 0)
				throw TM25Error("# of additional ray items must be <= 0: " + std::to_string(header_.n_addtl_items_4_7_1_14));
			// size of addtl text block
			section = "4.7.15";
			int addtl_textblock_size = f.Read<int>();
			if ((addtl_textblock_size < 0) || (addtl_textblock_size % 32) != 0)
				throw TM25Error("additional text block size must be nonnegative multiple of 32: " + std::to_string(addtl_textblock_size));
			// reserved for future use
			section = "4.7.16";
			f.ReadBytes<168>();
			// Known data flags block 4.7.2
			// Position flag
			section = "4.7.2.1";
			int position_flag = f.Read<int>();
			if (position_flag != 1) throw TM25Error("position flag must be 1: " + std::to_string(position_flag));
			// direction flag
			section = "4.7.2.2";
			int direction_flag = f.Read<int>();
			if (direction_flag != 1) throw TM25Error("direction flag must be 1: " + std::to_string(direction_flag));
			// radiant flux flag
			section = "4.7.2.3";
			int tmpi; // a temporary integer
			TestFlag01(tmpi = f.Read<int>(), "radiant flux flag");
			header_.rad_flux_flag_4_7_2_3 = static_cast<bool>(tmpi);
			// wavelength flag
			section = "4.7.2.4";
			TestFlag01(tmpi = f.Read<int>(), "wavelength flag");
			header_.lambda_flag_4_7_2_4 = static_cast<bool>(tmpi);
			if ((header_.spectrum_type_4_7_1_9 == 2) && (!header_.lambda_flag_4_7_2_4))
				throw TM25Error("wavelength flag must be 1 since spectrum type is 2: " + std::to_string(header_.lambda_flag_4_7_2_4));
			// luminous flux flag
			section = "4.7.2.5";
			TestFlag01(tmpi = f.Read<int>(), "luminous flux flag");
			header_.lum_flux_flag_4_7_2_5 = static_cast<bool>(tmpi);
			if (!(header_.rad_flux_flag_4_7_2_3 || header_.lum_flux_flag_4_7_2_5))
				throw TM25Error("luminous flux flag (" + std::to_string(header_.lum_flux_flag_4_7_2_5) +
					") and radiant flux flag (" + std::to_string(header_.rad_flux_flag_4_7_2_3) + ") cannot be both 0");
			if ((header_.spectrum_type_4_7_1_9 == 2) || (header_.spectrum_type_4_7_1_9 == 4))
				if (!(header_.rad_flux_flag_4_7_2_3 && (!header_.lum_flux_flag_4_7_2_5)))
					throw TM25Error("if spectrum type (" + std::to_string(header_.spectrum_type_4_7_1_9) + ") is 2 or 4, then luminous flux flag (" +
						std::to_string(header_.lum_flux_flag_4_7_2_5) + ") must be 0 and radiant flux flag (" + std::to_string(header_.rad_flux_flag_4_7_2_3) +
						") must be 1");
			// Stokes flag
			section = "4.7.2.6";
			TestFlag01(tmpi = f.Read<int>(), "Stokes flag");
			header_.stokes_flag_4_7_2_6 = static_cast<bool>(tmpi);
			if ((header_.stokes_flag_4_7_2_6 == 1) && (header_.rad_flux_flag_4_7_2_3 == 0))
				throw TM25Error("if Stokes flag is 1, then radiant flux flag (" + std::to_string(header_.rad_flux_flag_4_7_2_3)
					+ ") must be 1");
			// tristimulus flag
			section = "4.7.2.7";
			TestFlag01(tmpi = f.Read<int>(), "tristimulus flag");
			header_.tristimulus_flag_4_7_2_7 = static_cast<bool>(tmpi);
			if (header_.tristimulus_flag_4_7_2_7)
				{
				if (!header_.lum_flux_flag_4_7_2_5)
					throw TM25Error("if tristimulus flag is 1, then luminous flux flag (" + std::to_string(header_.lum_flux_flag_4_7_2_5)
						+ ") must be 1");
				if (header_.spectrum_type_4_7_1_9 != 0)
					throw TM25Error("if tristimulus flag is 1, then spectrum type (" + std::to_string(header_.spectrum_type_4_7_1_9)
						+ ") must be 0");
				}
			// spectrum index flag
			section = "4.7.2.8";
			TestFlag01(tmpi = f.Read<int>(), "spectrum index flag");
			header_.spectrum_index_flag_4_7_2_8 = static_cast<bool>(tmpi);
			if (header_.spectrum_index_flag_4_7_2_8 != (header_.spectrum_type_4_7_1_9 == 4))
				throw TM25Error("spectrum index flag (" + std::to_string(header_.spectrum_index_flag_4_7_2_8) +
					") set iff spectrum type == 4 (" + std::to_string(header_.spectrum_type_4_7_1_9) + ")");
			// description header block
			// name of light source
			section = "4.7.3.1";
			AssignChar32ArrayToString(f.Read<char32_t, 1000>(), header_.name_4_7_3_1);
			// manufacturer of light source
			section = "4.7.3.2";
			AssignChar32ArrayToString(f.Read<char32_t, 1000>(), header_.manufacturer_4_7_3_2);
			// creator of light source model
			section = "4.7.3.3";
			AssignChar32ArrayToString(f.Read<char32_t, 1000>(), header_.model_creator_4_7_3_3);
			// creator of ray file
			section = "4.7.3.4";
			AssignChar32ArrayToString(f.Read<char32_t, 1000>(), header_.rayfile_creator_4_7_3_4);
			// equipment / software
			section = "4.7.3.5";
			AssignChar32ArrayToString(f.Read<char32_t, 1000>(), header_.equipment_4_7_3_5);
			// camera
			section = "4.7.3.6";
			AssignChar32ArrayToString(f.Read<char32_t, 1000>(), header_.camera_4_7_3_6);
			// light source operation
			section = "4.7.3.7";
			AssignChar32ArrayToString(f.Read<char32_t, 1000>(), header_.lightsource_4_7_3_7);
			// additional information
			section = "4.7.3.8";
			AssignChar32ArrayToString(f.Read<char32_t, 1000>(), header_.additional_info_4_7_3_8);
			// data reference
			section = "4.7.3.9";
			AssignChar32ArrayToString(f.Read<char32_t, 1000>(), header_.data_reference_4_7_3_9);
			// spectral tables block
			section = "4.7.4";
			int M = header_.n_spectra_4_7_1_13;
			int bytecount = 0;
			for (int i = 1; i <= M; ++i)
				{
				int np = f.Read<int>();
				if (np <= 0)
					throw TM25Error("number of data pairs in spectral table " +
						std::to_string(i) + " (base 1) must be > 0: " + std::to_string(np));
				TSpectralTable st;
				st.idx_ = i;
				st.lambda_.reserve(np);
				st.weight_.reserve(np);
				for (int j = 0; j < np; ++j)
					{
					float tmp;
					st.lambda_.push_back(tmp = f.Read<float>());
					if (tmp <= 0.0f)
						throw TM25Error("wavelengths in spectral tables must be > 0, violated by " +
							std::to_string(tmp) + " at base 0 position " + std::to_string(j) +
							" in spectral table base 1 # " + std::to_string(i));
					st.weight_.push_back(tmp = f.Read<float>());
					if (tmp < 0.0f)
						throw TM25Error("weights in spectral tables must be >= 0, violated by " +
							std::to_string(tmp) + " at base 0 position " + std::to_string(j) +
							" in spectral table base 1 # " + std::to_string(i));
					}
				bytecount += (4 + np * 8);
				header_.spectra_4_7_4.push_back(std::move(st));
				}
			if ((bytecount % 32) > 0) // padding 
				f.ReadBytes(32 - (bytecount % 32));
			// addtl column header block
			section = "4.7.5";
			M = header_.n_addtl_items_4_7_1_14;
			for (int i = 0; i < M; ++i)
				{
				constexpr size_t N = 512 / sizeof(char32_t);
				std::u32string tmps;
				if (tmps.empty())
					throw TM25Error("each additional column must have a name with nonzero length");
				AssignChar32ArrayToString(f.Read<char32_t, N>(), tmps);
				header_.column_names_4_7_5.push_back(std::move(tmps));
				}
			// addtl text block
			section = "4.7.6";
			if (addtl_textblock_size > 0)
				{
				size_t N = addtl_textblock_size / sizeof(char32_t);
				AssignChar32ArrayToString(f.Read<char32_t>(N), header_.additional_text_4_7_6);
				}
			}
		catch (std::runtime_error err)
			{
			throw TM25Error("wrong file format in section " + section + ", " + err.what());
			}
		catch (...)
			{
			throw TM25Error("unknown exception in section " + section);
			}
		}

	template<typename TRayArray>
	void TBasicTM25RaySet<TRayArray>::ReadRayData(TReadFile& f)
		{
		size_t nRays = header_.n_rays_4_7_1_6;
		size_t nItems = items_.NTotalItems();
		ray_array_.Resize(nRays, nItems);
		std::vector<float> rayHolder(nItems);
		std::map<RayErrors, size_t> errMap;
		std::map<RayWarnings, size_t> warnMap;
		std::array<size_t, nStdItems> itemIndices = items_.ItemIndices();
		for (size_t i = 0; i < nRays; ++i)
			{
			f.Read<float>(rayHolder, nItems);
			CheckRay(errMap, warnMap, rayHolder, itemIndices, i);
			ray_array_.SetRay(i, rayHolder);
			}
		for (auto w : warnMap)
			{
			warnings_.push_back("TBasicTM25RaySet<TRayArray>::ReadRayData: " + ToString(w.first) +
				", first occurence at ray # (base 0) " + std::to_string(w.second));
			}
		std::stringstream s;
		for (auto e : errMap)
			{
			if (!s.str().empty()) s << "\n";
			s << ToString(e.first) << ", first occurence at ray # (base 0) " << std::to_string(e.second);
			}
		if (!s.str().empty())
			throw TM25Error("TBasicTM25RaySet<TRayArray>::ReadRayData:\n" + s.str());
		}

	template<typename TRayArray>
	std::string TBasicTM25RaySet<TRayArray>::ToString(RayWarnings w) const
		{
		switch (w)
			{
			case RayWarnings::kNotNormalized: return ("kNotNormalized");
			case RayWarnings::radFluxNotPositive: return ("radFluxNotPositive");
			case RayWarnings::wavelengthNotPositive: return ("wavelengthNotPositive");
			case RayWarnings::lumFluxNotPositive: return ("lumFluxNotPositive");
			case RayWarnings::lumFluxNegative: return ("lumFluxNegative");// if radiant flux is present 0 is allowed
			case RayWarnings::triXNegative: return ("triXNegative");
			case RayWarnings::triZNegative: return ("triZNegative");
			default: throw TM25Error
					 ("TBasicTM25RaySet<TRayArray>::ToString(RayWarnings w): unknown w ("
						 + std::to_string(static_cast<size_t>(w)) + "), this cannot happen");
			}
		}
	
	template<typename TRayArray>
	std::string TBasicTM25RaySet<TRayArray>::ToString(RayErrors e) const
		{
		switch (e)
			{
			case RayErrors::signalingNaN: return ("signalingNaN");
			case RayErrors::S1NotInPlusMinusOne: return ("S1NotInPlusMinusOne");
			case RayErrors::S2NotInPlusMinusOne: return ("S2NotInPlusMinusOne");
			case RayErrors::S3NotInPlusMinusOne: return ("S3NotInPlusMinusOne");
			case RayErrors::S123NotInPlusMinusOne: return ("S123NotInPlusMinusOne");
			default: throw TM25Error
					 ("TBasicTM25RaySet<TRayArray>::ToString(RayErrors e): unknown e ("
						 + std::to_string(static_cast<size_t>(e)) + "), this cannot happen");
			}
		}


	template<typename TRayArray>
	bool TBasicTM25RaySet<TRayArray>::CheckRay(
		std::map<RayErrors, size_t>& err,
		std::map<RayWarnings, size_t>& warn,
		const std::vector<float>& r,
		const std::array<size_t, nStdItems>& itemIndices,
		size_t i)
		{
		float eps = 10 * std::numeric_limits<float>::epsilon();
		bool rv = true;
		// some lambdas for more concise code below
		// set the error flag to the first line of occurrence
		auto setErr = [&err, i, &rv](RayErrors errtype)
			{
			if (err.find(errtype) == err.end())
				{
				err.insert({ errtype, i });
				rv = false;
				}
			};
		// set the warning flag to the first line of occurrence
		auto setWarn = [&warn, i, &rv](RayWarnings warntype)
			{
			if (warn.find(warntype) == warn.end())
				{
				warn.insert({ warntype, i });
				rv = false;
				}
			};
		// return the square of a float
		auto sqr = [](float f) {return f * f; };
		// convert RayItem from enum to size_t
		auto to_i = [](RayItem ri) {return static_cast<size_t>(ri); };
		// see if a RayItem is present
		auto present = [&itemIndices, &to_i](RayItem ri)
			{return itemIndices[to_i(ri)] != TRaySetItems::absent; };
		// get the value of a RayItem, precondition: present
		auto val = [&r, &to_i, &itemIndices](RayItem ri)
			{
			size_t idx = itemIndices[to_i(ri)];
			return r[idx];
			};
		//kNotNormalized,
		if (present(RayItem::kx) && present(RayItem::ky) && present(RayItem::kz))
			{
			float k = sqr(val(RayItem::kx)) + sqr(val(RayItem::ky)) + sqr(val(RayItem::kz));
			if (abs(k - 1.0f) > eps)
				setWarn(RayWarnings::kNotNormalized);
			}
		//radFluxNotPositive,
		if (present(RayItem::phi))
			if (val(RayItem::phi) <= 0)
				setWarn(RayWarnings::radFluxNotPositive);
		//wavelengthNotPositive,
		if (present(RayItem::lambda))
			if (val(RayItem::lambda) <= 0)
				setWarn(RayWarnings::wavelengthNotPositive);
		//lumFluxNotPositive,
		if (present(RayItem::Tri_Y))
			{
			if (!present(RayItem::phi) && val(RayItem::Tri_Y) <= 0)
				setWarn(RayWarnings::lumFluxNotPositive);
			//lumFluxNegative, // if radiant flux is present 0 is allowed
			if (present(RayItem::phi) && val(RayItem::Tri_Y) < 0)
				setWarn(RayWarnings::lumFluxNegative);
			}
		//triXNegative,
		if (present(RayItem::Tri_X))
			if (val(RayItem::Tri_X) < 0)
				setWarn(RayWarnings::triXNegative);
		//triZNegative
		if (present(RayItem::Tri_Z))
			if (val(RayItem::Tri_Z) < 0)
				setWarn(RayWarnings::triZNegative);
		// no NaNs
		for (float a : r)
			{
			if (std::isnan(a))
				setErr(RayErrors::signalingNaN);
			}
		//S1NotInPlusMinusOne,
		bool S1p = present(RayItem::S1);
		bool S2p = present(RayItem::S2);
		bool S3p = present(RayItem::S3);
		float S1v;
		float S2v;
		float S3v;
		if (S1p)
			if (abs(S1v = val(RayItem::S1)) > 1)
				setErr(RayErrors::S1NotInPlusMinusOne);
		//S2NotInPlusMinusOne,
		if (S2p)
			if (abs(S2v = val(RayItem::S2)) > 1)
				setErr(RayErrors::S2NotInPlusMinusOne);
		//S3NotInPlusMinusOne,
		if (S3p)
			if (abs(S3v = val(RayItem::S3)) > 1)
				setErr(RayErrors::S3NotInPlusMinusOne);
		//S123NotInPlusMinusOne,
		if (S1p && S2p && S3p)
			{
			float f = sqr(S1v) + sqr(S2v) + sqr(S3v);
			if (f > 1.0f)
				setErr(RayErrors::S123NotInPlusMinusOne);
			}
		return rv;
		}

	}