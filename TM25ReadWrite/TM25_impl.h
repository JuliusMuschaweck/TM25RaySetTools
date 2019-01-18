// template definitions
// to be included by TM25.h

#include "TM25Util.h"

namespace TM25
	{
	template<typename TRayArray> TBasicTM25RaySet<TRayArray>::TBasicTM25RaySet()
		{
		ray_array.resize(7);
		header.n_rays_4_7_1_6 = 1;
		ray_array.data[0] = 0;
		ray_array.data[1] = 0;
		ray_array.data[2] = 0;
		ray_array.data[3] = 0;
		ray_array.data[4] = 0;
		ray_array.data[5] = 1;
		ray_array.data[6] = 1;
		}

	template<typename TRayArray>
	void TBasicTM25RaySet<TRayArray>::Read(std::string filename)
		{
		TReadFile f(filename);
		ReadHeader(f);
		}

	template<typename TRayArray>
	void TBasicTM25RaySet<TRayArray>::ReadHeader(TReadFile& f)
		{
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
			header.version_4_7_1_2 = f.Read<int>();
			if (header.version_4_7_1_2 != 2013) throw TM25Error("file version should be 2013, is"
				+ std::to_string(header.version_4_7_1_2));
			// creation method, should be 0 or 1
			section = "4.7.1.3";
			header.creation_method_4_7_1_3 = f.Read<int>();
			if (!Is0or1(header.creation_method_4_7_1_3))
				warnings.push_back(section + ": creation method should be 0 or 1, is"
					+ std::to_string(header.creation_method_4_7_1_3));
			// total lum flux
			section = "4.7.1.4";
			header.phi_v_4_7_1_4 = f.Read<float>();
			if (!IsNonNegOrsNaN(header.phi_v_4_7_1_4))
				throw TM25Error("luminous flux is not nonnegative or NaN: " + std::to_string(header.phi_v_4_7_1_4));
			// total radiant flux
			section = "4.7.1.5";
			header.phi_4_7_1_5 = f.Read<float>();
			if (!IsNonNegOrsNaN(header.phi_4_7_1_5))
				throw TM25Error("radiant flux is not nonnegative or NaN: " + std::to_string(header.phi_4_7_1_5));
			// number of rays
			section = "4.7.1.6";
			header.n_rays_4_7_1_6 = f.Read<uint64_t>();
			if (header.n_rays_4_7_1_6 == 0)
				warnings.push_back(section + ": number of rays is zero");
			// file creation date and time
			section = "4.7.1.7";
			std::array<char, 28> tmp = f.Read<std::array<char, 28>>();
			AssignCharArrayToString<28>(tmp, header.file_date_time_str_4_7_1_7); // no error condition
			// ray start position
			section = "4.7.1.8";
			header.start_position_4_7_1_8 = f.Read<int>();
			if (header.start_position_4_7_1_8 < 0 || header.start_position_4_7_1_8 > 7)
				warnings.push_back(section + ": header.start_position not well defined (should be 0..7): "
				+ std::to_string(header.start_position_4_7_1_8));
			// spectral data identifier
			section = "4.7.1.9";
			header.spectrum_type_4_7_1_9 = f.Read<int>();
			if (!Is0to4( header.spectrum_type_4_7_1_9 ))
				throw TM25Error("Spectral data identifier must be in [0;4]: " + std::to_string(header.spectrum_type_4_7_1_9));
			// single wavelength
			section = "4.7.1.10";
			header.lambda_4_7_1_10 = f.Read<float>();
			if (header.spectrum_type_4_7_1_9 == 1 && !(header.lambda_4_7_1_10 > 0))
				throw TM25Error("Single wavelength must be positive number if spectral data identifier is 1: " + std::to_string(header.lambda_4_7_1_10));
			// minimum wavelength
			section = "4.7.1.11";
			header.lambda_min_4_7_1_11 = f.Read<float>();
			if (header.spectrum_type_4_7_1_9 >= 2 && !(header.lambda_min_4_7_1_11 > 0))
				warnings.push_back(section + ": minimum wavelength should be > 0 when spectral data is present: " + std::to_string(header.lambda_min_4_7_1_11));
			// maximum wavelength
			section = "4.7.1.12";
			header.lambda_max_4_7_1_12 = f.Read<float>();
			if (header.spectrum_type_4_7_1_9 >= 2 && !(header.lambda_max_4_7_1_12 > 0))
				warnings.push_back(section + ": maximum wavelength should be > 0 when spectral data is present: " + std::to_string(header.lambda_max_4_7_1_12));
			// no of spectral tables
			section = "4.7.13";
			header.n_spectra_4_7_1_13 = f.Read<int>();
			if ((header.spectrum_type_4_7_1_9 == 3) || (header.spectrum_type_4_7_1_9 == 4))
				if (header.n_spectra_4_7_1_13 <= 0)
					throw TM25Error("Spectrum identifier (" + std::to_string(header.spectrum_type_4_7_1_9) + ") requires spectral table, but there is none, # of tables is: "
						+ std::to_string(header.n_spectra_4_7_1_13));
			// no of addtl ray data items
			section = "4.7.14";
			header.n_addtl_items_4_7_1_14 = f.Read<int>();
			if (header.n_addtl_items_4_7_1_14 < 0)
				throw TM25Error("# of additional ray items must be <= 0: " + std::to_string(header.n_addtl_items_4_7_1_14));
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
			if (position_flag != 1) throw TM25Error("position flag must be 1: "+std::to_string(position_flag));
			// direction flag
			section = "4.7.2.2";
			int direction_flag = f.Read<int>();
			if (direction_flag != 1) throw TM25Error("direction flag must be 1: " + std::to_string(direction_flag));
			// radiant flux flag
			section = "4.7.2.3";
			int tmpi; // a temporary integer
			TestFlag01(tmpi = f.Read<int>(), "radiant flux flag");
			header.rad_flux_flag_4_7_2_3 = static_cast<bool>(tmpi);
			// wavelength flag
			section = "4.7.2.4";
			TestFlag01(tmpi = f.Read<int>(), "wavelength flag");
			header.lambda_flag_4_7_2_4 = static_cast<bool>(tmpi);
			if ((header.spectrum_type_4_7_1_9 == 2) && (! header.lambda_flag_4_7_2_4))
				throw TM25Error("wavelength flag must be 1 since spectrum type is 2: " + std::to_string(header.lambda_flag_4_7_2_4));
			// luminous flux flag
			section = "4.7.2.5";
			TestFlag01(tmpi = f.Read<int>(), "luminous flux flag");
			header.lum_flux_flag_4_7_2_5 = static_cast<bool>(tmpi);
			if (!(header.rad_flux_flag_4_7_2_3 || header.lum_flux_flag_4_7_2_5))
				throw TM25Error("luminous flux flag ("+ std::to_string(header.lum_flux_flag_4_7_2_5) + 
					") and radiant flux flag (" + std::to_string(header.rad_flux_flag_4_7_2_3)+") cannot be both 0");
			if ((header.spectrum_type_4_7_1_9 == 2) || (header.spectrum_type_4_7_1_9 == 4))
				if (! (header.rad_flux_flag_4_7_2_3 && (! header.lum_flux_flag_4_7_2_5)) )
					throw TM25Error("if spectrum type (" + std::to_string(header.spectrum_type_4_7_1_9) + ") is 2 or 4, then luminous flux flag (" +
						std::to_string(header.lum_flux_flag_4_7_2_5) + ") must be 0 and radiant flux flag (" + std::to_string(header.rad_flux_flag_4_7_2_3) +
						") must be 1");
			// Stokes flag
			section = "4.7.2.6";
			TestFlag01(tmpi = f.Read<int>(), "Stokes flag");
			header.stokes_flag_4_7_2_6 = static_cast<bool>(tmpi);
			if ((header.stokes_flag_4_7_2_6 == 1) && (header.rad_flux_flag_4_7_2_3 == 0))
				throw TM25Error("if Stokes flag is 1, then radiant flux flag (" + std::to_string(header.rad_flux_flag_4_7_2_3)
					+ ") must be 1");
			// tristimulus flag
			section = "4.7.2.7";
			TestFlag01(tmpi = f.Read<int>(), "tristimulus flag");
			header.tristimulus_flag_4_7_2_7 = static_cast<bool>(tmpi);
			if (header.tristimulus_flag_4_7_2_7)
				{
				if (! header.lum_flux_flag_4_7_2_5 )
					throw TM25Error("if tristimulus flag is 1, then luminous flux flag (" + std::to_string(header.lum_flux_flag_4_7_2_5)
						+ ") must be 1");
				if (header.spectrum_type_4_7_1_9 != 0)
					throw TM25Error("if tristimulus flag is 1, then spectrum type (" + std::to_string(header.spectrum_type_4_7_1_9)
						+ ") must be 0");
				}
			// spectrum index flag
			section = "4.7.2.8";
			TestFlag01(tmpi = f.Read<int>(), "spectrum index flag");
			header.spectrum_index_flag_4_7_2_8 = static_cast<bool>(tmpi);
			if (header.spectrum_index_flag_4_7_2_8 != (header.spectrum_type_4_7_1_9 == 4))
				throw TM25Error("spectrum index flag (" + std::to_string(header.spectrum_index_flag_4_7_2_8) +
					") set iff spectrum type == 4 (" + std::to_string(header.spectrum_type_4_7_1_9) + ")");
			// description header block
			// name of light source
			section = "4.7.3.1";
			AssignChar32ArrayToString(f.Read<char32_t, 1000>(), header.name_4_7_3_1);
			// manufacturer of light source
			section = "4.7.3.2";
			AssignChar32ArrayToString(f.Read<char32_t, 1000>(), header.manufacturer_4_7_3_2);
			// creator of light source model
			section = "4.7.3.3";
			AssignChar32ArrayToString(f.Read<char32_t, 1000>(), header.model_creator_4_7_3_3);
			// creator of ray file
			section = "4.7.3.4";
			AssignChar32ArrayToString(f.Read<char32_t, 1000>(), header.rayfile_creator_4_7_3_4);
			// equipment / software
			section = "4.7.3.5";
			AssignChar32ArrayToString(f.Read<char32_t, 1000>(), header.equipment_4_7_3_5);
			// camera
			section = "4.7.3.6";
			AssignChar32ArrayToString(f.Read<char32_t, 1000>(), header.camera_4_7_3_6);
			// light source operation
			section = "4.7.3.7";
			AssignChar32ArrayToString(f.Read<char32_t, 1000>(), header.lightsource_4_7_3_7);
			// additional information
			section = "4.7.3.8";
			AssignChar32ArrayToString(f.Read<char32_t, 1000>(), header.additional_info_4_7_3_8);
			// data reference
			section = "4.7.3.9";
			AssignChar32ArrayToString(f.Read<char32_t, 1000>(), header.data_reference_4_7_3_9);
			// spectral tables block
			section = "4.7.4";
			int M = header.n_spectra_4_7_1_13;
			int bytecount = 0;
			for (int i = 1; i <= M; ++i)
				{
				int np = f.Read<int>();
				if (np <= 0)
					throw TM25Error("number of data pairs in spectral table " +
						std::to_string(i) + " (base 1) must be > 0: " + std::to_string(np));
				TSpectralTable st;
				st.idx = i;
				st.lambda.reserve(np);
				st.weight.reserve(np);
				for (int j = 0; j < np; ++j)
					{
					float tmp;
					st.lambda.push_back(tmp = f.Read<float>());
					if (tmp <= 0.0f)
						throw TM25Error("wavelengths in spectral tables must be > 0, violated by " +
							std::to_string(tmp) + " at base 0 position " + std::to_string(j) +
							" in spectral table base 1 # " + std::to_string(i));
					st.weight.push_back(tmp = f.Read<float>());
					if (tmp < 0.0f)
						throw TM25Error("weights in spectral tables must be >= 0, violated by " +
							std::to_string(tmp) + " at base 0 position " + std::to_string(j) +
							" in spectral table base 1 # " + std::to_string(i));
					}
				bytecount += (4 + np * 8);
				header.spectra_4_7_4.push_back(std::move(st));
				}
			if ((bytecount % 32) > 0) // padding 
				f.ReadBytes(32 - (bytecount % 32));
			// addtl column header block
			section = "4.7.5";
			//
			section = "4.7.";
			//
			section = "4.7.";
			//
			section = "4.7.";
			//
			section = "4.7.";
			}
		catch (std::runtime_error err)
			{
			throw TM25Error("TM25 Error: wrong file format in section " + section + ", " + err.what());
			}
		catch (...)
			{
			throw TM25Error("TM25 Error: unknown exception in section " + section);
			}
		}
	}