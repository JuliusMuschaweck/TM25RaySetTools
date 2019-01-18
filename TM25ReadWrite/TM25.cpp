#include"TM25.h"
#include <stdexcept>
#include <limits>
namespace TM25 
	{

	TTM25Header::TTM25Header()
		{
		version_4_7_1_2 = 2013; // must be
		creation_method_4_7_1_3 = 0; // Simulation
		phi_v_4_7_1_4 = std::numeric_limits<float>::signaling_NaN(); // only radiant flux
		phi_4_7_1_5 = 1;
		n_rays_4_7_1_6 = 0;
		file_date_time_str_4_7_1_7 = "2013-09-04T08:30:29+01:00"; // example from TM25
		struct tm tm_time;
		tm_time.tm_year = 2013 - 1900; // years since 1900
		tm_time.tm_mon = 9 - 1; // Jan = 0
		tm_time.tm_mday = 4; // 1..31
		tm_time.tm_hour = 8;
		tm_time.tm_min = 30;
		tm_time.tm_sec = 29;
		tm_time.tm_isdst = 0;
		// file_date_time = std::chrono::system_clock::from_time_t(mktime(&tm_time));
		start_position_4_7_1_8 = 0; // unknown
		spectrum_type_4_7_1_9 = 0; // no wavelength information
		lambda_4_7_1_10 = std::numeric_limits<float>::signaling_NaN();
		lambda_min_4_7_1_11 = std::numeric_limits<float>::signaling_NaN();
		lambda_max_4_7_1_12 = std::numeric_limits<float>::signaling_NaN();
		n_spectra_4_7_1_13 = 0; // no spectra
		n_addtl_items_4_7_1_14 = 0; // no additional items
		rad_flux_flag_4_7_2_3 = true;
		lambda_flag_4_7_2_4 = false;
		lum_flux_flag_4_7_2_5 = false;
		stokes_flag_4_7_2_6 = false;
		tristimulus_flag_4_7_2_7 = false;
		spectrum_index_flag_4_7_2_8 = false;
		name_4_7_3_1 = U"default";
		manufacturer_4_7_3_2 = U"unknown";
		model_creator_4_7_3_3 = U"unknown";
		rayfile_creator_4_7_3_4 = U"unknown";
		equipment_4_7_3_5 = U"unknown";
		camera_4_7_3_6 = U"unknown";
		lightsource_4_7_3_7 = U"unknown";
		additional_info_4_7_3_8 = U"none";
		data_reference_4_7_3_9 = U"unknown";
		// std::vector<TSpectralTable>	spectra_4_7_4;
		// std::vector<std::u32string>	column_names_4_7_5;
		// std::u32string	additional_text_4_7_6;
		}

	float* TDefaultRayArray::get(const size_t pos, const  size_t num)
		{
		if ((pos + num) < data.size())
			return data.data() + pos;
		else
			throw std::runtime_error("TDefaultRayArray::get: pos + num out of bounds");
		}

	const float* TDefaultRayArray::get(const size_t pos, const size_t num) const
		{
		if ((pos + num) < data.size())
			return data.data() + pos;
		else
			throw std::runtime_error("TDefaultRayArray::get: pos + num out of bounds");
		}

	void TDefaultRayArray::resize(const size_t num)
		{
		data.resize(num);
		}

	}