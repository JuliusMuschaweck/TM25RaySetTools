#ifndef __TM25_H
#define __TM25_H
#include<string>
#include<vector>
#include<stdexcept>
#include<limits>
#include<cstdint>
// #include<chrono>
#include"ReadFile.h"

namespace TM25
	{

	class TM25Error : public std::runtime_error
		{
		public:
			explicit TM25Error(const std::string& msg) : std::runtime_error(msg) {};
		};

	struct TSpectralTable
		{
		int idx; // base 1, not 0, see 4.7.4
		std::vector<float> lambda;
		std::vector<float> weight;
		};

	using std::int32_t;
	using std::uint64_t;

	struct TTM25Header
		{
		int32_t	version_4_7_1_2;
		int32_t	creation_method_4_7_1_3;
		float	phi_v_4_7_1_4;
		float	phi_4_7_1_5;
		uint64_t n_rays_4_7_1_6;
		std::string	file_date_time_str_4_7_1_7;
//		std::chrono::time_point<std::chrono::system_clock> file_date_time;
		int32_t	start_position_4_7_1_8;
		int32_t	spectrum_type_4_7_1_9;
		float	lambda_4_7_1_10;
		float	lambda_min_4_7_1_11;
		float	lambda_max_4_7_1_12;
		int32_t	n_spectra_4_7_1_13;
		int32_t	n_addtl_items_4_7_1_14;
		bool	rad_flux_flag_4_7_2_3;
		bool	lambda_flag_4_7_2_4;
		bool	lum_flux_flag_4_7_2_5;
		bool	stokes_flag_4_7_2_6;
		bool	tristimulus_flag_4_7_2_7;
		bool	spectrum_index_flag_4_7_2_8;
		std::u32string	name_4_7_3_1;
		std::u32string	manufacturer_4_7_3_2;
		std::u32string	model_creator_4_7_3_3;
		std::u32string	rayfile_creator_4_7_3_4;
		std::u32string	equipment_4_7_3_5;
		std::u32string	camera_4_7_3_6;
		std::u32string	lightsource_4_7_3_7;
		std::u32string	additional_info_4_7_3_8;
		std::u32string	data_reference_4_7_3_9;
		std::vector<TSpectralTable>	spectra_4_7_4;
		std::vector<std::u32string>	column_names_4_7_5;
		std::u32string	additional_text_4_7_6;

		TTM25Header();
		};

	class TRaySet
		{
		public:
			enum class RayItem
				{
				x, y, z, kx, ky, kz, 
				phi, lambda, Y, S1, S2, S3, X, Z, SpectrumIdx
				};
			std::vector<RayItem> rayItems;
		};


	

	template<typename TRayArray> 
	class TBasicTM25RaySet
		{
		public:
			TBasicTM25RaySet();
			TBasicTM25RaySet(const TTM25Header& h, const TRayArray& r);
			TBasicTM25RaySet(const TTM25Header& h, TRayArray&& r);
			void Read(std::string filename);
			void Write(std::string filename) const;
			const TRayArray& RayArray() const;
			const TTM25Header& Header() const;
			const std::vector<std::string>& Warnings() const;
		private:
			void ReadHeader(TReadFile& f);
			TTM25Header header;
			TRayArray ray_array;
			std::vector<std::string> warnings;
		};
	
	class TDefaultRayArray
		{
		float* get(const size_t pos, const size_t num);
		const float* get(const size_t pos, const size_t num) const;
		void resize(const size_t num);
		private:
			std::vector<float> data;
			friend class TBasicTM25RaySet<TDefaultRayArray>;
		};

	using TTM25RaySet = TBasicTM25RaySet<TDefaultRayArray>;
	} // namespace TM25


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// template definitions
#include "TM25_impl.h"

	
#endif