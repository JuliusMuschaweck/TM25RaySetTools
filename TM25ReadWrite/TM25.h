/******************************************************************
Provided by Julius Muschaweck, JMO GmbH, Gauting to the public domain
under the Unlicense, see unlicense.txt in the repository
or http://unlicense.org/
2019-01-19
******************************************************************/
// TM25.h: The core functionality to read and write TM25 ray files.

#ifndef __TM25_H
#define __TM25_H
#include<string>
#include<vector>
#include<map>
#include<stdexcept>
#include<limits>
#include<cstdint>
// #include<chrono>
#include"ReadFile.h"

namespace TM25
	{

	enum class RayItem : size_t
		{
		x = 0, y = 1, z = 2, kx = 3, ky = 4, kz = 5,
		phi = 6, lambda = 7, Tri_Y = 8, S1 = 9, S2 = 10, S3 = 11,
		Tri_X = 12, Tri_Z = 13, spectrumIdx = 14, additional = 15
		};

	constexpr size_t nStdItems = 15;

	std::u32string RayItemToString(RayItem ri); 
	// returns U"x", U"y", .. U"additional" 

	RayItem StringToRayItem(const std::u32string& name); 
	// from name: StringToRayItem(RayItemToString(ri)) == ri
	// throws TM25Error if name does not match 

	struct TTM25Header; //forward

	class TRaySetItems
		{ // a sequence of standard items, 
		// followed by a sequence of user defined additional items
		public:
			TRaySetItems(); // all absent - empty
			explicit TRaySetItems(const TTM25Header& h);
			void MarkAsPresent(RayItem ri);
			void MarkAsAbsent(RayItem ri);
			void AddAdditionalItem(const std::u32string& name);
			bool IsPresent(RayItem ri) const;
			size_t NStdItems() const; // number of standard items present
			size_t NAdditionalItems() const; // number of additional items
			size_t NTotalItems() const; // convenience, sum of std+addtl
			RayItem ItemType(size_t i) const;
			// return the type of the i'th item in this ray set.
			// May be RayItem::additional; then you may query ItemName(i)
			// to find out what this item is.
			// throw if i >= NTotalItems()
			std::u32string ItemName(size_t i) const;
			// return the name of the i'th item in this ray set.
			// throw if i >= NTotalItems(), returns RayItemToString(ItemType(i))
			// for standard items, returns name given in AddAdditionalItem(name)
			// for additional items
			bool ContainsItems(const TRaySetItems& rhs) const;
			bool ContainsAdditionalItem(const std::u32string& name) const;
			// returns true if all items in rhs are present in this
			using TExtractionMap = std::vector<size_t>;
			TExtractionMap ExtractionMap(const TRaySetItems& rhs) const;
			// returns array of indices of this, throws TM25Error 
			// if Contains(rhs) == false
			// to be used in call to TRaySet::Extract
			static constexpr size_t absent = std::numeric_limits<size_t>::max();
			std::array<size_t, nStdItems> ItemIndices() const;
		private:
			size_t CheckRayItem(RayItem ri, const std::string& fn) const;
			std::array<bool, nStdItems> stdItems_;
			std::vector<std::u32string>	additionalItemNames_;
			size_t nTotal_;
		};

	class TM25Error : public std::runtime_error
		{
		public:
			explicit TM25Error(const std::string& msg) : std::runtime_error("TM25Error: "+msg) {};
		};

	struct TSpectralTable
		{
		int idx_; // base 1, not 0, see 4.7.4
		std::vector<float> lambda_;
		std::vector<float> weight_;
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
	

	template<typename TRayArray> 
	class TBasicTM25RaySet
		{
		public:
			TBasicTM25RaySet();
			TBasicTM25RaySet(const TTM25Header& h, const TRayArray& r);
			TBasicTM25RaySet(const TTM25Header& h, TRayArray&& r);
			
			void Read(std::string filename);
			void Write(std::string filename) const;
			
			const TRaySetItems& Items() const;
			const TRayArray& RayArray() const;
			const TTM25Header& Header() const;
			const std::vector<std::string>& Warnings() const;

			size_t NRays() const;
			std::vector<float> ExtractSingle(const TRaySetItems::TExtractionMap& em, size_t i) const;
			// returns the data of the i'th ray, containing the items in em
			// throw TM25Error if n >= NRays(), or if Items().Contains(em) == false
			TRayArray ExtractRange(const TRaySetItems::TExtractionMap& em,
				size_t i_begin, size_t i_end) const;
			// returns the data of the rays in range [i_begin, i_end[, containing the items in em.
			// note that i_end is one beyond the end of the range. The result will contain
			// i_end - i_begin rays, or no rays at all if i_end <= i_begin
			// throw TM25Error if i_end >= NRays(), or if Items().Contains(em) == false
			TRayArray ExtractAll(const TRaySetItems::TExtractionMap& em) const;
			// returns the data of the all rays, containing the items in em
			// throw TM25Error if Items().Contains(em) == false
			TRayArray ExtractSelection(const TRaySetItems::TExtractionMap& em,
				const std::vector<size_t>& idx) const;
			// returns the data of the rays according to the indices in idx, 
			// containing the items in em
			// throw TM25Error if any element of idx >= NRays(), or if Items().Contains(em) == false
		private:
			void ReadHeader(TReadFile& f);
			void ReadRayData(TReadFile& f);
			enum class RayWarnings
				{
				kNotNormalized,
				radFluxNotPositive,
				wavelengthNotPositive,
				lumFluxNotPositive,
				lumFluxNegative, // if radiant flux is present 0 is allowed
				triXNegative,
				triZNegative
				};
			enum class RayErrors
				{
				signalingNaN,
				S1NotInPlusMinusOne,
				S2NotInPlusMinusOne,
				S3NotInPlusMinusOne,
				S123NotInPlusMinusOne
				};
			std::string ToString(RayWarnings w) const;
			std::string ToString(RayErrors e) const;
			bool CheckRay(std::map<RayErrors, size_t>& err,
				std::map<RayWarnings, size_t>& warn,
				const std::vector<float>& r,
				const std::array<size_t, nStdItems>& itemIndices,
				size_t i); 
			// return true if ray is fully ok
			// sets entries in error and warning maps to i if not already set
			// thus, ray file can be fully read and first occurrence of any
			// error or warning is available
			TTM25Header header_;
			TRaySetItems items_;
			TRayArray ray_array_;
			std::vector<std::string> warnings_;
		};
	
	class TDefaultRayArray 
		// to be used as a template parameter for TBasicTM25RaySet, all public functions
		// must be provided
		// stored as a row major matrix, with rays = rows, items = columns
		// which makes it fast in C++ when rays are accessed sequentially
		{
		public:
			TDefaultRayArray();
			TDefaultRayArray(size_t nRays, size_t nItems);
			void Resize(size_t nRays, size_t nItems); // postcondition: filled with NaN
			void Clear();
			size_t NRays() const;
			size_t NItems() const;
			void SetRay(size_t i, const std::vector<float>& ray);
			// copies ray into the i'th row
			// throws TM25Error if ray.size() != nItems or i >= nRays
			void SetItem(size_t j, const std::vector<float>& item);
			// copies item into the j'th column
			// throws TM25Error if item.size() != nRays or i >= nRays
			void SetRayItem(size_t iray, size_t jitem, float r);
			// throws TM25Error if iray >= nRays(), or jitem >= nItems()
			std::vector<float> GetRay(size_t i) const;
			// returns i'th row as vector -- the safer but slower way
			// throws TM25Error if i >= nRays
			const float* GetRayDirect(size_t i) const;
			// returns i'th row as pointer to first element -- the dangerous, fast way
			// with rv = GetRayDirect(i), (rv, rv + NItems() ) is a valid range
			// throws TM25Error if i >= nRays
			const std::vector<float>& Data() const;
		private:
			size_t nRays_;
			size_t nItems_;
			std::vector<float> data_;
			friend class TBasicTM25RaySet<TDefaultRayArray>;
		};

	using TTM25RaySet = TBasicTM25RaySet<TDefaultRayArray>;
	} // namespace TM25


// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// template definitions
#include "TM25_impl.h"

	
#endif