/******************************************************************
Provided by Julius Muschaweck, JMO GmbH, Gauting to the public domain
under the Unlicense, see unlicense.txt in the repository
or http://unlicense.org/
2019-01-19
******************************************************************/

#include"TM25.h"
#include <stdexcept>
#include <limits>
#include <sstream>
#include <cassert>

namespace TM25 
	{

	class TRayItemNames
		{
		public:
			TRayItemNames()
				{
				itemStrings[RayItem::x] = U"x";
				itemStrings[RayItem::y] = U"y";
				itemStrings[RayItem::z] = U"z";
				itemStrings[RayItem::kx] = U"kx";
				itemStrings[RayItem::ky] = U"ky";
				itemStrings[RayItem::kz] = U"kz";
				itemStrings[RayItem::phi] = U"phi";
				itemStrings[RayItem::lambda] = U"lambda";
				itemStrings[RayItem::Tri_Y] = U"Tri_Y";
				itemStrings[RayItem::S1] = U"S1";
				itemStrings[RayItem::S2] = U"S2";
				itemStrings[RayItem::S3] = U"S3";
				itemStrings[RayItem::Tri_X] = U"Tri_X";
				itemStrings[RayItem::Tri_Z] = U"Tri_Z";
				itemStrings[RayItem::spectrumIdx] = U"spectrumIdx";
				itemStrings[RayItem::additional] = U"additional";
				stringItems[U"x"] = RayItem::x;
				stringItems[U"y"] = RayItem::y;
				stringItems[U"z"] = RayItem::z;
				stringItems[U"kx"] = RayItem::kx;
				stringItems[U"ky"] = RayItem::ky;
				stringItems[U"kz"] = RayItem::kz;
				stringItems[U"phi"] = RayItem::phi;
				stringItems[U"lambda"] = RayItem::lambda;
				stringItems[U"Tri_Y"] = RayItem::Tri_Y;
				stringItems[U"S1"] = RayItem::S1;
				stringItems[U"S2"] = RayItem::S2;
				stringItems[U"S3"] = RayItem::S3;
				stringItems[U"Tri_X"] = RayItem::Tri_X;
				stringItems[U"Tri_Z"] = RayItem::Tri_Z;
				stringItems[U"spectrumIdx"] = RayItem::spectrumIdx;
				stringItems[U"additional"] = RayItem::additional;

				};
			std::map<RayItem, std::u32string> itemStrings;
			std::map<std::u32string, RayItem> stringItems;
		};

	const TRayItemNames& RayItemNames()
		{ // Meyers' singleton
		static TRayItemNames names;
		return names;
		}

	// TRaySetItems implementation
	std::u32string RayItemToString(RayItem ri)
		{
		static auto to_i = [](RayItem rri) { return static_cast<size_t> (rri); };
		auto name = RayItemNames().itemStrings.find(ri);
		if (name == RayItemNames().itemStrings.end())
			{
			throw TM25Error("RayItemToString: illegal ray item "
				+ std::to_string(to_i(ri)));
			}
		return (*name).second;
		}
	// returns U"x", U"y", .. U"additional" 

	RayItem StringToRayItem(const std::u32string& name)
		{
		auto ri = RayItemNames().stringItems.find(name);
		if (ri == RayItemNames().stringItems.end())
			{
			throw TM25Error("StringToRayItem: item name " + ToString(name) + " is no standard item name");
			}
		return (*ri).second;
		}
	// from name: StringToRayItem(RayItemToString(ri)) == ri
	// throws TM25Error if name does not match 


	//	TRaySetItems implementation
	// a sequence of standard items, 
	// followed by a sequence of user defined additional items
			
	TRaySetItems::TRaySetItems()
		: nTotal_(0)
		{
		stdItems_.fill(false);
		}; 

	TRaySetItems::TRaySetItems(const TTM25Header& h)
		: nTotal_(0)
		{
		stdItems_.fill(false);
		MarkAsPresent(RayItem::x);
		MarkAsPresent(RayItem::y);
		MarkAsPresent(RayItem::z);
		MarkAsPresent(RayItem::kx);
		MarkAsPresent(RayItem::ky);
		MarkAsPresent(RayItem::kz);
		if (h.rad_flux_flag_4_7_2_3)
			MarkAsPresent(RayItem::phi);
		if (h.lambda_flag_4_7_2_4)
			MarkAsPresent(RayItem::lambda);
		if (h.lum_flux_flag_4_7_2_5)
			MarkAsPresent(RayItem::Tri_Y);
		if (h.stokes_flag_4_7_2_6)
			{
			MarkAsPresent(RayItem::S1);
			MarkAsPresent(RayItem::S2);
			MarkAsPresent(RayItem::S3);
			}
		if (h.tristimulus_flag_4_7_2_7)
			{
			MarkAsPresent(RayItem::Tri_X);
			MarkAsPresent(RayItem::Tri_Z);
			}
		if (h.spectrum_index_flag_4_7_2_8)
			MarkAsPresent(RayItem::spectrumIdx);
		for (auto s: h.column_names_4_7_5)
			AddAdditionalItem(s);
		}


	void TRaySetItems::MarkAsPresent(RayItem ri)
		{
		size_t i = CheckRayItem(ri, "MarkAsPresent");
		if (!stdItems_[i])
			{
			stdItems_[i] = true;
			++nTotal_;
			}
		}
			
	void TRaySetItems::MarkAsAbsent(RayItem ri)
		{
		size_t i = CheckRayItem(ri, "MarkAsAbsent");
		if (stdItems_[i])
			{
			stdItems_[i] = false;
			--nTotal_;
			}
		};
			
	void TRaySetItems::AddAdditionalItem(const std::u32string& name)
		{
		if (name.empty())
			throw TM25Error("TRaySetItems::AddAdditionalItem: empty name");
		if (ContainsAdditionalItem(name))
			throw TM25Error("TRaySetItems::AddAdditionalItem: duplicate name: "
				+ ToString(name));
		additionalItemNames_.push_back(name);
		++nTotal_;
		}
			
	bool TRaySetItems::IsPresent(RayItem ri) const
		{
		size_t i = CheckRayItem(ri, "IsPresent");
		return stdItems_[i];
		}
			
	size_t TRaySetItems::NStdItems() const // number of standard items present
		{
		return nTotal_ - additionalItemNames_.size();
		}
			
	size_t TRaySetItems::NAdditionalItems() const // number of additional items
		{
		return additionalItemNames_.size();
		}
			
	size_t TRaySetItems::NTotalItems() const // convenience, sum of std+addtl
		{
		return nTotal_;
		}
			
	RayItem TRaySetItems::ItemType(size_t i) const // throw if i >= NTotalItems()
		{
		if (i >= nTotal_)
			{
			std::stringstream s;
			s << "TRaySetItems::ItemType: i (" << i << ") >= NTotalItems() ("
				<< nTotal_ << ")";
			throw TM25Error(s.str());
			}
		if (i >= stdItems_.size())
			return RayItem::additional;
		for (size_t j = 0; j < nStdItems; ++j)
			{
			if (stdItems_[j])
				if (i == 0)
					return static_cast<RayItem>(j);
				else
					--i;
			}
		throw TM25Error("TRaySetItems::ItemType: logic error, this cannot happen");
		}
			
	std::u32string TRaySetItems::ItemName(size_t i) const
		{
		if (i >= nTotal_)
			{
			std::stringstream s;
			s << "TRaySetItems::ItemName: i (" << i << ") >= NTotalItems() ("
				<< nTotal_ << ")";
			throw TM25Error(s.str());
			}
		if (i >= stdItems_.size())
			return additionalItemNames_[i - stdItems_.size()];
		for (size_t j = 0; j < nStdItems; ++j)
			{
			if (stdItems_[j])
				if (i == 0)
					return RayItemToString(static_cast<RayItem>(j));
				else
					--i;
			}
		throw TM25Error("TRaySetItems::ItemName: logic error, this cannot happen");
		}
				
	bool TRaySetItems::ContainsItems(const TRaySetItems& rhs) const
		{
		for (size_t j = 0; j < nStdItems; ++j)
			if (rhs.stdItems_[j] && !(stdItems_[j]))
				return false;
		for (auto s : rhs.additionalItemNames_)
			{
			if (!ContainsAdditionalItem(s))
				return false;
			}
		return true;
		}
			
	bool TRaySetItems::ContainsAdditionalItem(const std::u32string& name) const
		{
		return
			std::find(additionalItemNames_.begin(), additionalItemNames_.end(), name)
			!= additionalItemNames_.end();
		}
			
	TRaySetItems::TExtractionMap TRaySetItems::ExtractionMap(
		const TRaySetItems& rhs) const
		{
		if (!ContainsItems(rhs))
			throw TM25Error("TRaySetItems::ExtractionMap: ContainsItems(rhs) == false");
		// using TExtractionMap = std::vector<size_t>;
		TExtractionMap rv;
		size_t thisIdx = 0;
		for (size_t j = 0; j < nStdItems; ++j)
			{
			if (stdItems_[j])
				{
				if (rhs.stdItems_[j])
					rv.push_back(thisIdx);
				++thisIdx;
				}
			}
		size_t nStd = nTotal_ - additionalItemNames_.size();
		if (thisIdx != nStd) 
			throw TM25Error("TRaySetItems::ExtractionMap: thisIdx != nStd, this cannot happen");
		for (auto s : rhs.additionalItemNames_)
			{
			for (size_t j = 0; j < additionalItemNames_.size(); ++j)
				{
				if (s == additionalItemNames_[j])
					rv.push_back(nStd + j);
				}
			}
		if (rv.size() != rhs.NTotalItems()) 
			throw TM25Error("TRaySetItems::ExtractionMap: rv.size() != rhs.NTotalItems(), this cannot happen");
		return rv;
		}
			// returns array of indices of this, throws TM25Error 
			// if Contains(rhs) == false
			// to be used in call to TRaySet::Extract
	std::array<size_t, nStdItems> TRaySetItems::ItemIndices() const
		{
		std::array<size_t, nStdItems> rv;
		rv.fill(absent);
		size_t idx = 0;
		for (size_t i = 0; i < nStdItems; ++i)
			{
			if (stdItems_[i])
				rv[i] = idx++;
			}
		return rv;
		}

	size_t TRaySetItems::CheckRayItem(RayItem ri, const std::string& fn) const
		{
		size_t rv = static_cast<size_t>(ri);
		if (rv >= stdItems_.size())
			throw TM25Error("TRaySetItems::" + fn + ": illegal RayItem: "
				+ std::to_string(rv));
		return rv;
		}

	// TTM25Header implementation
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

// ********************************************************************
	// TDefaultRayArray implementation
	TDefaultRayArray::TDefaultRayArray()
		: nRays_(0), nItems_(0) {};

	TDefaultRayArray::TDefaultRayArray(size_t nRays, size_t nItems)
		: nRays_(nRays), nItems_(nItems), 
		data_(nRays * nItems, std::numeric_limits<float>::signaling_NaN()) {};
	
	void TDefaultRayArray::Resize(size_t nRays, size_t nItems)
		{
		Clear();
		nRays_ = nRays;
		nItems_ = nItems;
		data_.resize(nRays * nItems);
		std::fill(data_.begin(), data_.end(), std::numeric_limits<float>::signaling_NaN());
		};

	void TDefaultRayArray::Clear()
		{
		nRays_ = nItems_ = 0;
		std::vector<float> tmp;
		data_.swap(tmp); // make sure previous storage is really deallocated
		}

	size_t TDefaultRayArray::NRays() const
		{
		return nRays_;
		}
	size_t TDefaultRayArray::NItems() const
		{
		return nItems_;
		}

	void TDefaultRayArray::SetRay(size_t i, const std::vector<float>& ray)
		{
		if (i >= nRays_)
			{
			std::stringstream s;
			s << "TDefaultRayArray::SetRay: i (" << i << ") >= NRays() (" << nRays_ << ")";
			throw TM25Error(s.str());
			}
		if (ray.size() != nItems_)
			{
			std::stringstream s;
			s << "TDefaultRayArray::SetRay: ray.size() (" << ray.size() << 
				") != NItems() (" << nItems_ << ")";
			throw TM25Error(s.str());
			}
		std::copy(ray.begin(), ray.end(), data_.begin() + i * NItems());
		}

	void TDefaultRayArray::SetItem(size_t j, const std::vector<float>& item)
		{
		if (j >= nItems_)
			{
			std::stringstream s;
			s << "TDefaultRayArray::SetItem: j (" << j << ") >= NItems() (" << nItems_ << ")";
			throw TM25Error(s.str());
			}
		if (item.size() != nRays_)
			{
			std::stringstream s;
			s << "TDefaultRayArray::SetItem: item.size() (" << item.size() <<
				") != NRays() (" << nRays_ << ")";
			throw TM25Error(s.str());
			}
		auto i_item = item.begin();
		for (size_t i = 0; i < nRays_; ++i)
			data_[i * nItems_ + j] = *(i_item++);
		}

	void TDefaultRayArray::SetRayItem(size_t iray, size_t jitem, float r)
		{
		if (iray >= nRays_)
			{
			std::stringstream s;
			s << "TDefaultRayArray::SetRayItem: iray (" << iray << ") >= NRays() (" 
				<< nRays_ << ")";
			throw TM25Error(s.str());
			}
		if (jitem >= nItems_)
			{
			std::stringstream s;
			s << "TDefaultRayArray::SetRayItem: jitem (" << jitem << ") >= NItems() (" 
				<< nItems_ << ")";
			throw TM25Error(s.str());
			}
		data_[iray * nItems_ + jitem] = r;
		}

	std::vector<float> TDefaultRayArray::GetRay(size_t i) const
		{
		if (i >= nRays_)
			{
			std::stringstream s;
			s << "TDefaultRayArray::GetRay: i (" << i << ") >= NRays() ("
				<< nRays_ << ")";
			throw TM25Error(s.str());
			}
		std::vector<float> rv(nItems_);
		std::copy(data_.begin() + i * nItems_, data_.begin() + (i + 1) * nItems_, rv.begin());
		return rv;
		}

	const float* TDefaultRayArray::GetRayDirect(size_t i) const
		{
		if (i >= nRays_)
			{
			std::stringstream s;
			s << "TDefaultRayArray::GetRayDirect: i (" << i << ") >= NRays() ("
				<< nRays_ << ")";
			throw TM25Error(s.str());
			}
		return data_.data() + i * nItems_;
		}

	const std::vector<float>& TDefaultRayArray::Data() const
		{
		return data_;
		}

	}