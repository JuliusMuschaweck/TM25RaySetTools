/******************************************************************
Provided by Julius Muschaweck, JMO GmbH, Gauting to the public domain
under the Unlicense, see unlicense.txt in the repository
or http://unlicense.org/
2019-01-19
******************************************************************/
#ifndef _CONSOLE
#include <stdafx.h>
#endif#undef min
#undef max

#include "TM25.h"
#include <stdexcept>
#include <limits>
#include <sstream>
#include <cassert>
#include <random>
#include "ZemaxBinary.h"
#include <chrono>
#include "Timer.h"
#include <limits>
#include <locale>
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
				itemStrings[RayItem::PolEllipseX] = U"PolEllipseX";
				itemStrings[RayItem::PolEllipseY] = U"PolEllipseY";
				itemStrings[RayItem::PolEllipseZ] = U"PolEllipseZ";
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
				stringItems[U"PolEllipseX"] = RayItem::PolEllipseX;
				stringItems[U"PolEllipseY"] = RayItem::PolEllipseY;
				stringItems[U"PolEllipseZ"] = RayItem::PolEllipseZ;
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
			MarkAsPresent(RayItem::PolEllipseX);
			MarkAsPresent(RayItem::PolEllipseY);
			MarkAsPresent(RayItem::PolEllipseZ);
			}
		if (h.tristimulus_flag_4_7_2_7)
			{
			MarkAsPresent(RayItem::Tri_X);
			MarkAsPresent(RayItem::Tri_Z);
			}
		if (h.spectrum_index_flag_4_7_2_8)
			MarkAsPresent(RayItem::spectrumIdx);
		for (auto s : h.column_names_4_7_5)
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

	TTM25Header::TSanityCheck::TSanityCheck()
		: msg(), nonfatalErrors(false), fatalErrors(false)
		{
		}

	void TTM25Header::TSanityCheck::Fatal(bool test, const std::string & s)
		{
		if (!test) return;
		msg += "fatal error: " + s + "\n";
		fatalErrors = true;
		}

	void TTM25Header::TSanityCheck::NonFatal(bool test, const std::string & s)
		{
		if (!test) return;
		msg += "nonfatal error: " + s + "\n";
		nonfatalErrors = true;
		}


	const TTM25Header::TSanityCheck TTM25Header::SanityCheck() const
		{
		TSanityCheck rv;
		// 4.7.1 file header block
		rv.Fatal(version_4_7_1_2 != 2013,
			"4.7.1.2: version is not 2013 -> " + std::to_string(version_4_7_1_2));
		auto goodFlux = [](float r)
			{return (r == 0.0f || std::isnan(r) || (std::isnormal(r) && r > 0.0f)); };
		rv.Fatal(!goodFlux(phi_v_4_7_1_4),
			"4.7.1.4: luminous flux is not positive normalized, zero or NaN -> " + std::to_string(phi_v_4_7_1_4));
		rv.Fatal(!goodFlux(phi_4_7_1_5),
			"4.7.1.5: radiant flux is not positive normalized, zero or NaN -> " + std::to_string(phi_4_7_1_5));
		rv.Fatal(n_rays_4_7_1_6 == 0, "4.7.1.6: # of rays is zero");
		rv.Fatal(spectrum_type_4_7_1_9 < 0 || spectrum_type_4_7_1_9 > 4,
			"4.7.1.9: spectrum type is not 0,1,2,3,4 -> " + std::to_string(spectrum_type_4_7_1_9));
		rv.Fatal(spectrum_type_4_7_1_9 == 1 && !(lambda_4_7_1_10 > 0),
			"4.7.1.10: spectrum type is 1 (single wavelength) but wavelength is not positive -> "
			+ std::to_string(lambda_4_7_1_10));
		bool minmaxlam = (spectrum_type_4_7_1_9 >= 2) && (spectrum_type_4_7_1_9 <= 4);
		rv.NonFatal(minmaxlam && !(lambda_min_4_7_1_11 > 0), "4.7.1.11: min. wavelength not positive -> "
			+ std::to_string(lambda_min_4_7_1_11));
		rv.NonFatal(minmaxlam && !(lambda_max_4_7_1_12 > 0), "4.7.1.12: max. wavelength not positive -> "
			+ std::to_string(lambda_max_4_7_1_12));
		rv.NonFatal(minmaxlam && !(lambda_max_4_7_1_12 >= lambda_min_4_7_1_11),
			"4.7.1.12: max wavelength not >= min wavelength -> " + std::to_string(lambda_max_4_7_1_12)
			+ " vs. " + std::to_string(lambda_min_4_7_1_11));
		rv.Fatal(n_spectra_4_7_1_13 < 0, "4.7.1.13: # of spectra < 0 -> " + std::to_string(n_spectra_4_7_1_13));
		rv.Fatal(n_addtl_items_4_7_1_14 < 0, "4.7.1.14: # of addtl items < 0 -> " + std::to_string(n_addtl_items_4_7_1_14));
		// 4.7.2 known data flags
		auto NotZeroOrOne = [](int32_t i) {return (i < 0) || (i > 1); };
		rv.Fatal(NotZeroOrOne(rad_flux_flag_4_7_2_3), "4.7.2.3: radiant flux flag not 0 or 1 -> "
			+ std::to_string(rad_flux_flag_4_7_2_3));
		rv.Fatal((spectrum_type_4_7_1_9 == 2 || spectrum_type_4_7_1_9 == 4) && !rad_flux_flag_4_7_2_3,
			"4.7.2.3: spectrum type is " + std::to_string(spectrum_type_4_7_1_9) + ", but radiant flux flag is not set");
		rv.Fatal(NotZeroOrOne(lambda_flag_4_7_2_4), "4.7.2.4: wavelength flag not 0 or 1 -> "
			+ std::to_string(lambda_flag_4_7_2_4));
		rv.Fatal((spectrum_type_4_7_1_9 == 2) && !lambda_flag_4_7_2_4, "4.7.2.4: spectrum type is 2, but wavelength flag not set");
		rv.Fatal(NotZeroOrOne(lum_flux_flag_4_7_2_5), "4.7.2.5: luminous flux flag not 0 or 1 -> "
			+ std::to_string(lum_flux_flag_4_7_2_5));
		rv.Fatal(!rad_flux_flag_4_7_2_3 && !lum_flux_flag_4_7_2_5, "4.7.2.5: both rad and lum flux flags are missing");
		rv.Fatal((spectrum_type_4_7_1_9 == 2 || spectrum_type_4_7_1_9 == 4) && lum_flux_flag_4_7_2_5,
			"4.7.2.5: spectrum type is " + std::to_string(spectrum_type_4_7_1_9) + ", but lum flux flag is set");
		rv.Fatal(NotZeroOrOne(stokes_flag_4_7_2_6), "4.7.2.6: stokes flag not 0 or 1 -> "
			+ std::to_string(stokes_flag_4_7_2_6));
		rv.Fatal(stokes_flag_4_7_2_6 && !rad_flux_flag_4_7_2_3, "4.7.2.3: stokes flag set but radiant flux flag not set");
		rv.Fatal(NotZeroOrOne(tristimulus_flag_4_7_2_7), "4.7.2.7: tristimulus flag not 0 or 1 -> "
			+ std::to_string(tristimulus_flag_4_7_2_7));
		rv.Fatal(tristimulus_flag_4_7_2_7 && !lum_flux_flag_4_7_2_5, "4.7.2.7: tristimulus flag set but lum flux flag not set");
		rv.Fatal(tristimulus_flag_4_7_2_7 && (spectrum_type_4_7_1_9 != 0),
			"4.7.2.7: tristimulus flag set but spectrum type " + std::to_string(spectrum_type_4_7_1_9) + " is not 0");
		rv.Fatal(NotZeroOrOne(spectrum_index_flag_4_7_2_8), "4.7.2.8: spectrum index flag not 0 or 1 -> "
			+ std::to_string(spectrum_index_flag_4_7_2_8));
		rv.Fatal(spectrum_index_flag_4_7_2_8 && (spectrum_type_4_7_1_9 != 4),
			"4.7.2.8: spectrum index flag set but spectrum type " + std::to_string(spectrum_type_4_7_1_9) + " is not 4");
		rv.Fatal(!spectrum_index_flag_4_7_2_8 && (spectrum_type_4_7_1_9 == 4),
			"4.7.2.8: spectrum index flag not set but spectrum type is 4");
		// 4.7.4 spectral tables 
		rv.Fatal(n_spectra_4_7_1_13 != spectra_4_7_4.size(),
			"4.7.4: # of spectra (4.7.1.13, " + std::to_string(n_spectra_4_7_1_13) + ") is not equal to number of spectral tables "
			+ std::to_string(spectra_4_7_4.size()));
		size_t i_sp = 0;
		for (const auto& sp : spectra_4_7_4)
			{
			++i_sp;
			rv.Fatal(sp.idx_ < 1, "4.7.4: index " + std::to_string(sp.idx_) + " must be > 0 in spectrum # " + std::to_string(i_sp));
			for (auto lam : sp.lambda_)
				rv.Fatal(!(lam > 0.0f), "4.7.4: non positive wavelength " + std::to_string(lam) + " in spectrum " + std::to_string(i_sp));
			for (auto wt : sp.weight_)
				rv.Fatal(!(wt >= 0.0f), "4.7.4: negative weight " + std::to_string(wt) + " in spectrum " + std::to_string(i_sp));
			}
		// 4.7.5 column labels
		rv.Fatal(n_addtl_items_4_7_1_14 != column_names_4_7_5.size(), "4.7.5: # of additional names: "
			+ std::to_string(column_names_4_7_5.size()) + " does not match 4.7.1.14: " + std::to_string(n_addtl_items_4_7_1_14));
		return rv;
		}

	// ********************************************************************
	// TDefaultRayArray implementation
	TDefaultRayArray::TDefaultRayArray()
		: nRays_(0), nItems_(0) {};

	TDefaultRayArray::TDefaultRayArray(size_t nRays, size_t nItems)
		: nRays_(nRays), nItems_(nItems),
		data_(nRays * nItems, std::numeric_limits<float>::signaling_NaN()) {};

	TDefaultRayArray::TDefaultRayArray(size_t nRays, size_t nItems, std::vector<float>&& data)
		: nRays_(nRays), nItems_(nItems)
		{
		size_t nflts = nRays * nItems;
		if (data.size() != nflts)
			throw std::runtime_error("TDefaultRayArray::TDefaultRayArray : data has wrong size");
		data_ = std::move(data);
		}

	TDefaultRayArray::TDefaultRayArray(size_t nRays, size_t nItems, const std::vector<float>& data)
		: nRays_(nRays), nItems_(nItems)
		{
		size_t nflts = nRays * nItems;
		if (data.size() != nflts)
			throw std::runtime_error("TDefaultRayArray::TDefaultRayArray : data has wrong size");
		data_.resize(nRays * nItems);
		std::copy(data.begin(), data.end(), data_.begin());
		}


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

	void TDefaultRayArray::SetRay(size_t i, const float* begin)
		{
		if (i >= nRays_)
			{
			std::stringstream s;
			s << "TDefaultRayArray::SetRay: i (" << i << ") >= NRays() (" << nRays_ << ")";
			throw TM25Error(s.str());
			}
		std::copy(begin, begin + NItems(), data_.begin() + i * NItems());
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

	float * TDefaultRayArray::GetRayDirect(size_t i)
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

	std::pair<TVec3f, TVec3f> TDefaultRayArray::BoundingBox() const
		{
		constexpr float minflt = std::numeric_limits<float>::lowest();
		constexpr float maxflt = std::numeric_limits<float>::max();

		std::pair<TVec3f, TVec3f> rv({ maxflt, maxflt, maxflt }, { minflt, minflt, minflt });

		auto check = [](float& lo, float& hi, float val)
			{
			if (lo > val) lo = val;
			if (hi < val) hi = val;
			};

		for (size_t i = 0; i < NRays(); ++i)
			{
			const float* iray = GetRayDirect(i);
			check(rv.first[0], rv.second[0], *iray);
			++iray;
			check(rv.first[1], rv.second[1], *iray);
			++iray;
			check(rv.first[2], rv.second[2], *iray);
			}
		return rv;
		}

	// uniform size_t random numbers
	class TUniformIntGen
		{
		public:
			TUniformIntGen() :gen(rd()), urd(0.0,1.0) {};

			// return random size_t rv in i0 <=rv < i1
			// precondition: i0 < i1
			size_t operator()(size_t i0, size_t i1)
				{
				double r = urd(gen);
				size_t rv = static_cast<size_t>(floor((1 - r)*i0 + r * i1));
				if (rv < i0)
					rv = i0;
				if (rv >= i1)
					rv = i1 - 1;
				return rv;
				};
		private:
			std::random_device rd;  //Will be used to obtain a seed for the random number engine
			std::mt19937 gen; //Standard mersenne_twister_engine seeded with rd()
			std::uniform_real_distribution<double> urd;
		};

	void TDefaultRayArray::Shuffle(size_t ibegin, size_t iend) // Fisher-Yates shuffle
		{
		if (iend <= (ibegin + 1))
			return;
		if (iend > nRays_)
			throw std::runtime_error("TDefaultRayArray::Shuffle: iend out of range");
		TUniformIntGen ig;
		for (size_t i = iend - 1; i > 0; --i)
			{
			size_t j = ig(0, i + 1);
			if (i != j)
				{
				float* ri = GetRayDirect(i);
				float* rj = GetRayDirect(j);
				std::swap_ranges(ri, ri + nItems_, rj);
				}
			}
		}
	void TDefaultRayArray::SetDataDirect(const std::vector<float>& rhs, size_t confirm_nRays, size_t confirm_nItems) // rhs.size() must be nRays * nItems. 
		{
		if (confirm_nRays != nRays_)
			throw std::runtime_error("TDefaultRayArray::SetDataDirect: nRays don't match");
		if (confirm_nItems != nItems_)
			throw std::runtime_error("TDefaultRayArray::SetDataDirect: nItems don't match");
		if (rhs.size() != nRays_ * nItems_)
			throw std::runtime_error("TDefaultRayArray::SetDataDirect: rhs has wrong size");
		data_ = rhs;
		}


	} // namespace TM25