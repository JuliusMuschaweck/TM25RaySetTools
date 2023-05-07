#include "ZemaxBinary.h"
#include <algorithm>
#include <ReadFile.h>
#include <WriteFile.h>

namespace TM25
	{

	TZemaxHeader::TZemaxHeader()
		{
		Identifier = 1010; // Format version ID, current value is 1010
		NbrRays = 0; // The number of rays in the file
		// must be <= 4,000,000,000
		char* c = Description;
		std::fill(c, c + sizeof(Description), char(0)); // A text description of the source
		std::string tmp = "Default Zemax Binary";
		std::copy(tmp.begin(), tmp.end(), Description);
		SourceFlux = 1; // The total flux in watts of this source
		RaySetFlux = 1; // The flux in watts represented by this Ray Set
		Wavelength = float(0.55); // The wavelength in micrometers, 0 if a composite
		InclinationBeg = 0;
		InclinationEnd = 0; // Angular range for ray set (Degrees)
		AzimuthBeg = 0;
		AzimuthEnd = 0; // Angular range for ray set (Degrees)
		DimensionUnits = 4; // METERS=0, IN=1, CM=2, FEET=3, MM=4
		LocX = LocY = LocZ = 0; // Coordinate Translation of the source
		RotX = RotY = RotZ = 0; // Source rotation (Radians) 
		ScaleX = ScaleY = ScaleZ = 1;// Currently unused
		unused1 = unused2 = unused3 = unused4 = 0;
		ray_format_type = 0; // 0 for flux only, 2 for spectral
		flux_type = 0; // if ray_format_type==0, then 0 for watts, 1 for lumens else 0 
		reserved1 = 0;
		reserved2 = 0;
		}


	TZemaxRaySet::TZemaxRaySet() // zero rays, flux only, 0.55 microns
		: // header_ default ctor is just fine
		format_type_{ TFormatType::flux_only },
		flux_type_{ TFluxType::radiometric },
		wavelength_{ float(0.55) }
		// data_ default ctor is just fine
		{
		};

	TZemaxRaySet::TZemaxRaySet(const std::string& filename)
		{
		this->Read(filename);
		};

	void TZemaxRaySet::SetDescription(const std::string& s)
		{
		char* c = header_.Description;
		std::fill(c, c + sizeof(header_.Description), char(0)); // A text description of the source
		if (s.length() < sizeof(header_.Description))
			std::copy(s.begin(), s.end(), header_.Description);
		else
			std::copy(s.begin(), s.begin() + sizeof(header_.Description) - 1, header_.Description);
		};

	std::string TZemaxRaySet::Description() const
		{
		return std::string(header_.Description);
		};

	void TZemaxRaySet::SetFormatType(TFormatType ft)
		{
		format_type_ = ft;
		if (ft == TFormatType::spectral)
			flux_type_ = TFluxType::radiometric;
		};

	TZemaxRaySet::TFormatType TZemaxRaySet::FormatType() const
		{
		return format_type_;
		};

	void TZemaxRaySet::SetFluxType(TFluxType ft)
		{
		flux_type_ = ft;
		if (ft == TFluxType::photometric)
			format_type_ = TFormatType::flux_only;
		};

	TZemaxRaySet::TFluxType TZemaxRaySet::FluxType() const
		{
		return flux_type_;
		};

	void TZemaxRaySet::SetWavelength(float lam_microns)
		{
		wavelength_ = lam_microns;
		header_.Wavelength = lam_microns;
		}

	float TZemaxRaySet::Wavelength() const
		{
		return wavelength_;
		}

	float TZemaxRaySet::MinWavelength() const
		{
		float rv = std::numeric_limits<float>::max();
		for (size_t i = 0; i < NRays(); ++i)
			{
			float tmp = data_[i * 8 + 7];
			if (tmp < rv)
				rv = tmp;
			}
		return rv;
		}

	float TZemaxRaySet::MaxWavelength() const
		{
		float rv = std::numeric_limits<float>::min();
		for (size_t i = 0; i < NRays(); ++i)
			{
			float tmp = data_[i * 8 + 7];
			if (tmp > rv)
				rv = tmp;
			}
		return rv;
		}


	TZemaxHeader TZemaxRaySet::Header() const
		{
		return header_;
		};

	void TZemaxRaySet::AddRay(float x, float y, float z, float kx, float ky, float kz, float flux, float lam)
		{
		data_.push_back(x);
		data_.push_back(y);
		data_.push_back(z);
		data_.push_back(kx);
		data_.push_back(ky);
		data_.push_back(kz);
		data_.push_back(flux);
		data_.push_back(lam);
		}

	void TZemaxRaySet::AddRay(const TRay_lam& ray)
		{
		AddRay(ray.x, ray.y, ray.z, ray.kx, ray.ky, ray.kz, ray.flux, ray.lam);
		}
	void TZemaxRaySet::AddRay(const TRay_fluxonly& ray)
		{
		AddRay(ray.x, ray.y, ray.z, ray.kx, ray.ky, ray.kz, ray.flux, header_.Wavelength);
		}

	std::size_t TZemaxRaySet::NRays() const
		{
		std::size_t rv = data_.size() / 8;
		return rv;
		}

	const std::vector<float>& TZemaxRaySet::Data() const
		{
		return data_;
		}

	void TZemaxRaySet::Read(const std::string& filename)
		{
		try {
			TZemaxHeader h;
			TM25::TReadFile f(filename);
			h = f.Read<TZemaxHeader>();
			if (h.Identifier != 1010)
				throw std::runtime_error("TZemaxRaySet::Read: Wrong format identifier in header of file " + filename);
			TFormatType ft;
			if (h.ray_format_type == 0)
				ft = TFormatType::flux_only;
			else if (h.ray_format_type == 2)
				ft = TFormatType::spectral;
			else
				throw std::runtime_error("TZemaxRaySet::Read: Unknown ray format type in header of file " + filename);
			header_ = h;
			size_t nrays = header_.NbrRays;
			for (size_t i = 0; i < nrays; ++i)
				if (ft == TFormatType::flux_only)
					AddRay(f.Read<TRay_fluxonly>());
				else
					AddRay(f.Read<TRay_lam>());
			if (!f.AtEof())
				throw std::runtime_error("TZemaxRaySet::Read: expected EOF after reading");
			}
		catch (std::runtime_error& err)
			{
			throw std::runtime_error("TZemaxRaySet::Read: Error reading file " + filename + ": " + err.what());
			}
		// test header consistency
		auto test = HeaderSanityCheck();
		if (!test.first)
			throw std::runtime_error("TZemaxRaySet::Read file " + filename + ": inconsistent header: " + test.second);
		// change back to mm if needed
		if (header_.DimensionUnits != 4)
			{
			float fac;
			switch (header_.DimensionUnits)
				{
					case 0: fac = 1000; break;				// meter
					case 1: fac = float(25.4); break;		// inch
					case 2: fac = 10; break;				// cm
					case 3: fac = float(25.4 * 12); break;	// feet
					default:
						throw std::runtime_error("TZemaxRaySet::Read: Unknown dimension unit flag");
				}
			size_t nrays = header_.NbrRays;
			for (size_t i = 0; i < nrays; ++i)
				{
				data_[i * 8] *= fac;		// x
				data_[i * 8 + 1] *= fac; // y
				data_[i * 8 + 2] *= fac; // z
				}
			header_.DimensionUnits = 4;
			}
		// set other internal values
		if (header_.ray_format_type == 0)
			format_type_ = TFormatType::flux_only;
		else
			format_type_ = TFormatType::spectral;
		if (header_.flux_type == 0)
			flux_type_ = TFluxType::radiometric;
		else
			flux_type_ = TFluxType::photometric;
		wavelength_ = header_.Wavelength;
		};

	void TZemaxRaySet::Write(const std::string& filename) const
		{
		try
			{
			TM25::TWriteFile f(filename);
			if (sizeof(TZemaxHeader) != 208)
				throw std::logic_error("TZemaxRaySet::Write: sizeof(TZemaxheader) != 208");
			f.Write<TZemaxHeader>(header_);
			size_t nrays = NRays();
			for (size_t i = 0; i < nrays; ++i)
				{
				const float* start = &(data_[i * 8]);
				if (FormatType() == TFormatType::flux_only)
					f.WriteRange(start, start + 7);
				else
					f.WriteRange(start, start + 8);
				}
			}
		catch (std::runtime_error& err)
			{
			throw std::runtime_error("TZemaxRaySet::Write: Error writing file " + filename + ": " + err.what());
			}

		};


	std::pair<bool, std::string> TZemaxRaySet::HeaderSanityCheck() const
		{
		std::pair<bool, std::string> rv = std::make_pair(true, std::string());
		auto check = [&rv](bool condition, std::string msg)
			{
			if (!condition)
				{
				rv.first = false;
				rv.second = rv.second + ", " + msg;
				}
			};
		check(header_.Identifier == 1010, "wrong format version ID");
		check(header_.NbrRays == NRays(), "wrong number of rays");
		check((header_.ray_format_type == 0) || (header_.ray_format_type == 2), "unknown ray format type");
		check((header_.flux_type == 0) || (header_.flux_type == 1), "unknown flux type");
		check(not((header_.ray_format_type == 2)
			&& (header_.flux_type == 1)), "flux type cannot be photometric with spectral format type");
		return rv;
		};


	} // namespace
