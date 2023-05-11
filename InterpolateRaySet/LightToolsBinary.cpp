#include "LightToolsBinary.h"
#include <algorithm>
#include <ReadFile.h>
#include <WriteFile.h>

namespace TM25 {

	TLightToolsHeader::TLightToolsHeader() // default empty header, zero rays, flux only
		// definition from Zemax help file.
		// quote: "OpticStudio only uses the NbrRays, DimensionUnits, ray_format_type, and flux_type parameters"
		: Signature { {'L', 'T', 'R', 'F'}}, 	// "LTRF"
		MajorVersion{ 2 },		// default 2, polarization 3
		MinorVersion{ 0 },		// 0
		DataType{ 0 },			// 0 radiometric, 1 photometric
		FarField{ 0 },			// 0 with xyz location, 1 without
		ColorInfo{ 0 },			// 0 - no color info, 1 - CIE tristimulus values written out 
								// 2 - Wavelength info is written out
		LengthUnits{ 1 },		// 0 - system, 1 - mm, 2 cm, 3 m, 4 inch, 5 feet, 6 microns, 7 nanometer
		Origin_X{ 0.0f },
		Origin_Y{ 0.0f },
		Origin_Z{ 0.0f },
		Flux{ 0.0f },				// total flux in watts or lumens according to DataType
		Polarization{ 0 }		// 0 none, 1 with Stokes and Direction vectors
		{};

		TLightToolsRaySet::TLightToolsRaySet()
			:
			// header_ just fine
			format_type_{ TFormatType::flux_only },
			flux_type_{ TFluxType::radiometric },
			wavelength_{ 550.0f }
			// data_ just fine
			{}; // zero rays, flux only, 0.55 microns, mm, at 0,0,0, flux 0

			TLightToolsRaySet::TLightToolsRaySet(const std::string& filename) // read from ray data file
				{
				this->Read(filename);
				}

			void TLightToolsRaySet::SetFormatType(TFormatType ft)
				{
				format_type_ = ft;
				if (ft == TFormatType::spectral)
					flux_type_ = TFluxType::radiometric;
				};

				// set the format type. If set to photometrix, flux type is set to flux only, since
				// photometric plus spectral is not allowed
			TLightToolsRaySet::TFormatType TLightToolsRaySet::FormatType() const // get the format type
				{
				return format_type_;
				};


			void TLightToolsRaySet::SetFluxType(TFluxType ft) // set the flux type. If photometric, format type is set to flux only, too
				{
				flux_type_ = ft;
				if (ft == TFluxType::photometric)
					format_type_ = TFormatType::flux_only;
				};

			TLightToolsRaySet::TFluxType TLightToolsRaySet::FluxType() const // get the flux type
				{
				return flux_type_;
				};

			void TLightToolsRaySet::SetWavelength(float lam_nanometer) // default is 550 nm
				{
				wavelength_ = lam_nanometer;
				}

			float TLightToolsRaySet::Wavelength() const // get the wavelength. 
				// If format type is spectral or tristimulus, this value is ignored.
				{
				return wavelength_;
				}


			float TLightToolsRaySet::MinWavelength() const
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
			
			float TLightToolsRaySet::MaxWavelength() const
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

				// rays are always stored internally with a separate wavelength field
			void TLightToolsRaySet::AddRay(float x, float y, float z, float kx, float ky, float kz, float flux, float lam)
				{
				data_.push_back(x);
				data_.push_back(y);
				data_.push_back(z);
				data_.push_back(kx);
				data_.push_back(ky);
				data_.push_back(kz);
				data_.push_back(flux);
				data_.push_back(lam);
				};

			void TLightToolsRaySet::AddRay(const TRay_lam& ray)
				{
				AddRay(ray.x, ray.y, ray.z, ray.kx, ray.ky, ray.kz, ray.flux, ray.lam);
				}

			void TLightToolsRaySet::AddRay(const TRay_fluxonly& ray) // 
				{
				AddRay(ray.x, ray.y, ray.z, ray.kx, ray.ky, ray.kz, ray.flux, wavelength_);
				}
			
			std::size_t TLightToolsRaySet::NRays() const
				{
				std::size_t rv = data_.size() / 8;
				return rv;
				}

			TLightToolsHeader TLightToolsRaySet::Header() const
				{
				return header_;
				}

			const std::vector<float>& TLightToolsRaySet::Data() const
				{
				return data_;
				}

			void TLightToolsRaySet::Read(const std::string& filename)
				{
				try
					{
					TLightToolsHeader h;
					TM25::TReadFile f(filename);
					h = f.Read<TLightToolsHeader>();
					
					if (h.Signature != std::array<char, 4> {'L', 'T', 'R', 'F'})
						{
						std::string s{ "Expected signature \"LTRF\", found " };
						s = (s + h.Signature[0]) + h.Signature[1] + h.Signature[2] + h.Signature[3];
						throw std::runtime_error(s);
						};
					if (h.MajorVersion == 3)
						throw std::runtime_error("cannot handle polarization");
					if (h.FarField != 0)
						throw std::runtime_error("cannot handle ray file without xyz location data");
					
					if (h.ColorInfo == 0)
						{
						format_type_ = TFormatType::flux_only;
						wavelength_ = 550;
						}
					else if (h.ColorInfo == 1)
						throw std::runtime_error("cannot handle tristimulus ray file");
					else if (h.ColorInfo == 2)
						format_type_ = TFormatType::spectral;
					else
						throw std::runtime_error("unknown ColorInfo " + std::to_string(h.ColorInfo));
					
					if (h.DataType == 0)
						flux_type_ = TFluxType::radiometric;
					else if (h.DataType == 1)
						flux_type_ = TFluxType::photometric;
					else 
						throw std::runtime_error("unknown DataType " + std::to_string(h.DataType));
					for (;;)
						{
						char c = f.Peek<char>();
						if (c == 'L')
							{
							std::array<char, 7> tmp = f.Peek< std::array<char, 7> >();
							if (tmp == std::array<char, 7>{'L', 'T', 'R', 'F', 'E', 'N', 'D'})
								break;
							}
						if (format_type_ == TFormatType::flux_only)
							AddRay(f.Read<TRay_fluxonly>());
						else
							AddRay(f.Read<TRay_lam>());
						}
					header_ = h;
					}
				catch (std::runtime_error& err)
					{
					throw std::runtime_error("TLightToolsRaySet::Read: Error reading file " + filename + ": " + err.what());
					}
				auto test = HeaderSanityCheck();
				if (!test.first)
					throw std::runtime_error("TLightToolsRaySet::Read file " + filename + ": inconsistent header: " + test.second);
				}
			//	void TLightToolsRaySet::Write(const std::string& filename) const;

			std::pair<bool, std::string>  TLightToolsRaySet::HeaderSanityCheck() const
				{
				return std::make_pair(true, "");
				}

			//private:
			//	TLightToolsHeader header_;
			//	// we use enums for format type and flux type for safe and explicit interface
			//	TFormatType format_type_;
			//	TFluxType flux_type_;
			//	float wavelength_;

			//	std::vector<float> data_;
			//};


	} // namespace TM25
