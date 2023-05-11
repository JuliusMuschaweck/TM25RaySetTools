/* Support reading and writing LightTools binary ray files.
* Documentation at 
* file:///C:/Program%20Files/Optical%20Research%20Associates/LightTools%202023.03/Help/index.html#page/illumination/liug_3.4.127.html##ww257161
header:

* 
*/

#ifndef __LIGHTTOOLSBINARY_H
#define __LIGHTTOOLSBINARY_H

#include <cstdint>
#include <vector>
#include <array>
#include <string>

namespace TM25
	{
	// The header can be created as a struct, following the documentation.
	struct TLightToolsHeader
		{
		TLightToolsHeader(); // default empty header, zero rays, flux only, 0.55 microns
		// definition from Zemax help file.
		// quote: "OpticStudio only uses the NbrRays, DimensionUnits, ray_format_type, and flux_type parameters"
		std::array<char, 4> Signature;	// "LTRF"
		uint32_t	MajorVersion;		// default 2, polarization 3
		uint32_t	MinorVersion;		// 0
		uint32_t	DataType;			// 0 radiometric, 1 photometric
		uint32_t	FarField;			// 0 with xyz location, 1 without
		uint32_t	ColorInfo;			// 0 - no color info, 
										// 1 - CIE tristimulus values written out -- we ignroe this
										// 2 - Wavelength info is written out
		uint32_t	LengthUnits;		// 0 - system, 1 - mm, 2 cm, 3 m, 4 inch, 5 feet, 6 microns, 7 nanometer
		float		Origin_X;
		float		Origin_Y;
		float		Origin_Z;
		float		Flux;				// total flux in watts or lumens according to DataType
		uint32_t	Polarization;		// 0 none, 1 with Stokes and Direction vectors
		};

	// The actual ray set class. We ignore polarization for now
	class TLightToolsRaySet
		{
		public:
			TLightToolsRaySet(); // zero rays, flux only, 0.55 microns, mm, at 0,0,0, flux 0

			explicit TLightToolsRaySet(const std::string& filename); // read from ray data file

			enum class TFormatType { flux_only, 
				// tristimulus, 
				spectral }; // safer and more explicit than 0 and 2
			// With TFormatType::flux_only, there will be x, y, z, kx, ky, kz, flux per ray
			// With TFormatType::spectral, there will be x, y, z, kx, ky, kz, flux, lambda per ray
			//   i.e. spectral means individual wavelength per ray.
			void SetFormatType(TFormatType ft);
			TFormatType FormatType() const; // get the format type

			enum class TFluxType { radiometric, photometric };  // ray flux means radiant or luminous flux, respectively
			void SetFluxType(TFluxType ft); // set the flux type. If photometric, format type is set to flux only, too
			TFluxType FluxType() const; // get the flux type

			void SetWavelength(float lam_nanometer); // default is 550 nm
			float Wavelength() const; // get the wavelength. 
			// If format type is spectral or tristimulus, this value is ignored.

			float MinWavelength() const;
			float MaxWavelength() const;

			struct TRay_lam  // single ray record for spectral
				{
				float x, y, z,
					kx, ky, kz,
					flux,
					lam;
				};

			struct TRay_fluxonly // single ray record for flux only
				{
				float x, y, z,
					kx, ky, kz,
					flux;
				};

			// rays are always stored internally with a separate wavelength field
			void AddRay(float x, float y, float z, float kx, float ky, float kz, float flux, float lam);
			void AddRay(const TRay_lam& ray);
			void AddRay(const TRay_fluxonly& ray); // 

			std::size_t NRays() const;

			TLightToolsHeader Header() const;

			const std::vector<float>& Data() const;

			void Read(const std::string& filename);
			void Write(const std::string& filename) const;

			std::pair<bool, std::string>  HeaderSanityCheck() const;

		private:
			TLightToolsHeader header_;
			// we use enums for format type and flux type for safe and explicit interface
			TFormatType format_type_;
			TFluxType flux_type_;
			float wavelength_;

			std::vector<float> data_;
		};
	} // namespace
#endif // header include


