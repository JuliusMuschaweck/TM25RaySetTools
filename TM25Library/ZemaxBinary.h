/* Support reading and writing Zemax binary ray files.
* Unfortunately, Zemax does not support TM25 ray files.
* The file format is documented in the Zemax help file, search for "Binary Source File format".
* In the Zemax OpticStudio 19.4 documentation, it is on page 657.
* The format itself is fairly simple: A 208 byte header, followed by rays, where each ray consists of
* seven floats (for flux only format) or eight floats (for spectral format).
*/

#ifndef __ZEMAXBINARY_H
#define __ZEMAXBINARY_H

#include <cstdint>
#include <vector>
#include <string>

namespace TM25
	{
	// The header can be created as a struct, following the documentation.
	struct TZemaxHeader
		{
		TZemaxHeader(); // default empty header, zero rays, flux only, 0.55 microns
		// definition from Zemax help file.
		// quote: "OpticStudio only uses the NbrRays, DimensionUnits, ray_format_type, and flux_type parameters"
		int32_t Identifier; // Format version ID, current value is 1010
		uint32_t NbrRays; // The number of rays in the file
		// must be <= 4,000,000,000
		char Description[100] = "Zemax binary header"; // A text description of the source
		float SourceFlux; // The total flux in watts of this source
		float RaySetFlux; // The flux in watts represented by this Ray Set
		float Wavelength; // The wavelength in micrometers, 0 if a composite
		float InclinationBeg, InclinationEnd; // Angular range for ray set (Degrees)
		float AzimuthBeg, AzimuthEnd; // Angular range for ray set (Degrees)
		uint32_t DimensionUnits; // METERS=0, IN=1, CM=2, FEET=3, MM=4
		float LocX, LocY, LocZ; // Coordinate Translation of the source
		float RotX, RotY, RotZ; // Source rotation (Radians) 
		float ScaleX, ScaleY, ScaleZ; // Currently unused
		float unused1, unused2, unused3, unused4;
		int32_t ray_format_type; // 0 for flux only, 2 for spectral
		int32_t flux_type; // if ray_format_type==0, then 0 for watts, 1 for lumens else 0 
		int32_t reserved1, reserved2;
		};

	// The actual ray set class
	class TZemaxRaySet
		{
		public:
			TZemaxRaySet(); // zero rays, flux only, 0.55 microns

			explicit TZemaxRaySet(const std::string& filename); // read from ray data file

			TZemaxRaySet(const TZemaxHeader& zh, std::vector<float>&& raydata); // ray set for this header
			void SetDescription(const std::string& s);	// set the description field of the header
			std::string Description() const;			// get the description field of the header

			enum class TFormatType { flux_only, spectral }; // safer and more explicit than 0 and 2
			// With TFormatType::flux_only, there will be x, y, z, kx, ky, kz, flux per ray
			// With TFormatType::spectral, there will be x, y, z, kx, ky, kz, flux, lambda per ray
			//   i.e. spectral means individual wavelength per ray.
			void SetFormatType(TFormatType ft);
			// set the format type. If set to photometrix, flux type is set to flux only, since
			// photometric plus spectral is not allowed
			TFormatType FormatType() const; // get the format type

			enum class TFluxType { radiometric, photometric };  // ray flux means radiant or luminous flux, respectively
			void SetFluxType(TFluxType ft); // set the flux type. If photometric, format type is set to flux only, too
			TFluxType FluxType() const; // get the flux type

			void SetWavelength(float lam_microns); // MICRONS!!!
			float Wavelength() const; // get the wavelength. 
			// If format type is spectral, this value is ignored.

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

			TZemaxHeader Header() const;

			const std::vector<float>& Data() const;

			void Read(const std::string& filename);
			void Write(const std::string& filename) const;

			std::pair<bool, std::string>  HeaderSanityCheck() const;

		private:
			TZemaxHeader header_;
			// we use enums for format type and flux type for safe and explicit interface
			TFormatType format_type_;
			TFluxType flux_type_;
			float wavelength_;

			std::vector<float> data_;

		};
	} // namespace
#endif // header include

