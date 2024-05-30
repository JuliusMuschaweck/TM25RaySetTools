#include "TranslateLightTools.h"
#include "Timer.h"

namespace TM25
	{
	TTM25RaySet LightToolsBinaryToTM25(const TLightToolsRaySet& ltRaySet)
		{
		TTM25Header header;
		TLightToolsHeader lth = ltRaySet.Header();
		// fields: (untouched are ok from default constructor)
		constexpr float nan = std::numeric_limits<float>::signaling_NaN();
		header.n_rays_4_7_1_6 = ltRaySet.NRays();
		header.file_date_time_str_4_7_1_7 = currentISO8601TimeUTC();
		if (ltRaySet.FormatType() == TLightToolsRaySet::TFormatType::spectral)
			{
			header.spectrum_type_4_7_1_9 = 2; // wavelength per ray
			header.lambda_4_7_1_10 = nan;
			header.lambda_min_4_7_1_11 = nan;
			header.lambda_max_4_7_1_12 = nan;
			header.lambda_flag_4_7_2_4 = true;
			}
		else
			{
			header.spectrum_type_4_7_1_9 = 1; // single wavelength for all rays
			header.lambda_4_7_1_10 = ltRaySet.Wavelength(); // Zemax: microns, TM25: nm
			header.lambda_min_4_7_1_11 = ltRaySet.Wavelength();
			header.lambda_max_4_7_1_12 = ltRaySet.Wavelength();
			header.lambda_flag_4_7_2_4 = false;
			}
		if (ltRaySet.FluxType() == TLightToolsRaySet::TFluxType::radiometric)
			{
			header.phi_v_4_7_1_4 = nan;
			header.phi_4_7_1_5 = lth.Flux;
			header.rad_flux_flag_4_7_2_3 = true;
			header.lum_flux_flag_4_7_2_5 = false;
			}
		else
			{
			header.phi_v_4_7_1_4 = lth.Flux;
			header.phi_4_7_1_5 = nan;
			header.rad_flux_flag_4_7_2_3 = false;
			header.lum_flux_flag_4_7_2_5 = true;
			}
		header.name_4_7_3_1 = ToU32String("LT ray file");
		header.rayfile_creator_4_7_3_4 = ToU32String("TM25 ray set tools, J. Muschaweck 2023");
		TTM25Header::TSanityCheck chk = header.SanityCheck();
		if (chk.nonfatalErrors)
			header.additional_text_4_7_6.append(ToU32String(chk.msg));
		if (chk.fatalErrors)
			throw std::runtime_error("LightToolsBinaryToTM25: fatal error in TM25 header sanity check: " + chk.msg);
		
		size_t nItems = 0;
		if (header.spectrum_type_4_7_1_9 == 1) // single wavelength
			nItems = 7;
		else
			nItems = 8;
		TDefaultRayArray rays(header.n_rays_4_7_1_6, nItems);
		auto zdatastart = ltRaySet.Data().begin();
		for (size_t i = 0; i < header.n_rays_4_7_1_6; ++i)
			{
			rays.SetRay(i, &(*zdatastart));
			zdatastart += nItems;
			}
		TTM25RaySet rv(header, rays);
		return rv;
		}
	
	TLightToolsRaySet TM25ToLightToolsBinary(const TTM25RaySet& rs)
		{
		TLightToolsHeader lth;
		TM25::TTM25Header tmh = rs.Header();
		size_t nRays = static_cast<uint32_t>(tmh.n_rays_4_7_1_6); // The number of rays in the file
		// luminous or radiant flux?
		if (tmh.rad_flux_flag_4_7_2_3)
			{
			lth.DataType = 0; // if ray_format_type==0, then 0 for watts, 1 for lumens else 0 
			lth.Flux = tmh.phi_4_7_1_5; // The total flux in watts of this source
			}
		else if (tmh.lum_flux_flag_4_7_2_5)
			{
			lth.DataType = 1; // if ray_format_type==0, then 0 for watts, 1 for lumens else 0 
			lth.Flux = tmh.phi_v_4_7_1_4;
			}
		else
			throw std::runtime_error("TM25ToZemaxBinary: one of rad_flux_flag_4_7_2_3 or lum_flux_flag_4_7_2_5 must be true");
		if (tmh.spectrum_type_4_7_1_9 == 2) // wavelength defined per ray
			{
			lth.ColorInfo = 2;
			}
		else
			{
			lth.ColorInfo = 0;
			}

		TLightToolsRaySet rv;


		TRaySetItems items = rs.Items();
		TRaySetItems needed;
		needed.MarkAsPresent(RayItem::x);
		needed.MarkAsPresent(RayItem::y);
		needed.MarkAsPresent(RayItem::z);
		needed.MarkAsPresent(RayItem::kx);
		needed.MarkAsPresent(RayItem::ky);
		needed.MarkAsPresent(RayItem::kz);
		needed.MarkAsPresent(RayItem::phi);
		bool withWavelength = lth.ColorInfo == 2;
		if (withWavelength) // wavelength per ray
			{
			needed.MarkAsPresent(RayItem::lambda);
			}
		if (items.ContainsItems(needed) == false)
			{
			throw std::runtime_error("TM25ToLightToolsBinary: Not all required ray items present");
			}
		TRaySetItems::TExtractionMap extractionMap = items.ExtractionMap(needed);
		TTM25RaySet::TRayArray rays = rs.ExtractAll(extractionMap);
		rv.SetRaysDirect(rays.Data(), rays.NRays(), withWavelength);
		rv.SetHeader(lth);
		return rv;
		}
	} // namespace