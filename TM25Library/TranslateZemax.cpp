#include "TranslateZemax.h"
#include "Timer.h"
namespace TM25
	{
	TTM25RaySet ZemaxBinaryToTM25(const TZemaxRaySet& zemaxRaySet)
		{
		TTM25Header header;
		TZemaxHeader zh = zemaxRaySet.Header();
		// fields: (untouched are ok from default constructor)
		constexpr float nan = std::numeric_limits<float>::signaling_NaN();
		header.n_rays_4_7_1_6 = zh.NbrRays;
		header.file_date_time_str_4_7_1_7 = currentISO8601TimeUTC();
		if (zemaxRaySet.FormatType() == TZemaxRaySet::TFormatType::spectral)
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
			header.lambda_4_7_1_10 = zh.Wavelength * 1000; // Zemax: microns, TM25: nm
			header.lambda_min_4_7_1_11 = zemaxRaySet.MinWavelength();
			header.lambda_max_4_7_1_12 = zemaxRaySet.MaxWavelength();
			header.lambda_flag_4_7_2_4 = false;
			}
		if (zemaxRaySet.FluxType() == TZemaxRaySet::TFluxType::radiometric)
			{
			header.phi_v_4_7_1_4 = nan;
			header.phi_4_7_1_5 = zh.SourceFlux;
			header.rad_flux_flag_4_7_2_3 = true;
			header.lum_flux_flag_4_7_2_5 = false;
			}
		else
			{
			header.phi_v_4_7_1_4 = zh.SourceFlux;
			header.phi_4_7_1_5 = nan;
			header.rad_flux_flag_4_7_2_3 = false;
			header.lum_flux_flag_4_7_2_5 = true;
			}
		header.name_4_7_3_1 = ToU32String(zh.Description);
		header.rayfile_creator_4_7_3_4 = ToU32String("TM25 ray set tools, J. Muschaweck 2023");
		TTM25Header::TSanityCheck chk = header.SanityCheck();
		if (chk.nonfatalErrors)
			header.additional_text_4_7_6.append(ToU32String(chk.msg));
		if (chk.fatalErrors)
			throw std::runtime_error("ZemaxBinaryToTM25: fatal error in TM25 header sanity check: " + chk.msg);
		
		size_t nItems = 0;
		if (header.spectrum_type_4_7_1_9 == 1) // single wavelength
			nItems = 7;
		else
			nItems = 8;
		TDefaultRayArray rays(header.n_rays_4_7_1_6, nItems);
		auto zdatastart = zemaxRaySet.Data().begin();
		for (size_t i = 0; i < header.n_rays_4_7_1_6; ++i)
			{
			rays.SetRay(i, &(*zdatastart));
			zdatastart += nItems;
			}
		rays.ChangeMicronsToNanometers();
		TTM25RaySet rv(header, rays);
		return rv;
		}
	
	TZemaxRaySet TM25ToZemaxBinary(const TTM25RaySet& rs)
		{
		TZemaxHeader zh;
		TM25::TTM25Header tmh = rs.Header();
		zh.NbrRays = static_cast<uint32_t>(tmh.n_rays_4_7_1_6); // The number of rays in the file
		// set description from name
		std::string tmp = ToString(tmh.name_4_7_3_1);
		if (tmp.length() > 99)
			tmp.resize(99);
		std::fill(zh.Description, zh.Description + 100, char(0));
		std::copy(tmp.begin(), tmp.end(), zh.Description);
		// luminous or radiant flux?
		if (tmh.rad_flux_flag_4_7_2_3)
			{
			zh.flux_type = 0; // if ray_format_type==0, then 0 for watts, 1 for lumens else 0 
			zh.SourceFlux = tmh.phi_4_7_1_5; // The total flux in watts of this source
			if (tmh.spectrum_type_4_7_1_9 == 2) // wavelength defined per ray
				{
				zh.ray_format_type = 2;
				zh.Wavelength = (tmh.lambda_min_4_7_1_11 + tmh.lambda_max_4_7_1_12) / 2;
				}
			else
				{
				zh.ray_format_type = 0;
				zh.Wavelength = tmh.lambda_4_7_1_10; // The wavelength in micrometers, 0 if a composite
				}
			}
		else if (tmh.lum_flux_flag_4_7_2_5)
			{
			zh.flux_type = 1; // if ray_format_type==0, then 0 for watts, 1 for lumens else 0 
			zh.SourceFlux = tmh.phi_v_4_7_1_4;
			zh.ray_format_type = 0; // Zemax cannot have individual wavelengths when flux is luminous
			if (tmh.spectrum_type_4_7_1_9 == 1) // single wavelength
				zh.Wavelength = tmh.lambda_4_7_1_10; // The wavelength in micrometers, 0 if a composite
			else
				zh.Wavelength = (tmh.lambda_min_4_7_1_11 + tmh.lambda_max_4_7_1_12) / 2;
			}
		else
			throw std::runtime_error("TM25ToZemaxBinary: one of rad_flux_flag_4_7_2_3 or lum_flux_flag_4_7_2_5 must be true");
		zh.RaySetFlux = zh.SourceFlux; // The flux in watts represented by this Ray Set
		zh.DimensionUnits = 4; // METERS=0, IN=1, CM=2, FEET=3, MM=4

		TRaySetItems items = rs.Items();
		TRaySetItems needed;
		needed.MarkAsPresent(RayItem::x);
		needed.MarkAsPresent(RayItem::y);
		needed.MarkAsPresent(RayItem::z);
		needed.MarkAsPresent(RayItem::kx);
		needed.MarkAsPresent(RayItem::ky);
		needed.MarkAsPresent(RayItem::kz);
		needed.MarkAsPresent(RayItem::phi);
		if (zh.ray_format_type == 2) // wavelength per ray
			{
			needed.MarkAsPresent(RayItem::lambda);
			}
		if (items.ContainsItems(needed) == false)
			{
			throw std::runtime_error("TM25ToZemaxBinary: Not all required ray items present");
			}
		TRaySetItems::TExtractionMap extractionMap = items.ExtractionMap(needed);
		TTM25RaySet::TRayArray rays = rs.ExtractAll(extractionMap);
		
		size_t n_rays = rays.NRays();
		size_t n_items = rays.NItems();
		float lam = zh.Wavelength;
		std::vector<float> raydata;
		for (size_t i = 0; i < n_rays; ++i)
			{
			const float* r = rays.GetRayDirect(i);
			if (zh.ray_format_type == 2) // wavelength per ray
				lam = *(r + 7);
			std::array<float, 8> oneray({ *r, *(r + 1), *(r + 2), *(r + 3), *(r + 4), *(r + 5), *(r + 6), lam * 0.001f });
			std::copy(oneray.begin(), oneray.end(), std::back_inserter(raydata));
			}
		TZemaxRaySet rv(zh, std::move(raydata));
		return rv;
		}
	} // namespace