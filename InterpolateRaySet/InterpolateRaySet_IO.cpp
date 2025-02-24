#include "InterpolateRaySet_IO.h"
#include <ASCIIRayFile.h>
#include <TranslateZemax.h>

using std::cout;
using std::endl;


TM25::TTM25RaySet ReadRaySet(const TInterpolateRaySetCfg& cfg, std::ostream& info)//const std::string& fn, std::ostream& info)
	{
	const TSection& rsc = cfg.Section("RaySetControl");
	TM25::TTM25RaySet rs;
	if (rsc.IsEmpty("inputRayFileName"))
		throw std::runtime_error("inputRayFileName must be supplied");
	const std::string format = rsc.String("inputRayFileFormat");
	const std::string fn = rsc.String("inputRayFileName");
	if (format.compare("TM25") == 0)
		{
		info << "Reading TM25 ray file " << fn << endl;
		rs.Read(fn);
		// rs.Make_k_unit();
		if (rs.Warnings().empty())
			info << "No warnings" << endl;
		for (auto w : rs.Warnings())
			{
			info << "Warning: " << w << "\n";
			}
		}
	else if (format.compare("ZemaxBinary") == 0)
		{
		info << "Reading Zemax binary ray file " << fn << endl;
		TM25::TZemaxRaySet zemaxRaySet(fn);
		rs = TM25::ZemaxBinaryToTM25(zemaxRaySet);
		}
	else if (format.compare("ZemaxText") == 0)
		{
		info << "Reading Zemax text ray file " << fn << endl;
		TM25::TASCIIRaySetOptions opts;
		opts.minHeaderLines = 1;
		opts.wavelengthUnit_ = TM25::TASCIIRaySetOptions::WLU::micrometer;
		rs = TM25::ReadGenericASCIIRaySet(fn, opts);
		}
	else if (format.compare("LightToolsBinary") == 0)
		{
		info << "Reading LightTools binary ray file " << fn << endl;
		TM25::TZemaxRaySet zemaxRaySet(fn);
		rs = TM25::ZemaxBinaryToTM25(zemaxRaySet);
		}

	else
		throw std::runtime_error("unknown input ray file format: " + format);
	return rs;
	}

void WriteRaySet(TM25::TTM25RaySet& rs, const TInterpolateRaySetCfg& cfg, std::ostream& info)
	{
	const TSection& rsc = cfg.Section("RaySetControl");
	const std::string format = rsc.String("outputRayFileFormat");
	const std::string fn = rsc.String("outputRayFileName");
	if (format.compare("TM25") == 0)
		{
		info << "Writing TM25 ray file " << fn << endl;
		rs.Write(fn);
		if (rs.Warnings().empty())
			info << "No warnings" << endl;
		for (auto w : rs.Warnings())
			{
			info << "Warning: " << w << "\n";
			}
		}
	else if (format.compare("ZemaxBinary") == 0)
		{
		info << "Writing  Zemax binary ray file " << fn << endl;
		TM25::TZemaxRaySet zemaxRaySet = TM25ToZemaxBinary(rs);
		zemaxRaySet.Write(fn);
		}
	}

