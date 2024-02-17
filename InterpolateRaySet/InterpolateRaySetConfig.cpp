#include "InterpolateRaySetConfig.h"

TRaySetControlSection::TRaySetControlSection() : TSection("RaySetControl") {};

void TRaySetControlSection::AddAllowedValues()
	{
	values_.insert({ "inputRayFileName",	MakeDefaultValueTokenSequence<Token::string>(std::string("inputRayFileName missing")) });
	values_.insert({ "inputRayFileFormat",	MakeDefaultValueTokenSequence<Token::string>(std::string("TM25")) });
	values_.insert({ "logFileName",			MakeDefaultValueTokenSequence<Token::string>(std::string("InterpolateRaySet.log")) });
	values_.insert({ "consoleOutput",		MakeDefaultValueTokenSequence<Token::boolean>(true) });
	values_.insert({ "scrambleInput",		MakeDefaultValueTokenSequence<Token::boolean>(true) });
	values_.insert({ "selectByMaxNumber",	MakeEmptyTokenSequence() });
	values_.insert({ "doDiagnostics",		MakeDefaultValueTokenSequence<Token::boolean>(true) });
	values_.insert({ "restrictToRelVirtualFocusDistance",MakeEmptyTokenSequence() });
	values_.insert({ "restrictToKz",		MakeEmptyTokenSequence() }); // backwards compatibility
	values_.insert({ "restrictToMinKz",		MakeEmptyTokenSequence() });
	values_.insert({ "restrictToXYBox",		MakeEmptyTokenSequence() });
	values_.insert({ "restrictToFirstNRays",		MakeEmptyTokenSequence() });
	values_.insert({ "restrictToEtendueThreshold", MakeDefaultValueTokenSequence<Token::real>(0.0) }); // ignore if 0.0
	values_.insert({ "nClip",				MakeDefaultValueTokenSequence<Token::integer>(0) });
	values_.insert({ "nOutputRays",			MakeDefaultValueTokenSequence<Token::integer>(0) });
	values_.insert({ "nNeighbors",			MakeDefaultValueTokenSequence<Token::integer>(10) });
	values_.insert({ "outputRayFileName",	MakeEmptyTokenSequence() });
	values_.insert({ "outputRayFileFormat",	MakeDefaultValueTokenSequence<Token::string>(std::string("TM25")) });
	values_.insert({ "phaseSpaceType",		MakeDefaultValueTokenSequence<Token::identifier>(std::string("VirtualFocusZSphere")) });
	values_.insert({ "ZPlane_z",			MakeDefaultValueTokenSequence<Token::real>(0.0) });
	values_.insert({ "ZCylinderRadius",		MakeDefaultValueTokenSequence<Token::real>(1.0) });
	values_.insert({ "characteristicCurveFileName",	MakeEmptyTokenSequence() });
	values_.insert({ "skewness_z_FileName",	MakeEmptyTokenSequence() });
	values_.insert({ "skewness_nBins",	MakeDefaultValueTokenSequence<Token::integer>(100) });
	values_.insert({ "skewness_binType",	MakeDefaultValueTokenSequence<Token::identifier>(std::string("sameSkewness")) });
	values_.insert({ "setTotalFlux",		MakeEmptyTokenSequence() });
	values_.insert({ "LuminanceLookupTable_FileName", MakeEmptyTokenSequence() });
	values_.insert({ "LuminanceLookupTable_xmin", MakeEmptyTokenSequence() });
	values_.insert({ "LuminanceLookupTable_xmax", MakeEmptyTokenSequence() });
	values_.insert({ "LuminanceLookupTable_ymin", MakeEmptyTokenSequence() });
	values_.insert({ "LuminanceLookupTable_ymax", MakeEmptyTokenSequence() });
	values_.insert({ "LuminanceLookupTable_kxmin", MakeEmptyTokenSequence() });
	values_.insert({ "LuminanceLookupTable_kxmax", MakeEmptyTokenSequence() });
	values_.insert({ "LuminanceLookupTable_kymin", MakeEmptyTokenSequence() });
	values_.insert({ "LuminanceLookupTable_kymax", MakeEmptyTokenSequence() });
	values_.insert({ "LuminanceLookupTable_xPoints", MakeEmptyTokenSequence() });
	values_.insert({ "LuminanceLookupTable_yPoints", MakeEmptyTokenSequence() });
	values_.insert({ "LuminanceLookupTable_kxPoints", MakeEmptyTokenSequence() });
	values_.insert({ "LuminanceLookupTable_kyPoints", MakeEmptyTokenSequence() });
	values_.insert({ "someRealVector", MakeEmptyTokenSequence() });
	values_.insert({ "someIntVector", MakeEmptyTokenSequence() });
	}

void TRaySetControlSection::AddAllowedKeywords()
	{
	keywords_.insert("VirtualFocusZSphere");
	keywords_.insert("ZPlane");
	keywords_.insert("ZCylinder");

	keywords_.insert("sameSkewness");
	keywords_.insert("sameEtendue");
	keywords_.insert("sameFlux");
	}; // default: empty


TInterpolateRaySetCfg::TInterpolateRaySetCfg()
	: TConfiguration{}
	{
	TRaySetControlSection rcs;
	rcs.InitAllowed();
	AddSection(std::move(rcs));
	}
