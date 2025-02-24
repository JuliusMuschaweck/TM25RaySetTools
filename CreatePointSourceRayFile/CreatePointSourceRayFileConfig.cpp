#include "CreatePointSourceRayFileConfig.h"

TPointSourceRayFileControlSection::TPointSourceRayFileControlSection() : TSection("PointSourceRayFileControl") {};

void TPointSourceRayFileControlSection::AddAllowedValues()
	{
	values_.insert({ "inputIntensityFileName",		MakeDefaultValueTokenSequence<Token::string>(std::string("inputRayFileName missing")) });
	values_.insert({ "delimiters",					MakeDefaultValueTokenSequence<Token::string>(std::string(" \t,;")) });
	values_.insert({ "equidistance_rel_epsilon",	MakeDefaultValueTokenSequence<Token::real>(1e-11) });
	values_.insert({ "logFileName",					MakeDefaultValueTokenSequence<Token::string>(std::string("CreatePointSourceRayFile.log")) });
	values_.insert({ "consoleOutput",				MakeDefaultValueTokenSequence<Token::boolean>(true) });
	values_.insert({ "doDiagnostics",				MakeDefaultValueTokenSequence<Token::boolean>(true) });
	values_.insert({ "nOutputRays",				MakeDefaultValueTokenSequence<Token::integer>(0) });
	values_.insert({ "rel_threshold",			MakeDefaultValueTokenSequence<Token::real>(0.1) });
	values_.insert({ "outputRayFileName",		MakeEmptyTokenSequence() });
	values_.insert({ "outputRayFileFormat",		MakeDefaultValueTokenSequence<Token::string>(std::string("TM25")) });
	values_.insert({ "setTotalFlux",			MakeEmptyTokenSequence() });
	}

void TPointSourceRayFileControlSection::AddAllowedKeywords()
	{
	}; // default: empty


TPointSourceRayFileCfg::TPointSourceRayFileCfg()
	: TConfiguration{}
	{
	TPointSourceRayFileControlSection rcs;
	rcs.InitAllowed();
	AddSection(std::move(rcs));
	}
