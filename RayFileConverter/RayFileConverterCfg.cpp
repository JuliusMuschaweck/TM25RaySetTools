#include"RayFileConverterCfg.h"

TRayFileConverterControlSection::TRayFileConverterControlSection() : TSection("RayFileConverterControl")
	{}

void TRayFileConverterControlSection::AddAllowedValues()
	{
	values_.insert({ "inputRayFileName",	MakeDefaultValueTokenSequence<Token::string>(std::string("inputRayFileName missing")) });
	values_.insert({ "inputRayFileFormat",	MakeDefaultValueTokenSequence<Token::string>(std::string("TM25")) });
	values_.insert({ "logFileName",			MakeDefaultValueTokenSequence<Token::string>(std::string("RayFileConverter.log")) });
	values_.insert({ "consoleOutput",		MakeDefaultValueTokenSequence<Token::boolean>(true) });
	values_.insert({ "outputRayFileName",	MakeEmptyTokenSequence() });
	values_.insert({ "outputRayFileFormat",	MakeDefaultValueTokenSequence<Token::string>(std::string("TM25")) });
	}

void TRayFileConverterControlSection::AddAllowedKeywords()
	{
	}

TRayFileConverterCfg::TRayFileConverterCfg()
	: TConfiguration{}
	{
	TRayFileConverterControlSection rcs;
	rcs.InitAllowed();
	AddSection(std::move(rcs));
	}


