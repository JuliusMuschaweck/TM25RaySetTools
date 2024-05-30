#ifndef __RAYFILECONVERTERCONFIG_H
#define __RAYFILECONVERTERCONFIG_H

#include <CfgFile.h>

class TRayFileConverterControlSection : public TSection
	{
	public:
		TRayFileConverterControlSection();
		virtual void AddAllowedValues();
		virtual void AddAllowedKeywords();
	};



class TRayFileConverterCfg : public TConfiguration
	{
	public:
		TRayFileConverterCfg();
	};

#endif