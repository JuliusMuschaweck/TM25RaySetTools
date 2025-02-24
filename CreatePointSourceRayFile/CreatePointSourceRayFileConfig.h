#ifndef __CREATEPOINTSOURCERAYFILECONFIG_H
#define __CREATEPOINTSOURCERAYFILECONFIG_H

#include <CfgFile.h>

class TPointSourceRayFileControlSection : public TSection
	{
	public:
		TPointSourceRayFileControlSection();
		virtual void AddAllowedValues();
		virtual void AddAllowedKeywords();
	};



class TPointSourceRayFileCfg : public TConfiguration
	{
	public:
		TPointSourceRayFileCfg();
	};

#endif

