#ifndef __INTERPOLATERAYSETCONFIG_H
#define __INTERPOLATERAYSETCONFIG_H

#include <CfgFile.h>

class TRaySetControlSection : public TSection
	{
	public:
		TRaySetControlSection();
		virtual void AddAllowedValues();
		virtual void AddAllowedKeywords();
	};



class TInterpolateRaySetCfg : public TConfiguration
	{
	public:
		TInterpolateRaySetCfg();
	};

#endif
