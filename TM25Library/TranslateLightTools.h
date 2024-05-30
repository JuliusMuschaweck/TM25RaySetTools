#ifndef __TRANSLATELIGHTTOOLS_H
#define __TRANSLATELIGHTTOOLS_H

#include <TM25.h>
#include <LightToolsBinary.h>

namespace TM25
	{
	TTM25RaySet LightToolsBinaryToTM25(const TLightToolsRaySet& zemaxRaySet);

	TLightToolsRaySet TM25ToLightToolsBinary(const TTM25RaySet& rs);
	}
#endif