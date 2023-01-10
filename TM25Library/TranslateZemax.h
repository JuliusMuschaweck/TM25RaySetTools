#ifndef __TRANSLATEZEMAX_H
#define __TRANSLATEZEMAX_H

#include "TM25.h"
#include "ZemaxBinary.h"

namespace TM25
	{
	TTM25RaySet ZemaxBinaryToTM25(const TZemaxRaySet& zemaxRaySet);

	TZemaxRaySet TM25ToZemaxBinary(const TTM25RaySet& rs);
	}
#endif