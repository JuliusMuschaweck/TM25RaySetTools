#ifndef __INTERPOLATERAYSET_IO_H
#define __INTERPOLATERAYSET_IO_H

#include <TM25.h>
#include "InterpolateRaySetConfig.h"
#include <iostream>
// #include <mutex>

TM25::TTM25RaySet ReadRaySet(const TInterpolateRaySetCfg& cfg, std::ostream& info);

void WriteRaySet(TM25::TTM25RaySet& rs, const TInterpolateRaySetCfg& cfg, std::ostream& info);



#endif
