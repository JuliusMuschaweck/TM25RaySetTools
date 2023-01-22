#include "PhaseSpaceGrid.h"
#include "InterpolateRaySetData.h"
#include <utility>
using namespace KDTree;


float AverageLuminance(const TInterpolateRaySetData& raySetData, const Def::TPointIdxArray& pts)
	{
	if (pts.empty())
		return 0;
	float volume = 0;
	float flux = 0;
	for (auto pi : pts)
		{
		volume += raySetData.volumes_[pi.pi_];
		flux += raySetData.cellFluxes_[pi.pi_];
		}
	float rv = flux / volume;
	return rv;
	}

float TotalFlux(const TInterpolateRaySetData& raySetData, const Def::TPointIdxArray& pts)
	{
	if (pts.empty())
		return 0;
	float flux = 0;
	for (auto pi : pts)
		{
		flux += raySetData.cellFluxes_[pi.pi_];
		}
	return flux;
	}

TPhaseSpaceGrid<float> FillDefaultPhaseSpaceGrid(const TInterpolateRaySetData& raySetData, const size_t nBins)
	{
	Def::TBox overallBB = BoundingBox(raySetData.kdtree_->Points());
	return FillPhaseSpaceGrid(overallBB, raySetData, {nBins, nBins, nBins, nBins});
	}


TPhaseSpaceGrid<float> FillPhaseSpaceGrid(const TPhaseSpaceGrid<float>::TBox& boundingBox,
	const TInterpolateRaySetData& raySetData, const TPhaseSpaceGrid<float>::T4DIndex& nBins)
	{
	// we create prod(nBins) boxes in the ray set's phase space, with prod(nBins+1) corner points.
	// we compute luminance averaged over each box
	// The phase space grid will consist of the nBins^4 center points of the boxes
	// we assign the corresponding luminance to each box

	// grid step sizes
 	float dx0 = (boundingBox.second[0] - boundingBox.first[0]) / nBins[0];
	float dx1 = (boundingBox.second[1] - boundingBox.first[1]) / nBins[1];
	float dk0 = (boundingBox.second[2] - boundingBox.first[2]) / nBins[2];
	float dk1 = (boundingBox.second[3] - boundingBox.first[3]) / nBins[3];

	// bounding box of phase space grid center points
	TPhaseSpaceGrid<float>::T4DPoint gridBB0
		{
		boundingBox.first[0] + dx0 / 2,
		boundingBox.first[1] + dx1 / 2,
		boundingBox.first[2] + dk0 / 2,
		boundingBox.first[3] + dk1 / 2
		};

	TPhaseSpaceGrid<float>::T4DPoint gridBB1
		{
		boundingBox.second[0] - dx0 / 2,
		boundingBox.second[1] - dx1 / 2,
		boundingBox.second[2] - dk0 / 2,
		boundingBox.second[3] - dk1 / 2
		};

	TPhaseSpaceGrid<float> rv(gridBB0, gridBB1, nBins);

	size_t maxpts = 0;
	size_t nNonzero = 0;
	for (size_t ix0 = 0; ix0 < nBins[0]; ++ix0)
		{
		std::cout << "FillPhaseSpaceGrid: ix0 = " << ix0 << " of " << nBins[0] << '\n';
		float thisx0 = boundingBox.first[0] + ix0 * dx0;
		for (size_t ix1 = 0; ix1 < nBins[1]; ++ix1)
			{
			float thisx1 = boundingBox.first[1] + ix1 * dx1;
			for (size_t ik0 = 0; ik0 < nBins[2]; ++ik0)
				{
				float thisk0 = boundingBox.first[2] + ik0 * dk0;
				for (size_t ik1 = 0; ik1 < nBins[3]; ++ik1)
					{
					float thisk1 = boundingBox.first[3] + ik1 * dk1;
					Def::TBox thisBB = std::make_pair(
						Def::TKDPoint{ thisx0, thisx1, thisk0, thisk1 },
						Def::TKDPoint{ thisx0 + dx0, thisx1 + dx1, thisk0 + dk0, thisk1 + dk1 }
					);
					//// Def::TPointIdxArray pts = raySetData.kdtree_->LocateOverlapping(thisBB);
					//Def::TPointIdxArray pts = raySetData.kdtree_->LocatePointsWithinBox(thisBB);
					//if (pts.size() > maxpts) maxpts = pts.size();
					//float L = AverageLuminance(raySetData, pts);
					//float U = dx0 * dx1 * dk0 * dk1; // the etendue of this box
					//// float L = TotalFlux(raySetData, pts) / U;
					const TKDTree& kdt = *(raySetData.kdtree_);
					Def::TNodeIdxArray overlappingLeafNodes = kdt.LocateOverlappingLeafNodes(thisBB);
					if (overlappingLeafNodes.empty())
						rv.Value({ ix0, ix1, ik0, ik1 }) = 0;
					else
						{
						if (maxpts < overlappingLeafNodes.size())
							maxpts = overlappingLeafNodes.size();
						Def::TPointIdxArray pts;
						for (auto ni : overlappingLeafNodes)
							{
							Def::TIdxIdx ii = kdt.Node(ni).beginpt_;
							pts.push_back(kdt.Index(ii));
							}
						float L = AverageLuminance(raySetData, pts);
						rv.Value({ ix0, ix1, ik0, ik1 }) = L;
						++nNonzero;
						}
					}
				}
			}
		}
	std::cout << "FillPhaseSpaceGrid: max # of points per box = " << maxpts << '\n';

	std::cout << "FillPhaseSpaceGrid: nonzero boxes = " << nNonzero << " out of " << rv.v_.size() <<'\n';

 	return rv;
	}