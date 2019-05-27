#include "PhaseSpace.h"

	void TestPhaseSpace()
		{
//		TPlanarPhaseSpace<float> ps;
		TZAxisStereographicSphericalPhaseSpace<double> ps2;
		ps2 = TZAxisStereographicSphericalPhaseSpace<double>({ 0,0,2 }, 2);
		TVec3d loc3{ 1,2,-0.5 };
		TVec3d dir3{ 0,0,-1 };
		TZAxisStereographicSphericalPhaseSpace<double>::TLoc loc2{ 4 / 3.,8 / 3. };
		auto dir2 = ps2.Dir(loc2, dir3);
		TVec3d dir3test = ps2.Dir_to_V3(loc2, dir2);
		};
