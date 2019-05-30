#include "PhaseSpace.h"

	void TestPhaseSpace()
		{
//		TPlanarPhaseSpace<float> ps;
		using TPS = TZAxisStereographicSphericalPhaseSpace<double>;
		TPS ps2;
		TVec3d center{ 0,0,2 };
		double radius = 2;
		ps2 = TPS(center, radius);
		TVec3d loc3{ 1,1.1,-0.5 };
		TVec3d dir3{ 0,0,-1 };
		bool b;
		TVec3d loc3Sphere;
		std::tie(b,loc3Sphere) = SphereIntersect<double>(center, radius, loc3, dir3);

		TPS::TLoc loc2 = ps2.Loc(loc3, dir3);
		auto loc3test = ps2.Loc_to_V3(loc2);
		auto dir2 = ps2.Dir(loc2, dir3);
		TVec3d dir3test = ps2.Dir_to_V3(loc2, dir2);


		};
