#include "playground.h"
#include"TM25.h"
#include <iostream>
#include <string>


void playground()
{
	TM25::TTM25RaySet rs;
	// this small ray file is on GitHub
	std::string fn;
	std::cout << "enter path and file name: ";
	std::cin >> fn;

	rs.Read(fn);
	for (auto w : rs.Warnings())
	{
		std::cout << w << "\n";
	}
	auto bb = rs.RayArray().BoundingBox();
	std::cout << "Bounding Box: x in [" << bb.first[0] << ',' << bb.second[0] << "], y in ["
		<< bb.first[1] << ',' << bb.second[1] << "], z in [" << bb.first[2] << ',' << bb.second[2] << "]" << std::endl;
	auto vf = rs.VirtualFocus();
	std::cout << "Virtual focus: (" << vf[0] << ',' << vf[1] << ',' << vf[2] << ")\n";
}