/******************************************************************
Provided by Julius Muschaweck, JMO GmbH, Gauting to the public domain
under the Unlicense, see unlicense.txt in the repository
or http://unlicense.org/
2019-01-19
******************************************************************/

// TM25ReadWrite.cpp : Defines the entry point for the console application.
// and provides a typical usage example.

#include <TM25.h>
#include <iostream>
int main()
{	
	try
		{
 		TestLinAlg3();
		TM25::TTM25RaySet rs;
		// this small ray file is on GitHub
		rs.Read("../rayfile_LERTDUW_S2WP_green_100k_20161013_IES_TM25.TM25RAY");
		// this 500 MB ray file with 20 million rays can be downloaded from www.osram-os.com
		//rs.Read("../rayfile_LERTDUW_S2WP_green_20M_20161013_IES_TM25.TM25RAY");
		for (auto w : rs.Warnings())
			{
			std::cout << w << "\n";
			}
		auto bb = rs.RayArray().BoundingBox();
		std::cout << "Bounding Box: x in [" << bb.first[0] << ',' << bb.second[0] << "], y in ["
			<< bb.first[1] << ',' << bb.second[1] << "], z in [" << bb.first[2] << ',' << bb.second[2] << "]"<< std::endl;
		
		TM25::TTM25RaySet rs2;
		// this small ray file is on GitHub
		rs2.Read("../CreatePointSourceRayFile/farfield.TM25RAY");
		for (auto w : rs2.Warnings())
			{
			std::cout << w << "\n";
			}


		
		}
	catch (TM25::TM25Error e)
		{
		std::cout << e.what() << std::endl;
		}
	catch (std::runtime_error e)
		{
		std::cout << "std::runtime_error: " << e.what() << std::endl;
		}
	catch (...)
		{
		std::cout << "unknown error" << std::endl;
		}
	return 0;
}

