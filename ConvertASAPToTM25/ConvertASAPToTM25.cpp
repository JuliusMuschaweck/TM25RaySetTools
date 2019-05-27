// ConvertASAPToTM25.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "ReadASAPRaySet.h"

int main(int argc, char *argv[]) 
	{
	if (argc != 3)
		{
		std::cout << "ConvertASAPToTM25: call me ConvertASAPToTM25 <ASAP DIS file> <TM25 file>";
		return 1;
		}
	try
		{
		auto rs = ReadASAPRaySet(argv[1]);
		rs.Write(argv[2]);

		TM25::TTM25RaySet rs2;
		rs2.Read(argv[2]);
		auto sc = rs2.Header().SanityCheck();
		std::cout << rs2.NRays() << '\n';
		}
	catch (std::runtime_error err)
		{
		std::cout << "Error: " << err.what() << "\n";
		}
	}

