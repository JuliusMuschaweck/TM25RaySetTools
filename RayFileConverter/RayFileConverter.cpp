// RayFileConverter.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <TM25.h>
#include <RaySet_IO.h>
#include "RayFileConverterCfg.h"


int main(int argc, char* argv[])
    {
	using std::cout;
	using std::endl;
	cout << "RayFileConverter, by Julius Muschaweck, 2024" << endl;
	if (argc < 2)
		{
		cout << "Call me 'RayFileConverter cfg_fn', where cfgFn is the name of the configuration file" << endl;
		return 0;
		}
	try
		{
		std::string cfgFn{ argv[1] };
		cout << "parsing configuration file " << cfgFn << endl;
		TRayFileConverterCfg cfg;
		cfg.ParseCfgFile(cfgFn);
		const TSection& rsc = cfg.Section("RayFileConverterControl");
		TLogPlusCout logS(rsc.Bool("consoleOutput"), rsc.String("logFileName"));
		// read the raw input ray file
		logS << "reading ray file " << rsc.String("inputRayFileFormat") << endl;
		TM25::TTM25RaySet rs = ReadRaySet(rsc, logS);
		logS << "writing ray file " << rsc.String("outputRayFileFormat") << endl;
		WriteRaySet(rs, rsc, logS);
		}
	catch (TM25::TM25Error e)
		{
		std::cout << e.what() << std::endl;
		}
	catch (std::runtime_error e)
		{
		std::cout << "std::runtime_error: " << e.what() << std::endl;
		}
	//catch (...)
	//	{
	//	std::cout << "unknown error" << std::endl;
	//	}
	return 0;

     }

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
