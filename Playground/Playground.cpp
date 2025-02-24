//
// Playground.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "pch.h"
#include <string>
#include <iostream>
#include <fstream>
#include <array>
#include <vector> 
#include "../InterpolateRaySet/KDTree.h"
#include "ParseString.h"
#include <ReadFile.h>
#include <ASCIIRayFile.h>
#include <Tokenize.h>
#include <ReadASCIIMatrix.h>


void TestBuf()
	{
	//std::ifstream f("..\\rayfile_LERTDUW_S2WP_20161017_IES_TM25.zip", std::ios::binary);
	std::ifstream f("..\\rayfile_LERTDUW_S2WP_blue_100k_20161013_IES_TM25.TM25RAY", std::ios::binary);
	std::vector<char> buf(1024 * 1024 * 10);
	size_t n = buf.size();
	size_t n_read = f.rdbuf()->sgetn(buf.data(), n);
	n_read = f.rdbuf()->sgetn(buf.data(), n);
	bool b = f.eof(); // false! buffer knows nothing about file which owns it
	}

int main()
{
	std::string fn2 = "..\\TM25Library\\TestASCIIMatrix.txt";
	TM25::TNumberList nl;
	nl.Read(fn2);

	TM25::TVectorMatrixList vml;
	vml.Read(fn2);

	TM25::TReadASCIIMatrix am;
	am.Read(fn2, " \t,;", 1e-10);

	double test = am.Interpolate(1, 10); 
	test = am.Interpolate(1.1, 10);
	test = am.Interpolate(1, 20.1);
	test = am.Interpolate(2, 10);


	std::vector<float> dum;
	std::cout << dum.capacity() << '\n';
	dum.push_back(2.0f);
	std::cout << dum.capacity() << '\n';


	std::string s;
	size_t tmp2 = sizeof(s); // 40 bytes, buffer[16]
	std::cout << tmp2 << ' ' << (s.begin() == s.end()) << ' ' << s.capacity() <<
		' ' << s.empty() << ' ' << sizeof(s.begin()) << ' ' << reinterpret_cast<uint64_t>(s.data()+0)
		<< ' ' << static_cast<uint8_t>(*(s.data())) << 'x' << std::endl;
	char* d1 = s.data();
	const char* d3 = s.c_str();
	std::cout << (d1 == d3) << std::endl;
	s = "c";
	char* d2 = s.data();
	std::cout << (d1 == d2) << std::endl;
	d3 = s.c_str();
	std::cout << (d1 == d3) << std::endl;

	KDTree::TestKDTree4D();
	TestTokenize();
	
	std::string fn = "../TestRayFiles/ZemaxText.txt";
	TM25::TASCIIRaySetOptions opts;
	opts.minHeaderLines = 1;
	opts.wavelengthUnit_ = TM25::TASCIIRaySetOptions::WLU::micrometer;
	TM25::TTM25RaySet rs = TM25::ReadGenericASCIIRaySet(fn, opts);
	fn = "../TestRayFiles/ZemaxWavelengthText.txt";
	TM25::TTM25RaySet rs2 = TM25::ReadGenericASCIIRaySet(fn, opts);
	fn = "../TestRayFiles/deleteme.txt";
	TM25::TTM25RaySet rs3 = TM25::ReadGenericASCIIRaySet(fn, opts);

	auto items = SplitString(""," \t");
	items = SplitString(" abc ", " \t");
	items = SplitString("abc ", " \t");
	items = SplitString(" abc", " \t");
	items = SplitString(" abc def\tg\t h", " \t");

	double d = std::stod("123.45");
	std::vector<double> tmp;
	std::cout << sizeof(std::vector<double>) << '\n';
	std::cout << sizeof(std::string) << '\n';

	TM25::TReadFile rf("../TM25Library/Timer.h");
	while (!rf.AtEof())
		{
		std::string s = rf.ReadLine(true);
		std::cout << s;
		}

	std::cout << "Hello World!\n";
	if (KDTree::Def::dim == 2)
		KDTree::TestKDTree2D("TestKDTree.m");
	if (KDTree::Def::dim == 4)
		KDTree::TestKDTree4D();
	TestBuf();
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
