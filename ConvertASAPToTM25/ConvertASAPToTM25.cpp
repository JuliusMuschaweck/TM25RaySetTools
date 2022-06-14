// ConvertASAPToTM25.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "ReadASAPRaySet.h"
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>

std::string ReadFN()
{
	wchar_t filename[MAX_PATH];

	OPENFILENAME ofn;
	ZeroMemory(&filename, sizeof(filename));
	ZeroMemory(&ofn, sizeof(ofn));
	ofn.lStructSize = sizeof(ofn);
	ofn.hwndOwner = NULL;  // If you have a window to center over, put its HANDLE here
	ofn.lpstrFilter = L"Text Files\0*.txt\0Any File\0*.*\0";
	ofn.lpstrFile = filename;
	ofn.nMaxFile = MAX_PATH;
	ofn.lpstrTitle = L"Select a File, yo!";
	ofn.Flags = OFN_DONTADDTORECENT | OFN_FILEMUSTEXIST;

	if (GetOpenFileNameW(&ofn))
	{
		return std::string(filename, filename + MAX_PATH);
	}
	else
	{
		return "dum";
	}
}

int main(int argc, char *argv[]) 
	{
	std::string fn = ReadFN();
	std::cout << fn << '\n';
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

