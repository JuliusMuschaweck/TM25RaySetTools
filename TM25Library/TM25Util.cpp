/******************************************************************
Provided by Julius Muschaweck, JMO GmbH, Gauting to the public domain
under the Unlicense, see unlicense.txt in the repository
or http://unlicense.org/
2019-01-19
******************************************************************/

// implementation of TM25Util.h

#ifndef _CONSOLE
#include <stdafx.h>
#endif
#include "TM25Util.h"
#include <vector>
#include <locale>
namespace TM25
	{

	std::string ToString(const std::u32string& s, char notranslation)
		{
		std::vector<char> rv(s.length());
		std::locale loc = std::locale();
		std::use_facet<std::ctype<std::u32string::value_type>>(loc).
			narrow(s.data(), s.data() + s.length(), notranslation, rv.data());
		return std::string(rv.data(), rv.size());
		}

	std::u32string ToU32String(const std::string& s)
		{
		std::vector<char32_t> rv(s.length());
		std::locale loc = std::locale();
		std::use_facet<std::ctype<std::u32string::value_type>>(loc).
			widen(s.data(), s.data() + s.length(), rv.data());
		return std::u32string(rv.data(), rv.size());
		}

	void AssignChar32ArrayToString(const std::vector<char32_t>& source, std::u32string& target)
		{
		target.clear();
		size_t N = source.size();
		for (int i = 0; i < N; ++i)
			{
			char32_t c = source[i];
			if (c != 0)
				target.push_back(c);
			else
				break;
			}
		}
	}
