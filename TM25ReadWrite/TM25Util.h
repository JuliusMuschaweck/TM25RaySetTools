/******************************************************************
Provided by Julius Muschaweck, JMO GmbH, Gauting to the public domain
under the Unlicense, see unlicense.txt in the repository
or http://unlicense.org/
2019-01-19
******************************************************************/

// TM25Util.h: Utility functions



#ifndef __TM25Util_H
#define __TM25Util_H
#include <string>
#include <array>
#include <vector>

namespace TM25
	{
	// assign string from fixed length char array, taking 0 termination into account
	template< int N>
	void AssignCharArrayToString(const std::array<char, N>& source, std::string& target);

	template< int N>
	void AssignChar32ArrayToString(const std::array<char32_t, N>& source, std::u32string& target);

	void AssignChar32ArrayToString(const std::vector<char32_t>& source, std::u32string& target);


	std::string ToString(const std::u32string& s, char notranslation = '?');

	std::u32string ToU32String(const std::string& s);



// template definitions

	template< int N>
	void AssignCharArrayToString(const std::array<char, N>& source, std::string& target)
		{
		target.clear();
		for (int i = 0; i < N; ++i)
			{
			char c = source[i];
			if (c != 0)
				target.push_back(c);
			else
				break;
			}
		}

	template< int N>
	void AssignChar32ArrayToString(const std::array<char32_t, N>& source, std::u32string& target)
		{
		target.clear();
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


#endif
