#ifndef __TM25Util_H
#define __TM25Util_H
#include <string>
#include <array>

namespace TM25
	{
	// assign string from fixed length char array, taking 0 termination into account
	template< int N>
	void AssignCharArrayToString(const std::array<char, N>& source, std::string& target);

	template< int N>
	void AssignChar32ArrayToString(const std::array<char32_t, N>& source, std::u32string& target);


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
