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
#include <numeric>
#include <algorithm>

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

	// sort v, but return index of permutation instead of modifying v
	// postcondition: i < j => v[rv[i]] <= v[rv[j]]
	template<typename Compare = decltype(std::less())>
	std::vector<size_t> IndexSort(const std::vector<float>& v, Compare comp = std::less());

	std::string_view TrimWhiteSpace(std::string_view s); // trims ' ' and '\t' from left and right

	std::vector<std::string_view> Tokenize(std::string_view s, std::string_view delims);

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

	// sort v, but return index of permutation instead of modifying v
	// postcondition: i < j => v[rv[i]] <= v[rv[j]]
	template<typename Compare>
	std::vector<size_t> IndexSort(const std::vector<float>& v, Compare comp)
		{
		std::vector<size_t> rv(v.size());
		std::iota(rv.begin(), rv.end(), 0);
		auto pred = [&v, comp](size_t lhs, size_t rhs) -> bool {return comp(v[lhs], v[rhs]); };
		std::sort(rv.begin(), rv.end(), pred);
		return rv;
		}

	}


#endif
