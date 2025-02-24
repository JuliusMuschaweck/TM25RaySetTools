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

	std::string_view TrimWhiteSpace(std::string_view s)
		{
		std::string_view rv(s);
		rv.remove_prefix(std::min(rv.find_first_not_of(" \t"), rv.size()));
		auto trim_pos = rv.find_last_not_of(" \t");
		rv.remove_suffix(rv.size() - trim_pos - 1);
		return rv;
		}

	std::vector<std::string_view> Tokenize(std::string_view s, std::string_view delims)
		{
		std::vector<std::string_view> rv;
		size_t pos = s.find_first_not_of(delims);
		size_t end = s.size();
		while (pos < end)
			{
			size_t count = s.find_first_of(delims, pos) - pos;
			rv.push_back(s.substr(pos, count));
			pos = s.find_first_not_of(delims, pos+count);
			}
		return rv;
		}

	} // namespace TM25

