#include "ParseString.h"

namespace 
	{
	inline bool IsDelim(char c, const std::string& delims)
		{
		for (auto p : delims)
			{
			if (c == p)
				return true;
			}
		return false;
		}
	}


std::vector<std::pair<size_t, size_t>> SplitString(const std::string& s, const std::string& delims)
	{
	std::vector<std::pair<size_t, size_t>> rv;
	const size_t start = 0;
	const size_t end = s.length();
	size_t cur = 0;
	auto icur = s.begin();
	while (cur != end)
		{
		// skip delims
		while (cur != end && IsDelim(*icur, delims))
			{
			++cur; ++icur;
			}
		if (cur != end) // found a non delim
			{
			size_t itemstart = cur;
			while (cur != end && !(IsDelim(*icur, delims)))
				{
				++cur; ++icur;
				}
			size_t itemend = cur;
			rv.push_back(std::make_pair(itemstart, itemend));
			}
		}
	return rv;
	}

