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
	SplitString(rv, s, delims);
	return rv;
	}

void SplitString(std::vector<std::pair<size_t, size_t>>& rv, const std::string& s, const std::string& delims)
	{
	rv.resize(0); // in typical C++ library, size will be zero but capacity remains unchanged.
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
	}

std::vector<double> StringToVector(const std::string& s, const std::string& delims)
// s: sequence of real numbers separated by delims
// returns vector of these real numbers, possibly empty
// throws std::logic_error if format is wrong
	{
	std::vector<std::pair<size_t, size_t>> pos = SplitString(s, delims);
	std::vector<double> rv;
	const char* const s0 = s.c_str();
	for (std::pair<size_t, size_t> ipos : pos)
		{
		rv.push_back(CharRangeToDouble(s0 + ipos.first, s0 + ipos.second));
		}
	return rv;
	}

double CharRangeToDouble(const char* const  begin, const char* const end)
// thin wrapper around std::stod
	{
	return std::stod(std::string(begin, end));
	}

