#ifndef __PARSESTRING_H
#define __PARSESTRING_H

#include<vector>
#include<string>
#include<stdexcept>

// SplitString:
// Input: string s, string delim.
// A delimiter is any character in delim.
// s is a sequence of substrings which do not contain a delimiter
// separated by any number of delimiters
// starting and ending delimiters are trimmed.
// Output: vector of start/end position pairs. End position is one after last, as usual in C++
// rv may be empty. 
std::vector<std::pair<size_t, size_t>> SplitString(const std::string& s, const std::string& delims = " \t");

void SplitString(std::vector<std::pair<size_t, size_t>>& rv, const std::string& s, const std::string& delims = " \t");
// use this overload to reuse rv, saving the allocate/deallocate effort

std::vector<double> StringToVector(const std::string& s, const std::string& delims = " \t");
// s: sequence of real numbers separated by delims
// returns vector of these real numbers, possibly empty
// throws std::runtime_error if format is wrong

template<size_t N>
std::array<double, N> StringToArray(std::vector<std::pair<size_t, size_t>>& pos, const std::string& s, const std::string& delims = " \t");
// Convert string containing known-length-sequence of numbers without memory allocation overhead
// s: sequence of real numbers separated by delims
// returns vector of these real numbers
// throws std::logic_error if format is wrong or if the number of reals is not N

double CharRangeToDouble(const char* const  begin, const char* const end);
// thin wrapper around std::stod
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


template<size_t N>
std::array<double, N> StringToArray(std::vector<std::pair<size_t, size_t>>& pos, const std::string& s, const std::string& delims)
    {
	SplitString(pos, s, delims);
	if (pos.size() != N)
		throw std::runtime_error("StringToArray: expected " + std::to_string(N) + " substrings");
	const char* const s0 = s.c_str();
	size_t i = 0;
	std::array<double, N> rv;
	for (std::pair<size_t, size_t> ipos : pos)
		{
		rv[i++] = CharRangeToDouble(s0 + ipos.first, s0 + ipos.second);
		}
	return rv;
	}



#endif
