#ifndef __PARSESTRING_H
#define __PARSESTRING_H

#include<vector>
#include<string>

// SplitString:
// Input: string s, string delim.
// A delimiter is any character in delim.
// s is a sequence of substrings which do not contain a delimiter
// separated by any number of delimiters
// starting and ending delimiters are trimmed.
// Output: vector of start/end position pairs. End position is one after last, as usual in C++
// rv may be empty. 
std::vector<std::pair<size_t, size_t>> SplitString(const std::string& s, const std::string& delims = " \t");



#endif
