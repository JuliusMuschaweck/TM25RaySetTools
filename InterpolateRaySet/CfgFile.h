#ifndef __CfgFile_h
#define __CfgFile_h
// general configuration file syntax
// line oriented
// Syntax
/*
% % this is a comment
% [Section] % anything after a % outside a string is a comment, too
% [xx] starts a section named [xx].
% followed by name=value pairs
% names are identifiers made of a..z, A..Z, 0..9, _ i.e. alnum + '_', starting with no digit
% blanks, i.e. chr(32) and tab='\t'=chr(9) characters are skipped
% values can be strings, integers, reals, integer and real vectors, the unquoted strings 'false' and 'true'
% and any program defined keywords (which are names)
% strings: "airhgai d ad a"
% integers: [+|-] 0 | (1..9 {0..9}) , see c++ stoi documentation
% reals: [+|-] 0..9 {0..9} [ . [0..9 {0..9}]] [ e|E [+|-] 0..9 {0..9}], see c++ stod documentation
% vectors: '{' real {; real}'}' for real vectors, '{' integer {; integer} '}' for integer vectors
% there is a special section [Variables] where variables are defined as name=value pairs
% these variables are then replaced by their values in later sections, on the right side of name=value assignments
% variable names must be unequal to program defined keywords
% later: escape sequences in strings, if then ...
% Example:
% (see implementation of void TestConfiguration() in CfgFile.cpp)

[Variables]
hopp="topp"
answerToQuestionOfLifeTheUniverseAndEverything=42
flag=true
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[SomeSectionYouExpectInYourProgram]
hopp=hopp % legal: variable replacement only on right side of assignment
irrational=sqrt(2)
doSomething=flag
*/

#include <set>
#include <string>
#include <map>
#include <vector>
#include <variant>
#include <iosfwd>
#include <type_traits>
#include "Tokenize.h"
#include <sstream>

using TKeywordSet = std::set<std::string>;

//struct TConfigValueType {
//	static const size_t keyword = 0;
//	static const size_t boolean = 1;
//	static const size_t string = 2;
//	static const size_t integer = 3;
//	static const size_t real = 4;
//	static const size_t intVector = 5;
//	static const size_t realVector = 6;
//	};
//
//using TConfigValue = std::variant<
//	std::string,
//	bool,
//	std::string,
//	int,
//	double,
//	std::vector<int>,
//	std::vector<double> >;
//
//using TConfigValueSet = std::map<std::string, TConfigValue>;

class TLineFile;

class TSection
	{
	public:
			// construction
		TSection() = default;
		explicit TSection(const std::string& name, bool unknownValuesAllowed = false, bool unknownKeywordsAllowed = false);
		void InitAllowed(); // calls virtual AddAllowed... functions
		std::string SetName(const std::string& name); // return previous name
		std::string Name() const { return name_; }; // return previous name

		bool SetUnknownValuesAllowed(bool rhs); // return previous
		bool SetUnknownKeywordsAllowed(bool rhs); // dito
			// read content from file. Can be called multiple times
		void ReadContent(TLineFile& f, const TSection& variables);
			// keywords: the right hand side identifiers
			// override in derived class to define which keywords are allowed
		virtual void AddAllowedKeywords() {}; // default: empty
		const TKeywordSet Keywords() const { return keywords_; };
		
			// values: the allowed value names for valuename=value entries
			// override in derived class to define which value names are allowed
			// sets allowed values to TTokenSequence with only endofline token
		virtual void AddAllowedValues() {};
		void AddAllowedValue(const std::string& n);
		bool ContainsValue(const std::string& n) const;
			// low level value access
			// throw if not found
		const TTokenSequence& Value(const std::string& n) const;
			// throw if unknownValuesAllowed== false and not in list of allowed values
		void SetValue(const std::string& name, TokenNS::TSIt begin, TokenNS::TSIt end);

			// high level value access
		bool IsEmpty(const std::string& n) const;
		std::string Keyword(const std::string& n) const;
		bool Bool(const std::string& n) const;
		std::string String(const std::string& n) const;
		int Int(const std::string& n) const;
		double Real(const std::string& n) const;
		std::vector<int> IntVector(const std::string& n) const;
		std::vector<double> RealVector(const std::string& n) const;

			// multivalues: the allowed value names for valuename{i,j,...}=something entries
		using TMultiValue = std::tuple<std::vector<int>, TTokenSequence>;
		virtual void AddAllowedMultiValues(){}; // default none allowed
		void AddAllowedMultiValue(const std::string& name);
		bool ContainsMultiValue(const std::string& name) const;
		std::vector<TMultiValue> MultiValue(const std::string& n) const;
		void SetMultiValue(const std::string& name, const TMultiValue& val);

		virtual ~TSection() = default;
	protected:
		template<typename TRV, Token Tok>
		TRV GetSingleTokenValue(const std::string& valueName) const;

		template<typename Num>
		std::pair<std::vector<Num>, TokenNS::TSIt> ParseVector(TokenNS::TSIt begin, TokenNS::TSIt end,
			const std::string& entryname) const;
		void ParseLine(const TTokenSequence& t, size_t currentLine, const std::string& fn, const TSection& variables);
		bool unknownValuesAllowed_ = false;
		bool unknownKeywordsAllowed_ = false;
		TKeywordSet keywords_;
		std::string name_;
		using TValueSet = std::map<std::string, TTokenSequence>;
		TValueSet values_;
		using TMultiValueSet = std::map<std::string, std::vector<TMultiValue>>;
		TMultiValueSet multiValues_;
	};

class TConfiguration
	{
	public:
		TConfiguration();

		using CStrRef = const std::string&;
		void ParseCfgFile(CStrRef fn);
		std::string Content() const;

		bool ContainsSection(CStrRef name) const;
		const TSection& Section(CStrRef name) const;		
		// use for special TSection derived classes
		void AddSection(TSection&& s); //
		// use for standard TSection derived classes which provide virtual AddAllowed... funcitons
		template<typename SectionType> 
		void AddSection();
		
		bool UnknownSectionsAllowed() const; // default: false
		bool SetUnknownSectionsAllowed(bool yesno); // return previous value

		// const TConfigValue& Value(CStrRef section, CStrRef name) const; // convenience

		virtual ~TConfiguration() = default;
	protected:
		std::map<std::string, TSection> sections_; // default: empty Variables section
		bool unknownSectionsAllowed_;
		std::stringstream cfgContent_;

	};

void TestConfiguration();

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// template definitions
template<>
double TSection::GetSingleTokenValue<double, Token::real>(const std::string& valueName) const;

template<>
int TSection::GetSingleTokenValue<int, Token::integer>(const std::string& valueName) const;

template<typename TRV, Token Tok>
TRV TSection::GetSingleTokenValue(const std::string& valueName) const
	{
	const TTokenSequence& v = Value(valueName);
	if (v.size() != 2)
		throw std::runtime_error("TSection::GetSingleTokenValue: token sequence too long for value " + valueName
		+ " in section " + name_);
	Token t = std::get<Token>(v.front());
	if (t != Tok)
		throw std::runtime_error("TSection::GetSingleTokenValue: Wrong token type for value " + valueName
		+ ", expected " + TokenNS::TokenToString(Tok) + ", found " + TokenNS::TokenToString(t) + " in section " + name_);
	using iTRV = TRV;
	return std::get<iTRV>(std::get<TTokenValue>(v.front()));
	}

template<typename Num>
std::pair<std::vector<Num>, TokenNS::TSIt> TSection::ParseVector(TokenNS::TSIt begin, TokenNS::TSIt end,
	const std::string& entryname) const
	{
	using TNum = Num;
	static_assert(std::is_same<TNum, double>() || std::is_same<TNum, int>());
	auto ConvertToTNum = [&entryname, this](TToken t) -> TNum
		{ //precond. t is double or integer
		if (std::is_same<TNum, double>())
			{
			if (std::get<Token>(t) == Token::integer)
				{
				bool ok;
				t = MakeRealIf(t, ok);
				}
			}
		else // TNum == int
			{
			if (std::get<Token>(t) == Token::real)
				{
				bool ok;
				t = MakeIntegerIf(t, ok);
				if (!ok)
					throw std::runtime_error("TSection::ParseVector: cannot convert real "
					+ std::to_string(std::get<double>(std::get<TTokenValue>(t)))
					+ " to integer for value " + entryname + " in section " + this->name_);
				}
			}
		return std::get<TNum>(std::get<TTokenValue>(t));
		};
	auto TNumStr = []()->std::string
		{
		if (std::is_same<TNum, double>())
			return "real";
		else
			return "integer";
		};
	auto AssertNumToken = [&entryname, this, TNumStr](TokenNS::TSIt found)
		{
		if ((Token::integer != std::get<Token>(*found)) && (Token::real != std::get<Token>(*found)))
			throw std::runtime_error("TSection::ParseVector: expected " + TNumStr() +
			", found " + TokenNS::TokenToString(std::get<Token>(*found)) + " for entry " + entryname +
			" at position " + std::to_string(std::get<size_t>(*found)) + " in section " + this->name_);
		};
	auto AssertToken = [&entryname, this](Token expected, TokenNS::TSIt found)
		{
		if (expected != std::get<Token>(*found))
			throw std::runtime_error("TSection::IntVector: expected " + TokenNS::TokenToString(expected) +
			", found " + TokenNS::TokenToString(std::get<Token>(*found)) + " for entry " + entryname +
			" at position " + std::to_string(std::get<size_t>(*found)) + " in section " + this->name_);
		};
	std::vector<TNum> rv;
	TokenNS::TSIt cur = begin;	
	AssertToken(Token::lBrace, cur);
	cur = AdvanceToken(cur, end);
	AssertNumToken(cur);
	rv.push_back(ConvertToTNum(*cur));
	cur = AdvanceToken(cur, end);
		while (std::get<Token>(*cur) == Token::comma)
		{
		cur = AdvanceToken(cur, end);
		AssertNumToken(cur);
		rv.push_back(ConvertToTNum(*cur));
		cur = AdvanceToken(cur, end);
		}
	AssertToken(Token::rBrace, cur);
	cur = AdvanceToken(cur, end);
	return std::make_pair(rv, cur);
	}


template<typename SectionType>
void TConfiguration::AddSection()
	{
	using MySection = SectionType;
	MySection s;
	s.InitAllowed();
	AddSection(std::move(s));
	}

// include guard
#endif 