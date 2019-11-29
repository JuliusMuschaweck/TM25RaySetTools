#include "CfgFile.h"

#include<fstream>
#include<chrono>

// line oriented file input helper class

class TLineFile
	{
	public:
		explicit TLineFile(const std::string& fn);
		std::string  Buf() const { return buf_; };
		size_t LineNo() const { return lineno_; };
		bool Eof() const { return (!pushback_flag_ && f_.eof()); };
		void NextLine();
		void PushBackLine();
		std::string Filename()const { return fn_; };
	private:
		std::string fn_;
		std::ifstream f_;
		size_t lineno_;
		static const size_t bufsize = 10000;
		char buf_[bufsize];
		bool pushback_flag_ = false;
	};

TLineFile::TLineFile(const std::string& fn)
	: fn_(fn), f_(fn), lineno_(0)
	{
	if (!f_.is_open())
		throw std::runtime_error("TLineFile::TLineFile: cannot open " + fn);
	}

void TLineFile::NextLine()
	{
	if (Eof())
		throw std::runtime_error("TLineFile::NextLine: at eof at line " + std::to_string(lineno_)
		+ " when reading " + fn_);
	if (pushback_flag_)
		{
		pushback_flag_ = false;
		return;
		}
	++lineno_;
	f_.getline(buf_, bufsize);
	if (!f_.good() && !Eof())
		throw std::runtime_error("TLineFile::NextLine: read failure (line too long?) when reading " + fn_
		+ " at line no. " + std::to_string(lineno_));
	}

void TLineFile::PushBackLine()
	{
	if (pushback_flag_)
		throw std::runtime_error("TLineFile::PushBackLine: cannot push back twice");
	pushback_flag_ = true;
	}



// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// TSection implementation

// construction

TSection::TSection(const std::string& name, bool unknownValuesAllowed, bool unknownKeywordsAllowed)
	:name_(name), unknownValuesAllowed_(unknownValuesAllowed), unknownKeywordsAllowed_(unknownKeywordsAllowed)
	{
	}

void TSection::InitAllowed()
	{
	AddAllowedKeywords();
	AddAllowedValues();
	AddAllowedMultiValues();
	};

std::string TSection::SetName(const std::string & name)
	{
	std::string rv = name_;
	name_ = name;
	return rv;
	}

bool TSection::SetUnknownValuesAllowed(bool rhs)
	{
	bool rv = unknownValuesAllowed_;
	unknownKeywordsAllowed_ = rhs;
	return rv;
	}

bool TSection::SetUnknownKeywordsAllowed(bool rhs)
	{
	bool rv = unknownKeywordsAllowed_;
	unknownKeywordsAllowed_ = rhs;
	return rv;
	}

void TSection::ReadContent(TLineFile& f, const TSection& variables)
	{
	while (!f.Eof())
		{
		f.NextLine();
		TTokenSequence t = Tokenize(f.Buf()); // t contains at least endofline -- cannot be 		TokenNS::TSIt cur = t.begin();
		TokenNS::TSIt cur = t.begin();
		TokenNS::TSIt end = t.end();
		if (std::get<Token>(*cur) == Token::endofline) // empty line, only space and comments
			{
			continue;
			}
		if (std::get<Token>(*cur) == Token::lBracket) // next section starts
			{
			f.PushBackLine();
			break;
			}
		ParseLine(std::move(t), f.LineNo(), f.Filename(), variables);
		}
	}

// values

void TSection::AddAllowedValue(const std::string & n)
	{
	values_[n] = MakeEmptyTokenSequence();
	}

bool TSection::ContainsValue(const std::string & n) const
	{
	return values_.find(n) != values_.end();
	}

const TTokenSequence & TSection::Value(const std::string & n) const
	{
	auto it = values_.find(n);
	if (it != values_.end())
		return it->second;
	else
		{
		throw std::runtime_error("TSection::Value: entry not found: " + n + " in section " + name_);
		}
	}
void TSection::SetValue(const std::string & name, TokenNS::TSIt begin, TokenNS::TSIt end)
	{
	if (unknownValuesAllowed_)
		values_[name] = TTokenSequence(begin, end);
	else
		{
		auto it = values_.find(name);
		if (it != values_.end())
			it->second = TTokenSequence(begin, end);
		else
			throw std::runtime_error("TSection::SetValue: value not allowed: " + name + " in section " + name_);
		}
	}

// high level value access

// helper function

bool TSection::IsEmpty(const std::string& n) const
	{
	const TTokenSequence& v = Value(n);
	if (v.size() != 1)
		return false;
	Token t = std::get<Token>(v.front());
	return (t == Token::endofline);
	}

std::string TSection::Keyword(const std::string & n) const
	{
	return GetSingleTokenValue< std::string, Token::identifier>(n);
	}

bool TSection::Bool(const std::string & n) const
	{
	return GetSingleTokenValue< bool, Token::boolean>(n);
	}

std::string TSection::String(const std::string & n) const
	{
	return GetSingleTokenValue< std::string, Token::string>(n);
	}

int TSection::Int(const std::string & n) const
	{
	return GetSingleTokenValue< int, Token::integer>(n);
	}

double TSection::Real(const std::string & n) const
	{
	return GetSingleTokenValue< double, Token::real>(n);
	}

std::string Vec2Str(const std::vector<int> vi)
	{
	std::string rv = "(";
	bool start = true;
	for (auto i : vi)
		{
		if (!start)
			rv.append(";");
		rv.append(std::to_string(i));
		}
	rv.append(")");
	return rv;
	};
std::string Vec2Str(const std::vector<double> vd)
	{
	std::string rv = "(";
	bool start = true;
	for (auto d : vd)
		{
		if (!start)
			rv.append(";");
		rv.append(std::to_string(d));
		}
	rv.append(")");
	return rv;
	};

std::vector<int> TSection::IntVector(const std::string & n) const
	{
	TTokenSequence v = Value(n);
	return ParseVector<int>(v.begin(), v.end(), n).first;
	}

std::vector<double> TSection::RealVector(const std::string & n) const
	{
	TTokenSequence v = Value(n);
	return ParseVector<double>(v.begin(), v.end(), n).first;
	}

void TSection::AddAllowedMultiValue(const std::string & name)
	{
	multiValues_[name] = std::vector<TMultiValue>();
	}

bool TSection::ContainsMultiValue(const std::string & name) const
	{
	return multiValues_.find(name) != multiValues_.end();
	}

std::vector<TSection::TMultiValue> TSection::MultiValue(const std::string & n) const
	{
	auto it = multiValues_.find(n);
	bool found = (it != multiValues_.end());
	if (!found )
		throw std::runtime_error("TSection::MultiValue: unknown value " + n + " in section " + name_);
	return it->second;
	}

void TSection::SetMultiValue(const std::string & name, const TMultiValue & val)
	{
	auto it = multiValues_.find(name);
	bool found = (it != multiValues_.end());
	if (!found && !unknownValuesAllowed_)
		throw std::runtime_error("TSection::SetMultiValue: unknown value "+name+" in section "+name_);
	if (!found)
		std::tie(it,found) = multiValues_.insert({ name, std::vector<TMultiValue>() });
	it->second.push_back(val);
	}

void TSection::AddAllowedCommand(const std::string & name)
	{
	allowedCommands_.insert(name);
	}

bool TSection::ContainsCommand(const std::string & name) const
	{
	return allowedCommands_.find(name) != allowedCommands_.end();
	}

const TSection::TCommandList & TSection::CommandList() const
	{
	return commandList_;
	}


void TSection::ParseLine(const TTokenSequence& t, size_t currentLine, const std::string& fn, const TSection& variables)
	{
	TTokenSequence ts = t;
	for (const auto& v : variables.values_)
		ts = ReplaceVariable(ts, v.first, v.second);
	ts = ReplaceExpressions(ts);
	auto AssertToken = [currentLine, fn, this](Token expected, TokenNS::TSIt found)
		{
		if (expected != std::get<Token>(*found))
			throw std::runtime_error("TSection::ParseLine: expected " + TokenNS::TokenToString(expected) +
			", found " + TokenNS::TokenToString(std::get<Token>(*found)) + " in line # " + std::to_string(currentLine)
			+ " at position " + std::to_string(std::get<size_t>(*found)) + " in section " + this->name_
			+ " of file " + fn);
		};
	TokenNS::TSIt cur = ts.begin();
	TokenNS::TSIt end = ts.end();
	AssertToken(Token::identifier, cur);
	std::string n = std::get<std::string>(std::get<TTokenValue>(*cur));

	auto AssertNumToken = [&currentLine, &n, this](TokenNS::TSIt found)
		{
		if ((Token::integer != std::get<Token>(*found)) && (Token::real != std::get<Token>(*found)))
			throw std::runtime_error("TSection::ParseLine: expected real or integer"
			", found " + TokenNS::TokenToString(std::get<Token>(*found)) + " for entry " + n +
			" at position " + std::to_string(std::get<size_t>(*found)) + "of line "
			+ std::to_string(currentLine) + " in section " + this->name_);
		};

	cur = AdvanceToken(cur, end);
	if (std::get<Token>(*cur) == Token::equals) // value assignment
		{
		cur = AdvanceToken(cur, end);
		if (ContainsValue(n) || unknownValuesAllowed_)
			SetValue(n, cur, end);
		else
			throw std::runtime_error("TSection::ParseLine: unknown value name " + n + "in line # "
			+ std::to_string(currentLine) + " in section " + this->name_
			+ " of file " + fn);
		return;
		}
	if (std::get<Token>(*cur) == Token::lBrace) // multivalue assignment
		{
		if (ContainsMultiValue(n) || unknownValuesAllowed_)
			{
			std::vector<int> iv;
			std::tie(iv, cur) = ParseVector<int>(cur, end, n);
			AssertToken(Token::equals, cur);
			cur = AdvanceToken(cur, end);
//			AssertNumToken(cur);
//			bool ok;
//			TToken rhs = MakeRealIf(*cur, ok); // cannot fail
//			double drhs = std::get<double>(std::get<TTokenValue>(rhs));
//			SetMultiValue(n, {iv,drhs});
			SetMultiValue(n, {iv,TTokenSequence(cur,end)});
			}
		else
			throw std::runtime_error("TSection::ParseLine: unknown value name " + n + "in line # "
			+ std::to_string(currentLine) + " in section " + this->name_
			+ " of file " + fn);
		return;
		}
	if (std::get<Token>(*cur) == Token::lParen) // command
		{
		if (ContainsCommand(n) || unknownValuesAllowed_)
			{
			TCommandArgs args;
			cur = AdvanceToken(cur, end);
			while (!IsOneOf(*cur, { Token::rParen, Token::endofline }))
				{
				TTokenSequence arg;
				while (!IsOneOf(*cur, { Token::comma, Token::rParen, Token::endofline }))
					{
					arg.push_back(*cur);
					cur = AdvanceToken(cur, end);
					}
				args.push_back(arg);
				if (std::get<Token>(*cur) == Token::comma)
					{
					cur = AdvanceToken(cur, end);
					continue;
					}
				break;
				}
			AssertToken(Token::rParen, cur);
			cur = AdvanceToken(cur, end);
			commandList_.push_back(std::make_pair(n, args));
			}
		else
			{
			throw std::runtime_error("TSection::ParseLine: unknown command name " + n + "in line # "
				+ std::to_string(currentLine) + " in section " + this->name_
				+ " of file " + fn);
			}
		return;
		}
	throw std::runtime_error("TSection::ParseLine: syntax error after " + n + "in line # "
		+ std::to_string(currentLine) + " in section " + this->name_
		+ " of file " + fn);
	};



// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %% TConfiguration


TConfiguration::TConfiguration()
	: unknownSectionsAllowed_(false) 
	{
	AddSection(TSection("Variables", true, true));
	}

void TConfiguration::ParseCfgFile(CStrRef fn)
	{
	// std::chrono::system_clock cl;
	cfgContent_ << "% "
		// << cl.now() 
		<< " \n"
		<< "% called with configuration file " << fn
		<< "\n% which has the following content:" << std::endl;
	{
	TLineFile f(fn);
	while (!f.Eof())
		{
		f.NextLine();
		cfgContent_ << f.Buf() << std::endl;
		}
	}
	TLineFile f(fn);
	while (!f.Eof())
		{
		f.NextLine();
		TTokenSequence t = Tokenize(f.Buf()); // t contains at least endofline -- cannot be empty
		TokenNS::TSIt cur = t.begin();
		TokenNS::TSIt end = t.end();
		if (t.empty())
			throw std::runtime_error("ParseCfgFile: Token list empty when reading " + fn
			+ " at line no. " + std::to_string(f.LineNo()) + ":this cannot happen");
		if (std::get<Token>(*cur) == Token::endofline) // comment line
			{
			continue;
			}
		if (std::get<Token>(*cur) != Token::lBracket)
			throw std::runtime_error("ParseCfgFile: Expect section start ('[') when reading " + fn
			+ " at line no. " + std::to_string(f.LineNo()) + "found " + f.Buf());
		cur = AdvanceToken(cur, end);
		if (std::get<Token>(*cur) != Token::identifier)
			throw std::runtime_error("ParseCfgFile: Expect section name when reading " + fn
			+ " at line no. " + std::to_string(f.LineNo()) + " found " + f.Buf());
		std::string sectionname = std::get<std::string>(std::get<TTokenValue>(*cur));
		cur = AdvanceToken(cur, t.end());
		if (std::get<Token>(*cur) != Token::rBracket)
			throw std::runtime_error("ParseCfgFile: Expect section end (']') when reading " + fn
			+ " at line no. " + std::to_string(f.LineNo()) + "found " + f.Buf());
		cur = AdvanceToken(cur, t.end());
		if (std::get<Token>(*cur) != Token::endofline)
			throw std::runtime_error("ParseCfgFile: Expect only section ('[xxx]') when reading " + fn
			+ " at line no. " + std::to_string(f.LineNo()) + "found " + f.Buf());
		bool contains = ContainsSection(sectionname);
		if (contains)
			{
			TSection& s = sections_.find(sectionname)->second;
			s.ReadContent(f, Section("Variables"));
			}
		else
			{
			if (UnknownSectionsAllowed())
				{
				TSection s(sectionname, true);
				s.ReadContent(f, Section("Variables"));
				sections_[sectionname] =  std::move(s);
				}
			else
				throw std::runtime_error("ParseCfgFile: Unknown section ["+sectionname+"] when reading " + fn
					+ " at line no. " + std::to_string(f.LineNo()));
			}
		}
	}

std::string TConfiguration::Content() const
	{
	return cfgContent_.str();
	}

bool TConfiguration::ContainsSection(CStrRef name) const
	{
	auto it = sections_.find(name);
	return (it != sections_.end());
	}

const TSection & TConfiguration::Section(CStrRef name) const
	{
	auto it = sections_.find(name);
	if (it != sections_.end())
		return it->second;
	else
		throw std::runtime_error("TConfiguration::Section: not found: "+name);
	}

void TConfiguration::AddSection(TSection&& s)
	{
	sections_[s.Name()] = std::move(s);
	}

bool TConfiguration::UnknownSectionsAllowed() const
	{
	return unknownSectionsAllowed_;
	}

bool TConfiguration::SetUnknownSectionsAllowed(bool yesno)
	{
	bool rv = unknownSectionsAllowed_;
	unknownSectionsAllowed_ = yesno;
	return rv;
	}

// fully specialized template functions
template<>
int TSection::GetSingleTokenValue<int, Token::integer>(const std::string& valueName) const
	{
	const TTokenSequence& v = Value(valueName);
	if (v.size() != 2)
		throw std::runtime_error("TSection::GetSingleTokenValue: token sequence too long for value " + valueName
		+ " in section " + name_);
	Token t = std::get<Token>(v.front());
	if (t == Token::real)
		{
		bool ok;
		TToken tt = MakeIntegerIf(v.front(), ok);
		if (!ok)
			throw std::runtime_error("TSection::GetSingleTokenValue: MakeInteger failed for value " + valueName
			+ " in section " + name_);
		return std::get<int>(std::get<TTokenValue>(tt));
		}
	if (t != Token::integer)
		throw std::runtime_error("TSection::GetSingleTokenValue: Wrong token type for value " + valueName
		+ ", expected integer, found " + TokenNS::TokenToString(t) + " in section " + name_);
	return std::get<int>(std::get<TTokenValue>(v.front()));
	}

template<>
double TSection::GetSingleTokenValue<double, Token::real>(const std::string& valueName) const
	{
	const TTokenSequence& v = Value(valueName);
	if (v.size() != 2)
		throw std::runtime_error("TSection::GetSingleTokenValue: token sequence too long for value " + valueName
		+ " in section " + name_);
	Token t = std::get<Token>(v.front());
	if (t == Token::integer)
		{
		bool ok;
		TToken tt = MakeRealIf(v.front(), ok);
		if (!ok)
			throw std::runtime_error("TSection::GetSingleTokenValue: MakeRealIf failed: this cannot happen");
		return std::get<double>(std::get<TTokenValue>(tt));
		}
	if (t != Token::real)
		throw std::runtime_error("TSection::GetSingleTokenValue: Wrong token type for value " + valueName
		+ ", expected real, found " + TokenNS::TokenToString(t) + " in section " + name_);
	return std::get<double>(std::get<TTokenValue>(v.front()));
	}



// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// test

class TSomeSectionYouExpectInYourProgramSection: public TSection
	{
	public:
		TSomeSectionYouExpectInYourProgramSection()
			: TSection("SomeSectionYouExpectInYourProgram") {};
		virtual void AddAllowedValues() 
			{
			values_.insert({ "topp",MakeEmptyTokenSequence() });
			values_.insert({ "irrational",MakeEmptyTokenSequence() });
			values_.insert({ "doSomething",MakeEmptyTokenSequence() });
			values_.insert({ "hasDefault42", TTokenSequence{TToken{Token::integer, TTokenValue{42}, 0}, TToken{ Token::endofline, TTokenValue(), 0 } } });
			};
	};

class TTestCfg : public TConfiguration
	{
	public:
		TTestCfg()
			: TConfiguration()
			{
			TSomeSectionYouExpectInYourProgramSection ds;
			ds.InitAllowed();
			AddSection(std::move(ds));
			};
	};


void TestConfiguration()
	{
	TTestCfg cfg;
	cfg.ParseCfgFile("Test.cfg");
	/*
	[Variables]
	hopp="topp"
	answerToQuestionOfLifeTheUniverseAndEverything=42
	flag=true
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	[SomeSectionYouExpectInYourProgram]
	topp=hopp % topp will be a string with value "topp"
	irrational=sqrt(2)
	doSomething=flag
	*/
	std::string topp = cfg.Section("SomeSectionYouExpectInYourProgram").String("topp");
	double sqrt2 = cfg.Section("SomeSectionYouExpectInYourProgram").Real("irrational");
	bool doIt = cfg.Section("SomeSectionYouExpectInYourProgram").Bool("doSomething");
	int _42 = cfg.Section("SomeSectionYouExpectInYourProgram").Int("hasDefault42");
	}
