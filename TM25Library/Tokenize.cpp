#include "Tokenize.h"
#include <sstream>
#include <functional>
#include <map>
#include <algorithm>

TTokenSequence MakeEmptyTokenSequence()
	{
	TTokenSequence rv;
	rv.push_back({ Token::endofline, TTokenValue(), 0 });
	return rv;
	}

TToken MakeRealIf(const TToken & rhs, bool & ok)
	{
	ok = true;
	if (std::get<Token>(rhs) == Token::real)
		return rhs;
	else if (std::get<Token>(rhs) == Token::integer)
		{
		TTokenValue val = static_cast<double>(std::get<int>(std::get<TTokenValue>(rhs)));
		return {Token::real, val, std::get<size_t>(rhs) };
		}
	ok = false;
	return rhs;
	}

TToken MakeIntegerIf(const TToken & rhs, bool & ok)
	{
	ok = true;
	if (std::get<Token>(rhs) == Token::integer)
		return rhs;
	else if (std::get<Token>(rhs) == Token::real)
		{
		double d = std::get<double>(std::get<TTokenValue>(rhs));
		double test = std::round(d);
		bool fits =
			(d <= static_cast<double>(std::numeric_limits<int>::max()))
			&& (d >= static_cast<double>(std::numeric_limits<int>::min()));
		if ((d == test) && fits)
			{
			TTokenValue val = static_cast<int>(d);
			return { Token::integer, val, std::get<size_t>(rhs) };
			}
		}
	ok = false;
	return rhs;
	}

TTokenSequence Tokenize(const std::string& line, bool skipSpace, bool skipComment)
	{
	using StrIt = std::string::const_iterator;
	StrIt c  = line.begin();
	TTokenSequence rv;
	auto AddVal = [&rv,&line](Token t, TTokenValue&& v,  StrIt cc)
		{rv.push_back(std::make_tuple(t, std::move(v), cc - line.begin())); };
	auto Add = [&rv, &AddVal](Token t, StrIt cc)
		{AddVal(t, TTokenValue(), cc); };
	while (c != line.end())
		{
		switch (*c)
			{
				case '[': Add(Token::lBracket, c); ++c; break;
				case ']': Add(Token::rBracket, c); ++c; break;
				case '(': Add(Token::lParen, c);  ++c; break;
				case ')': Add(Token::rParen, c);  ++c; break;
				case '{': Add(Token::lBrace, c);  ++c; break;
				case '}': Add(Token::rBrace, c);  ++c; break;
				case '%': if (!skipComment) AddVal(Token::comment, std::string(c+1,line.end()), c);
					c = line.end();
					break;
				case ';': Add(Token::semicolon, c);  ++c; break;
				case ',': Add(Token::comma, c); ++c; break;
				//case '\'': Add(Token::quote, c);  ++c; break;
				case '"': 
					{
					StrIt c0 = c;
					while (++c != line.end() && (*c) != '"');
					if (c == line.end())
						{
						std::stringstream s;
						s << "Tokenize: string has no ending \" after starting at position "
							<< std::to_string(c0 - line.begin())
							<< " (" << *c0 << ") of " << line;
						throw std::runtime_error(s.str());
						}
					else
						++c;
					// ok since c-c0 >= 2 (at least two ++)
					AddVal(Token::string, std::string(c0+1, c-1), c0);
					}
					break;
				case '\'':
					{
					StrIt c0 = c;
					while (++c != line.end() && (*c) != '\'');
					if (c == line.end())
						{
						std::stringstream s;
						s << "Tokenize: string has no ending \' after starting at position "
							<< std::to_string(c0 - line.begin())
							<< " (" << *c0 << ") of " << line;
						throw std::runtime_error(s.str());
						}
					else
						++c;
					// ok since c-c0 >= 2 (at least two ++)
					AddVal(Token::string, std::string(c0+1, c-1), c0);
					}
					break;
				case '=': Add(Token::equals, c);  ++c; break;
				case '+': Add(Token::plus, c);  ++c; break;
				case '-': Add(Token::minus, c);  ++c; break;
				case '*': Add(Token::times, c);  ++c; break;
				case '/': Add(Token::divide, c);  ++c; break;
				default: // identifier, space, integer, real
					if (isalpha(*c) || (*c == '_')) // identifier
						{
						StrIt c0 = c;
						while (++c != line.end() && (isalnum(*c) || (*c == '_')));
						std::string ident(c0, c);
						if (ident == "true")
							AddVal(Token::boolean, true, c0);
						else if (ident == "false")
							AddVal(Token::boolean, false, c0);
						else
							AddVal(Token::identifier, ident, c0);
						}
					else if (isspace(*c)) // space
						{
						StrIt c0 = c;
						while (++c != line.end() && isspace(*c));
						if (!skipSpace) Add(Token::space, c0);
						}
					else if (isdigit(*c)) // integer or real
						{
						StrIt c0 = c;
						size_t ipos, rpos;
						int i = std::stoi(line.substr(c - line.begin()), &ipos);
						double d = std::stod(line.substr(c - line.begin()), &rpos);
						if (rpos == ipos)
							AddVal(Token::integer, i, c0);
						else
							AddVal(Token::real, d, c0);
						c += rpos;
						}
					else
						{
						std::stringstream s;
						s << "Tokenize: cannot process token at position "
							<< std::to_string(c - line.begin())
							<< " (" << *c << ") of " << line;
						throw std::runtime_error(s.str());
						}
			} // switch
		} // while
	Add(Token::endofline, c);
	return rv;
	}

TTokenSequence ReplaceVariable(const TTokenSequence & seq, const std::string & name, const TTokenSequence & value)
	{
	TTokenSequence rv;
	auto t = seq.begin();
	while (t != seq.end())
		{
		if (std::get<Token>(*t) == Token::identifier)
			{
			TTokenValue val = std::get<TTokenValue>(*t);
			std::string sval = std::get<std::string>(val);
			if (sval == name)
				{
				for (auto v : value)
					{
					std::get<size_t>(v) = std::get<size_t>(*t); // assign position in original string
					if (std::get<Token>(v) != Token::endofline)
						rv.push_back(v);
					}
				}
			else
				rv.push_back(*t);
			}
		else
			rv.push_back(*t);
		++t;
		}
	return rv;
	}

bool IsOneOf(const TToken & tok, std::vector<Token> toks)
	{
	return (std::find(toks.begin(), toks.end(), std::get<Token>(tok)) != toks.end());
	}

std::string TestTokenize()
	{
	std::string s1 = "[](){};, \"stringwith'quote'\" 'stringwith\"doublequotes\"' =+-*/ horst _horst H0 10 10. 10.0 1e+1 1e+a";
	std::string s2 = "+- true false %comment";
	std::string s3 = "%comment";
	std::string s4 = "\"noclosingdq";
	std::string s5 = "\'noclosingq";
	std::string s6 = "error #";

	std::string somevar = "'some string' horst";
	std::string s7 = "somevar 10 somevar 20";
	std::string s8 = "-1,-sin(0);-2+3.0*sin(pi()/6); atan2(1,0)/pi(); 4.0-1.0, if(true,1,0); min(1,4,-exp(1)), max(1,4.0)";
	try
		{
		
		auto t = Tokenize(s1);
		t = Tokenize(s2);
		t = Tokenize(s3);
		//t = Tokenize(s4);
		//t = Tokenize(s5);
		//t = Tokenize(s6);
		auto repl = Tokenize(somevar);
		t = Tokenize(s7);
		t = ReplaceVariable(t, "somevar", repl);
		t = Tokenize(s8);
		t = ReplaceExpressions(t);
		}
	catch (std::runtime_error e)
		{
		return e.what();
		}
	return "ok";
	}

namespace TokenNS{
	using TSIt = TTokenSequence::const_iterator;

	std::tuple<TToken, TSIt> Expression(TSIt begin, TSIt end);
	bool IsFunction(const std::string& ident);
	}

TTokenSequence ReplaceExpressions(const TTokenSequence & seq)
	{
	TTokenSequence rv;
	auto t = seq.begin();
	while (t != seq.end())
		{
		Token tt = std::get<0>(*t);
		if ((tt == Token::plus) || (tt == Token::minus) || (tt == Token::real) 
			|| (tt == Token::integer) || (tt == Token::lParen)
			|| (tt == Token::identifier 
				&& TokenNS::IsFunction(std::get<std::string>(std::get<TTokenValue>(*t))))
			) 
			{
			auto res = TokenNS::Expression(t, seq.end());
			rv.push_back(std::get<0>(res));
			t = std::get<1>(res);
			}
		else
			{
			rv.push_back(*t);
			++t;
			}
		}
	return rv;
	}

TokenNS::TSIt AdvanceToken(TokenNS::TSIt begin, TokenNS::TSIt end)
	{
	if (begin == end)
		throw std::runtime_error("AdvanceToken: begin == end: this cannot happen");
	if (std::get<Token>(*begin) == Token::endofline)
		return begin;
	++begin;
	if (begin == end)
		throw TTokenError("AdvanceToken: Unexpected end of expression", --begin);
	return begin;
	}


namespace TokenNS{


	Token Which(TSIt it) { return std::get<Token>(*it); }
	Token Which(std::tuple<TToken, TSIt> res) { return std::get<Token>(std::get<TToken>(res)); }

	std::string TokenToString(Token t)
		{
		switch (t)
			{
				case (Token::lBracket): return "[";
				case (Token::rBracket): return "]";
				case (Token::lParen): return "(";
				case (Token::rParen): return ")";
				case (Token::lBrace): return "{";
				case (Token::rBrace): return "}";
				case (Token::comment): return "comment (%)";
				case (Token::semicolon): return ";";
				case (Token::comma): return ",";
				case (Token::string): return "string";
				case (Token::equals): return "=";
				case (Token::plus): return "+";
				case (Token::minus): return "-";
				case (Token::times): return "*";
				case (Token::divide): return "/";
				case (Token::endofline): return "end of line";
				case (Token::identifier): return "identifier";
				case (Token::space): return "space";
				case (Token::integer): return "integer";
				case (Token::real): return "real";
				case (Token::boolean): return "boolean (true or false)";
				default: throw std::runtime_error("TokenToString: unknown token: this cannot happen");
			}
		}

	std::string TTokenToString(TToken t)
		{
		std::stringstream rv;
		rv << TokenToString(std::get<Token>(t)) << ": ";
		TTokenValue v = std::get<TTokenValue>(t);
		// using TTokenValue = std::variant<std::monostate, std::string, int, double, bool>;
		auto ps = std::get_if<std::string>(&v);
		if (ps) rv << *ps;
		auto pi = std::get_if<int>(&v);
		if (pi) rv << *pi;
		auto pb = std::get_if<bool>(&v);
		if (pb) rv << ( (*pb)? "true" : "false");
		rv << " at position " << std::get<size_t>(t);
		return rv.str();
		}


	void AssertToken(TSIt it, Token tt, const std::string& where)
		{
		if (Which(it) != tt)
			throw TTokenError(where + ": expected " + TokenToString(tt)
			+ ", found " + TokenToString(Which(it)), it);
		}

	// return values are tuple of result (integer or real, value) and iterator where to continue
	std::tuple<TToken, TSIt> Expression(TSIt begin, TSIt end);
	std::tuple<TToken, TSIt> Expression2(TSIt begin, TSIt end, bool minusPrefix = false);
	std::tuple<TToken, TSIt> Term(TSIt begin, TSIt end);
	std::tuple<TToken, TSIt> Factor(TSIt begin, TSIt end);
	std::tuple<TToken, TSIt> Function(TSIt begin, TSIt end);

	std::tuple<TToken, TSIt> Function0(TSIt begin, TSIt end);
	std::tuple<TToken, TSIt> Function1(TSIt begin, TSIt end);
	std::tuple<TToken, TSIt> Function2(TSIt begin, TSIt end);
	std::tuple<TToken, TSIt> Function3(TSIt begin, TSIt end);
	std::tuple<TToken, TSIt> Function_ge2(TSIt begin, TSIt end);

	double GetExpressionResult(const std::tuple<TToken, TSIt>& rhs)
		{
		const TTokenValue& v = std::get<TTokenValue>(std::get<TToken>(rhs));
		if (std::holds_alternative<int>(v))
			return static_cast<double>(std::get<int>(v));
		if (std::holds_alternative<double>(v))
			return std::get<double>(v);
		throw TTokenError("GetExpressionResult: not int or double: this cannot happen", 
			std::get<TSIt>(rhs));
		}

	/*std::tuple<bool, int> GetExpressionIntResultIf(const std::tuple<TToken, TSIt>& rhs)
		{
		const TTokenValue& v = std::get<TTokenValue>(std::get<TToken>(rhs));
		if (std::holds_alternative<int>(v))
			return { true, std::get<int>(v) };
		if (std::holds_alternative<double>(v))
			{
			double d = std::get<double>(v);
			double test = std::round(d);
			bool fits =
				(d <= static_cast<double>(std::numeric_limits<int>::max()))
				&& (d >= static_cast<double>(std::numeric_limits<int>::min()));
			if ((d == test) && fits)
				return { true, static_cast<int>(d) };
			else
				return { false, 0 };
			}
		throw TTokenError("GetExpressionIntResultIf: not int or double: this cannot happen", std::get<size_t>(v));
		}

	int GetExpressionIntResult(const std::tuple<TToken, TSIt>& rhs)
		{
		const TTokenValue& v = std::get<TTokenValue>(std::get<TToken>(rhs));
		if (std::holds_alternative<int>(v))
			return std::get<int>(v);
		throw TTokenError("GetExpressionIntResult: not int: this cannot happen", std::get<size_t>(v));
		}*/


	// nullary functions
	using TF0 = std::function<double()>;

	class TF0Map
		{
		public:
			bool IsTF0(const std::string& ident) const {return (fmap_.find(ident) != fmap_.end());};
			TF0 Func(const std::string& ident) const
				{
				auto irv = fmap_.find(ident);
				if (irv == fmap_.end())
					throw std::runtime_error("TF0Map::Func: function not found: " + ident + ", this cannot happen");
				return (*irv).second;
				};
		private:
			TF0Map();
			friend const TF0Map& F0Map();
			std::map<std::string, TF0> fmap_;
		};		

	TF0Map::TF0Map()
		{
		fmap_["nan"] = TF0([]() {return NAN; });
		fmap_["pi"] = TF0([]() {return abs(atan2(0,-1)); });
		}

	const TF0Map& F0Map()
		{
		static TF0Map f0m;
		return f0m;
		}


	// unary functions
	using TF1 = std::function<double(double)>;

	class TF1Map 
		{
		public:
			bool IsTF1(const std::string& ident) const {return (fmap_.find(ident) != fmap_.end());};
			TF1 Func(const std::string& ident) const
				{
				auto irv = fmap_.find(ident);
				if (irv == fmap_.end())
					throw std::runtime_error("TF1Map::Func: function not found: "+ident+", this cannot happen");
				return (*irv).second;
				};
		private:
			TF1Map();
			friend const TF1Map& F1Map();
			std::map<std::string, TF1> fmap_;
		};

	TF1Map::TF1Map()
		{
		fmap_["abs"] = TF1([](double x) {return abs(x); });
		fmap_["exp"] = TF1([](double x) {return exp(x); });
		fmap_["log"] = TF1([](double x) {return log(x); });
		fmap_["log10"] = TF1([](double x) {return log10(x); });
		fmap_["sqrt"] = TF1([](double x) {return sqrt(x); });
		fmap_["sin"] = TF1([](double x) {return sin(x); });
		fmap_["asin"] = TF1([](double x) {return asin(x); });
		fmap_["cos"] = TF1([](double x) {return cos(x); });
		fmap_["acos"] = TF1([](double x) {return acos(x); });
		fmap_["tan"] = TF1([](double x) {return tan(x); });
		fmap_["atan"] = TF1([](double x) {return atan(x); });
		fmap_["ceil"] = TF1([](double x) {return ceil(x); });
		fmap_["floor"] = TF1([](double x) {return floor(x); });
		fmap_["trunc"] = TF1([](double x) {return trunc(x); });
		fmap_["round"] = TF1([](double x) {return round(x); });
		}	

	const TF1Map& F1Map() 
		{
		static TF1Map f1m;
		return f1m;
		}

	// binary functions
	using TF2 = std::function<double(double, double)>;

	class TF2Map
		{
		public:
			bool IsTF2(const std::string& ident) const{return (fmap_.find(ident) != fmap_.end());};
			TF2 Func(const std::string& ident) const
				{
				auto irv = fmap_.find(ident);
				if (irv == fmap_.end())
					throw std::runtime_error("TF2Map::Func: function not found: " + ident + ", this cannot happen");
				return (*irv).second;
				};
		private:
			TF2Map();
			friend const TF2Map& F2Map();
			std::map<std::string, TF2> fmap_;
		};		

	TF2Map::TF2Map()
		{
		fmap_["atan2"] = TF2([](double lhs, double rhs) {return atan2(lhs, rhs); });
		fmap_["mod"] = TF2([](double lhs, double rhs) {return fmod(lhs, rhs); });
		fmap_["pow"] = TF2([](double lhs, double rhs) {return pow(lhs, rhs); });
		};

	const TF2Map& F2Map()
		{
		static TF2Map f2m;
		return f2m;
		};

	// <function>

	bool IsFunction(const std::string& ident)
		{
		bool rv = F0Map().IsTF0(ident) || F1Map().IsTF1(ident) || F2Map().IsTF2(ident)
			|| (ident == "min") || (ident == "max") 
			|| (ident == "if");
		return rv;
		}

	std::tuple<TToken, TSIt> Function(TSIt begin, TSIt end)
		{
		if (Which(begin) != Token::identifier)
			throw std::runtime_error("TokenNS::Function: no identifier: this cannot happen");
		const std::string& ident = std::get<std::string>(std::get<TTokenValue>(*begin));
		if (F0Map().IsTF0(ident))
			return Function0(begin, end);
		if (F1Map().IsTF1(ident))
			return Function1(begin, end);
		if (F2Map().IsTF2(ident))
			return Function2(begin, end);
		if ((ident == "min") || (ident == "max"))
			return Function_ge2(begin, end);
		if ((ident == "if"))
			return Function3(begin, end);
		throw std::runtime_error("TokenNS::Function: ínvalid function: " + ident + ": this cannot happen");
		}

	// <function0>
	std::tuple<TToken, TSIt> Function0(TSIt begin, TSIt end)
		{
		const std::string& ident = std::get<std::string>(std::get<TTokenValue>(*begin));
		std::function<double()> f = F0Map().Func(ident);
		TSIt cur = AdvanceToken(begin, end);
		AssertToken(cur, Token::lParen, "Function0");
		cur = AdvanceToken(cur, end);
		AssertToken(cur, Token::rParen, "Function0");
		cur = AdvanceToken(cur, end);
		TToken rv = std::make_tuple(Token::real, TTokenValue(f()), std::get<size_t>(*begin));
		return std::make_tuple(rv, cur);
		}

	// <function1>
	std::tuple<TToken, TSIt> Function1(TSIt begin, TSIt end)
		{
		const std::string& ident = std::get<std::string>(std::get<TTokenValue>(*begin));
		std::function<double(double)> f = F1Map().Func(ident);
		TSIt cur = AdvanceToken(begin, end);
		AssertToken(cur, Token::lParen, "Function1");
		cur = AdvanceToken(cur, end);
		std::tuple<TToken, TSIt> iarg = Expression(cur, end);
		double arg = GetExpressionResult(iarg);
		cur = std::get<TSIt>(iarg);
		AssertToken(cur, Token::rParen, "Function1");
		cur = AdvanceToken(cur, end);
		TToken rv = std::make_tuple(Token::real, TTokenValue(f(arg)), std::get<size_t>(*begin));
		return std::make_tuple(rv, cur);
		}

	// <function2>
	std::tuple<TToken, TSIt> Function2(TSIt begin, TSIt end)
		{
		const std::string& ident = std::get<std::string>(std::get<TTokenValue>(*begin));
		std::function<double(double, double)> f = F2Map().Func(ident);
		TSIt cur = AdvanceToken(begin, end);
		AssertToken(cur, Token::lParen, "Function2");
		cur = AdvanceToken(cur, end);
		std::tuple<TToken, TSIt> iarg = Expression(cur, end);
		double arg1 = GetExpressionResult(iarg);
		cur = std::get<TSIt>(iarg);
		AssertToken(cur, Token::comma, "Function2");
		cur = AdvanceToken(cur, end);
		iarg = Expression(cur, end);
		double arg2 = GetExpressionResult(iarg);
		cur = std::get<TSIt>(iarg);
		AssertToken(cur, Token::rParen, "Function2");
		cur = AdvanceToken(cur, end);
		TToken rv = std::make_tuple(Token::real, TTokenValue(f(arg1, arg2)), std::get<size_t>(*begin));
		return std::make_tuple(rv, cur);
		}

	// <function3>
	std::tuple<TToken, TSIt> Function3(TSIt begin, TSIt end)
		{
		const std::string& ident = std::get<std::string>(std::get<TTokenValue>(*begin));
		if (!(ident == "if"))
			throw TTokenError("Function3: expected if, found " + ident + ", this cannot happen", begin);
		TSIt cur = AdvanceToken(begin, end);
		AssertToken(cur, Token::lParen, "Function3");
		cur = AdvanceToken(cur, end);
		AssertToken(cur, Token::boolean, "Function3");
		bool yesno = std::get<bool>(std::get<TTokenValue>(*cur));
		cur = AdvanceToken(cur, end);
		AssertToken(cur, Token::comma, "Function3");
		cur = AdvanceToken(cur, end);
		std::tuple<TToken, TSIt> iarg = Expression(cur, end);
		double arg1 = GetExpressionResult(iarg);
		cur = std::get<TSIt>(iarg);
		AssertToken(cur, Token::comma, "Function3");
		cur = AdvanceToken(cur, end);
		iarg = Expression(cur, end);
		double arg2 = GetExpressionResult(iarg);
		cur = std::get<TSIt>(iarg);
		AssertToken(cur, Token::rParen, "Function3");
		cur = AdvanceToken(cur, end);
		double res = yesno ? arg1 : arg2;
		TToken rv = std::make_tuple(Token::real, TTokenValue(res), std::get<size_t>(*begin));
		return std::make_tuple(rv, cur);
		}

	// <function_ge2>
	std::tuple<TToken, TSIt> Function_ge2(TSIt begin, TSIt end)
		{
		const std::string& ident = std::get<std::string>(std::get<TTokenValue>(*begin));
		std::vector<double> args;
		TSIt cur = AdvanceToken(begin, end);
		AssertToken(cur, Token::lParen, "Function_ge2");
		cur = AdvanceToken(cur, end);
		std::tuple<TToken, TSIt> iarg = Expression(cur, end);
		args.push_back(GetExpressionResult(iarg));
		cur = std::get<TSIt>(iarg);
		while (std::get<Token>(*cur) == Token::comma)
			{
			cur = AdvanceToken(cur, end);
			iarg = Expression(cur, end);
			args.push_back(GetExpressionResult(iarg));
			cur = std::get<TSIt>(iarg);
			}
		AssertToken(cur, Token::rParen, "Function_ge2");
		cur = AdvanceToken(cur, end);
		double res;
		if (ident == "min")
			res = *(std::min_element(args.begin(), args.end()));
		else if (ident == "max")
			res = *(std::max_element(args.begin(), args.end()));
		else
			throw TTokenError("Function3: expected min or max, found " + ident + ", this cannot happen", begin);
		TToken rv = std::make_tuple(Token::real, TTokenValue(res), std::get<size_t>(*begin));
		return std::make_tuple(rv, cur);
		}


	std::tuple<TToken, TSIt> Factor(TSIt begin, TSIt end)
		{
		switch (std::get<Token>(*begin))
			{
				case Token::real: 
					return std::make_tuple(*begin, AdvanceToken(begin, end));
				case Token::integer:
					return std::make_tuple(*begin, AdvanceToken(begin, end));
				case Token::lParen:
					{
					TSIt cur = AdvanceToken(begin, end);
					auto rv = Expression(cur, end);
					cur = std::get<TSIt>(rv);
					AssertToken(cur, Token::rParen, "Function3");
					cur = AdvanceToken(cur, end);
					std::get<TSIt>(rv) = cur;
					return rv;
					}
				case Token::identifier:
					{
					const std::string& ident = std::get<std::string>(std::get<TTokenValue>(*begin));
					if (IsFunction(ident))
						return Function(begin, end);
					else
						throw TTokenError("Factor: expected function, found "+ident, begin);
					}
				default:
					throw TTokenError("Factor: expected real, integer, ( or function", begin);
			}
		}

	std::tuple<TToken, TSIt> Term(TSIt begin, TSIt end)
		{
		std::tuple<TToken, TSIt> lhs = Factor(begin, end);
		TSIt cur = std::get<TSIt>(lhs);
		if (Which(cur) == Token::times || Which(cur) == Token::divide)
			{
			bool times = (Which(cur) == Token::times);
			cur = AdvanceToken(cur, end);
			std::tuple<TToken, TSIt> rhs = Term(cur, end);
			size_t pos = std::get<size_t>(std::get<TToken>(rhs));
			double drv;
			double dlhs = GetExpressionResult(lhs);
			double drhs = GetExpressionResult(rhs);
			if (times)
				drv = dlhs * drhs;
			else // divide
				drv = dlhs / drhs;
			return { {Token::real, TTokenValue(drv), pos}, std::get<TSIt>(rhs) };
			} 
		else
			return lhs;
		}

	std::tuple<TToken, TSIt> Expression2(TSIt begin, TSIt end, bool minusPrefix)
		{
		std::tuple<TToken, TSIt> lhs = Term(begin, end);
		if (minusPrefix)
			{
			if (Which(lhs) == Token::real)
				std::get<double>(std::get<TTokenValue>(std::get<TToken>(lhs))) *= (-1.0);
			else
				std::get<int>(std::get<TTokenValue>(std::get<TToken>(lhs))) *= (-1);
			}
		TSIt cur = std::get<TSIt>(lhs);
		if (Which(cur) == Token::plus || Which(cur) == Token::minus)
			{
			bool plus = (Which(cur) == Token::plus);
			cur = AdvanceToken(cur, end);
			std::tuple<TToken, TSIt> rhs = Expression2(cur, end);
			size_t pos = std::get<size_t>(std::get<TToken>(rhs));
			double drv;
			double dlhs = GetExpressionResult(lhs);
			double drhs = GetExpressionResult(rhs);
			if (plus)
				drv = dlhs + drhs;
			else // minus
				drv = dlhs - drhs;
			return { {Token::real, TTokenValue(drv), pos}, std::get<TSIt>(rhs) };
			}
		else
			return lhs;
		}

	std::tuple<TToken, TSIt> Expression(TSIt begin, TSIt end)
		{
		if (Which(begin) == Token::plus)
			return Expression2(AdvanceToken(begin,end), end);
		if (Which(begin) == Token::minus)
			return Expression2(AdvanceToken(begin, end), end, true);
		return Expression2(begin, end);
		}


	}

	/**/