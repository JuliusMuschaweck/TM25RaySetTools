#ifndef __Tokenize_h
#define __Tokenize_h

#include <variant>
#include <string>
#include<tuple>
#include<vector>

enum class Token
	{
	lBracket,
	rBracket,
	lParen,
	rParen,
	lBrace,
	rBrace,
	comment, // <comment> ::= '%' {any char} until end of line
	semicolon,
	comma,
	string, // <string> ::= '"' {any char except '"'} '"' | ''' {any char except '''} '''
	equals,
	plus,
	minus,
	times,
	divide,
	endofline,
	identifier, // ('_' | alpha) { '_' | alnum}
	space,
	integer, // see stoi
	real, // see stof
	boolean // <boolean> ::= "true" | "false"
	};

using TTokenValue = std::variant<std::monostate, std::string, int, double, bool>;

using TToken = std::tuple<Token, TTokenValue, size_t>; // size_t = position in string

using TTokenSequence = std::vector<TToken>;
TTokenSequence MakeEmptyTokenSequence(); // contains only Token::endofline

template<Token tok, typename TValue>
TTokenSequence MakeDefaultValueTokenSequence(TValue val); 


TToken MakeRealIf(const TToken& rhs, bool& ok); // ok if token is real or integer
TToken MakeIntegerIf(const TToken& rhs, bool& ok); // ok if token is integer or (real and has integer value within range)

// Tokenize: token sequence always ends with Token::endofline. May throw std::runtime_error
TTokenSequence Tokenize(const std::string& line, bool skipSpace = true, bool skipComment = true);


// ReplaceVariable: any Token::identifier in "seq" which matches "name" is replaced by "value"
//TODO:  open problem: string position of TToken is then unclear
TTokenSequence ReplaceVariable(const TTokenSequence& seq,
	const std::string& name, const TTokenSequence& value);


bool IsOneOf(const TToken& tok, std::vector<Token> toks);

namespace TokenNS
	{
	using TSIt = TTokenSequence::const_iterator;
	std::string TokenToString(Token t);
	std::string TTokenToString(TToken t);
	}

class TTokenError : public std::runtime_error
	{
	public:
		TTokenError(const std::string& what, TokenNS::TSIt at)
			: std::runtime_error(what), at_(at) {};
		TokenNS::TSIt at_;
	};

// increments iterator except if begin is Token::endofline: then does nothing
// Thus, return value always points to a valid entry as long as token sequence
// ends with Token::endofline.
// May throw runtime_error if begin==end and TTokenError if ++begin==end
TokenNS::TSIt AdvanceToken(TokenNS::TSIt begin, TokenNS::TSIt end);

/*
Numeric expressions: 

<expression> ::= ['+'|'-'] <expression2>
<expression2>::= <term> | <term> '+' <expression2> |   <term> "-" <expression2>
<term>       ::= <factor> | <factor> '*' <term>  |   <factor> "/" <term>
<factor>     ::= <real> | <int> | '(' <expression> ')' | <function>
<function>   ::= <function0> | <function1> | <function2> | <function3> | <function_ge2>
<function0>  ::= "nan" | "pi" '(' ')'
<function1>  ::= "abs" | "exp" | "log" | "log10"  | "sqrt"
					| "sin" | "asin" | "cos" | "acos" | "tan" | "atan" | "ceil" | "floor"
					| "trunc" | "round"
					'(' <expression> ')'
<function2>  ::= "atan2" | "mod" | "pow" '(' <expression> ',' <expression> ')'
<function3>  ::= "if" '(' ("true"|"false") ',' <expression> ',' <expression> ')'
<function_ge2>::= "min" | "max" '(' <expression> ',' <expression> {',' expression} ')'

expressions start with a real, integer, '+', '-', '(' or an identifier which is a valid function
*/

// ReplaceExpressions: any expression in "seq" replaced by a single token with its result.
TTokenSequence ReplaceExpressions(const TTokenSequence& seq);

std::string TestTokenize();


// template definition
template<Token tok, typename TValue>
TTokenSequence MakeDefaultValueTokenSequence(TValue val)
	{
	return TTokenSequence{ TToken{tok, TTokenValue{val}, 0}, TToken{ Token::endofline, TTokenValue(), 0 } };
	}

// include guard
#endif