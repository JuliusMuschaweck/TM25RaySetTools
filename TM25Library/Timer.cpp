#define _CRT_SECURE_NO_WARNINGS
#include"timer.h"
#include<sstream>
#include <iomanip> // for put_time
#include <time.h>
TTimer::TTimer( )
	: tic_( TClock::now( ) )
	{
	toc_ = tic_;
	};

void TTimer::tic( )
	{
	tic_ = TClock::now( );
	};

double TTimer::toc( )
	{
	toc_ = TClock::now( );
	return elapsed( );
	};

double TTimer::elapsed( ) const
	{
	TD rv = toc_ - tic_;
	return rv.count( );
	};

std::string currentISO8601TimeUTC() // see https://techoverflow.net/2018/03/30/iso8601-utc-time-as-stdstring-using-c11-chrono/
	{
	auto now = std::chrono::system_clock::now();
	auto itt = std::chrono::system_clock::to_time_t(now);
	std::ostringstream ss;
	ss << std::put_time(gmtime(&itt), "%FT%TZ");
	return ss.str();
	}
