#include"timer.h"

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
