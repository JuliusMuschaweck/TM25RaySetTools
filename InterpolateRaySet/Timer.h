#include<chrono>
class TTimer
	{
	public:
		using TClock = std::chrono::high_resolution_clock;
		using TP = std::chrono::time_point < TClock >;
		using TD = std::chrono::duration < double >;
		TTimer( );
		void tic( ); // start timer
		double toc( ); // stop timer, return elapsed seconds
		double elapsed( ) const;
	private:
		TP tic_;
		TP toc_;
	};

