#include <Sobol.h>
#include <fstream>
void TestSobol()
	{
	TSobol<2> sob;
	TSobol<6> sob6;
	// fails as it should: TSobol<7> sob7;
	std::ofstream f("sobol.txt");
	for (size_t i = 0; i < 1024; ++i)
		{
		auto v = sob();
		f << v[0] << '\t' << v[1] << '\n';
		}
	}