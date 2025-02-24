#pragma once
#include<array>
#include<stdexcept>

template<size_t dim>
class TSobol
	{
	public:
		TSobol();
		std::array<double,dim> operator()(); // returns Sobol numbers within ]0;1[
	private:
		static const size_t maxbit_ = 30;
		static const size_t maxdim_ = 6;
		static const std::array<int, maxdim_> mdeg_;
		size_t in_;
		std::array<size_t, maxdim_> ix_;
		std::array<size_t*, maxbit_> iu_;
		static const std::array<size_t,maxdim_> ip_;
		std::array<size_t, maxbit_ * maxdim_> iv_;
		static const std::array<size_t, maxbit_ * maxdim_> iv0_;
		const double fac_;
	};

void TestSobol();

// %%%%%%%% template definition %%%%%%%%%%%%%%%%

template<size_t dim>
const std::array<int, TSobol<dim>::maxdim_> TSobol<dim>::mdeg_ = {1,2,3,3,4,4};

template<size_t dim>
const std::array<size_t, TSobol<dim>::maxdim_> TSobol<dim>::ip_ = {0,1,1,2,1,4};

template<size_t dim>
const std::array<size_t, TSobol<dim>::maxbit_ * TSobol<dim>::maxdim_> 
	TSobol<dim>::iv0_ = {1,1,1,1,1,1,3,1,3,3,1,1,5,7,7,3,3,5,15,11,5,15,13,9};

template<size_t dim>
TSobol<dim>::TSobol()
		: fac_ (1.0 / (1<<maxbit_))
	{
	size_t dimtest[maxdim_+1-dim]; // static assert dim <= maxdim
	dimtest[0] = 0;
	int j,k,l;
	size_t i,ipp;
	iv_ = iv0_;
	for (k=0; k<maxdim_; ++k)
		ix_[k] = 0;
	in_ = 0;
	for (j=0,k=0; j<maxbit_; ++j,k+= maxdim_)
		iu_[j] = &iv_[k];
	for (k=0; k< maxdim_; ++k)
		{
		for (j=0; j<mdeg_[k]; ++j)
			iu_[j][k] <<= (maxbit_-1-j);
		for (j=mdeg_[k]; j<maxbit_; ++j)
			{
			ipp = ip_[k];
			i = iu_[j-mdeg_[k]][k];
			i ^= (i >> mdeg_[k]);
			for (l=mdeg_[k]-1; l>=1; --l)
				{
				if (ipp&1)
					i ^= iu_[j-l][k];
				ipp >>= 1;
				}
			iu_[j][k] = i;
			}
		}
	};

template<size_t dim>
std::array<double,dim> TSobol<dim>::operator()()
	{
	size_t im = in_++;
	int j,k;
	for (j=0; j < maxbit_; ++j)
		{
		if (!(im&1))
			break;
		im >>= 1;
		}
	if (j >= maxbit_)
		throw std::runtime_error("TSobol::operator()(): maxbit_ too small");
	im = j * maxdim_;
	std::array<double,dim> rv;
	for (k=0; k < dim; ++k)
		{
		ix_[k] ^= iv_[im+k];
		rv[k] = ix_[k] * fac_;
		}
	return rv;
	}