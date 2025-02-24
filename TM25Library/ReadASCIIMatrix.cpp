#include "ReadASCIIMatrix.h"

#include<ranges>
#include<span>
#include<string_view>
#include "ReadFile.h"
#include "TM25Util.h"

namespace rng = std::ranges;

namespace TM25
	{

	// **********************************************************************************
	void TNumberList::Read(const std::string& fn, const std::string& delims)
		{
		TReadFile f(fn);
		Read(f);
		}
	void TNumberList::Read(TReadFile& f, const std::string& delims)
		{
		Clear();
		size_t ln = 0;
		std::string numstart{ "+-0123456789" };
		while (!f.AtEof())
			{
			std::string l = f.ReadLine();
			++ln;
			auto tokens = Tokenize(TrimWhiteSpace(l), delims);
			if (!tokens.empty())
				{
				auto test = tokens.front();
				if (numstart.find(test.front()) != std::string::npos) // else comment line
					{
					size_t num = 0;
					for (auto tok : tokens)
						{
						double d;
						const char* first = tok.data();
						const char* last = first + tok.size();
						std::from_chars_result res = std::from_chars(first, last, d);
						if (res.ptr == last)
							{
							++num;
							num_.push_back(d);
							}
						else
							{
							throw std::invalid_argument(std::format("TNumberList::Read: invalid number in line {}: {} ", ln, tok));
							}
						}
					num_per_line_.push_back(num);
					}
				}
			} // while
		// now num_ and num_per_line_ are filled and won't be mutated. lineStarts contains the indices of num_entries
		size_t start = 0;
		for (auto i : num_per_line_)
			{
			lines_.push_back(std::span<double>(num_.begin() + start, i));
			start += i;
			}
		}

	const std::vector<std::span<double>>& TNumberList::Lines() const // the list of number lines
		{
		return lines_;
		}

	// low level access functions
	const std::vector<double>& TNumberList::Num() const
		{
		return num_;
		};
	const std::vector<size_t>& TNumberList::Num_per_line() const
		{
		return num_per_line_;
		}

	void TNumberList::Clear()
		{
		num_.clear();
		num_per_line_.clear();
		lines_.clear();
		}


	// **************************************************************************************
	void TVectorMatrixList::Read(TReadFile& f, const std::string& delims)
		{
		TNumberList nl;
		nl.Read(f, delims);
		CreateFromNumberList(nl);
		}

	void TVectorMatrixList::Read(const std::string& fn, const std::string& delims)
		{
		TReadFile f(fn);
		Read(f, delims);
		}
	
	size_t To_size_t(double d)
		{
		if (d > 0 
			&& d < std::numeric_limits<size_t>::max()
			&& d == std::trunc(d))
			{
			return static_cast<size_t>(d);
			}
		else
			throw std::runtime_error(std::format("ReadASCIIMatrix To_size_t: {} is no size_t",d));
		}
	
	double GetScalar(const std::span<double>& spd)
		{
		if (spd.size() != 1)
			throw std::runtime_error(std::format("ReadASCIIMatrix GetScalar: Not a scalar line"));
		return spd.front();
		}

	std::vector<double> GetVector(const std::span<double>& spd, size_t n)
		{
		if (spd.size() != n)
			throw std::runtime_error(std::format("ReadASCIIMatrix GetScalar: Not a length {} vector line",n));
		return std::vector<double> (spd.begin(), spd.end());
		}

	
	void TVectorMatrixList::CreateFromNumberList(const TNumberList& nl)
		{
		auto iline = nl.Lines().begin();
		auto end = nl.Lines().end();
		while (iline != end)
			{
			// see if vector or matrix
			if (iline->empty() || iline->size() > 2)
				throw std::invalid_argument(std::format(
					"TVectorMatrixList::CreateFromNumberList: vector or matrix size expected"));
			if (iline->size() == 1) // vector
				{
				size_t nx = To_size_t(iline->front());
				TVec v;
				v.reserve(nx);
				++iline;
				for (size_t i = 0; i < nx; ++i)
					{
					if (iline == end)
						throw std::runtime_error(std::format(
							"TVectorMatrixList::CreateFromNumberList: not enough scalars"));
					v.push_back(GetScalar(*iline));
					++iline;
					}
				vec_mat_.push_back(std::move(v));
				}
			else // must be 2 -> Matrix
				{
				size_t nx = To_size_t((*iline)[0]);
				size_t ny = To_size_t((*iline)[1]);


				TMat m;
				m.reserve(nx);
				++iline;
				for (size_t i = 0; i < nx; ++i)
					{
					if (iline == end)
						throw std::runtime_error(std::format(
							"TVectorMatrixList::CreateFromNumberList: not enough scalars"));
					m.push_back(GetVector(*iline, ny));
					++iline;
					}
				vec_mat_.push_back(std::move(m));
				}
			}
		}

	const std::vector<TVectorMatrixList::TVecMat>& TVectorMatrixList::VecMat() const
		{
		return vec_mat_;
		}

	// ********************************************************************************************

	const TReadASCIIMatrix::TVec& TReadASCIIMatrix::x() const
		{
		return x_;
		}
	const TReadASCIIMatrix::TVec& TReadASCIIMatrix::y() const
		{
		return y_;
		}
	const TReadASCIIMatrix::TMat& TReadASCIIMatrix::z() const
		{
		return z_;
		}
	
	static bool IsEquidistant(const std::span<const double>& spd, double rel_eps)
		{
		if (spd.size() < 2)
			return false;
		// strictly ascending : back - front > 0
		double d_expected = (spd.back() - spd.front()) / (spd.size() - 1);
		double absmax = std::max({ spd.front(), spd.back() });
		double maxerr = 2 * rel_eps * absmax / d_expected;
		double curr = spd[0];
		for (auto it : spd.subspan(1))
			{
			if (std::abs((it - curr) - d_expected) > maxerr)
				return false;
			curr = it;
			}
		return true;
		}
	
	void TReadASCIIMatrix::Read(const std::string& fn, std::string delims, double rel_eps)
		{
		TVectorMatrixList vml;
		vml.Read(fn, delims);
		auto vm = vml.VecMat(); // vector of TVecMat
		if (vm.size() != 3)
			throw std::runtime_error(std::format(
				"TReadASCIIMatrix::Read: require 3 vector/matrix items, not {}", vm.size()));
		if (!std::holds_alternative<TVec>(vm[0]))
			throw std::runtime_error(std::format(
				"TReadASCIIMatrix::Read: require vector item for x, not matrix"));
		if (!std::holds_alternative<TVec>(vm[1]))
			throw std::runtime_error(std::format(
				"TReadASCIIMatrix::Read: require vector item for y, not matrix"));
		if (!std::holds_alternative<TMat>(vm[2]))
			throw std::runtime_error(std::format(
				"TReadASCIIMatrix::Read: require matrix item for z, not vector"));
		TVec& x = std::get<TVec>(vm[0]);
		TVec& y = std::get<TVec>(vm[1]);
		TMat& z = std::get<TMat>(vm[2]);
		Set(std::move(x), std::move(y), std::move(z), rel_eps);
		}


	void TReadASCIIMatrix::Set(const TVec & x, const TVec & y, const TMat & z, double rel_eps)
		{
		size_t nx = x.size();
		size_t ny = y.size();
		if (z.size() != nx || z[0].size() != ny)
			throw std::runtime_error(std::format(
				"TReadASCIIMatrix::Set: require {} by {} matrix, not {} by {} ",nx, ny, z.size(), z[0].size()));
		is_x_equidistant_ = IsEquidistant(x, rel_eps);
		is_y_equidistant_ = IsEquidistant(y, rel_eps);
		x_ = x;
		y_ = y;
		z_ = z;
		}

	void TReadASCIIMatrix::Set(TVec&& x, TVec&& y, TMat&& z, double rel_eps)
		{
		size_t nx = x.size();
		size_t ny = y.size();
		if (z.size() != nx || z[0].size() != ny)
			throw std::runtime_error(std::format(
				"TReadASCIIMatrix::Set: require {} by {} matrix, not {} by {} ", nx, ny, z.size(), z[0].size()));
		is_x_equidistant_ = IsEquidistant(x, rel_eps);
		is_y_equidistant_ = IsEquidistant(y, rel_eps);
		x_ = std::move(x);
		y_ = std::move(y);
		z_ = std::move(z);
		}

	static void Get_i_u(const std::vector<double> v, double x, bool is_equidistant, size_t& ix, double& ux)
		// precond: v not empty, x in [front,back]
		{
		if (x == v.front()) // also handles x_.size() = 1
			{
			ix = 0;
			ux = 0;
			return;
			}
		// now x_.size() > 1 
		if (x == v.back())
			{
			ix = v.size() - 2;
			ux = 1;
			return;
			}
		// now front < x < back
		if (is_equidistant)
			{
			double xstep = (v.back() - v.front()) / (v.size() - 1);
			double steps = (x - v.front()) / xstep;
			ix = static_cast<size_t>(std::floor(steps));
			ux = steps - ix;
			}
		else // bisect
			{
			size_t i0 = 0;
			size_t i1 = v.size();
			while (i1 > i0 + 1)
				{
				size_t imid = (i0 + i1)/2;
				double test = v[imid];
				if (test > x) // strictly in lower half
					i1 = imid;
				else
					i0 = imid;
				}
			ix = i0;
			double x0 = v[i0];
			double x1 = v[i1];
			ux = (x - x0) / (x1 - x0);
			}
		if (ix < 0 || ix >= (v.size() - 1) || ux < 0 || ux > 1)
			throw std::logic_error("ReadASCIIMatrix: Get_i_u: this cannot happen");
		}

	double TReadASCIIMatrix::Interpolate(double x, double y, double ooB) const
		{
		if (x < x_.front() || x > x_.back() || y < y_.front() || y > y_.back())
			return ooB;
		size_t ix;
		double ux;
		Get_i_u(x_, x, is_x_equidistant_, ix, ux);
		size_t iy;
		double uy;
		Get_i_u(y_, y, is_y_equidistant_, iy, uy);
		double z00 = z_[ix][iy];
		double z01 = z_[ix][iy+1];
		double z10 = z_[ix+1][iy];
		double z11 = z_[ix+1][iy+1];
		double rv = z00 * (1 - ux) * (1 - uy) + z01 * (1 - ux) * uy
			+ z10 * ux * (1 - uy) + z11 * ux * uy;
		return rv;
		}


	static TReadASCIIMatrix::TVec ReadNumberLine(TReadFile& f, const std::string& delims)
		{
		std::string_view ln = TrimWhiteSpace(f.ReadLine());
		TReadASCIIMatrix::TVec rv;
		size_t pos = ln.find_first_not_of(delims);
		while (pos != ln.size())
			{
			size_t tokend = ln.find_first_of(delims, pos);
			auto tok = ln.substr(pos, tokend - pos);

			}
		// TODO
		return rv;
		}

	double TReadASCIIMatrix::zMax() const
		{
		double t = std::numeric_limits<double>::lowest();
		for (const auto& iz : z_)
			{
			t = std::max(t, std::ranges::max(iz));
			}
		return t;
		}
	
	double TReadASCIIMatrix::zMin() const
		{
		double t = std::numeric_limits<double>::max();
		for (const auto& iz : z_)
			{
			t = std::min(t, std::ranges::min(iz));
			}
		return t;
		}

	
	}