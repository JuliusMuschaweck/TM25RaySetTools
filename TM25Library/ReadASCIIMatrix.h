/* Read x - y matrix from ASCII file
Format:
Line 1: nx
Line 2 .. nx+1: x values
line nx+2: ny
line nx+3 .. nx+2+ny: y values
line nx+ny+3: nx ny (repeated)
line nx+ny+4 .. nx+ny+(nx*ny)+3: z(x,y) values 

whitespace at beginning and end is 

*/
#ifndef __ReadASCIIMatrix_H
#define __ReadASCIIMatrix_H
// #include "ReadFile.h"
#include <span>
#include <variant>
#include <string>
#include <vector>
#include <charconv>
#include <cmath>
#endif __ReadASCIIMatrix_H

namespace TM25
	{

	/*
	First, we define a generic text file format for structured numbers -> TNumberList

	We read an ASCII file that contain real numbers, line by line
	1. White space, " \t", is trimmed
	2. Lines that are empty or do not start with "-0123456789" are discarded as comments.
	3. The line is tokenized into substrings, separated by delims, default " \t,;"
		Delims can be mixed, "0 1\t,-2.3e-5;,3" -> {"0","1","-2.3e-5","3"}
		The line must contain only tokens that can be converted to real numbers	
	4. Each token is converted to a real number. 
	5. These numbers are push_back'ed to an array of real
	6. For each line, the # of numbers is counted and push_back'ed to an array of size_t

	Second, based on above we define a text file format for a sequence of vectors and matrices:

	nx
	x_1
	...
	x_nx

	nrow ncol
	z_1_1 ... z_1_ncol
	...
	z_nrow_1 ... z_nrow_ncol

	Third, based on above, we define a matrix of z-values on a rectangular x-y grid
	1. nx x-values, strictly ascending
	2. ny y-values, strictly ascending
	3. nx by ny z values 

	Example:

		A comment line 

		Empty lines are comments, too
		nx: # of x values
		2
		followed by nx single number lines
		1.1
		2.2
		ny: # of y values
		3
		followed by ny single number lines
		10
		11.1
		15
		followed by nx by ny numbers
		1 2 3E2
		4 5 -6.0e-1

	Second, we define a text file format for a matrix of values on a rectangular x/y grid
	x, y are vectors that can have different length, nx and ny, and need only be strictly ascending, not equidistant
	The format starts with the x and y values, followed by an nx by ny matrix of z values

	*/

	class TReadFile; // forward decl from ReadFile.h

	class TNumberList
		{
		public:
			void Read(const std::string& fn, const std::string& delims = " \t,;");
			void Read(TReadFile& f, const std::string& delims = " \t,;");
			const std::vector<std::span<double>>& Lines() const; // the list of number lines

			// low level access functions
			const std::vector<double>& Num() const;
			const std::vector<size_t>& Num_per_line() const;
		private:
			void Clear();
			std::vector<double> num_; // the array of all numbers in the file
			std::vector<size_t> num_per_line_;
			std::vector<std::span<double> > lines_;
		};

	class TVectorMatrixList
		{
		public:
			void Read(TReadFile& f, const std::string& delims = " \t,;");
			void Read(const std::string& fn, const std::string& delims = " \t,;");
			void CreateFromNumberList(const TNumberList& nl);
			using TVec = std::vector<double>;
			using TMat = std::vector<TVec>;
			using TVecMat = std::variant<TVec, TMat>;
			const std::vector<TVecMat>& VecMat() const;
		private:
			std::vector<TVecMat> vec_mat_;
		};

	class TReadASCIIMatrix
		{
		public:
			using TVec = std::vector<double>;
			using TMat = std::vector<TVec>;
			const TVec& x() const;
			const TVec& y() const;
			const TMat& z() const;
			
			void Read(const std::string& fn, std::string delims, double rel_eps);

			void Set(const TVec& x, const TVec& y, const TMat& z, double rel_eps);

			void Set(TVec&& x, TVec&& y, TMat&& z, double rel_eps);

			double Interpolate(double x, double y, double ooB = NAN) const;

			double zMax() const;
			double zMin() const;

		private:
			bool is_x_equidistant_ = false;
			bool is_y_equidistant_ = false;
			TVec x_;
			TVec y_;
			TMat z_;
		};

	}
