/******************************************************************
Provided by Julius Muschaweck, JMO GmbH, Gauting to the public domain
under the Unlicense, see unlicense.txt in the repository
or http://unlicense.org/ 
2019-01-19
******************************************************************/
// ReadFile.h 
// class TReadFile provides read access to large binary files
// template functions allow type safe reading of nearly any type
// TReadFile uses an internal buffer, reading data from file in chunks
// of static const size_t buf_size, currently 10 MB.
// On a standard desktop PC, it took about 2 seconds to read 500 MB
// of floats, including the sanity checking of the content.

// Error handling: Throws TReadFileError : public std::runtime_error
// if file cannot be opend, is not open when reading, or too short.


#ifndef __ReadFile_H
#define __ReadFile_H
#include <string>
#include <array>
#include <vector>
#include <streambuf>
#include <stdexcept>
#include <fstream>

namespace TM25
	{
	class TReadFile
		{
		public:
			explicit TReadFile(const std::string& filename);
			// opens file for reading, determines file size

			template<typename T> 
			T Read();
			// read and return a T from its binary representation
			// advance file position
			// may throw if file not open or too short

			template<typename T, size_t N>
			std::array<T, N>  Read();
			// read and return an array<N> of T from its binary representation
			// advance file position
			// use when N is known at compile time
			// may throw if file not open or too short

			template<typename T>
			std::vector<T> Read(const size_t n);
			// read and return a vector size n of T from its binary representation
			// advance file position
			// use when n is known at run time only
			// may throw if file not open or too short

			template<typename T>
			void Read(std::vector<T>& v, const size_t n);
			// read and assign a vector size n of T from its binary representation
			// advance file position
			// use when n is known at run time only
			// may throw if file not open or too short

			template<typename T>
			T Peek();
			// read and return a T from its binary representation
			// do not advance file position
			// may throw if file not open or too short
		
			template<size_t N>
			std::array<char, N> ReadBytes();
			// read and return an array<N> of char from its binary representation
			// advance file position
			// use when N is known at compile time
			// may throw if file not open or too short

			std::vector<char> ReadBytes(const size_t n);
			// read and return a vector size n of char from its binary representation
			// advance file position
			// use when n is known at run time only
			// may throw if file not open or too short

			template<typename T>
			bool ReadIf(T& t);
			// read a T from its binary representation
			// advance file position
			// returns true on success, false if file not open or too short

			template<size_t N>
			bool ReadBytesIf(std::array<char, N>& t);
			// read an array<N> of char from its binary representation
			// advance file position
			// use when N is known at compile time
			// returns true on success, false if file not open or too short

			bool ReadBytesIf(std::vector<char>& t, const size_t n);
			// read a vector size n of char from its binary representation
			// advance file position
			// use when n is known at run time only
			// returns true on success, false if file not open or too short

			bool AtEof() const;
			// returns true if all bytes in the file have been read

		private:
			bool ReadDirectIf(char* t, size_t n);
			// read n char's from binary representation
			// advance file position
			// t must point to sufficient memory
			// returns true on success, throws TM25Error if file not open or too short
			bool PeekDirectIf(char* t, size_t n);
			
			size_t BufAvail() const; // bytes remaining in buffer
			void Advance(size_t n);
			void Overflow(); // read new data from file and fill buffer

			size_t ReadDirectFromFile(char* t, size_t n); // return bytes read

			size_t pos; // position in the file.
			// size_t maxpos; // max. position in the file == file size
			size_t buf_pos; // position in the buffer buf
			size_t buf_end; // max. position in the buffer + 1, < bufsize iff eof has been loaded
			bool all_read_from_f; // true if reading from f has reached eof;
			static const size_t buf_size = 10 * 1024 * 1024; // read chunks of 10 MB;
			// static constexpr size_t overlap_size = 1024 * 1024; // 1 MB overlap
			std::vector<char> buf;
			std::ifstream f;
		};

	class TReadFileError : public std::runtime_error
		{
		// the exception thrown by the member functions
		public:
			explicit TReadFileError(const std::string& msg) : std::runtime_error(msg) {};
		};
	}

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// template definitions

namespace TM25
	{


	template<typename T>
	T TReadFile::Read()
		{
		T rv;
		char* p = reinterpret_cast<char*>(&rv);
		if (!ReadDirectIf(p, sizeof(T)))
			throw TReadFileError("TReadFile::Read: failed");
		return rv;
		}

	template<typename T, size_t N>
	std::array<T, N>  TReadFile::Read()
		{
		std::array<T, N> rv;
		char* p = reinterpret_cast<char*>(&(rv[0]));
		if (!ReadDirectIf(p, N * sizeof(T)))
			throw TReadFileError("TReadFile::Read: failed");
		return rv;
		}

	template<typename T>
	std::vector<T> TReadFile::Read(const size_t n)
		{
		std::vector<T> rv(n);
		Read<T>(rv, n);
		return rv;
		}

	template<typename T>
	void TReadFile::Read(std::vector<T>& v, const size_t n)
		{
		v.resize(n);
		char* p = reinterpret_cast<char*>(&(v[0]));
		if (!ReadDirectIf(p, n * sizeof(T)))
			throw TReadFileError("TReadFile::Read: failed");
		}

	template<typename T>
	T TReadFile::Peek()
		{
		T rv;
		char* p = reinterpret_cast<char*>(&rv);
		if (!PeekDirectIf(p, sizeof(T)))
			throw TReadFileError("TReadFile::Peek: failed");
		return rv;
		}

	template<size_t N>
	std::array<char, N> TReadFile::ReadBytes()
		{
		return Read<std::array<char, N>>();
		}

	template<typename T>
	bool TReadFile::ReadIf(T& t)
		{
		T rv;
		char* p = reinterpret_cast<char*>(&rv);
		return ReadDirectIf(p, sizeof(T));
		}

	template<size_t N>
	bool TReadFile::ReadBytesIf(std::array<char, N>& t)
		{
		return ReadDirectIf(&(t[0]), N);
		}
	}
#endif