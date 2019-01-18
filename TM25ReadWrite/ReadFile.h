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
			T Peek();
			// read and return a T from its binary representation
			// do not advance file position
			// may throw if file not open or too short

			
			template<typename T>
			bool PeekOk();
			// check if file is long enough
			// do not advance file position
			// does not throw

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
			// returns true on success, false if file not open or too short
			bool PeekDirectIf(char* t, size_t n);
			size_t pos;
			size_t maxpos;
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
		char* p = reinterpret_cast<char*>(&(rv[0]));
		if (!ReadDirectIf(p, n * sizeof(T)))
			throw TReadFileError("TReadFile::Read: failed");
		return rv;
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

	template<typename T>
	bool TReadFile::PeekOk()
		{
		return (maxpos - pos) >= sizeof(T);
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