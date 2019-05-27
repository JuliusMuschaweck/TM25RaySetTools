#ifndef __WriteFile_H
#define __WriteFile_H

#include <string>
#include <array>
#include <vector>
#include <fstream>

namespace TM25
	{
	class TWriteFile
		{
		public:
			explicit TWriteFile(const std::string& fn); 
			// opens file for binary writing
			
			template<typename T>
			size_t Write(const T& t);
			// writes the binary representation of t to file
			// returns # of bytes written == sizeof(T)

			size_t WriteCString(const std::string& s);
			size_t WriteCString(const char* s);
			// wrítes null terminated character string, returns bytes written = length + 1

			size_t WriteUTF32TextBlock(const std::u32string& s, size_t nBytes);
			// writes at most nBytes % 4 characters from s, then writes 0 until bBytes written

			template<typename Iterator> // for u32string fields
			size_t WriteRange(Iterator begin, Iterator end);

			size_t WriteZeroBytes(size_t n);

			template<typename T, size_t N>
			size_t WriteArray(const std::array<T, N>& t);
			
			template<typename T>
			size_t WriteVector(const std::vector<T>& t);

			void WriteBytes(const void* p, size_t n);

			size_t BytesWritten() const;

			~TWriteFile();
		private:
			std::vector<char> buf_;
			static const size_t buf_size_ = 10 * 1024 * 1024; // write chunks of 10 MB;
			size_t bufpos_;
			size_t written_;
			std::ofstream f_;
			
			void WriteBuf();
			void WriteBytesDirect(const void* p, size_t n);

		};


	} // namespace

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// %%%% template definitions %%%%%%%%%%%%%%%%%%%

namespace TM25
	{
	template<typename T>
	inline size_t TWriteFile::Write(const T & t)
		{
		WriteBytes(&t, sizeof(T));
		return sizeof(T);
		}

	template<typename Iterator>
	inline size_t TWriteFile::WriteRange(Iterator begin, Iterator end)
		{
		size_t rv = 0;
		using it = Iterator;
		for (it i = begin; i != end; ++i)
			rv += Write(*i);
		return rv;
		}

	template<typename T, size_t N>
	inline size_t TWriteFile::WriteArray(const std::array<T, N>& t)
		{
		size_t n = N * sizeof(T);
		WriteBytes(t.data(), n);
		return n;
		}

	template<typename T>
	inline size_t TWriteFile::WriteVector(const std::vector<T>& t)
		{
		size_t n = t.size() * sizeof(T);
		WriteBytes(t.data(), n);
		return n;
		}

	} // namespace


#endif // include guard
