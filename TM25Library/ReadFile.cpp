/******************************************************************
Provided by Julius Muschaweck, JMO GmbH, Gauting to the public domain
under the Unlicense, see unlicense.txt in the repository
or http://unlicense.org/
2019-01-19
******************************************************************/
// implementation of non-template functions of ReadFile.h

#ifndef _CONSOLE
	#include <stdafx.h>
#endif
#include "ReadFile.h"
#include <cstring> // for memcpy
#include <sys/stat.h>
#include "ParseString.h"

namespace TM25
	{

	size_t FileSize(const std::string& fn)
		{
#ifdef _MSC_VER
		struct _stat64 ss;
		int test = _stat64(fn.c_str(), &ss);
#else	
		// we could use SFINAE, see https://stackoverflow.com/questions/257288/templated-check-for-the-existence-of-a-class-member-function
		// to determine if stat64 is available 
		struct stat ss;
		int test = stat(fn.c_str(), &ss);
#endif
		static_assert(sizeof(ss.st_size) == 8, "FileSize: stat.st_size is not 8 byte");
		if (test == 0)
			return ss.st_size;
		else
			throw TReadFileError("FileSize: cannot open file " + fn);
		}

	TReadFile::TReadFile(const std::string& filename)
		{
		try
			{
			buf.resize(buf_size);
			f.open(filename, std::ios_base::in | std::ios_base::binary);
			if (!f.good())
				throw std::runtime_error("file not found");
			pos = 0;
			// maxpos = FileSize(filename);
			buf_pos = 0;
			buf_end = 0;
			all_read_from_f = false;
			Overflow();
			}
		catch (std::runtime_error e)
			{
			throw TReadFileError("TReadFile::TReadFile: cannot open file " + filename + ", error: " + e.what());
			}
		}

	std::vector<char> TReadFile::ReadBytes(const size_t n)
		{
		std::vector<char> rv(n);
		if (!ReadBytesIf(rv, n))
			throw TReadFileError("TReadFile::ReadBytes: cannot read");
		else
			return rv;
		}

	bool TReadFile::ReadBytesIf(std::vector<char>& t, const size_t n)
		{
		t.resize(n);
		return ReadDirectIf(&(t[0]), n);
		}

	std::string TReadFile::ReadLine(bool includeEOL)
		{
		std::string rv;
		ReadLine(rv, includeEOL);
		return rv;
		}

	const std::string& TReadFile::ReadLine(std::string& rv, bool includeEOL)
		{
		rv.resize(0); // does not free the underlying memory! 
		while (!AtEof()) // see what comes next
			{
			char c;
			if (ReadIf<char>(c))
				{ // next char is available
				if (c != 10 && c != 13)
					{ // a non-eol character
					rv.push_back(c);
					}
				else // a newline
					{
					if (includeEOL)
						rv.push_back(c);
					if (c == 13) // see if next is 10, Windows standard
						{
						if (PeekDirectIf(&c, 1) && c == 10) // yes
							{
							Advance(1);
							if (includeEOL)
								rv.push_back(c);
							}
						}
					else // c == 10, see if next is 13, nonstandard
						{
						if (PeekDirectIf(&c, 1) && c == 13) // yes
							{
							Advance(1);
							if (includeEOL)
								rv.push_back(c);
							}
						}
					return rv;
					}
				}
			else // ReadIf failed for whatever reason -- done
				return rv;
			}
		return rv;
		}

	bool TReadFile::AtEof() const
		{
		return all_read_from_f && (buf_pos == buf_end);
		}

	bool TReadFile::ReadDirectIf(char* t, size_t n)
		{
		if (n <= BufAvail())
			// the easy case: just copy data from buf and advance pos
			{
			std::memcpy(t, buf.data() + buf_pos, n);
			Advance(n);
			}
		else
			// insufficient data in buffer
			{
			if (n <= buf_size)
				// if n <= bufsize fill buffer, copy and advance
				{
				Overflow();
				if (n > BufAvail())
					throw TReadFileError("TReadFile::ReadDirectIf: file too short");
				std::memcpy(t, buf.data() + buf_pos, n);
				Advance(n);
				}
			else
				// n very large, > bufsize, loop through buffers and fill large output
				{
				while (n > buf_size)
					{
					if (!ReadDirectIf(t, buf_size))
						throw TReadFileError("TReadFile::ReadDirectIf: file too short");
					t += buf_size;
					n -= buf_size;
					}
				if (!ReadDirectIf(t, n))
					throw TReadFileError("TReadFile::ReadDirectIf: file too short");
				}
			}
		return true;
		}

	bool TReadFile::PeekDirectIf(char* t, size_t n)
		{
		if (n <= BufAvail())
			// the easy case: just copy data from buf
			{
			std::memcpy(t, buf.data() + buf_pos, n);
			}
		else
			// insufficient data in buffer
			{
			if (n <= buf_size)
				// if n <= bufsize fill buffer, copy and advance
				{
				Overflow();
				if (n > BufAvail())
					// throw TReadFileError("TReadFile::PeekDirectIf: file too short");
					return false;
				std::memcpy(t, buf.data() + buf_pos, n);
				}
			else
				{
				throw TReadFileError("TReadFile::PeekDirectIf: cannot peek variables larger than "
					+ std::to_string(buf_size) + " bytes");
				}
			}
		return true;
		}

	size_t TReadFile::BufAvail() const // bytes remaining in buffer
		{
#ifndef NDEBUG
		if (buf_pos > buf_end)
			throw TReadFileError("TReadFile::BufAvail(): buf_pos > buf_end, this cannot happen");
#endif		 
		return buf_end - buf_pos;
		}

	void TReadFile::Advance(size_t n)
		{
#ifndef NDEBUG
		if ((buf_pos + n) > buf_end)
			throw TReadFileError("TReadFile::Advance(): buf_pos + n > buf_end, this cannot happen");
#endif		 
		buf_pos += n;
		pos += n;
		}

	void TReadFile::Overflow() // read new data from file and fill buffer
		{
#ifndef NDEBUG
		if (buf_pos > buf_end)
			throw TReadFileError("TReadFile::BufAvail(): buf_pos > buf_end, this cannot happen");
#endif		 
		size_t n_remain = buf_end - buf_pos;
		std::memmove(buf.data(), buf.data() + buf_pos, n_remain);
		buf_pos = 0;
		buf_end = n_remain;
		size_t n_req = buf_size - buf_end;
		size_t n_read = ReadDirectFromFile(buf.data() + buf_end, n_req);
		buf_end += n_read;
#ifndef NDEBUG
		if (buf_end > buf_size)
			throw TReadFileError("TReadFile::BufAvail(): buf_end > buf_size, this cannot happen");
#endif		 

		}

	size_t TReadFile::ReadDirectFromFile(char* t, size_t n)
		{
		if ((n > 0) && all_read_from_f)
			// throw TReadFileError("TReadFile::ReadDirectFromFile(): trying to read beyond eof, this cannot happen");
			return 0;
		size_t rv;
		try
			{
			rv = f.rdbuf()->sgetn(t, n);
			if (rv < n)
				all_read_from_f = true;
			}
		catch (std::runtime_error e)
			{
			throw TReadFileError("TReadFile::ReadDirectFromFile(): runtime error " + std::string(e.what()));
			}
		return rv;
		}

	}