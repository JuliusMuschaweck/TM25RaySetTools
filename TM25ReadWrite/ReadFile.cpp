#include "ReadFile.h"

#include <sys/stat.h>

namespace TM25
	{
	size_t FileSize(const std::string& fn)
		{
		#ifdef _MSC_VER
 		struct _stat64 ss;
		int test = _stat64(fn.c_str(), &ss);
		#else
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
			f.open(filename,std::ios_base::in);
			if (!f.good())
				throw std::runtime_error("file not found");
			pos = 0;
			maxpos = FileSize(filename);
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
		return ReadDirectIf(&(t[0]),n);
		}

	bool TReadFile::AtEof() const
		{
		return f.eof();
		}

	bool TReadFile::ReadDirectIf(char* t, size_t n)
		{
		std::streambuf* buf = f.rdbuf();
		if (buf && f.good())
			{
			size_t avail = buf->in_avail();
			if (maxpos - pos >= n)
				{
				size_t nread = buf->sgetn(t, n);
				avail = buf->in_avail();
				if (nread != n)
					{
					buf->pubseekpos(pos);
					throw TReadFileError("TReadFile::ReadDirectIf: file too short");
					}
				}
			else
				throw TReadFileError("TReadFile::ReadDirectIf: file too short");
			}
		else
			{
			throw TReadFileError("TReadFile::ReadDirectIf: file not open");
			}
		pos += n;
		return true;
		}

	bool TReadFile::PeekDirectIf(char* t, size_t n)
		{
		std::streambuf* buf = f.rdbuf();
		if (buf && f.good())
			{
			if (maxpos - pos >= n)
				{
				size_t nread = buf->sgetn(t, n);
				std::streambuf::off_type signed_n = static_cast<std::streamoff>(nread);
				buf->pubseekoff( - signed_n, std::ios_base::cur);
				if (nread != n)
					{
					throw TReadFileError("TReadFile::ReadDirectIf: file too short");
					}
				}
			else
				throw TReadFileError("TReadFile::ReadDirectIf: file too short");
			}
		else
			{
			throw TReadFileError("TReadFile::ReadDirectIf: file not open");
			}
		return true;
		}

	}

