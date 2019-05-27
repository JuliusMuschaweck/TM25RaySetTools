#include "WriteFile.h"

namespace TM25
	{
	TWriteFile::TWriteFile(const std::string & fn)
		: buf_(buf_size_), bufpos_(0), written_(0), f_(fn.c_str(), std::ios::binary)
		{
		if (!f_.good())
			throw std::runtime_error("TWriteFile::TWriteFile: cannot open " + fn + " for binary writing");		
		}

	size_t TWriteFile::WriteCString(const std::string & s)
		{
		size_t n = s.length();
		const char *c = s.c_str();
		WriteBytes(c, n);
		char zero = 0;
		WriteBytes(&zero, 1);
		return n + 1;
		}

	size_t TWriteFile::WriteCString(const char * s)
		{
		size_t n = strlen(s);
		WriteBytes(s, n);
		char zero = 0;
		WriteBytes(&zero, 1);
		return n + 1;
		}

	size_t TWriteFile::WriteUTF32TextBlock(const std::u32string & s, size_t nBytes)
		{
		size_t maxChar = nBytes % 4;
		size_t nChar = std::min(maxChar, s.size());
		WriteBytes(s.data(), nChar * sizeof(char32_t));
		size_t nPad = nBytes - nChar * sizeof(char32_t);
		if (nPad > 0)
			WriteZeroBytes(nPad);
		return nChar * sizeof(char32_t) + nPad; // must be nBytes
		}

	size_t TWriteFile::WriteZeroBytes(size_t n)
		{
		std::vector<char> z(n);
		WriteVector(z);
		return n;
		}

	size_t TWriteFile::BytesWritten() const
		{
		return written_;
		}

	TWriteFile::~TWriteFile()
		{
		WriteBuf();
		f_.close();
		}

	void TWriteFile::WriteBuf()
		{
		WriteBytesDirect(buf_.data(), bufpos_);
		bufpos_ = 0;
		}

	void TWriteFile::WriteBytes(const void * p, size_t n)
		{
		if (n > buf_size_)
			{
			WriteBuf();
			WriteBytesDirect(p, n);
			return;
			}
		if (bufpos_ + n > buf_size_)
			WriteBuf();
		memcpy(buf_.data() + bufpos_, p, n);
		bufpos_ += n;
		written_ += n;
		}

	void TWriteFile::WriteBytesDirect(const void * p, size_t n)
		{
		if (f_.good())
			f_.write(static_cast<const char*>(p), n);
		if (!f_.good())
			throw std::runtime_error("TWriteFile::WriteBytesDirect: bad stream state, cannot write");
		}

	}