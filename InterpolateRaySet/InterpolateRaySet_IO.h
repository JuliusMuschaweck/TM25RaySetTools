#ifndef __INTERPOLATERAYSET_IO_H
#define __INTERPOLATERAYSET_IO_H

#include <TM25.h>
#include "InterpolateRaySetConfig.h"
#include <iostream>
#include <mutex>

TM25::TTM25RaySet ReadRaySet(const TInterpolateRaySetCfg& cfg, std::ostream& info);

void WriteRaySet(TM25::TTM25RaySet& rs, const TInterpolateRaySetCfg& cfg, std::ostream& info);

class ComposeStream : public std::ostream
	{
	public:
		ComposeStream();
		void LinkStream(std::ostream& out);

	private:
		struct ComposeBuffer : public std::streambuf
			{
			void addBuffer(std::streambuf* buf);
			virtual int overflow(int c);
			private:
				std::vector<std::streambuf*>    bufs;
			};
		ComposeBuffer myBuffer;
	};

class TThreadSafe_stdcout
	{
	public:
		template<typename T>
		void Write(const T& t);
	private:
		std::mutex m_;
	};

class TLogPlusCout : public ComposeStream
	{
	public:
		TLogPlusCout(bool doCout, const std::string& logfn);
	private:
		std::unique_ptr<std::ofstream> log_;
	};

// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

template<typename T>
void TThreadSafe_stdcout::Write(const T& t)
	{
	std::lock_guard<std::mutex> lock(m_);
	std::cout << t;
	}



#endif