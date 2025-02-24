#include <LogPlusCout.h>
#include <fstream>

using std::cout;
using std::endl;

// create an ostream that sends output to many other ostreams
void ComposeStream::ComposeBuffer::addBuffer(std::streambuf* buf)
	{
	bufs.push_back(buf);
	}
int ComposeStream::ComposeBuffer::overflow(int c)
	{
	auto do_sputc = [](std::streambuf* buf, int c)
		{
		buf->sputc(c);
		};
	for (auto buf : bufs)
		{
		do_sputc(buf, c);
		}
	return c;
	}

ComposeStream::ComposeStream()
	:std::ostream(nullptr)
	{
	std::ostream::rdbuf(&myBuffer);
	}
void ComposeStream::LinkStream(std::ostream& out)
	{
	out.flush();
	myBuffer.addBuffer(out.rdbuf());
	}


TLogPlusCout::TLogPlusCout(bool doCout, const std::string& logfn)
	{
	if (doCout)
		LinkStream(std::cout);
	if (!(logfn.empty()))
		{
		log_.reset(new std::ofstream(logfn));
		if (!(log_->good()))
			{
			cout << "warning: cannot open log file " << logfn << " for writing" << endl;
			log_.reset(nullptr);
			}
		}
	if (log_)
		LinkStream(*log_);
	}
