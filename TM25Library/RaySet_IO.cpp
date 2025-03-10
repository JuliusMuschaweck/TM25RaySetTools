#include <RaySet_IO.h>
#include <ASCIIRayFile.h>
#include <TranslateZemax.h>
#include <TranslateLightTools.h>

using std::cout;
using std::endl;


TM25::TTM25RaySet ReadRaySet(const TSection& rsc, std::ostream& info)
	{
	TM25::TTM25RaySet rs;
	if (rsc.IsEmpty("inputRayFileName"))
		throw std::runtime_error("inputRayFileName must be supplied");
	const std::string format = rsc.String("inputRayFileFormat");
	const std::string fn = rsc.String("inputRayFileName");
	if (format.compare("TM25") == 0)
		{
		info << "Reading TM25 ray file " << fn << endl;
		rs.Read(fn);
		// rs.Make_k_unit();
		if (rs.Warnings().empty())
			info << "No warnings" << endl;
		for (auto w : rs.Warnings())
			{
			info << "Warning: " << w << "\n";
			}
		}
	else if (format.compare("ZemaxBinary") == 0)
		{
		info << "Reading Zemax binary ray file " << fn << endl;
		TM25::TZemaxRaySet zemaxRaySet(fn);
		rs = TM25::ZemaxBinaryToTM25(zemaxRaySet);
		}
	else if (format.compare("ZemaxText") == 0)
		{
		info << "Reading Zemax text ray file " << fn << endl;
		TM25::TASCIIRaySetOptions opts;
		opts.minHeaderLines = 1;
		opts.wavelengthUnit_ = TM25::TASCIIRaySetOptions::WLU::micrometer;
		rs = TM25::ReadGenericASCIIRaySet(fn, opts);
		}
	else if (format.compare("LightToolsBinary") == 0)
		{
		info << "Reading LightTools binary ray file " << fn << endl;
		TM25::TLightToolsRaySet lightToolsRaySet(fn);
		rs = TM25::LightToolsBinaryToTM25(lightToolsRaySet);
		}

	else
		throw std::runtime_error("unknown input ray file format: " + format);
	return rs;
	}

void WriteRaySet(TM25::TTM25RaySet& rs, const TSection& rsc, std::ostream& info)
	{
	const std::string format = rsc.String("outputRayFileFormat");
	const std::string fn = rsc.String("outputRayFileName");
	if (format.compare("TM25") == 0)
		{
		info << "Writing TM25 ray file " << fn << endl;
		rs.Write(fn);
		if (rs.Warnings().empty())
			info << "No warnings" << endl;
		for (auto w : rs.Warnings())
			{
			info << "Warning: " << w << "\n";
			}
		}
	else if (format.compare("ZemaxBinary") == 0)
		{
		info << "Writing  Zemax binary ray file " << fn << endl;
		TM25::TZemaxRaySet zemaxRaySet = TM25ToZemaxBinary(rs);
		zemaxRaySet.Write(fn);
		}
	else if (format.compare("LightToolsBinary") == 0)
		{
		info << "Writing  LightTools binary ray file " << fn << endl;
		TM25::TLightToolsRaySet ltRaySet = TM25ToLightToolsBinary(rs);
		ltRaySet.Write(fn);
		}
	else
		throw TM25::TM25Error("WriteRaySet: unknown output format");
	}

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
