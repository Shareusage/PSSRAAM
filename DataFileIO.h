#ifndef DATAFILEIO_H
#define DATAFILEIO_H

#include <windows.h>
#include <vector>
#include <string>
using namespace std;

//structure header file signal
namespace dataIO
{
	struct SSignal {
		char ftitle[8];         // code header, indicating the binary file filling signal
		char NameElement[100];  // name of the element from which the signal is removed
		char NameIO[100];       // the name of the output signal from which is removed
		int ne;                 // Element ID is removed from the signal
		int nio;                // ID output signal from which is removed
		int TypeSignal;         // signal type 0 - integer, 1 - real, 2 - complex
		int TypeReg;            // flag that the signal was taken uniformly with a given period
		double start;           // the start time point (used to write a regular situation)
		double period;          // the sampling period of the signal
		int NumPoint;           // many points
		int NumDataName;        // the number of elements in the structure of a signal which contains a data structure
		int res[50];
		//after the structure Signal, in the file is a list of the names of the structural elements of the signal. size lines - 40 characters
	};

	enum DataFileFormat { plain = 0, rlezip = 1, unknown = 255 };

	class DataFileIO
	{
	public:

		// Constructor - opening a file named
		// 1. reading the header and data extraction of series
		// 2. calculation steps and extreme values of the arguments
		DataFileIO(const wchar_t *fname);

		// Destructor - closing the file and release of resources
		~DataFileIO(void);

		// Checking the status of the error
		int GetError(void) const
		{
			return Error;
		}

		// Obtaining information about the type of data in the file -
		// does not support the file format
		int GetDataInfo(void);

		// Get the length of the sequence in the file (how many points in the series)
		__int64 GetNumPoint(void) const
		{
			return NumPoint;
		}

		// Getting the number of series in the file (ie pairs TIME - VALUE)
		int GetSeriesCount(void) const
		{
			return int(SeriesCount);
		}

		// Getting the min and max values in the argument
		double GetMinArg(void) const
		{
			return xmin;
		}
		double GetMaxArg(void) const
		{
			return xmax;
		}


		// 03.04.2008
		const SSignal *GetHeader(void) const
		{
			return &header;
		}
		const char *GetElementName(int i) const
		{
			return NameBuf[i];
		}

		// Getting the value of the step argument
		double GetArgStep(void) const
		{
			return ArgStep;
		}

		// Series Index number is the minimum and maximum values
		// on the interval from pos to pos + length and written in MinY and MaxY
		int GetMinMaxByIndex(unsigned int Index,
			__int64 pos, unsigned int length,
			double &MinY, double &MaxY);

		// Download the min and max values for all series in the interval from pos
		// to pos + length arrays and MinY MaxY
		int LoadMinMaxSignal(__int64 pos, unsigned int length,
			double *MinY, double *MaxY);

		// Just download the signal samples in the buffer (not extreme)
		int LoadSignal(unsigned int Index, __int64 pos, unsigned int length,
			double *Y, double *x = nullptr);

		// Returns FILETIME structure over time the file was created
		FILETIME GetCreationTime(void);

	private:
		// Search in compressed values
		__int64 BinSearch(double val);

		// Reading of a file argument in the case where more than 50 million marks
		double ReadArg(int index);

		// Reading reports in the buffer with the pointer
		int ReadMarks(__int64 pos, int length, double *buf);

		// Reading marks the current position of the file pointer
		//int ReadNextMarks(int length, double *buf);

		// File descriptor data
		HANDLE hFile;

		// error status
		int Error;

		// The file header with service information
		SSignal header;

		// The position of the pointer to the beginning of the data in the file
		LONG fDataPtr;

		// Container term with the names of series
		char **NameBuf;

		// Signal parameters file
		__int64 NumPoint;           // the number of marks in the file
		unsigned int SeriesCount;   // the number of series in the file

		// The limits of variation of the argument
		double xmin, xmax;
		double ArgStep;

		// Variable format file
		int FileFormat;

		// For a compressed file format - stored values of the argument
		double *argument;

		// A sign of a large file - work with the file instead of in memory
		bool bLargeFile;

		// The number of recorded marks in a compressed file
		int zipPoint;
	};

	///////////////////////////////////////////////////////////////////////////////
	class DataFileIOWriter
	{
	public:
		// Opening and header record
		DataFileIOWriter(const wchar_t *fname, unsigned int nSeries);
		DataFileIOWriter(const wchar_t *fname, vector<string> nSeries);
		// close
		~DataFileIOWriter();

		// Record count of all series
		int WriteMark(double t, double *series);

		// error status
		int GetError(void) const
		{
			return Error;
		}
		// close files
		void Close();

	private:
		// Descriptor file for writing
		HANDLE hFile;

		// The value of the last recorded argument
		double lastarg{0};

		// The number of series in the file
		unsigned int SeriesCount;

		// error
		int Error;
	};
}
#endif // DATAFILEIO_H
