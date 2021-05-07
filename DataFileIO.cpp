//#include "stdafx.h"
#include "DataFileIO.h"
#include <string.h>
#include <stdio.h>

#include <iostream>
#include <fstream>
using namespace std;
using namespace dataIO;
#ifdef QT_CREATOR
#pragma comment(lib, "user32")
#endif

//#define _USE_READ_FILE_DUMP_

void ShowErrorF(LPTSTR lpszFunction)
{
    MessageBox(nullptr, lpszFunction, L"ShowErrorF", MB_OK | MB_ICONWARNING);
}

///////////////////////////////////////////////////////////////////////////////
// Search series value by argument
// Searches for the index of the argument in the array index of the compressed file
// The result - the index values
// Example: arg = {2, 3, 5, 7, 11}
// Find t = 6 - answer 2 -> arg[2] = 5
// see series[index][2] to result series with t = 6
__int64 DataFileIO::BinSearch(double val)
{
    __int64 l, r, m;
    l = 0;
    r = zipPoint - 1;

    do {
        m = (r + l) / 2;

        if (ReadArg(m) > val + 1e-16) {
            r = m;
        } else {
            l = m + 1;
        }
    } while (l < r);

    if (ReadArg(l) <= val + 1e-16)
        l++;

    //if( l == 0 || l == zipPoint) return -1;
    //return l-1;

    /// It is based on the correctness of the input argument
    if (l == 0) l = 1;

#ifdef _USE_READ_FILE_DUMP_
    // dumps
    if (1) {
        ofstream tmpf("dump.txt", ios_base::app);
        tmpf << "BinSearch() for " << val << " returned " << l << " (arg[l] = " << ReadArg(l) << endl;
        tmpf.close();
    }
#endif

    return l - 1;
}

///////////////////////////////////////////////////////////////////////////////
// Access to the values of the argument depending on the file size
double DataFileIO::ReadArg(int index)
{
    // When the memory is not enough for the buffer - working with files
    if (bLargeFile == true) {
        double varg = 0;
        DWORD br = 0;
        // Setting a pointer to the first count of the required data in the file
        __int64 offset = (index * (SeriesCount + 1)) * sizeof(double);
        LONG lDistLow = static_cast<LONG>(fDataPtr + offset);
        LONG lDistHigh = static_cast<LONG>((fDataPtr + offset) >> 32);
        SetFilePointer(hFile, lDistLow, &lDistHigh, FILE_BEGIN);
        ReadFile(hFile, static_cast<PVOID>(&varg), sizeof(double), &br, nullptr);
        return varg;
    }
    // Otherwise, use a buffer in memory
    else {
        return argument[index];
    }
}

///////////////////////////////////////////////////////////////////////////////
// Function to read from a file contiguous array of samples in buffer
int DataFileIO::ReadMarks(__int64 pos, int length, double *buf)
{
    if (FileFormat == plain) {
        // How many bytes are occupied by one mark - size double * the number of records
        unsigned int marklenbt = sizeof(double) * (SeriesCount + 1);

        // Offset label pos number of bytes from the beginning of the data
        __int64 offset = pos * marklenbt;

        // Setting a pointer to the first count of the required data in the file
        LONG lDistLow = static_cast<LONG>(fDataPtr + offset);
        LONG lDistHigh = static_cast<LONG>((fDataPtr + offset) >> 32);
        SetFilePointer(hFile, lDistLow, &lDistHigh, FILE_BEGIN);

        //How to read - block - blocklenbt - block length in bytes
        unsigned int blocklenbt = length * marklenbt;

        // The number of bytes read
        DWORD br = 0;

        // reading data
        int r = ReadFile(hFile, static_cast<PVOID>(buf), blocklenbt, &br, 0);

        // Processing read errors
        if (br < blocklenbt) {
            ShowErrorF(L"DataFileIO::ReadMarks() - Error ReadFile()");
            return -1;
        }
        return 0;
    }
    if (FileFormat == rlezip) {
        double t1 = xmin + pos * ArgStep;
        double t2 = xmin + (pos + length - 1) * ArgStep;

        // The index of the first mark in the file
        int i1 = BinSearch(t1);
        // The index of the last mark in the file
        int i2 = BinSearch(t2);

#ifdef _USE_READ_FILE_DUMP_
        // dumps
        if (1) {
            ofstream tmpf("dump.txt", ios_base::app);
            tmpf << "ReadMarks() for " << pos << ", " << length << endl;
            tmpf << "t1 = " << t1 << ", i1 = " << i1 << endl;
            tmpf << "t2 = " << t2 << ", i2 = " << i2 << endl;
            tmpf.close();
        }
#endif
        /// Between them, read from a file into a buffer

        //How many bytes are occupied by one mark - size double * the number of records
        unsigned int marklenbt = sizeof(double) * (SeriesCount + 1);

        // Offset label pos number of bytes from the beginning of the data
        __int64 offset = i1 * marklenbt;

        // Setting a pointer to the first count of the required data in the file
        LONG lDistLow = static_cast<LONG>(fDataPtr + offset);
        LONG lDistHigh = static_cast<LONG>((fDataPtr + offset) >> 32);
        SetFilePointer(hFile, lDistLow, &lDistHigh, FILE_BEGIN);

        // How to read - block - blocklenbt - block length in bytes
        unsigned int blocklenbt = (i2 - i1 + 1) * marklenbt;

        // Temporary buffer for marks
        char *tempbuf = new char[blocklenbt];

        // The number of bytes read
        DWORD br = 0;

        // Reading data into the buffer marks
        int r = ReadFile(hFile, static_cast<PVOID>(tempbuf), blocklenbt, &br, nullptr);

        // Processing read errors
        if (br < blocklenbt) {
            ShowErrorF(L"DataFileIO::ReadMarks() - Error ReadFile()");
            delete[] tempbuf;
            return -1;
        }

        // Pointer to the data
        double *data = (double *)(tempbuf);

        // Sampling and filling the array output marks
        for (int i = 0; i < length; i++) {
            double t = xmin + (pos + i) * ArgStep;
            int index = BinSearch(t);
            buf[i * (SeriesCount + 1)] = ReadArg(index);
            for (int j = 0; j < SeriesCount; j++) {
                buf[i * (SeriesCount + 1) + j + 1] = data[(index - i1) * (SeriesCount + 1) + j + 1];
            }
        }
        delete[] tempbuf;
        return 0;
    }

    // Returns an error flag
    return -1;
}

///////////////////////////////////////////////////////////////////////////////
// Designer - opens the file and read the title, along with the necessary information
//-----------------------------------------------------------------------------
DataFileIO::DataFileIO(const wchar_t *fname)
{
    // The initial setting values
    xmin = xmax = -1;
    Error = 0;
    NameBuf = 0;
    fDataPtr = 0;
    NumPoint = -1;
    SeriesCount = -1;
    bLargeFile = false;

    hFile = CreateFileW(fname, GENERIC_READ, FILE_SHARE_READ, 0, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, 0);
    if (hFile != INVALID_HANDLE_VALUE) {
        DWORD dwFileSizeHigh;
        __int64 qwFileSize = GetFileSize(hFile, &dwFileSizeHigh);
        qwFileSize += (((__int64)dwFileSizeHigh) << 32);

        // Variable to store the number of bytes read
        DWORD br = 0;
        // Buffer to store the signature file
        char signature[8] = {0};
        // status variable operation
        int st = 0;
        // Reading the signature file
        st = ReadFile(hFile, static_cast<PVOID>(signature), 8, &br, nullptr);
        if (st == 0) {
            // Error reading from file signatures
            CloseHandle(hFile);
            ShowErrorF(L"DataFileIO::DataFileIO - Error read signal signature!");
            Error = -1;
            return;
        }

        // Setting the file pointer at the beginning
        SetFilePointer(hFile, 0, nullptr, FILE_BEGIN);

        // Reading the header file
        st = ReadFile(hFile, static_cast<PVOID>(&header), sizeof(SSignal), &br, nullptr);

        if (header.NumDataName > 0) {
            SeriesCount = header.NumDataName;

            NameBuf = new char *[SeriesCount];
            for (unsigned int i = 0; i < SeriesCount; i++) {
                NameBuf[i] = new char[40];
                ReadFile(hFile, NameBuf[i], 40, &br, nullptr);
            }
        } else { // ( header.NumDataName == 0 )
            SeriesCount = 1;
            NameBuf = nullptr;
        }

        fDataPtr = sizeof(SSignal) + header.NumDataName * 40;

        NumPoint = (qwFileSize - fDataPtr) / (SeriesCount + 1) / sizeof(double);

        if (memcmp(signature, "Digizip", 7) == 0) {
            // Large file - this is when a buffer for arguments greater,
            // than 8 * N ~ 400Mb => N ~ 50 million.
            if (NumPoint > 20000000) {
                bLargeFile = true;
                argument = new double[0];
            } else {
                bLargeFile = false;

                // Allocate memory for the buffer argument
                try {
                    argument = new double[NumPoint];
                } catch (exception *e) {
                    //throw std::invalid_argument( "Error allocating memory for arguments" );
                }

                // Download values of the argument in the buffer memory
                for (int k = 0; k < NumPoint; k++) {
                    SetFilePointer(hFile, fDataPtr + (k * (SeriesCount + 1))* sizeof(double), 0, FILE_BEGIN);
                    ReadFile(hFile, static_cast<PVOID>(&(argument[k])), sizeof(double), &br, nullptr);
                }
            }

            // Defining the limits of the argument
            if (NumPoint > 0) {
                xmin = ReadArg(0);
                xmax = ReadArg(NumPoint - 1);
            } else {
                xmin = xmax = 0;
            }

            // Step of the argument
            ArgStep = header.period;

            // The number of points in a compressed file
            zipPoint = NumPoint;

            // The true number of points (like after uncompressing)
            NumPoint = static_cast<int>((xmax - xmin) / ArgStep);

            // Data file format
            FileFormat = rlezip;
        } else if (memcmp(signature, "Digital", 7) == 0) {
            // Reading the minimum value of the argument
            SetFilePointer(hFile, fDataPtr, nullptr, FILE_BEGIN);
            ReadFile(hFile, static_cast<PVOID>(&xmin), sizeof(double), &br, nullptr);

            // Getting resampling
            // Reading 2nd discrete value to determine the hours you-sampling
            SetFilePointer(hFile, fDataPtr + (SeriesCount + 1)* sizeof(double), nullptr, FILE_BEGIN);
            ReadFile(hFile, static_cast<PVOID>(&ArgStep), sizeof(double), &br, nullptr);

            if (xmin >= 0)
                ArgStep -= xmin;

            if (ArgStep <= 0.0) {
                //MessageBoxA( nullptr, "The error in the calculation of resampling dt", "DataFileIO::LoadInfo()", MB_OK | MB_ICONERROR );
                //throw std::invalid_argument( "The error in the calculation of resampling dt" );
            }

            // Determination of the maximum value of the argument
            xmax = ArgStep * NumPoint;

            // Data file format
            FileFormat = plain;
        } else {
            MessageBoxA(nullptr, "Unknown file format", "format error", 0x30);
            Error = -1;
        }
    } else {
        // Error opening file
        //ShowErrorF(L"Error opening file.");
        Error = -1;
    }
}

///////////////////////////////////////////////////////////////////////////////

DataFileIO::~DataFileIO(void)
{
    // Close the file
    CloseHandle(hFile);

    // Exemption dynamic structures - line titles series
    if (NameBuf != nullptr) {
        for (int i = 0; i < header.NumDataName; i++) {
            delete[] NameBuf[i];
        }
        delete[] NameBuf;
    }

    // Deleting an array of arguments (if any) for the compressed format
    if (FileFormat == rlezip) {
        delete[] argument;
    }
}

///////////////////////////////////////////////////////////////////////////////
// Getting information about the signals in the file - not supported
int DataFileIO::GetDataInfo(void)
{
    return header.TypeSignal;
}

///////////////////////////////////////////////////////////////////////////////
// For all series is the minimum and maximum values
// on the interval from pos to pos + length and written in MinY and MaxY
int DataFileIO::LoadMinMaxSignal(__int64 pos, unsigned int length,
                                 double *MinY, double *MaxY)
{
    // The input is the initial reference index, the length of the scanning sequence
    // and arrays of minima and maxima for each series

    // The empirical index of the last frame
    __int64 lastn = pos + length;

    // Checking for a valid argument
    // If the area is completely right or completely to the right domain
    if (((pos + length) <= 0) || (pos >= NumPoint)) {
        // Outside the domain of definition - no signal
        for (unsigned int j = 0; j < SeriesCount; j++) {
            // Resetting the extreme
            MinY[j] = MaxY[j] = 0;
        }
        // Shutdown function
        return 0;
    }

    // Correction parameters if the output of the beginning of the file
    if (pos < 0) {
        // Search extremum will be a zero element
        pos = 0;
        // And the length of the scanning sequence will be less
        length = static_cast<unsigned int>(length + pos);
    }

    // Correction parameters if the output of the end of the file
    else if (lastn > NumPoint) {
        // Shave points beyond file
        length = static_cast<unsigned int>(NumPoint - pos);
    }

    //-------------------------------------------------------------------------
    if (FileFormat == plain) {
        // Allocate memory for the buffer marks
        double *localbuf = new double[length * (SeriesCount + 1)];

        // Memory allocation error - not my fault if you ask too much
        if (localbuf == nullptr) {
            ShowErrorF(L"DataFileIO::LoadMinMaxSignal() - Error create localbuf!");
            return -1;
        }

        // Reading marks packs a local buffer
        if (ReadMarks(pos, length, localbuf) != 0) {
            ShowErrorF(L"DataFileIO::LoadMinMaxSignal() - Error in ReadMarks()");
            delete[] localbuf;
            return -1;
        }

        // Initialization of the highs and lows of the first value of the buffer
        for (unsigned int j = 0; j < SeriesCount; j++) {
            MinY[j] = MaxY[j] = localbuf[j + 1];
        }

        // Cycle view of values and determining the extrema
        for (unsigned int i = 0; i < length; i++) {
            for (unsigned int j = 0; j < SeriesCount; j++) {
                double val = localbuf[i * (SeriesCount + 1) + j + 1];
                if (val > MaxY[j]) MaxY[j] = val;
                if (val < MinY[j]) MinY[j] = val;
            }
        }

        // Exemption taken memory
        delete[] localbuf;

        return 0;
    }

    if (FileFormat == rlezip) {
        /// Enough to read the marks on t1 (pos) to t2 (pos + length)
        /// and among them seek min and max
        double t1 = xmin + pos * ArgStep;
        double t2 = xmin + (pos + length - 1) * ArgStep;

        __int64 index1 = BinSearch(t1);
        __int64 index2 = BinSearch(t2);

#ifdef _USE_READ_FILE_DUMP_
        // damps
        if (1) {
            ofstream tmpf("dump.txt", ios_base::app);
            tmpf << "LoadMinMaxSignal() for " << pos << " .. " << pos + length - 1 << endl;
            tmpf << "t1 = " << t1 << ", i1 = " << index1 << endl;
            tmpf << "t2 = " << t2 << ", i2 = " << index2 << endl;
            tmpf.close();
        }
#endif

        /// Between them, read from a file into a buffer

        // How many bytes are occupied by one mark - size double * the number of records
        unsigned int marklenbt = sizeof(double) * (SeriesCount + 1);

        // Offset label pos number of bytes from the beginning of the data
        __int64 offset = index1 * marklenbt;

        // Setting a pointer to the first count of the required data in the file
        LONG lDistLow = static_cast<LONG>(fDataPtr + offset);
        LONG lDistHigh = static_cast<LONG>((fDataPtr + offset) >> 32);
        SetFilePointer(hFile, lDistLow, &lDistHigh, FILE_BEGIN);

        // How to read - block - blocklenbt - block length in bytes
        unsigned int blocklenbt = (index2 - index1 + 1) * marklenbt;

        // Temporary buffer for marks
        char *tempbuf = new char[blocklenbt];

        // The number of bytes read
        DWORD br = 0;

        // reading data
        int r = ReadFile(hFile, static_cast<PVOID>(tempbuf), blocklenbt, &br, nullptr);

        // Processing read errors
        if (br < blocklenbt) {
            ShowErrorF(L"DataFileIO::ReadMarks() - Error ReadFile()");
            delete[] tempbuf;
            return -1;
        }

        double *data = (double *)(tempbuf);

        // Initialization of the highs and lows of the first value of the buffer
        for (unsigned int j = 0; j < SeriesCount; j++) {
            MinY[j] = MaxY[j] = data[j + 1];
        }

        for (int i = 0; i < (index2 - index1 + 1); i++) {
            for (unsigned int j = 0; j < SeriesCount; j++) {
                double val = data[i * (SeriesCount + 1) + j + 1];
                if (val > MaxY[j]) MaxY[j] = val;
                if (val < MinY[j]) MinY[j] = val;
            }
        }
        // Freeing memory
        delete[] tempbuf;
        return 0;
    }

    return -1;
}

///////////////////////////////////////////////////////////////////////////////
// Series Index number is the minimum and maximum values
// n the interval from pos to pos + length and written in MinY and MaxY
int DataFileIO::GetMinMaxByIndex(unsigned int Index,
                                 __int64 pos, unsigned int length,
                                 double &MinY, double &MaxY)
{
    // Check from the fool
    if (Index >= SeriesCount) {
        ShowErrorF(L"DataFileIO::GetMinMaxByIndex() - Error bad series index!");
        return -1;
    }

    // Memory buffers for extreme
    double *miny = new double[SeriesCount];
    double *maxy = new double[SeriesCount];

    LoadMinMaxSignal(pos, length, miny, maxy);

    // Filling needed for the selected series
    MinY = miny[Index];
    MaxY = maxy[Index];

    // Freeing memory
    delete[] miny;
    delete[] maxy;

    return 0;
}

///////////////////////////////////////////////////////////////////////////////
// Download a ring into the buffer - all given samples
// Reasonably necessary to specify the number of samples to read - memory is limited
int DataFileIO::LoadSignal(unsigned int Index,
                           __int64 pos, unsigned int length,
                           double *Y, double *x)
{
    // Check from the fool
    if (Index >= SeriesCount) {
        ShowErrorF(L"DataFileIO::LoadSignal() - Error bad series index!");
        return -1;
    }

    // Checking for a valid argument
    if (pos >= NumPoint) {
        // If you are asked to download the file value is
        //ShowErrorF(L"DataFileIO::LoadSignal() - Error bad start position!");
        return 0;
    }

    // Actually read long - what to take from the file
    unsigned int readlen = length;

    // Correction parameters if the output of the end of the file
    if ((length + pos) > NumPoint) {
        readlen = static_cast<unsigned int>(NumPoint - pos);
    }

    // Allocate memory for the buffer for one mark
    double *localbuf = new double[SeriesCount + 1];

    // Erase Cycles values
    for (unsigned int i = 0; i < readlen; i++) {
        // reading marks
        ReadMarks(pos + i, 1, localbuf);
        // A copy of the desired value
        Y[i] = localbuf[Index + 1];
        // And the argument if ordered
        if (x != nullptr) {
            x[i] = localbuf[0];
        }
    }

    // release the memory taken
    delete[] localbuf;

    // Do not fall for the values are reset
    for (unsigned int i = readlen; i < length; i++) {
        Y[i] = 0;
        if (x != nullptr) {
            x[i] = 0;
        }
    }

    return 0;
}

///////////////////////////////////////////////////////////////////////////////
// Preparation time file - creation, modification, access
FILETIME DataFileIO::GetCreationTime(void)
{
    FILETIME ft;
    GetFileTime(hFile, &ft, nullptr, nullptr);
    return FILETIME(ft);
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// Constructor record in the file
DataFileIOWriter::DataFileIOWriter(const wchar_t *fname, unsigned int nSeries)
{
    Error = 0;
    // Creating a File
    hFile = CreateFileW(fname,
                        GENERIC_WRITE,
                        0,
                        nullptr,
                        CREATE_ALWAYS,
                        FILE_ATTRIBUTE_NORMAL,
                        nullptr);

    // Processing errors opening
    if (hFile == INVALID_HANDLE_VALUE) {
        ShowErrorF(L"DataFileIOWriter::DataFileIOWriter() - Error CreateFile()!");
        Error = -1;
        return;
    }

    SSignal header;
    DWORD br = 0;

    ZeroMemory(&header, sizeof(SSignal));
    strcpy((char *)&header, "Digital");
    // The fact that the series does not mean there is one thing they
    header.NumDataName = nSeries;
    if (nSeries == 0) {
        SeriesCount = 1;
    } else {
        SeriesCount = nSeries;
    }

    // header record
    WriteFile(hFile, &header, sizeof(SSignal), &br, nullptr);

    // Record series names (yet empty)
    char ze[40] = {0};
    for (unsigned int i = 0; i < SeriesCount; i++) {
        wsprintfA(ze, "Preview_Series_%2u", i);
        WriteFile(hFile, ze, 40, &br, nullptr);
    }
}
// Constructor record in the file
DataFileIOWriter::DataFileIOWriter(const wchar_t *fname, vector<string> nSeries)
{
    Error = 0;
    // Creating a File
    hFile = CreateFileW(fname,
                        GENERIC_WRITE,
                        0,
                        NULL,
                        CREATE_ALWAYS,
                        FILE_ATTRIBUTE_NORMAL,
                        NULL);

    // Processing errors opening
    if (hFile == INVALID_HANDLE_VALUE) {
        //      ShowErrorF("DataFileIOWriter::DataFileIOWriter() - Error CreateFile()!");
        Error = -1;
        return;
    }

    SSignal header;
    DWORD br = 0;

    ZeroMemory(&header, sizeof(SSignal));
    strcpy((char *)&header, "Digital");
    // The fact that the series does not mean there is one thing they
    header.NumDataName = nSeries.size();
    if (nSeries.size() == 0) {
        SeriesCount = 1;
    } else {
        SeriesCount = nSeries.size();
    }

    // header record
    WriteFile(hFile, &header, sizeof(SSignal), &br, NULL);

    // Record series names (yet empty)
    char ze[40] = { 0 };
    for (unsigned int i = 0; i < SeriesCount; i++) {
        //wsprintfA(ze, "Preview_Series_%2u", i);
        wsprintfA(ze, nSeries[i].c_str(), i);
        WriteFile(hFile, ze, 40, &br, NULL);
    }
}
///////////////////////////////////////////////////////////////////////////////
// Exemption object entries in the file
DataFileIOWriter::~DataFileIOWriter()
{
    CloseHandle(hFile);
}

///////////////////////////////////////////////////////////////////////////////
// Record mark in the file
int DataFileIOWriter::WriteMark(double t, double *series)
{
    if (Error != 0) return -1;
    DWORD br = 0;
    WriteFile(hFile, &t, sizeof(double), &br, nullptr);
    WriteFile(hFile, series, SeriesCount * sizeof(double), &br, nullptr);
    return 0;
}

void dataIO::DataFileIOWriter::Close()
{
    CloseHandle(hFile);
}
