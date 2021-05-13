#include <QCoreApplication>
#include <iostream>
#include "common.h"
#include "pssraam.h"

#include <QTime>                                                                // УДАЛИТЬ ПО ОТЛАДКЕ
#include <QFile>                                                                // УДАЛИТЬ ПО ОТЛАДКЕ
#include <QTextStream>                                                          // УДАЛИТЬ ПО ОТЛАДКЕ
#include <QTextStream>

//Для автозапуска программы ананлиза
#include <QProcess>                                                             // УДАЛИТЬ ПО ОТЛАДКЕ
//Для программы ананлиза
#include "DataFileIO.h"                                                         // УДАЛИТЬ ПО ОТЛАДКЕ

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    // Удаление предыдущей версии файла анализа - PSSRAAM.sig2d                 // УДАЛИТЬ ПО ОТЛАДКЕ
    QFile fileSig2d("D:/Egor/NAN/QT/build-PSSRAAM-Desktop_Qt_5_12_10_MinGW_32_bit-Debug/PSSRAAM.sig2d");
    if(fileSig2d.exists())                                                      // УДАЛИТЬ ПО ОТЛАДКЕ
    {
        fileSig2d.remove();                                                     // УДАЛИТЬ ПО ОТЛАДКЕ
        cout << "PSSRAAM.sig2d - renewed." << endl;                              // УДАЛИТЬ ПО ОТЛАДКЕ
    }
    else
        cout << "PSSRAAM.sig2d - not exists. Created." << endl;                           // УДАЛИТЬ ПО ОТЛАДКЕ

    // Организация вывода для программы анализа
    QString qstr("PSSRAAM.sig2d");                                              // УДАЛИТЬ ПО ОТЛАДКЕ
    const wchar_t *fileNameDyn = reinterpret_cast<const wchar_t *>(qstr.utf16());   // УДАЛИТЬ ПО ОТЛАДКЕ
    vector<string> series_name;                                             // УДАЛИТЬ ПО ОТЛАДКЕ
    series_name.push_back("R");                                             // УДАЛИТЬ ПО ОТЛАДКЕ
    series_name.push_back("I");                                             // УДАЛИТЬ ПО ОТЛАДКЕ
    dataIO::DataFileIOWriter *_writers = new dataIO::DataFileIOWriter(fileNameDyn, series_name);    // УДАЛИТЬ ПО ОТЛАДКЕ


    QTime time;  //вспом. контроль времени вычисления                           УДАЛИТЬ ПО ОТЛАДКЕ
    time.start();//вспом. контроль времени вычисления                           УДАЛИТЬ ПО ОТЛАДКЕ

    STRCT str;
    Pssraam cl(str);

    // вывод в файл
    QString strX = "D:/Egor/NAN/QT/PSSRAAM/Quadratures_PSSRAAM.txt";
    QFile fileX(strX);
    if (fileX.open(QIODevice::WriteOnly))
    {
        QTextStream writeStream(&fileX);

        complex <double> q;
        for(int i = 0; i <= 2850; ++i)
        {
            str.takt = i;
            q = cl.Model(str);

            // вывод в файл
            if(q.real() != q.imag())
            {
                QString rez;
                rez.sprintf("%14.10f\t%14.10f", q.real(), q.imag());
                writeStream << rez << endl;
            }

            // вывод для программы анализа
            double series[2] = { q.real(), q.imag() };                          //УДАЛИТЬ ПО ОТЛАДКЕ
            _writers->WriteMark(i, series);                                     //УДАЛИТЬ ПО ОТЛАДКЕ
        }
        // вывод в файл
        _writers->Close();                                                      //УДАЛИТЬ ПО ОТЛАДКЕ
        fileX.close();
    }
    else
        cout << "Can't read file." << endl;

    //вспом. контроль времени вычисления                                        УДАЛИТЬ ПО ОТЛАДКЕ
    //cout << "\nModeling time: " << time.elapsed() << " mls." << endl;           //УДАЛИТЬ ПО ОТЛАДКЕ
    printf("%s %07d %s", "\nModeling time: ", time.elapsed(), "mls.");

    //Автозапуск программы ананлиза
    // QProcess process;                                                           //УДАЛИТЬ ПО ОТЛАДКЕ
    // process.start("D:/Egor/NAN/QT/build-PSSRAAM-Desktop_Qt_5_12_10_MinGW_32_bit-Debug/PSSRAAM.bat");
    return a.exec();
}
