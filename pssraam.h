// PSSRAAM - Passive Signal Surface Reflected Air-Air Modelling

#ifndef PSSRAAM_H
#define PSSRAAM_H

#include "common.h"
#include <complex.h>
#include <QVector3D>
#include <random>
#include <QRandomGenerator>

using namespace std;

// Constants defining
#define _c      299792458                          // light speed definition
#define toRad   0.01745329251994329576923690768489 // from degree to radians convert (M_PI / 180)
#define div_2_c 6.6712819039630409915115342894984e-9 // 2 / light speed definition (2/с) - оптимизационная константа для увеличени быстродействия

// Own used type to unbound from the name of input data external structure
typedef STRCT inData;

// Own used input data names to unbound from field's names of input data external structure
#define f_0         param1     // frequency of the probing signal
#define T_0         param2     // duration of radar probe pulses
#define T_r         param3     // radar probing pulses period in the pack
#define T_s         param4     // sequence period of radar probing pulses pack
#define delt        param5     // sampling interval of the modelling signal in time
#define T_sig       param6     // reflected signal total modelling time
#define delTet_a_ed param7     // radar's main lobe width in the angular plane
#define delTet_a_bd param8     // radar's main lobe width in the azimuthal plane
#define G_T         param9     // transmission gain  factor of radar's main lobe
#define G_R         param10    // receive gain  factor of radar's main lobe
#define q_T         param11    // transmit attenuation coefficient of radar's side lobe
#define q_R         param12    // receive attenuation coefficient of radar's side lobe
#define P_0         param13    // radar transmitter power
#define deleps_sld  param14    // radar's side lobe width in the angular plane
#define V_f_m       param15    // radar's carrier velocity module
#define b_f_Vd      param16    // azimuth of radar's carrier velocity vector
#define eps_f_Vd    param17    // angle of radar's carrier velocity vector
#define D_f         param18    // radar's carrier range (distance)
#define b_f_d       param19    // radar's carrier azimuth
#define eps_f_d     param20    // radar's carrier angle
#define b_ad        param21    // radar's main lobe azimuth
#define eps_ad      param22    // radar's main lobe angle
#define a_S         param23    // earth's surface specific EPR
#define V_W_mean    param24    // current wind speed 10 m. earth's surface above
#define isSurfaceMove param25    // index of wind movable elements presence (vegetation, etc.)
#define isDataRenewed param26    // input data update flag (ture - update, false - unchanged)
#define takt        param27    // index of trajectory's point for PSSRAA-modeling(external counter)

// Template for forming filter coefficients
typedef struct
{
    double a0 = 0;
    double a1 = 0;
    double b1 = 0;
    double b2 = 0;
} Filtr;

// Template for quadrature components modelling process
typedef struct
{
    double nx = 0;
    double nx_1 = 0;
    double X = 0;
    double X_1 = 0;
    double X_2 = 0;
    double ny = 0;
    double ny_1 = 0;
    double Y = 0;
    double Y_1 = 0;
    double Y_2 = 0;
} Process;

//  PSSRAAM.cs - Passive Signal Surfase Reflected Air-Air Modeling class
class Pssraam
{
public:
    // Class constructor !!! НУЖНО ПЕРЕДАВАТЬ СТРУКТУРУ, ЧТО БЫ ИЗБЕЖАТЬ ПЕРЕСЧЕТА КОНСТАНТ В МЕТОДЕ МОДЕЛИРОВАНИЯ
    Pssraam(const inData &iData);

    // Passive Signal Surfase Reflected Air-Air Modeling method (quadrature components computing)
    complex <double> calculate (inData &iData);

private:
    // Member's fields
    // Containers for Mersenne Twister pseudorandom number generators
    mt19937 engine_0_1;
    mt19937 engine_0_05;
    mt19937 engine_0_D_nML_XY;
    mt19937 engine_0_D_nSL_XY;

    // Fields - the random sequence generators components with required expectation and deviation
    normal_distribution<double> distribution_0_1;
    normal_distribution<double> distribution_0_05;
    normal_distribution<double> distribution_0_D_nML_XY;
    normal_distribution<double> distribution_0_D_nSL_XY;

    // Сoordinates Unit Vectors (UV)
    QVector3D   e_f_V,      // UV along the velocity vector (OcXcYcZc)                                      
                e_1_l,      // UV of main lobe edge in angular plane(OсXсYсZс)
                e_1_h,      // UV of main lobe edge in angular plane(OсXсYсZс)
                e_1_sl;     // UV of side lobe edge in angular plane(OсXсYсZс)

    // Tick parameters for calculating components for SL and ML
    int N_s;
    int N_sig;
    int N_sl_min;
    int N_sl_max;
    int N_ml_min;
    int N_ml_max;

    //Signal power parameters for calculating components for SL and ML
    double A_s_SL;      // константа - оптим. для исключения из циклов пересчета
    double A_s_ML;      // не деленные на f_0^2, т.к. она может меняться ???
    double a_nML;
    double a_nSL;
    double a_DML;
    double a_DSL;
    double D_nML_XY;
    double D_nSL_XY;

    // Auxiliary variables for optimization of calculations (preliminary and general
    double delTet_a_e_div2; //  division operations)
    double deleps_sl;

    // Forming filtrs coefficients fields
    Filtr filtr_nML;
    Filtr filtr_nSL;
    Filtr filtr_DML;
    Filtr filtr_DSL;

    // Quadrature components modelling processes
    Process proc_nML;
    Process proc_DML;
    Process proc_nSL;
    Process proc_DSL;    

    // Quadrature components with Main lobe received signal
    double X_ML;
    double Y_ML;

    // Quadrature components with Side lobe received signal
    double X_SL;
    double Y_SL;

    // Quadrature components total received signal
    double X_A;
    double Y_A;

    // METHODS
    // The forming filtr's coefficients calculation method
    inline void setFiltr_ (Filtr &filtr, double a_xxL, double deltaT, double D);

    // Set quadrature components modelling processes
    inline void setProcess_ (Process &proc, normal_distribution<double> distribution,
                             mt19937 &engine);

    // Renew quadrature components modelling processes
    inline void renewProcess_ (Process &proc, normal_distribution<double> distribution,
                               mt19937 &engine, Filtr filtr);

    // Recalculation input and internal data
    void preCalculate_ (const inData &iData);
};

#endif // PSSRAAM_H
