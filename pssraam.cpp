#include "pssraam.h"
#include <chrono>
#include <QtMath>

#include <iostream>

using namespace std;

////////////////////////////////////////////////////////////////////////////////
// Class constructor  !!! НУЖНО ПЕРЕДАВАТЬ СТРУКТУРУ, ЧТО БЫ:
//  1.  ИЗБЕЖАТЬ ЦИКЛИЧЕСКОГО ПЕРЕСЧЕТА КОНСТАНТ В МЕТОДЕ МОДЕЛИРОВАНИЯ
//  2.  ПРОИНИЦИАЛИЗИРОВАТЬСЯ И КАК СЛЕДСТВИЕ ... см. п. 3
//  3.  НЕ ВЫЛЕТЕТЬ, ЕСЛИ ПЕРВЫЙ ВЫЗОВ НЕ СОВПАДЕТ С НАЧАЛОМ ПАЧКИ И ПОЛЯ НЕ БУДУТ ПРОИНИЦИАЛИЗИРОВАНЫ
////////////////////////////////////////////////////////////////////////////////

Pssraam::Pssraam(const inData &iData)
{
    // Auxiliary variables for preliminary calculations initialization
    deleps_sl = iData.deleps_sld * toRad;
    double delTet_a_e = iData.delTet_a_ed * toRad;
    delTet_a_e_div2 = delTet_a_e / 2;
    N_sig = qRound64(iData.T_sig / iData.delt);

    //Constant signal power parameters calculation and initialization
    // A_s_SL         константы - оптимизация для исключения из циклов пересчета
    // A_s_ML;        не деленные на f_0^2, т.к. f_0 может изменяться ??????
    double delTet_a_b = iData.delTet_a_bd * toRad;
    A_s_ML = (iData.P_0 * iData.G_T * iData.G_R * iData.a_S * delTet_a_b * _c *
              _c) / (64 * M_PI * M_PI * M_PI);
    A_s_SL = qPow(10, -iData.q_T / 10.0) * qPow(10, -iData.q_R / 10.0) *
             deleps_sl * A_s_ML;

    //Random seeding contaners for a not reproducible random sequence
    // unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

    //Constant seeding contaners for a reproducible random sequence
    unsigned seed = 2;

    engine_0_1.seed(seed);
    engine_0_05.seed(seed);
    engine_0_D_nML_XY.seed(seed);
    engine_0_D_nSL_XY.seed(seed);

    // Initialization random sequence generators components with required expectation and deviation
    distribution_0_1 = normal_distribution<double>(0, 1);
    distribution_0_05 = normal_distribution<double>(0, 0.5);
    //distribution_0_D_nML_XY = normal_distribution<double>(0, D_nML_XY); (пере-)инициализируются в recalculate_() каждый раз, в зависи-
    //distribution_0_D_nSL_XY = normal_distribution<double>(0, D_nSL_XY); мости от текущих значений D_nML_XY и D_nSL_XY

    // Class's fields initialization in accordance current input data
    preCalculate_ (iData);
//    cout << "N_sl_min = " << N_sl_min << endl;
//    cout << "N_sl_max = " << N_sl_max << endl;
//    cout << "N_ml_min = " << N_ml_min << endl;
//    cout << "N_ml_max = " << N_ml_max << endl;
}

// Passive Signal Surfase Reflected Air-Air Modeling method (quadrature components computing)
complex <double> Pssraam::calculate(inData &iData)
{
    if(iData.isDataRenewed == true)
    {
        preCalculate_ (iData);
        iData.isDataRenewed = false; // Reset flag of renewed input data
                                     // after reading and recalculation
                                     // Лучше, если это будет делать программа
                                     // сценария самостоятельно (устанавливать
                                     // и сбрасывать этот флаг), что бы исключить
                                     // случайный доступ по записи к элементам
                                     // структуры вх. данных из класса Pssraam(const inData &iData)
    }

    // Modelling quadrature components initialization
    X_A = 0;
    Y_A = 0;
    X_ML = 0;
    Y_ML = 0;
    X_SL = 0;
    Y_SL = 0;

    int m = iData.takt % N_s;
    if( m == 0)
    {
        setProcess_ (proc_DML, distribution_0_05, engine_0_05);
        setProcess_ (proc_DSL, distribution_0_05, engine_0_05);
        if (iData.isSurfaceMove && iData.V_W_mean)
        {
            setProcess_ (proc_nML, distribution_0_D_nML_XY, engine_0_D_nML_XY);
            setProcess_ (proc_nSL, distribution_0_D_nSL_XY, engine_0_D_nSL_XY);
        }
    }
    else
    {
        // Signal imitation for Side Lobe
        if(m >= N_sl_min && m < N_sl_max)
        {
            renewProcess_ (proc_DSL, distribution_0_1, engine_0_1,  filtr_DSL);
            if (iData.isSurfaceMove && iData.V_W_mean)
            {
                renewProcess_ (proc_nSL, distribution_0_1, engine_0_1,  filtr_nSL);
            }
            else
            {
                proc_nSL.X = qSqrt(D_nSL_XY);
                proc_nSL.Y = proc_nSL.X;
            }
            X_SL = proc_DSL.X * proc_nSL.X - proc_DSL.Y * proc_nSL.Y;
            Y_SL = proc_DSL.X * proc_nSL.Y + proc_DSL.Y * proc_nSL.X;
        }
        else
        {
            // Signal imitation for Main Lobe
            if(m >= N_ml_min && m <= N_ml_max)
            {
                renewProcess_ (proc_DML, distribution_0_1, engine_0_1,  filtr_DML);
                if (iData.isSurfaceMove && iData.V_W_mean)
                {
                    renewProcess_ (proc_nML, distribution_0_1, engine_0_1,  filtr_nML);
                }
                else
                {
                    proc_nML.X = qSqrt(D_nML_XY);
                    proc_nML.Y = proc_nML.X;
                }
                X_ML = proc_DML.X * proc_nML.X - proc_DML.Y * proc_nML.Y;
                Y_ML = proc_DML.X * proc_nML.Y + proc_DML.Y * proc_nML.X;
            }
            else    // Out of receive range, no imitation
            {
                X_SL = 0;
                Y_SL = 0;
                X_ML = 0;
                Y_ML = 0;
            }
        }
    }
    X_A = X_SL + X_ML;
    Y_A = Y_SL + Y_ML;

    complex <double> q;
    q.real(X_A);
    q.imag(Y_A);

    return q;
}

// The forming filtr's coefficients calculation method
inline void Pssraam::setFiltr_ (Filtr &filtr, double a_xxL, double deltaT, double D)
{
    double gamma = a_xxL * deltaT;
    double p = qExp( -gamma );
    double a_0 = p * p * p * (1 + gamma) - p * (1 - gamma);
    double a_1 = 1.0 - 4 * p * p * gamma - p * p * p * p;
    filtr.a0 = qSqrt(D) * qSqrt((a_1 * a_1 + qSqrt(a_1 * a_1 - 4 * a_0 * a_0)) / 2);
    filtr.a1 = a_0 * qSqrt(D) / a_1;
    filtr.b1 = 2 * p;
    filtr.b2 = - p * p;
}

// Recalculation input and internal data
void Pssraam::preCalculate_ (const inData &iData)
{
    double b_f_V    = iData.b_f_Vd * toRad;
    double eps_f_V  = iData.eps_f_Vd * toRad;
    e_f_V.setX(qCos(eps_f_V) * qCos(b_f_V));
    e_f_V.setY(qSin(eps_f_V));
    e_f_V.setZ(qCos(eps_f_V) * qSin(b_f_V));

    // From degrees to radians
    double eps_f = iData.eps_f_d * toRad;

    // Flight altitude
    int H_f = qRound(iData.D_f * qSin(eps_f));

    // Recalculation coordinates Unit Vectors
    double b_a    = iData.b_ad * toRad;
    double eps_a  = iData.eps_ad * toRad;
    double delTet_a_e = iData.delTet_a_ed * toRad;
    delTet_a_e_div2 = delTet_a_e / 2;  // optimizing subsequent division
    QVector3D e_l;
    e_l.setX(qCos(eps_a - delTet_a_e_div2) * qCos(b_a));
    e_l.setY(qSin(eps_a - delTet_a_e_div2));
    e_l.setZ(qCos(eps_a - delTet_a_e_div2) * qSin(b_a));
    QVector3D e_h;
    e_h.setX(qCos(eps_a + delTet_a_e_div2) * qCos(b_a));
    e_h.setY(qSin(eps_a + delTet_a_e_div2));
    e_h.setZ(qCos(eps_a + delTet_a_e_div2) * qSin(b_a));

    e_1_l.setX(e_l.x() * qCos(b_f_V) + e_l.z() * qSin(b_f_V));
    e_1_l.setY(e_l.y());
    e_1_l.setZ(-e_l.x() * qSin(b_f_V) + e_l.z() * qCos(b_f_V));

    e_1_h.setX(e_h.x() * qCos(b_f_V) + e_h.z() * qSin(b_f_V));
    e_1_h.setY(e_h.y());
    e_1_h.setZ(-e_h.x() * qSin(b_f_V) + e_h.z() * qCos(b_f_V));


    QVector3D e_sl;
    e_sl.setX(qCos(eps_a - delTet_a_e_div2 - deleps_sl) * qCos(b_a));
    e_sl.setY(qSin(eps_a - delTet_a_e_div2 - deleps_sl));
    e_sl.setZ(qCos(eps_a - delTet_a_e_div2 - deleps_sl) * qSin(b_a));

    e_1_sl.setX(e_sl.x() * qCos(b_f_V) + e_sl.z() * qSin(b_f_V));
    e_1_sl.setY(e_sl.y());
    e_1_sl.setZ(-e_sl.x() * qSin(b_f_V) + e_sl.z() * qCos(b_f_V));

    // Dopler's frequency offsets calculation
    double f_h_D;       //  for main lobe upper edge
    f_h_D = div_2_c * iData.f_0 * iData.V_f_m * (e_f_V.x()*e_1_h.x() + e_f_V.y() *
            e_1_h.y() + e_f_V.z() * e_1_h.z());
    double f_l_D;       //  for main lobe lower edge
    f_l_D = div_2_c * iData.f_0 * iData.V_f_m * (e_f_V.x()*e_1_l.x() + e_f_V.y() *
            e_1_l.y() + e_f_V.z() * e_1_l.z());
    double f_sl_D;      //  for side lobe left edge
    f_sl_D = div_2_c * iData.f_0 * iData.V_f_m * (e_f_V.x()*e_1_sl.x() + e_f_V.y() *
             e_1_sl.y() + e_f_V.z() * e_1_sl.z());

    // Range of Dopler's frequenciey calculation
    double delf_ML_D;   //  with main lobe received signal
    delf_ML_D = qAbs(f_h_D - f_l_D);
    double delf_SL_D;   //  with side lobe received signal
    delf_SL_D = qAbs(f_l_D - f_sl_D);

    // Main lobe receive time within one period of probing signal
    double t_ml_min;    // start receiving
    t_ml_min = div_2_c * H_f / qSin(qAbs(eps_a) + delTet_a_e_div2);
    double t_ml_max;    // end receiving
    t_ml_max = div_2_c * H_f / qSin(qAbs(eps_a) - delTet_a_e_div2);
    if(t_ml_max > iData.T_r)
        t_ml_max = iData.T_r;

    // Side lobe receive time within one period of probing signal
    double t_sl_min;    // start receiving
    t_sl_min = div_2_c * H_f;
    double t_sl_x;
    t_sl_x = t_sl_min / qSin(deleps_sl + qAbs(eps_a) + delTet_a_e_div2);
    double t_sl_max;    // end receiving
    t_sl_max = qMax(t_sl_x, t_ml_min);
    if(t_sl_max > iData.T_r)      // !!!!!!!!!?????? САМОДЕЯТЕЛЬНОСТЬ см. ср. 157
        t_sl_max = iData.T_r;     // !!!!!!!!!?????? САМОДЕЯТЕЛЬНОСТЬ см. ср. 157

    // Tick parameters for SL and ML calculation
    N_sl_min = qRound(t_sl_min / iData.delt);
    N_sl_max = qRound(t_sl_max / iData.delt);
    N_ml_min = qRound(t_ml_min / iData.delt);
    N_ml_max = qRound(t_ml_max / iData.delt);

    // Main lobe receive signal power parameters calculation
    double r_Gmin;
    r_Gmin = H_f / qTan(qAbs(eps_a) + delTet_a_e_div2);
    double r_Gmax;
    r_Gmax = H_f / qTan(qAbs(eps_a) - delTet_a_e_div2);
    double sig_2_nML;
    sig_2_nML = iData.T_0 * A_s_ML * (1 / (r_Gmin * r_Gmin) - 1 / (r_Gmax * r_Gmax)) /
                (iData.T_r * 2 * iData.f_0 * iData.f_0);

    // Computing the signal quadrature components dispersion
    // reflected from Earth surface and received with antenna's main lobes
    D_nML_XY = sig_2_nML / 2;

    double delf_nML;
    double delf_nSL;
    delf_nML = 3.936E-2 * iData.f_0 * qPow(iData.V_W_mean, 1.3) / _c;
    delf_nSL = delf_nML;
    double tay_nML;
    if(delf_nML != 0)
    {
        tay_nML = 1 / (2.0 * delf_nML);
        a_nSL = 2.0 / tay_nML;
        a_nML = a_nSL;
    }
    double tay_DML;

    if(delf_ML_D != 0)      // В ПРОТИВНОМ СЛУЧАЕ КОФФ. БУДЕТ НЕ РАССЧИТАН!!!!!!!! ПРЕДУСМОТРЕТЬ!!!!!
    {
        tay_DML = 1 / (2 * delf_ML_D);
        a_DML = 2 / tay_DML;
    }

    // Side lobe receive signal power parameters calculation
    double sig_2_nSL;
    sig_2_nSL = iData.T_0 * A_s_SL / (iData.T_r * 2 * H_f * H_f * iData.f_0 * iData.f_0);

    // Computing the signal quadrature components dispersion
    // reflected from Earth surface and received with antenna's side lobes
    D_nSL_XY = sig_2_nSL / 2;

    double tay_DSL;
    if(delf_SL_D != 0)      // В ПРОТИВНОМ СЛУЧАЕ КОФФ. БУДЕТ НЕ РАССЧИТАН!!!!!!!! ПРЕДУСМОТРЕТЬ!!!!!
    {
        tay_DSL = 1.0 / (2.0 * delf_SL_D);
        a_DSL = 2.0 / tay_DSL;
    }

    // Tick parameters for SL and ML calculation
    N_s = qRound64(iData.T_r / iData.delt);

    // Initialization random sequence generators components with required
    // expectation and deviation
    distribution_0_D_nML_XY = normal_distribution<double>(0, D_nML_XY);
    distribution_0_D_nSL_XY = normal_distribution<double>(0, D_nSL_XY);

    if(delf_nML != 0)
    {
        setFiltr_(filtr_nML, a_nML, iData.delt, D_nML_XY);
        setFiltr_(filtr_nSL, a_nSL, iData.delt, D_nSL_XY);
    }
    setFiltr_(filtr_DML, a_DML, iData.delt, 0.5);
    setFiltr_(filtr_DSL, a_DSL, iData.delt, 0.5);
}

// Set quadrature components modelling processes
inline void Pssraam::setProcess_ (Process &proc, normal_distribution<double>
                                                    distribution, mt19937 &engine)
{
    proc.nx = distribution_0_1(engine_0_1);
    proc.nx_1 = distribution_0_1(engine_0_1);
    proc.ny = distribution_0_1(engine_0_1);
    proc.ny_1 = distribution_0_1(engine_0_1);
    proc.X = 0;
    proc.X_1 = distribution(engine);
    proc.X_2 = distribution(engine);
    proc.Y = 0;
    proc.Y_1 = distribution(engine);
    proc.Y_2 = distribution(engine);   
}

// Renew quadrature components modelling processes
inline void Pssraam::renewProcess_ (Process &proc, normal_distribution<double>
                                     distribution, mt19937 &engine,  Filtr filtr)
{
    proc.X = filtr.a0 * proc.nx + filtr.a1 * proc.nx_1 + filtr.b1 * proc.X_1 +
             filtr.b2 * proc.X_2;
    proc.Y = filtr.a0 * proc.ny + filtr.a1 * proc.ny_1 + filtr.b1 * proc.Y_1 +
             filtr.b2 * proc.Y_2;
    proc.nx_1 = proc.nx;
    proc.ny_1 = proc.ny;
    proc.nx = distribution(engine);
    proc.ny = distribution(engine);
    proc.X_2 = proc.X_1;
    proc.X_1 = proc.X;
    proc.Y_2 = proc.Y_1;
    proc.Y_1 = proc.Y;
}
