/*
Migrated from FORTRAN IV for CDC 6600 system
Orig. author: K. R. Brown et al.
---
Author: Topaz(L99)
Start at 2021/02/17 19:15 BJT
First successful run at ---
Release at --- BJT
*/

#include <iostream>
#include <fstream>
#include <cmath>
#include <complex>
using namespace std;

bool kwargs_verbose = false;

// COMMON/CCPINJ
namespace CCPINJ
{
    double uk;
    int leg;
    double *atp1; 
    short *atp2;
    double am, bmax, cexv;
}

// Subroutine OUT
void OUT(double (&x)[6], double (&q)[6], double (&xf)[6], double (&qf)[6], int &legmax, int npath);
namespace WOUT
{
    double dc1[50][12], drf1[50], dvf1[50], kc1[50], ck1[50], wtime[50][7];
    int kc = 0;
}

// Subroutine SWITCH
void SWITCH(double &amass, int legmax, int imax, double (&c)[12], double (&qoe)[6], double (&xoe)[6], int hmax);
// COMMON/WSWIT
namespace WSWITCH
{
    double e[12][13], dc[12];
    double drf, dvf, ck, evt;
    int kcount;
    double burnt;
}

// Subroutine COAST
void COAST(double (&xo)[6], double (&qo)[6], double (&xf)[6], double (&qf)[6], double (&phi)[6][6], double (&dphi)[6][6], int no);
// COMMON/WCOAST
namespace WCOAST
{
    double psy, alpha, ft;
}

// Subroutine AMULT
void AMULT();
// COMMON/CCOAST
namespace WAMULT
{
    double ann[2][2], bnn[3][3], xxo[6][6], ro[3], vo[3], r[3], v[3];
    double dro[3], dvo[3], dr[3], dv[3], dan[2][2], dbnn[3][3], dxxo[6][6];
}

// Subroutine RKGO31
void RKGO31(double (&xo)[6], double (&qo)[6], double (&xf)[6], double (&qf)[6], double (&z)[12][12], double &evt, int &hmax, int &legmax, int &no);

// Subroutine RKSTEP
// This program advances Yn and Ydn by a step of size H to Yn1 and Ydn1
// using a 4th-order R-K numerical integration scheme.
// If N is positive, all elements of the matrices Yn and Ydn are advanced.
// Otherwise only the first column of each matrix is updated.
void RKSTEP(double (&yn)[6][13], double (&ydn)[6][13], double tn, double (&yn1)[6][13], double (&ydn1)[6][13], double h, int n, int &leg);

// Subroutine YDDRHS
// Compute basic quantities common to many components of Ydd.
void YDDRHS(double (&y)[6][13], double (&ydd)[6][13], double h, double t, int n);

// Subroutine BVEVAL, a.k.a OPDK3A
// This version deals with a 5-constraint right-end boundary-value problem,
// where the five constrained functions are the three components of the
// orbital angular velocity vector (R cross V) and the first two components of the vector
// whose direction is the direction of pericenter and whose magnitude is the orbital eccentricity.
// Thus, in effect, all of the six classical orbital elements are constrained
// except the mean anomaly, which is free.
void BVEVAL(double (&xf)[6], double (&qf)[6], double (&z)[12][12], double (&c)[12], double (&e)[12][13], double (&dc)[12]);

// Subroutine ADJUST
void ADJUST(double (&qo)[6], double (&e)[12][13], double (&z)[12][12], double &drf, double &dvf, double &ck, int legmax);
// COMMON/WADJ
namespace WADJUST
{
    double a[8], dxf[6];
}

// Subroutine SOLVE
void SOLVE(double (&a)[12][13], int legmax);


void OUT(double (&x)[6], double (&q)[6], double (&xf)[6], double (&qf)[6], int &legmax, int npath)
{
    using namespace CCPINJ;
    using namespace WSWITCH;
    using namespace WADJUST;
    using namespace WCOAST;
    using namespace WOUT;
    double h[3], e1[3], dc1[50][12], drf1[50];
    if (npath < 0)
    {   
        if (leg < 3)
            cout << "\n\nIteration number: " << kcount << endl;
        if (leg != 1)
        {
            cout << "\nCoast arc, leg = " << leg - 1 << endl;
            cout << "State at end:    " << x[0] << "    " << x[1] << "    " << x[2] << "    " << x[3] << "    " << x[4] << "    " << x[5];
            cout << "\nCostate at end:    " << q[0] << "    " << q[1] << "    " << q[2] << "    " << q[3] << "    " << q[4] << "    " << q[5] << endl;
            cout << "Psy = " << psy << "    Alpha = " << alpha << "\nCalculated coast time = " << ft << endl;
        }
    }
    if (npath <= 0)
    {
        h[0] = x[1] * x[5] - x[2] * x[4];
        h[1] = x[2] * x[3] - x[0] * x[5];
        h[2] = x[0] * x[4] - x[1] * x[3];
        double hm = sqrt(h[0] * h[0] + h[1] * h[1] + h[2] * h[2]);
        double rm = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
        e1[0] = -(x[0] / rm + (h[1] * x[5] - h[2] * x[4])/ uk);
        e1[1] = -(x[1] / rm + (h[2] * x[3] - h[0] * x[5])/ uk);
        e1[2] = -(x[2] / rm + (h[0] * x[4] - h[1] * x[3])/ uk);
        double em = sqrt(e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2]);
        double energy = -uk / rm + 0.5 * (x[3] * x[3] + x[4] * x[4] + x[5] * x[5]);
        double aaxis = -uk / (2 * energy);
        double rmin = aaxis * (1 - em);
        double rmax = aaxis * (1 + em);
        double period = 6.2831853 * sqrt(fabs(aaxis * aaxis * aaxis / uk));
        cout << "Semimajor axis = " << aaxis << "    Rmin = " << rmin << "    Rmax = " << rmax << "    Energy = " << energy;
        cout << "\nPeriod = " << period << "    |H| = " << hm << "    H vector:    " << h[0] << "    " << h[1] << "    " << h[2] << endl;
        cout << "|R| = " << rm << "    |E| = " << em << "    E vector:    " << e1[0] << "    " << e1[1] << "    " << e1[2] << endl;
        if (npath < 0)
        {
            cout << "\nBurn arc, leg = " << leg << "    Mass at end of leg = " << am;
            cout << "\nState at end:    " << xf[0] << "    " << xf[1] << "    " << xf[2] << "    " << xf[3] << "    " << xf[4] << "    " << xf[5];
            cout << "\nCostate an end:    " << qf[0] << "    " << qf[1] << "    " << qf[2] << "    " << qf[3] << "    " << qf[4] << "    " << qf[5] << endl;
            if (leg < legmax)
                return;
            kc = (kcount - 1) % 50;
            for (int i = 0; i < leg; ++i)
                wtime[kc][i] = atp1[i + 1] - atp1[i];
            return;
        }
        double dete = 1;
        for (int i = 0; i < leg + 6; ++i)
            dete *= e[i][i];
        wtime[kc][6] = burnt;
        cout << "Total burn time = " << wtime[kc][6] << "    Arc times: ";
        for (int i = 0; i < leg; ++i)
            cout << wtime[kc][i] << "    ";
        cout << "\nDC:    " << dc[0] << "    " << dc[1] << "    " << dc[2] << "    " << dc[3] << "    " << dc[4] << "    " << dc[5] << "\n       " << dc[6] << "    " << dc[7] << "    " << dc[8] << "    " << dc[9] << "    " << dc[10] << "    " << dc[11] << endl;
        cout << "Determinant of E = " << dete;
        cout << "\nDiagonal of E:    ";
        for (int i = 0; i < leg + 6; ++i)
            cout << e[i][i] << "    ";
        cout << "\nDXF:    ";
        for (int i = 0; i < 6; ++i)
            cout << dxf[i] << "    "; 
        cout << "\nDRF = " << drf << "    DVF = " << dvf << "    CK = " << ck << "    EVT = " << evt << endl;
        cout << "CK = min. of:    ";
        for (int i = 0; i < leg + 2; ++i)
            cout << a[i] << "    ";
        cout << "\nChange requisited in initial costate, switching times and final time:" << endl;
        for (int i = 0; i < leg + 7; ++i)
            cout << e[i][leg + 6] << "    ";
        cout << "\n\nEnd of iteration " << kcount << endl;
        cout << "New Qo:    " << qf[0] << "    " << qf[1] << "    " << qf[2] << "    " << qf[3] << "    " << qf[4] << "    " << qf[5] << endl;
        cout << "New switch times:    ";
        for (int i = 1; i < leg + 1; ++i)
            cout << atp1[i] << "    ";
        cout << endl;
        for (int i = 1; i < 12; ++i)
            dc1[kc][i] = dc[i];
        ck1[kc] = ck;
        drf1[kc] = drf;
        dvf1[kc] = dvf;
        kc1[kc] = kcount;
    }
    else
    {
        if (kcount > 50)
            kcount = 50;
        cout << "\n\nSummary Tables:\n" << "Iteration    Total burn time    Length of burn and coast arcs" << endl;
        for (int i = 0; i < kcount; ++i)
        {
            cout << kc1[i] << "        " << wtime[i][6] << "        ";
            for (int j = 0; j < leg; ++j)
                cout << wtime[i][j] << "    ";
            cout << endl;
        }
    }
    return;
}


void SWITCH(double &amass, int legmax, int imax, double (&c)[12], double (&qoe)[6], double (&xoe)[6], int hmax)
{
    using namespace CCPINJ;
    using namespace WSWITCH;
    double xo[6], qo[6], xf[6], qf[6], z[12][12], qo1[6], dummy[13], phi[6][6], dphi[6][6], dz[12][12];
    kcount = 0;
    int no = 0;
    int kmax = imax;
    evt = 1E-8;
    for (int i = 0; i < 6; ++i)
        qo1[i] = qoe[i];
    cout << "Constants:\nGravity constant = " << uk << "\nInitial mass = " << amass << "\nMass rate = " << bmax
         << "\nExhaust velocity = " << cexv << "\nHMax = " << hmax;
    cout << "    IMax = " << imax << "    LegMax = " << legmax;
    cout << "\nInitial arcs:\nTime(s)       Type of arc" << endl;
    for (int i = 0; i < legmax + 1; ++i)
        cout << atp1[i] << "        " << (atp2[i] < 0 ? "Burn" : "Coast") << endl;
    cout << "Initial state Xo: " << xoe[0] << "  " << xoe[1] << "  " << xoe[2] << "  " << xoe[3] << "  " << xoe[4] << "  " << xoe[5] << endl;
    cout << "Estimated costate Qo: " << qo1[0] << "  " << qo1[1] << "  " << qo1[2] << "  " << qo1[3] << "  " << qo1[4] << "  " << qo1[5] << endl;
    cout << "Desired final C: " << c[0] << "  " << c[1] << "  " << c[2] << "  " << c[3] << "  " << c[4] << "  " << c[5] << "  " << c[6] << endl;
    do
    {
        ++kcount;
        leg = 1;
        am = amass;
        burnt = 0.0;
        for (int i = 0; i < 6; ++i)
        {
            qo[i] = qo1[i];
            xo[i] = xoe[i];
        }
        for (int i = 0; i < 12; ++i)
        {
            e[i][12] = 0;
            for (int j = 0; j < 12; ++j)
            {
                e[i][j] = 0;
                z[i][j] = 0;
            }
        }
        for (int i = 0; i < 6; ++i)
            z[i + 6][i] = 1;
        /// Dangerous Zone ------------------------------
        double um, cbmu, r2, rs, c3, c4, ump;
        if (atp2[0] < 0)
            goto LSWITCH_6;

LSWITCH_7:
        COAST(xo, qo, xf, qf, phi, dphi, no);
        for (int i = 0; i < 6; ++i)
            for (int j = 0; j < leg + 5; ++j)
            {
                dz[i][j] = 0;
                dz[i + 6][j] = 0;
                for (int k = 0; k < 6; ++k)
                {
                    dz[i][j] += phi[i][k] * z[k][j];
                    dz[i + 6][j] += phi[i][k] * z[k + 6][j] + dphi[i][k] * z[k][j];
                }
            }
        for (int i = 0; i < 12; ++i)
            for (int j = 0; j < leg + 5; ++j)
                z[i][j] = dz[i][j];
        um = sqrt(qf[0] * qf[0] + qf[1] * qf[1] + qf[2] * qf[2]);
        cbmu = cexv * bmax / am / um;
        for (int i = 0; i < 3; ++i)
            z[i + 3][leg + 5] = -cbmu * qf[i];
        if (leg <= 1)
            goto LSWITCH_12;
        e[leg + 5][legmax + 6] = ump - um;
        for (int i = 0; i < leg + 6; ++i)
        {
            e[leg + 5][i] = -e[leg + 5][i];
            for (int j = 0; j < 3; ++j)
                e[leg + 5][i] += qf[j] / um * z[j + 6][i];
        }
        for (int i = 0; i < 6; ++i)
        {
            qo[i] = qf[i];
            xo[i] = xf[i];   
        }
        ++leg;
            
LSWITCH_6:  
        RKGO31(xo, qo, xf, qf, z, evt, hmax, legmax, no);
        burnt += atp1[leg] - atp1[leg - 1];
        OUT(xo, qo, xf, qf, legmax, -1);
        if (leg >= legmax)
            goto LSWITCH_18;
        um = sqrt(qf[0] * qf[0] + qf[1] * qf[1] + qf[2] * qf[2]);
        ump = um;
        cbmu = cexv * bmax / am / um;
        for (int i = 0; i < 3; ++i)
            z[i + 3][leg + 5] = cbmu * qf[i];
        for (int i = 0; i < leg + 6; ++i)
            for (int j = 0; j < 3; ++j)
                e[leg + 6][i] += qf[j] / um * z[j + 6][i];

LSWITCH_12:
        r2 = xf[0] * xf[0] + xf[1] * xf[1] + xf[2] * xf[2];
        rs = xf[0] * qf[0] + xf[1] * qf[1] + xf[2] * qf[2];
        c3 = - uk / r2 / sqrt(r2);
        c4 = -3 * c3 * rs / r2;
        e[leg + 5][legmax + 6] = xf[3] * qf[3] + xf[4] * qf[4] + xf[5] * qf[5] - c3 * rs;
        for (int i = 0; i < 3; ++i)
        {
            dummy[i] = -c3 * qf[i] - c4 * xf[i];
            dummy[i + 3] = qf[i + 3];
            dummy[i + 6] = -c3 * xf[i];
            dummy[i + 9] = xf[i + 3];
        }
        for (int i = 0; i < leg + 6; ++i)
            for (int j = 0; j < 12; ++j)
                e[leg + 5][i] -= dummy[j] * z[j][i];
        for (int i = 0; i < 6; ++i)
        {
            qo[i] = qf[i];
            xo[i] = xf[i];
        }
        ++leg;
        if (atp2[leg - 1] > 0)
            goto LSWITCH_7;
        else
        {
            e[6][6] = 0;
            goto LSWITCH_6;
        }
        

     
        /// Dangerous Zone End --------------------------
LSWITCH_18:
        BVEVAL(xf, qf, z, c, e, dc);
        for (int i = 0; i < 6; ++i)
        {
            e[i][legmax + 6] = dc[i];
            dc[i + 6] = e[i + 6][legmax + 6];
            e[legmax + 5][i] = qo1[i];
        }
        ADJUST(qo1, e, z, drf, dvf, ck, legmax);
        evt = max(1E-14, 1E-8 * dvf * dvf);
        evt = min(evt, 1E-8);
        cout << endl;
        OUT(xf, qf, xo, qo1, legmax, 0);
        if (kcount >= kmax)
            break;
    } 
    while (drf >= 0.1 || dvf >= 0.005);
    OUT(xo, qoe, xf, qf, legmax, 1);
    return;
}


void COAST(double (&xo)[6], double (&qo)[6], double (&xf)[6], double (&qf)[6], double (&phi)[6][6], double (&dphi)[6][6], int no)
{
    using namespace CCPINJ;
    using namespace WCOAST;
    using namespace WAMULT;
    complex<double> calph, capsy;
    double ho[3];
    // Time of coast
    double t = atp1[leg] - atp1[leg - 1];
    // Qo (Q at start of coast)
    for (int i = 0; i < 3; ++i)
    {
        ro[i] = xo[i];
        vo[i] = xo[i + 3];
        dro[i] = qo[i];
        dvo[i] = qo[i + 3];
    }
    // Ro, Vo: state at start of coast
    int jump = 0;
    double rmo = sqrt(ro[0] * ro[0] + ro[1] * ro[1] + ro[2] * ro[2]);
    double drmo = (ro[0] * dro[0] + ro[1] * dro[1] + ro[2] * dro[2]) / rmo;
    double sigo = ro[0] * vo[0] + ro[1] * vo[1] + ro[2] * vo[2];
    double dsigo = vo[0] * dro[0] + vo[1] * dro[1] + vo[2] * dro[2]
                    + ro[0] * dvo[0] + ro[1] * dvo[1] + ro[2] * dvo[2];
    alpha = vo[0] * vo[0] + vo[1] * vo[1] + vo[2] * vo[2] - 2 * uk / rmo;
    ho[0] = ro[1] * vo[2] - ro[2] * vo[1];
    ho[1] = ro[2] * vo[0] - ro[0] * vo[2];
    ho[2] = ro[0] * vo[1] - ro[1] * vo[0];
    double po = (ho[0] * ho[0] + ho[1] * ho[1] + ho[2] * ho[2]) / uk;
    psy = t / po;
    if (alpha < 0)
    {
        calph.real(sqrt(-alpha));
        calph.imag(0);
    }
    else
    {
        calph.real(0);
        calph.imag(sqrt(alpha));       
    }
    double s0, s1, s2, s3, rm;
    do
    {
        capsy = psy * calph;
        s0 = cos(capsy).real();
        s1 = (sin(capsy) / calph).real();
        s2 = (s0 - 1) / alpha;
        s3 = (s1 - psy) / alpha;
        ft = rmo * s1 + sigo * s2 + uk * s3;
        rm = rmo * s0 + sigo * s1 + uk * s2;
        if (jump == 1)
            break;
        psy += (t - ft) / rm;
        if (fabs(t - ft) < 0.0001)
            jump = 1;
    } 
    while (1);
    double fm1 = - uk * s2 / rmo; 
    double f = 1 + fm1;
    double fd = -uk * s1 / rm / rmo;
    double g = ft - uk * s3;
    double gdm1 = -uk * s2 / rm;
    double gd = 1 + gdm1;
    double ukr3 = uk / rm / rm / rm;
    double ukro3 = uk / rmo / rmo / rmo;
    double dalph = 2 * (vo[0] * dvo[0] + vo[1] * dvo[1] + vo[2] * dvo[2] 
                        + ukro3 * (ro[0] * dro[0] + ro[1] * dro[1] + ro[2] * dro[2]));
    double dapa = dalph / alpha;
    double dapa2 = dapa / alpha;
    double dpsy = - (drmo * s1 + dsigo * s2 + rmo * (psy * s0 - s1) * dapa * 0.5
                     + sigo * (psy * s1 * 0.5 - s2) * dapa + uk * (psy - 1.5 * s1 + psy * s0 * 0.5) * dapa2) / rm;
    double ds0 = (alpha * dpsy + 0.5 * psy * dalph) * s1;
    double ds1 = s0 * dpsy + (psy * s0 - s1) * dapa * 0.5;
    double ds2 = s1 * dpsy + (0.5 * psy * s1 - s2) * dapa;
    double ds3 = s2 * dpsy + (psy - 1.5 * s1 + 0.5 * psy * s0) * dapa2;
    double s4 = (s2 - psy * psy * 0.5) / alpha;
    double ds4 = s3 * dpsy + (psy * psy * 0.5 - 2 * s2 + 0.5 * psy * s1) * dapa2;
    double s5 = (s3 - psy * psy * psy / 6) / alpha;
    double ds5 = s4 * dpsy + (psy * psy * psy / 6 + (2 * psy - 2.5 * s1 + 0.5 * psy * s0) / alpha) * dapa2;
    double u = s2 * ft + uk * (psy * s4 - 3 * s5);
    double du = ds2 * ft + uk * (dpsy * s4 + psy * ds4 - 3 * ds5);
    double drm = s0 * drmo + ds0 * rmo + s1 * dsigo + ds1 * sigo + uk * ds2;
    double df = (-uk * ds2 - fm1 * drmo) / rmo;
    double dg = -uk * ds3;
    double dgd = (-uk * ds2 - gdm1 * drm) / rm;
    double r01 = rmo * rm;
    double dro1 = rm * drmo + drm * rmo;
    double dfd = (-uk * ds1 - fd * dro1) / r01;
    for (int i = 0; i < 3; ++i)
    {
        dr[i] = ro[i] * df + vo[i] * dg + dro[i] * f + dvo[i] * g;
        qf[i] = dr[i];
        r[i] = ro[i] * f + vo[i] * g;
        dv[i] = ro[i] * dfd + vo[i] * dgd + dro[i] * fd + dvo[i] * gd;  
        qf[i + 3] = dv[i];
        v[i] = ro[i] * fd + vo[i] * gd;
    }
    // R, V: state at end of coast
    if (no != 1)
    {
        double dukr3 = -3 * ukr3 * drm / rm;
        double duro3 = -3 * ukro3 * drmo / rmo;
        double s1ro = s1 / rmo;
        double ds1ro = (ds1 - s1ro * drmo) / rmo;
        double s1r = s1 / rm;
        double ds1r = (ds1 - s1r * drm) / rm;
        double r02 = 1 / rmo / rmo;
        double r2 = 1 / rm / rm;
        double uuko3 = -u * ukro3;
        double duuk3 = -du * ukro3 - u * duro3;
        ann[0][0] = -fd * s1ro - fm1 * r02; 
        ann[0][1] = -fd * s2;
        ann[1][0] = fm1 * s1ro + uuko3;
        ann[1][1] = fm1 * s2;
        double dumm1 = ann[0][0];
        dan[0][0] = -fd * ds1ro - dfd * s1ro + ukro3 * ds2 + duro3 * s2;
        double dummy = dan[0][0];
        dan[0][1] = -dfd * s2 - fd * ds2;
        dan[1][0] = fm1 * ds1ro + df * s1ro + duuk3;
        dan[1][1] = fm1 * ds2 + df * s2;
        AMULT();
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
            {
                dxxo[i][j] = dbnn[i][j];
                xxo[i][j] = bnn[i][j];
            }
        for (int i = 0; i < 3; ++i)
        {
            dxxo[i][i] += df;
            xxo[i][i] += f;
        }
        ann[0][0] = ann[0][1];
        ann[1][0] = ann[1][1];
        ann[0][1] = -gdm1 * s2;
        ann[1][1] = g * s2 - u;
        dan[0][0] = dan[0][1];
        dan[1][0] = dan[1][1];
        dan[0][1] = -gdm1 * ds2 - dgd * s2;
        dan[1][1] = -du + dg * s2 + g * ds2;
        AMULT();
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
            {
                dxxo[i][j + 3] = dbnn[i][j];
                xxo[i][j + 3] = bnn[i][j];
            }
        for (int i = 0; i < 3; ++i)
        {
            dxxo[i][i + 3] += dg;
            xxo[i][i + 3] += g;
        }
        ann[1][0] = -ann[0][0];
        ann[1][1] = -ann[0][1];
        ann[0][0] = -fd * s1r - gdm1 * r2;
        ann[0][1] = u * ukr3 - gdm1 * s1r;
        dan[1][0] = -dan[0][0];
        dan[1][1] = -dan[0][1];
        dan[0][0] = -dfd * s1r - fd * ds1r + ukr3 * ds2 + dukr3 * s2;
        dan[0][1] = -gdm1 * ds1r - dgd * s1r + du * ukr3 + u * dukr3;
        AMULT();
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
            {
                dxxo[i + 3][j + 3] = dbnn[i][j];
                xxo[i + 3][j + 3] = bnn[i][j];
            }
        for (int i = 3; i < 6; ++i)
        {
            dxxo[i][i] += dgd;
            xxo[i][i] += gd;
        }
        ann[0][1] = ann[0][0];
        ann[1][1] = ann[1][0];
        ann[1][0] = -dumm1;
        ann[0][0] = -fd * (s0 / r01 + r2 + r02) - uuko3 * ukr3;
        dan[0][1] = dan[0][0];
        dan[1][1] = dan[1][0];
        dan[1][0] = -dummy;
        dan[0][0] = -uuko3 * dukr3 - duuk3 * ukr3 - dfd * (s0 / r01 + r2 + r02)
                    - fd * ((ds0 - s0 * dro1 / r01) / r01 - 2 * (r2 * drm / rm + r02 * drmo / rmo));
        AMULT();
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
            {
                dxxo[i + 3][j] = dbnn[i][j];
                xxo[i + 3][j] = bnn[i][j];
            }
        for (int i = 0; i < 3; ++i)
        {
            dxxo[i + 3][i] += dfd;
            xxo[i + 3][i] += fd;
        }
        // Xxo: partial of Xf respect to Xo
        // DXxo: partial of Qf respect to Xo
        for (int i = 0; i < 6; ++i)
            for (int j = 0; j < 6; ++j)
            {
                phi[i][j] = xxo[i][j];
                dphi[i][j] = dxxo[i][j];
            }
    }
    for (int i = 0; i < 3; ++i)
    {
        xf[i] = r[i];
        xf[i + 3] = v[i];
    }
    return;
}


void AMULT()
{
    using namespace WAMULT;
    double da[3][2], a[3][2];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j< 2; ++j)
        {
            da[i][j] = dr[i] * ann[0][j] + dv[i] * ann[1][j] + r[i] * dan[0][j] + v[i] * dan[1][j];
            a[i][j] = ann[0][j] * r[j] + ann[1][j] * v[j];
        }
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
        {
            dbnn[i][j] = a[i][0] * dro[j] + a[i][1] * dvo[j] + da[i][0] * ro[j] + da[i][1] * vo[j];
            bnn[i][j] = a[i][0] * ro[j] + a[i][1] * vo[j];
        }
    return;
}


void RKGO31(double (&xo)[6], double (&qo)[6], double (&xf)[6], double (&qf)[6], double (&z)[12][12], double &evt, int &hmax, int &legmax, int &no)
{
    using namespace CCPINJ;
    double yn[6][13], ydn[6][13], y3h[6][13], yd3h[6][13], ym[6][13], ydm[6][13], ey[6], eyd[6];
    for (int i = 0; i < 3; ++i)
    {
        yn[i][0] = xo[i];
        yn[i + 3][0] = qo[i];
        ydn[i][0] = xo[i + 3];
        ydn[i + 3][0] = qo[i + 3];
        for (int j = 1; j < leg + 6; ++j)
        {
            yn[i][j] = z[i][j - 1];
            yn[i + 3][j] = z[i + 6][j - 1];
            ydn[i][j] = z[i + 3][j - 1];
            ydn[i + 3][j] = z[i + 9][j - 1];
        }
    }
    double tf = atp1[leg];
    double t0 = atp1[leg - 1];
    double h = copysign(hmax, tf - t0);
    double tn = t0;
    while (1)
    {
        if (fabs(h) > fabs(tf - tn))
            h = tf - tn;
        if (fabs(h) > hmax)
            h = copysign(hmax, h);
        RKSTEP(yn, ydn, tn, y3h, yd3h, h, 1, leg);
        RKSTEP(yn, ydn, tn, ym, ydm, h / 3, 0, leg);
        RKSTEP(ym, ydm, tn + h / 3, ym, ydm, h / 3, 0, leg);
        RKSTEP(ym, ydm, tn + h / 3 * 2, ym, ydm, h / 3, 0, leg);
        for (int i = 0; i < 6; ++i)
        {
            ey[i] = 0.125E-1 * (ym[i][0] - y3h[i][0]);
            eyd[i] = 0.125E-1 * (ydm[i][0] - yd3h[i][0]);
        }
        double ev2min = 1E-13 * (ydm[0][0] * ydm[0][0] + ydm[1][0] * ydm[1][0] + ydm[2][0] * ydm[2][0]);
        double evl2 = max(evt, ev2min);
        double r = (eyd[0] * eyd[0] + eyd[1] * eyd[1] + eyd[2] * eyd[2]) / evl2;
        for (int i = 0; i < 6; ++i)
        {
            yn[i][0] = ym[i][0] + ey[i];
            ydn[i][0] = ydm[i][0] + eyd[i];
        }
        for (int i = 0; i < 6; ++i)
            for (int j = 1; j < leg + 6; ++j)
            {
                yn[i][j] = y3h[i][j];
                ydn[i][j] = yd3h[i][j];
            }
        if (fabs(h) >= fabs(tf - tn))
            break;
        tn += h;
        r = (r < 0.04 ? 0.04 : r);
        h = h / pow(r, 0.125);
    }
    if (leg <= legmax)
    {
        YDDRHS(yn, yd3h, 1, atp1[leg], 0);
        for (int i = 0; i < 3; ++i)
        {
            z[i][legmax + 5] = ydn[i][0];
            z[i + 3][legmax + 5] = yd3h[i][0];
            z[i + 6][legmax + 5] = ydn[i + 3][0];
            z[i + 9][legmax + 5] = yd3h[i + 3][0];
        }
    }
    for (int i = 0; i < 3; ++i)
    {
        xf[i] = yn[i][0];
        xf[i + 3] = ydn[i][0];
        qf[i] = yn[i + 3][0];
        qf[i + 3] = ydn[i + 3][0];
        for (int j = 1; j < leg + 6; ++j)
        {
            z[i][j - 1] = yn[i][j];
            z[i + 3][j - 1] = ydn[i][j];
            z[i + 6][j - 1] = yn[i + 3][j];
            z[i + 9][j - 1] = ydn[i + 3][j];
        }
    }
    am -= bmax * (atp1[leg] - atp1[leg - 1]);
    return;
}


void RKSTEP(double (&yn)[6][13], double (&ydn)[6][13], double tn, double (&yn1)[6][13], double (&ydn1)[6][13], double h, int n, int &leg)
{
    double d1[6][13], d2[6][13], d3[6][13], y[6][13];
    int jmax = 1;
    if (n > 0)
        jmax = leg + 6;
    double h2 = 0.5 * h;
    YDDRHS(yn, d1, h, tn, n);
    for (int j = 0; j < jmax; ++j)
        for (int i = 0; i < 6; ++i)
            y[i][j] = yn[i][j] + h2 * (ydn[i][j] + 0.25 * d1[i][j]);
    YDDRHS(y, d2, h, tn + h2, n);
    for (int j = 0; j < jmax; ++j)
        for (int i = 0; i < 6; ++i)
            y[i][j] = yn[i][j] + h * (ydn[i][j] + 0.5 * d2[i][j]);
    YDDRHS(y, d3, h, tn + h, n);
    for (int j = 0; j < jmax; ++j)
        for (int i = 0; i < 6; ++i)
        {
            yn1[i][j] = yn[i][j] + h * (ydn[i][j] + (d1[i][j] + 2 * d2[i][j]) / 6);
            ydn1[i][j] = ydn[i][j] + (d1[i][j] + 4 * d2[i][j] + d3[i][j]) / 6;
        }    
    return;
}


void YDDRHS(double (&y)[6][13], double (&ydd)[6][13], double h, double t, int n)
{
    using namespace CCPINJ;
    double r[3], u[3], b[6][6];
    // Compute basic quantities common to many components of Ydd
    for (int i = 0; i < 3; ++i)
    {
        r[i] = y[i][0];
        u[i] = y[i + 3][0];
    }
    double am1 = am - (t - atp1[leg - 1]) * bmax;
    double r2 = 1 / (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
    double u2 = 1 / (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
    double ru = r[0] * u[0] + r[1] * u[1] + r[2] * u[2];
    double alpha = -h * uk * r2 * sqrt(r2);
    double beta = h * sqrt(u2) * cexv * bmax / am1;
    double gamma = -3 * alpha * r2 * ru;
    // Compute Rdd and Udd
    for (int i = 0; i < 3; ++i)
    {
        ydd[i][0] = r[i] * alpha + u[i] * beta;  
        ydd[i + 3][0] = r[i] * gamma + u[i] * alpha;
    }
    // Decide whether Wdd is required at this time
    if (n <= 0)
        return;
    // Compute the additional quantities common to many components of Wdd
    double delta = -3 * alpha * r2;
    double epsil = -beta * u2;
    double zeta = -5 * gamma * r2;
    // Compute the matrix B needed in the matrix equation Wdd = B * W
    for (int j = 0; j < 3; ++j)
    {
        double rrj = delta * r[j];
        double ruj = epsil * u[j];
        double urj = zeta * r[j] + delta * u[j];
        double uuj = rrj;
        for (int i = 0; i < 3; ++i)
        {
            double q = 0;
            if (i == j) 
                q = 1;
            b[i][j] = r[i] * rrj + q * alpha;
            b[i][j + 3] = u[i] * ruj + q * beta;
            b[i + 3][j] = r[i] * urj + u[i] * uuj + q * gamma;
            b[i + 3][j + 3] = b[i][j];
        }
    }
    // Perform the matrix multiplication B * W to get Wdd
    for (int i = 0; i < 6; ++i)
        for (int j = 1; j < leg + 6; ++j)
        {
            double sum = 0;
            for (int k = 0; k < 6; ++k)
                sum += b[i][k] * y[k][j];
            ydd[i][j] = sum;
        }
    if (leg <= 1)
        return;
    double rddmb = -beta * bmax / am1;
    for (int j = 7; j < leg + 6; ++j)
        for (int i = 0; i < 3; ++i)
        {
            ydd[i][j] -= copysign(rddmb, atp2[j - 7]) * u[i];
        }
    return;
}


void BVEVAL(double (&xf)[6], double (&qf)[6], double (&z)[12][12], double (&c)[12], double (&e)[12][13], double (&dc)[12])
{
    using namespace CCPINJ;
    double r[3], v[3], g[7][6] = {0}, dummy[12];
    for (int i = 0; i < 3; ++i)
    {
        r[i] = xf[i];
        v[i] = xf[i + 3];
    }
    double r2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
    g[0][1] = v[2];
    g[0][2] = -v[1];
    g[1][2] = v[0];
    g[0][4] = -r[2];
    g[0][5] = r[1];
    g[1][5] = -r[0];
    double rm = sqrt(r2);
    double r3 = rm * r2;
    double c1 = -1 / rm + (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) / uk;
    double c2 = -(r[0] * v[0] + r[1] * v[1] + r[2] * v[2]) / uk;
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            if (i > j)
            {
                g[i][j] = -g[j][i];
                g[i][j + 3] = -g[j][i + 3];
            }
            g[i + 3][j] = r[i] * r[j] / r3 - v[i] * v[j] / uk;
            g[i + 3][j + 3] = (r[i] * v[j] * 2 - v[i] * r[j]) / uk;
        }
        g[i][i] = 0;
        g[i][i + 3] = 0;
        g[i + 3][i] += c1;
        g[i + 3][i + 3] += c2;
    }
    for (int i = 0; i < 3; ++i)
    {
        dc[i] = 0;
        for (int j = 0; j < 3; ++j)
            dc[i] += g[i][j] * r[j];  
    }
    for (int i = 0; i < 2; ++i)
    {
        double sum = 0;
        for (int j = 0; j < 3; ++j)
            sum += g[i][j] * dc[j];
        dc[i + 3] = c[i + 3] + r[i] / rm + sum / uk;
    }
    for (int i = 0; i < 3; ++i)
        dc[i] = c[i] - dc[i];
    for (int i = 0; i < 5; ++i)
        for (int j = 0; j < leg + 6; ++j)
        {
            e[i][j] = 0;
            for (int k = 0; k < 6; ++k)
                e[i][j] += g[i][k] * z[k][j];
        }
    double rs = xf[0] * qf[0] + xf[1] * qf[1] + xf[2] * qf[2];
    double c3 = -uk / r3;
    double c4 = -3 * c3 * rs / r2;
    dc[5] = xf[3] * qf[3] + xf[4] * qf[4] + xf[5] * qf[5] - c3 * rs;
    for (int i = 0; i < 3; ++i)
    {
        dummy[i] = -c3 * qf[i] - c4 * xf[i];
        dummy[i + 3] = qf[i + 3];
        dummy[i + 6] = -c3 * xf[i];
        dummy[i + 9] = xf[i + 3];
    }
    for (int j = 0; j < leg + 6; ++j)
        for (int k = 0; k < 12; ++k)
            e[5][j] -= dummy[k] * z[k][j];
    return;
}


void ADJUST(double (&qo)[6], double (&e)[12][13], double (&z)[12][12], double &drf, double &dvf, double &ck, int legmax)
{
    using namespace CCPINJ;
    using namespace WADJUST;
    SOLVE(e, legmax);
    for (int i = 0; i < 6; ++i)
    {
        dxf[i] = 0;
        for (int k = 0; k < legmax + 6; ++k)
            dxf[i] += z[i][k] * e[k][legmax + 6];
    }
    drf = sqrt(dxf[0] * dxf[0] + dxf[1] * dxf[1] + dxf[2] * dxf[2]);
    dvf = sqrt(dxf[3] * dxf[3] + dxf[4] * dxf[4] + dxf[5] * dxf[5]);
    double du2 = sqrt(e[0][legmax + 6] * e[0][legmax + 6] + e[1][legmax + 6] * e[1][legmax + 6] + e[2][legmax + 6] * e[2][legmax + 6]);
    double dud2 = sqrt(e[3][legmax + 6] * e[3][legmax + 6] + e[4][legmax + 6] * e[4][legmax + 6] + e[5][legmax + 6] * e[5][legmax + 6]);
    a[0] = 0.2 / du2;
    a[1] = 0.0003 / dud2;
    a[2] = 1;
    if (atp2[0] > 0)
        a[2] = 0.5 * (atp1[1] - atp1[0]) / fabs(e[6][legmax + 6]);
    ck = min(min(1.0, a[0]), min(a[1], a[2]));
    for (int i = 1; i < legmax; ++i)
    {
        a[i + 2] = 0.5 * (atp1[i + 1] - atp1[i]) / fabs(e[i + 6][legmax + 6] - e[i + 5][legmax + 6]);
        ck = min(ck, a[i + 2]);
    }
    for (int i = 0; i < 6; ++i)
    {
        atp1[i + 1] += ck * e[i + 6][legmax + 6];
        qo[i] += ck * e[i][legmax + 6];
    }
    double um = sqrt(qo[0] * qo[0] + qo[1] * qo[1] + qo[2] * qo[2]);
    for (int i = 0; i < 6; ++i)
        qo[i] /= um;
    return;
}


void SOLVE(double (&a)[12][13], int legmax)
{
    for (int n = 0; n < legmax + 6; ++n)
    {
        int ibig = n;
        for (int i = n; i < legmax + 6; ++i)
            if (fabs(a[i][n]) > fabs(a[ibig][n]))
                ibig = i;
        if (ibig != n)
            for (int j = n; j < legmax + 7; ++j)
            {
                double q = a[n][j];
                a[n][j] = a[ibig][j];
                a[ibig][j] = q;
            }
        for (int i = 0; i < legmax + 6; ++i)
            if (i != n)
            {
                double q = a[i][n] / a[n][n];
                for (int k = n + 1; k < legmax + 7; ++k)  
                    a[i][k] -= q * a[n][k];
            }
    }
    for (int i = 0; i < legmax + 6; ++i)
        a[i][legmax + 6] /= a[i][i];
    return;
}



int main(int argc, char* argv[])
{
    using namespace CCPINJ;
    int ncase;
    double amass;
    int hmax;
    double xoe[6], qoe[6];
    double c[12] = {0};
    int legmax, imax;
    
    cout << "OPGUID Solver v0.1\n";
    if (argc == 2)
    {
        if (argv[1] == "--verbose")
        {
            kwargs_verbose = true;
            cout << "kwargs_verbose: True";
        }
    }

    // 读入文件，原始文件中 5 号穿孔卡
    ifstream pc5("case.txt", ios::in);
    pc5 >> ncase;
    int nc = 0;
    while (nc < ncase)
    {
        ++nc;
        cout << "========================\n\nCase number: " << nc << endl;
        pc5 >> uk >> amass >> bmax >> cexv >> hmax;
        pc5 >> xoe[0] >> xoe[1] >> xoe[2] >> xoe[3] >> xoe[4] >> xoe[5];
        pc5 >> qoe[0] >> qoe[1] >> qoe[2] >> qoe[3] >> qoe[4] >> qoe[5];
        pc5 >> c[0] >> c[1] >> c[2] >> c[3] >> c[4] >> c[5] >> c[6];
        pc5 >> legmax >> imax;
        atp1 = new double[legmax];
        atp2 = new short[legmax];  
        for (int i = 0; i < legmax + 1; ++i) /// 注意: 5 个断点对应 4 个弧段
            pc5 >> atp1[i] >> atp2[i];
        SWITCH(amass, legmax, imax, c, qoe, xoe, hmax);
        delete []atp1;
        delete []atp2;
    }
    pc5.close();
    cout << "\n\n" << endl;
    return 0;
}
