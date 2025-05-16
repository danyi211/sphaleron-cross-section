/**
 * NNPDF tutorial: How to compute the gluon PDF central value 
 * and 1-sigma error band in C++, using LHAPDF interface.
 *
 * Compile with:
 * g++ tutorial.cxx -o tutorial_cxx -I $(shell lhapdf-config --incdir)
 * -L $(shell lhapdf-config --libdir) -lLHAPDF
 * 
 * Author: The NNPDF Collaboration - 2012
 */

#include <iostream>
#include <cmath>
#include <vector>
#include "LHAPDF/LHAPDF.h"
using namespace std;
using namespace LHAPDF;

int f1 = 0;
int f2 = 0;

double get_factor(double x, double m){
  double res = 0.;
  double a = -0.005;
  if( x < m ) res = x * ( 1./m + a*(x - m) ) - 1;
  return res;
}

double get_amp(double rs, double m){

  double E = rs/1000.;
  if(E < m){
    return 0;
  }else{
    return 1;
  }

}

double get_amp_TW(double rs, double c, double m){

  double E = rs/1000.;
  double aW = 1./30.;
  double pi = 4.*atan(1.);

  double amp;
  amp = exp( c * 4.*pi / aW * get_factor(E, m) );

  return amp;
}

//======================================//

int main()
{ 

  // Loading the NNPDF set
  //LHAPDF::initPDFSet("NNPDF21_100.LHgrid"); 
  //LHAPDF::initPDFSet("NNPDF30_nlo_as_0119_atlas.LHgrid"); 
  //LHAPDF::initPDFSet("NNPDF23_nlo_as_0119.LHgrid"); 
  //initPDFSet("MSTW2008nlo90cl_nf4");   
  //initPDFSet("CT14nnlo");  
  initPDFSet("CT10");    
  //const int NumberReplicas = LHAPDF::numberPDF();

  double mw = 80.4;
  double GeVsqinv_to_nb = 3.89 * pow(10., 5); 
  double GeVsqinv_to_fb = 3.89 * pow(10., 11); 
  double sig_EW = 1./(mw * mw) * GeVsqinv_to_fb;

  // int NEsph = 50;
  // vector<double> m;
  // for(int i=0; i<NEsph; i++){
  //   double dm = 8. + 0.5 * i;
  //   m.push_back(dm);
  // }

  int NEsph = 30;

  double E = 13.; // TeV
  double mmin = 8.;
  double mmax = 11.;

  // double E = 27.; // TeV
  // double mmin = 8.;
  // double mmax = 22.;

  // double E = 100.; // TeV
  // double mmin = 8.;
  // double mmax = 60.;

  double mstep = (mmax - mmin)/NEsph;

  vector<double> m;

  vector<double> sig;  
  for(int i=0; i<NEsph+1; i++){
    double dm = mmin + mstep * i;
    m.push_back(dm);
  }

  double rootS = E * 1000.;
  double S = pow(rootS, 2);
    
  int NtPoints = 10000;  

  for(int j=0; j<NEsph+1; j++){

    double t_min = m[j]*m[j]/S;
    double t_max = 1.;
    //double del = (t_max - t_min)/NtPoints;

    double u_max = 0.;
    double u_min = log(t_min);
    double del = (u_max - u_min)/NtPoints;

    sig.push_back(0.);

    for (int it = 0; it < NtPoints; it++){
      LHAPDF::initPDF(0); // <- Replica 0 == Central Value (AVG)

      // double ss = exp(slog_min + is * ds_log);
      // double rs = sqrt(ss);
      // double yac = ss;
      // double rt = rs / rootS;

      double u = u_min + it * del;
      double tau = exp(u);
      double rt = sqrt(tau);
      double yac = tau;
      double rs = rt * rootS;

      // double tt = t_min + it * del;
      // double rt = sqrt(tt);
      // double yac = 1.;
      // double rs = rt * rootS;

      int NyPoints = 500;
      double ymin = -log(1./rt);
      double ymax = +log(1./rt);
      double dy = (ymax - ymin)/NyPoints;
      double centrise = dy / 2.;

      // uu //
      double integ_uu = 0.;
      double integ_ud = 0.;
      double integ_us = 0.;    
      double integ_uc = 0.;    
      double integ_dd = 0.;
      double integ_ds = 0.;
      double integ_dc = 0.;
      double integ_ss = 0.;
      double integ_sc = 0.;    
      double integ_cc = 0.;
      for (int iy = 0; iy < NyPoints; iy++){
        double y = ymin + iy * dy + centrise;
        double x1 = rt * exp( y);
        double x2 = rt * exp(-y);
        double Q = rs;
        double D_1, D_2;

        //cout << it <<" "<< rt <<" "<< dt <<" "<< rs <<" "<< tau <<" "<< y <<" "<< x1 <<" "<< x2 << endl; 


        //-- u 
        D_1 = xfx(x1, Q, 2) / x1;      
        D_2 = xfx(x2, Q, 2) / x2;
        integ_uu += D_1 * D_2 * dy;

        D_2 = xfx(x2, Q, 1) / x2;
        integ_ud += 2. * D_1 * D_2 * dy;

        D_2 = xfx(x2, Q, 3) / x2;
        integ_us += 2. * D_1 * D_2 * dy;

        D_2 = xfx(x2, Q, 4) / x2;
        integ_uc += 2. * D_1 * D_2 * dy;

        //-- d
        D_1 = xfx(x1, Q, 1) / x1;      
        D_2 = xfx(x2, Q, 1) / x2;
        integ_dd += D_1 * D_2 * dy;

        D_2 = xfx(x2, Q, 3) / x2;
        integ_ds += 2. * D_1 * D_2 * dy;

        D_2 = xfx(x2, Q, 4) / x2;
        integ_dc += 2. * D_1 * D_2 * dy;

        //-- s
        D_1 = xfx(x1, Q, 3) / x1;      
        D_2 = xfx(x2, Q, 3) / x2;
        integ_ss += D_1 * D_2 * dy;

        D_2 = xfx(x2, Q, 4) / x2;
        integ_sc += 2. * D_1 * D_2 * dy;

        //-- c
        D_1 = xfx(x1, Q, 4) / x1;      
        D_2 = xfx(x2, Q, 4) / x2;
        integ_cc += D_1 * D_2 * dy;

      }
      double Lumi = 0;
      Lumi = integ_uu + integ_ud + integ_us + integ_uc
           + integ_dd + integ_ds + integ_dc 
           + integ_ss + integ_sc + integ_cc; 

        double amp = get_amp(rs, m[j]);        
        //cout << rs <<"  "<< amp <<"  "<< sig_EW <<"  "<< Lumi << endl;
        double d_sig = sig_EW * Lumi * amp / 12.;
        sig[j] += d_sig * yac * del;

      }

      // cout << sig[j] << "  " << m[j] <<"  "<< E << endl;
      cout << sig[j] << "  " << m[j] << endl;

    }
  return 0;

  }
