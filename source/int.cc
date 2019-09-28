#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <assert.h>
#include <stdio.h>
#include <vector>
#include <complex>
#define PI 3.14159265358979323846
#define delta_k 0.0007853981633974483

int indexing(int k_mag){
  if (k_mag <= 1){
    return - k_mag/2;
  }
  else{
    return 1280 - k_mag/2;
  }
}


double matter_power_lin(std::vector<double> &xData, std::vector<double> &yData, double x, bool extrapolate)
{
   int size = xData.size();

   int iii = 0;                                                                  // find left end of interval for interpolation
   if ( x >= xData[size - 2] )                                                 // special case: beyond right end
   {
      iii = size - 2;
   }
   else
   {
      while ( x > xData[iii+1] ) iii++;
   }
   double xL = xData[iii], yL = yData[iii], xR = xData[iii+1], yR = yData[iii+1];      // points on either side (unless beyond ends)
   if ( !extrapolate )                                                         // if beyond ends of array and not extrapolating
   {
      if ( x < xL ) yR = yL;
      if ( x > xR ) yL = yR;
   }

   double dydx = ( yR - yL ) / ( xR - xL );                                    // gradient

   return yL + dydx * ( x - xL );                                              // linear interpolation
}

int
main(int argc, char* argv[])
{ 
  const std::complex<double> rr(1.0,0.0);  
  const std::complex<double> ii(0.0,1.0);  

  double temp_read;
  
  std::ifstream file3("data/cfield_real.txt");
  std::vector<double> cfield_real;
  while(!file3.eof())
  {
    file3 >> temp_read;
    cfield_real.push_back(temp_read);
  }

  std::cout << "Finished Reading cr" << "\n";

  cfield_real.pop_back();

  std::vector<std::vector<std::vector<std::complex<double> > > > cfield(1280,std::vector<std::vector<std::complex<double> > >(1280,std::vector <std::complex<double> >(641,(0.,0.))));

  int ll;
  int mm;
  int nn;
  for(ll=0;ll<1280;ll++){
    std::cout << "COMPLEXITY" << ll << "\n";
    for(mm=0;mm<1280;mm++){
        for(nn=0;nn<641;nn++){
          cfield[ll][mm][nn] = rr*cfield_real[ll*1280*641+mm*641+nn];
      }
    }
  }

  cfield_real.clear();
  cfield_real.shrink_to_fit();

  std::ifstream file4("data/cfield_imag.txt");
  std::vector<double> cfield_imag;
  while(!file4.eof())
  {
    file4 >> temp_read;
    cfield_imag.push_back(temp_read);
  }

  std::cout << "Finished Reading ci" << "\n";

  cfield_imag.pop_back();

  for(ll=0;ll<1280;ll++){
    std::cout << "COMPLEXITY" << ll << "\n";
    for(mm=0;mm<1280;mm++){
        for(nn=0;nn<641;nn++){
          cfield[ll][mm][nn] = cfield[ll][mm][nn]+ii*cfield_imag[ll*1280*641+mm*641+nn];
      }
    }
  }

  cfield_imag.clear();
  cfield_imag.shrink_to_fit();

  std::vector<std::vector<std::vector<std::complex<double> > > > cfield2(1280,std::vector<std::vector<std::complex<double> > >(1280,std::vector <std::complex<double> >(1280,(0.,0.))));

  for(ll=0;ll<1280;ll++){
    std::cout << "COPYING" << ll << "\n";
    for(mm=0;mm<1280;mm++){
      for(nn=0;nn<641;nn++){
        cfield2[ll][mm][nn] = cfield[ll][mm][nn];
      }
    }
  }

  cfield.clear();
  cfield.shrink_to_fit();

  for(ll=1;ll<1280;ll++){
    std::cout << "READING " << ll << "\n";
    for(mm=1;mm<1280;mm++){
      for(nn=641;nn<1280;nn++){
        cfield2[ll][mm][nn] = std::conj(cfield2[1280-ll][1280-mm][1280-nn]);
      }
    }
  }

  for(ll=1;ll<1280;ll++){
    std::cout << "READING1" << ll << "\n";
    for(nn=641;nn<1280;nn++){
      cfield2[ll][0][nn] = std::conj(cfield2[1280-ll][0][1280-nn]);
    }
  }

  for(mm=1;mm<1280;mm++){
    std::cout << "READING2" << mm << "\n";
    for(nn=641;nn<1280;nn++){
      cfield2[0][mm][nn] = std::conj(cfield2[0][1280-mm][1280-nn]);
    }
  }

  for(nn=641;nn<1280;nn++){
    std::cout << "READING3" << nn << "\n";
    cfield2[0][0][nn] = std::conj(cfield2[0][0][1280-nn]);
  }

  std::ifstream file0("data/xi_mag.txt");
  std::vector<int> xi_mag;

  while(!file0.eof())
  {
  file0 >> temp_read;
  xi_mag.push_back(temp_read);
  }

  std::cout << "Finished Reading xi" << "\n";

  std::ifstream file1("data/xj_mag.txt");
  std::vector<int> xj_mag;

  while(!file1.eof())
  {
  file1 >> temp_read;
  xj_mag.push_back(temp_read);
  }

  std::cout << "Finished Reading xj" << "\n";

  std::ifstream file2("data/xk_mag.txt");
  std::vector<int> xk_mag;

  while(!file2.eof())
  {
  file2 >> temp_read;
  xk_mag.push_back(temp_read);
  }

  std::cout << "Finished Reading xk" << "\n";

  std::ifstream file5("data/k_h.txt");
  std::vector<double> k_h;

  while(!file5.eof())
  {
  file5 >> temp_read;
  k_h.push_back(temp_read);
  }

  std::cout << "Finished Reading kh" << "\n";

  std::ifstream file6("data/p_h.txt");
  std::vector<double> p_h;

  while(!file6.eof())
  {
  file6 >> temp_read;
  p_h.push_back(temp_read);
  }

  std::cout << "Finished Reading ph" << "\n";

  xi_mag.pop_back();
  xj_mag.pop_back();
  xk_mag.pop_back();
  k_h.pop_back();
  p_h.pop_back();

  std::cout << "Finished Reading" << "\n";

  int N_kl_x;
  int N_kl_y;
  int N_kl_z;

  std::complex<double> Int;

  int sum_index = 54330094;
  int i;

  double h = 0.6777;
  double upperlimit = 470.*PI/(4000./h);
  double lowerlimit = 38.*PI/(4000./h);

  double M_kl;
  double M_ks;
  double M_ks_counterpart;
  double dot_ls;
  double dot_ls_counterpart;
  double F2;
  double F2_counterpart;
  // double f;
  double g;

  double xi_magg;
  double xj_magg;
  double xk_magg;

  std::ofstream myfile1 ("data/Integral_real.txt");
  std::ofstream myfile2 ("data/Integral_imag.txt");
  std::ofstream myfile3 ("data/delta_real.txt");
  std::ofstream myfile4 ("data/delta_imag.txt");

  for(N_kl_x=-10;N_kl_x<=10;N_kl_x=N_kl_x+2){
      for(N_kl_y=-10;N_kl_y<=10;N_kl_y=N_kl_y+2){
      	for(N_kl_z=-10;N_kl_z<=10;N_kl_z=N_kl_z+2){
          M_kl = std::sqrt(N_kl_x*N_kl_x+N_kl_y*N_kl_y+N_kl_z*N_kl_z);
      		Int = 0.;

      		for(i=0;i<sum_index;i++){
            xi_magg = xi_mag[i];
            xj_magg = xj_mag[i];
            xk_magg = xk_mag[i];
            M_ks = std::sqrt(xi_magg*xi_magg+xj_magg*xj_magg+xk_magg*xk_magg);
            dot_ls = N_kl_x*xi_magg+N_kl_y*xj_magg+N_kl_z*xk_magg;
            M_ks_counterpart = std::sqrt((xi_magg-N_kl_x)*(xi_magg-N_kl_x)+(xj_magg-N_kl_y)*(xj_magg-N_kl_y)+(xk_magg-N_kl_z)*(xk_magg-N_kl_z));
            dot_ls_counterpart = -N_kl_x*(xi_magg-N_kl_x)-N_kl_y*(xj_magg-N_kl_y)-N_kl_z*(xk_magg-N_kl_z);
            F2 = 5./7.-0.5*dot_ls*(M_kl/M_ks+M_ks/M_kl)/(M_kl*M_ks)+2.*dot_ls*dot_ls/(7.*M_kl*M_kl*M_ks*M_ks);
            F2_counterpart = 5./7.-0.5*dot_ls_counterpart*(M_kl/M_ks_counterpart+M_ks_counterpart/M_kl)/(M_kl*M_ks_counterpart)+2.*dot_ls_counterpart*dot_ls_counterpart/(7.*M_kl*M_kl*M_ks_counterpart*M_ks_counterpart);
            g = F2/(2.*matter_power_lin(k_h,p_h,M_ks_counterpart*delta_k,true))+F2_counterpart/(2.*matter_power_lin(k_h,p_h,M_ks*delta_k,true));
      			Int = Int + g*cfield2[indexing(int(xi_magg))][indexing(int(xj_magg))][indexing(int(xk_magg))]*cfield2[indexing(-int(xi_magg)+N_kl_x)][indexing(-int(xj_magg)+N_kl_y)][indexing(-int(xk_magg)+N_kl_z)];
            // f=F2*matter_power_lin(k_h,p_h,M_ks*delta_k,true)+F2_counterpart*matter_power_lin(k_h,p_h,M_ks_counterpart*delta_k,true);
            // std::cout << "f*d=" << f*cfield2[640+N_kl_x/2][640+N_kl_y/2][640+N_kl_z/2] << " dot=" << cfield2[xi[i]][xj[i]][xk[i]]*cfield2[1280+N_kl_x/2-xi[i]][1280+N_kl_y/2-xj[i]][1280+N_kl_z/2-xk[i]];
      		}

      		std::cout << 4000.*4000.*4000.*Int*4.*PI*((upperlimit/h)*(upperlimit/h)*(upperlimit/h)-(lowerlimit/h)*(lowerlimit/h)*(lowerlimit/h))/(3.*sum_index*(2.*PI)*(2.*PI)*(2.*PI)) << " " << N_kl_x << " " << N_kl_y << " " << N_kl_z << " " << cfield2[indexing(N_kl_x)][indexing(N_kl_y)][indexing(N_kl_z)] << "\n";
          myfile1 << real(4000.*4000.*4000.*Int*4.*PI*((upperlimit/h)*(upperlimit/h)*(upperlimit/h)-(lowerlimit/h)*(lowerlimit/h)*(lowerlimit/h))/(3.*sum_index*(2.*PI)*(2.*PI)*(2.*PI))) << "\n";
          myfile2 << imag(4000.*4000.*4000.*Int*4.*PI*((upperlimit/h)*(upperlimit/h)*(upperlimit/h)-(lowerlimit/h)*(lowerlimit/h)*(lowerlimit/h))/(3.*sum_index*(2.*PI)*(2.*PI)*(2.*PI))) << "\n";
          myfile3 << real(cfield2[indexing(N_kl_x)][indexing(N_kl_y)][indexing(N_kl_z)]) << "\n";
          myfile4 << imag(cfield2[indexing(N_kl_x)][indexing(N_kl_y)][indexing(N_kl_z)]) << "\n";

      	}
      }
  }

  myfile1.close();
  myfile2.close();
  myfile3.close();
  myfile4.close();

  cfield2.clear();
  cfield2.shrink_to_fit();

  xi_mag.clear();
  xi_mag.shrink_to_fit();
  xj_mag.clear();
  xj_mag.shrink_to_fit();
  xk_mag.clear();
  xk_mag.shrink_to_fit();
  k_h.clear();
  k_h.shrink_to_fit();
  p_h.clear();
  p_h.shrink_to_fit();


  return 0.;
}
