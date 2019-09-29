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
#define delta_k 0.0012566370614359172
#define linear_growth 3.06854

int indexing(int k_mag){
  if (k_mag <= 1){
    return - k_mag/2;
  }
  else{
    return 800 - k_mag/2;
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
  const std::complex<double> rr_lg(3.06854,0.0);  
  const std::complex<double> ii_lg(0.0,3.06854);  

  double temp_read;
  
  // std::ifstream file3("data/cfield_real_t.txt");
  // std::vector<double> cfield_real;
  // while(!file3.eof())
  // {
  //   file3 >> temp_read;
  //   cfield_real.push_back(temp_read);
  // }

  // std::cout << "Finished Reading cr" << "\n";

  // cfield_real.pop_back();

  // std::vector<std::vector<std::vector<std::complex<double> > > > cfield(800,std::vector<std::vector<std::complex<double> > >(800,std::vector <std::complex<double> >(401,(0.,0.))));

  int ll;
  int mm;
  int nn;
  // for(ll=0;ll<800;ll++){
  //   std::cout << "COMPLEXITY" << ll << "\n";
  //   for(mm=0;mm<800;mm++){
  //       for(nn=0;nn<401;nn++){
  //         cfield[ll][mm][nn] = rr*cfield_real[ll*800*401+mm*401+nn];
  //     }
  //   }
  // }

  // cfield_real.clear();
  // cfield_real.shrink_to_fit();

  // std::ifstream file4("data/cfield_imag_t.txt");
  // std::vector<double> cfield_imag;
  // while(!file4.eof())
  // {
  //   file4 >> temp_read;
  //   cfield_imag.push_back(temp_read);
  // }

  // std::cout << "Finished Reading ci" << "\n";

  // cfield_imag.pop_back();

  // for(ll=0;ll<800;ll++){
  //   std::cout << "COMPLEXITY" << ll << "\n";
  //   for(mm=0;mm<800;mm++){
  //       for(nn=0;nn<401;nn++){
  //         cfield[ll][mm][nn] = cfield[ll][mm][nn]+ii*cfield_imag[ll*800*401+mm*401+nn];
  //     }
  //   }
  // }

  // cfield_imag.clear();
  // cfield_imag.shrink_to_fit();

  // std::vector<std::vector<std::vector<std::complex<double> > > > cfield2(800,std::vector<std::vector<std::complex<double> > >(800,std::vector <std::complex<double> >(800,(0.,0.))));

  // for(ll=0;ll<800;ll++){
  //   std::cout << "COPYING" << ll << "\n";
  //   for(mm=0;mm<800;mm++){
  //     for(nn=0;nn<401;nn++){
  //       cfield2[ll][mm][nn] = cfield[ll][mm][nn];
  //     }
  //   }
  // }

  // cfield.clear();
  // cfield.shrink_to_fit();

  // for(ll=1;ll<800;ll++){
  //   std::cout << "READING " << ll << "\n";
  //   for(mm=1;mm<800;mm++){
  //     for(nn=401;nn<800;nn++){
  //       cfield2[ll][mm][nn] = std::conj(cfield2[800-ll][800-mm][800-nn]);
  //     }
  //   }
  // }

  // for(ll=1;ll<800;ll++){
  //   std::cout << "READING1" << ll << "\n";
  //   for(nn=401;nn<800;nn++){
  //     cfield2[ll][0][nn] = std::conj(cfield2[800-ll][0][800-nn]);
  //   }
  // }

  // for(mm=1;mm<800;mm++){
  //   std::cout << "READING2" << mm << "\n";
  //   for(nn=401;nn<800;nn++){
  //     cfield2[0][mm][nn] = std::conj(cfield2[0][800-mm][800-nn]);
  //   }
  // }

  // for(nn=401;nn<800;nn++){
  //   std::cout << "READING3" << nn << "\n";
  //   cfield2[0][0][nn] = std::conj(cfield2[0][0][800-nn]);
  // }

  // /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  std::ifstream file33("data/cfield_real_t.txt");
  std::vector<double> cfield_real_ini;
  while(!file33.eof())
  {
    file33 >> temp_read;
    cfield_real_ini.push_back(temp_read);
  }

  std::cout << "Finished Reading cr ini" << "\n";

  cfield_real_ini.pop_back();

  std::vector<std::vector<std::vector<std::complex<double> > > > cfield_ini(800,std::vector<std::vector<std::complex<double> > >(800,std::vector <std::complex<double> >(401,(0.,0.))));

  for(ll=0;ll<800;ll++){
    std::cout << "COMPLEXITY ini" << ll << "\n";
    for(mm=0;mm<800;mm++){
        for(nn=0;nn<401;nn++){
          cfield_ini[ll][mm][nn] = rr*cfield_real_ini[ll*800*401+mm*401+nn];
      }
    }
  }

  cfield_real_ini.clear();
  cfield_real_ini.shrink_to_fit();

  std::ifstream file44("data/cfield_imag_t.txt");
  std::vector<double> cfield_imag_ini;
  while(!file44.eof())
  {
    file44 >> temp_read;
    cfield_imag_ini.push_back(temp_read);
  }

  std::cout << "Finished Reading ci ini" << "\n";

  cfield_imag_ini.pop_back();

  for(ll=0;ll<800;ll++){
    std::cout << "COMPLEXITY ini" << ll << "\n";
    for(mm=0;mm<800;mm++){
        for(nn=0;nn<401;nn++){
          cfield_ini[ll][mm][nn] = cfield_ini[ll][mm][nn]+ii*cfield_imag_ini[ll*800*401+mm*401+nn];
      }
    }
  }

  cfield_imag_ini.clear();
  cfield_imag_ini.shrink_to_fit();

  std::vector<std::vector<std::vector<std::complex<double> > > > cfield_ini2(800,std::vector<std::vector<std::complex<double> > >(800,std::vector <std::complex<double> >(800,(0.,0.))));

  for(ll=0;ll<800;ll++){
    std::cout << "COPYING ini" << ll << "\n";
    for(mm=0;mm<800;mm++){
      for(nn=0;nn<401;nn++){
        cfield_ini2[ll][mm][nn] = cfield_ini[ll][mm][nn];
      }
    }
  }

  cfield_ini.clear();
  cfield_ini.shrink_to_fit();

  for(ll=1;ll<800;ll++){
    std::cout << "READING ini" << ll << "\n";
    for(mm=1;mm<800;mm++){
      for(nn=401;nn<800;nn++){
        cfield_ini2[ll][mm][nn] = std::conj(cfield_ini2[800-ll][800-mm][800-nn]);
      }
    }
  }

  for(ll=1;ll<800;ll++){
    std::cout << "READING1 ini" << ll << "\n";
    for(nn=401;nn<800;nn++){
      cfield_ini2[ll][0][nn] = std::conj(cfield_ini2[800-ll][0][800-nn]);
    }
  }

  for(mm=1;mm<800;mm++){
    std::cout << "READING2 ini" << mm << "\n";
    for(nn=401;nn<800;nn++){
      cfield_ini2[0][mm][nn] = std::conj(cfield_ini2[0][800-mm][800-nn]);
    }
  }

  for(nn=401;nn<800;nn++){
    std::cout << "READING3 ini" << nn << "\n";
    cfield_ini2[0][0][nn] = std::conj(cfield_ini2[0][0][800-nn]);
  }

    std::ifstream file0("data/xi_mag1.txt");
  std::vector<int> xi_mag;

  while(!file0.eof())
  {
  file0 >> temp_read;
  xi_mag.push_back(temp_read);
  }

  std::cout << "Finished Reading xi" << "\n";

  std::ifstream file1("data/xj_mag1.txt");
  std::vector<int> xj_mag;

  while(!file1.eof())
  {
  file1 >> temp_read;
  xj_mag.push_back(temp_read);
  }

  std::cout << "Finished Reading xj" << "\n";

  std::ifstream file2("data/xk_mag1.txt");
  std::vector<int> xk_mag;

  while(!file2.eof())
  {
  file2 >> temp_read;
  xk_mag.push_back(temp_read);
  }

  std::cout << "Finished Reading xk" << "\n";

  std::ifstream file5("../Example/data/k_h_lin.txt");
  std::vector<double> k_h;

  while(!file5.eof())
  {
  file5 >> temp_read;
  k_h.push_back(temp_read);
  }

  std::cout << "Finished Reading kh" << "\n";

  std::ifstream file6("../Example/data/p_h_lin.txt");
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

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  std::cout << "Finished Reading" << "\n";

  int N_kl_x;
  int N_kl_y;
  int N_kl_z;

  // int N_ks_x;
  // int N_ks_y;
  // int N_ks_z;

  // std::complex<double> Int;

  // int i;

  // double h = 0.6777;

  // double M_kl;
  // double M_ks_counterpart;
  // double dot_ls;
  // double F2;
  // double comb;

  // int upper_int = 400;
  // double upperlimit = delta_k * upper_int;

  // int sum_index;

  // std::ofstream myfile1 ("final_test/rhs_real.txt");
  // std::ofstream myfile2 ("final_test/rhs_imag.txt");
  // std::ofstream myfile3 ("final_test/lhs_real.txt");
  // std::ofstream myfile4 ("final_test/lhs_imag.txt");

  // for(N_ks_x=92;N_ks_x<=200;N_ks_x=N_ks_x+12){
  //     for(N_ks_y=92;N_ks_y<=200;N_ks_y=N_ks_y+12){
  //      for(N_ks_z=-92;N_ks_z>=-200;N_ks_z=N_ks_z-12){
  //        Int = 0.;
  //        sum_index = 0;

  //         for(N_kl_x=-upper_int;N_kl_x<=upper_int;N_kl_x=N_kl_x+2){
  //          for(N_kl_y=-upper_int;N_kl_y<=upper_int;N_kl_y=N_kl_y+2){
  //            for(N_kl_z=-upper_int;N_kl_z<=upper_int;N_kl_z=N_kl_z+2){
  //              if(((N_kl_x*N_kl_x+N_kl_y*N_kl_y+N_kl_z*N_kl_z)<=upper_int*upper_int)&&((N_kl_x*N_kl_x+N_kl_y*N_kl_y+N_kl_z*N_kl_z)!=0)&&(((N_ks_x-N_kl_x)*(N_ks_x-N_kl_x)+(N_ks_y-N_kl_y)*(N_ks_y-N_kl_y)+(N_ks_z-N_kl_z)*(N_ks_z-N_kl_z))!=0)){
  //                sum_index = sum_index+1;
  //                M_kl = std::sqrt(N_kl_x*N_kl_x+N_kl_y*N_kl_y+N_kl_z*N_kl_z);
  //                dot_ls = N_kl_x*N_ks_x+N_kl_y*N_ks_y+N_kl_z*N_ks_z;
  //                M_ks_counterpart = std::sqrt((N_ks_x-N_kl_x)*(N_ks_x-N_kl_x)+(N_ks_y-N_kl_y)*(N_ks_y-N_kl_y)+(N_ks_z-N_kl_z)*(N_ks_z-N_kl_z));
  //                comb = (dot_ls-M_kl*M_kl)/(M_ks_counterpart*M_kl);
  //                F2 = 5./7.+0.5*comb*(M_kl/M_ks_counterpart+M_ks_counterpart/M_kl)+2.*comb*comb/7.;
  //                Int = Int + F2*cfield_ini2[indexing(N_kl_x)][indexing(N_kl_y)][indexing(N_kl_z)]*cfield_ini2[indexing(-N_kl_x+N_ks_x)][indexing(-N_kl_y+N_ks_y)][indexing(-N_kl_z+N_ks_z)];
  //              }
  //            }
  //           }
  //        // }

  //        std::cout << 2500.*2500.*2500.*Int*4.*PI*upperlimit*upperlimit*upperlimit/(3.*sum_index*(2.*PI)*(2.*PI)*(2.*PI)) << " " << N_ks_x << " " << N_ks_y << " " << N_ks_z << " " << cfield2[indexing(N_ks_x)][indexing(N_ks_y)][indexing(N_ks_z)]-cfield_ini2[indexing(N_ks_x)][indexing(N_ks_y)][indexing(N_ks_z)] << "\n";
  //        myfile1 << real(2500.*2500.*2500.*Int*4.*PI*upperlimit*upperlimit*upperlimit/(3.*sum_index*(2.*PI)*(2.*PI)*(2.*PI))) << "\n";
  //         myfile2 << imag(2500.*2500.*2500.*Int*4.*PI*upperlimit*upperlimit*upperlimit/(3.*sum_index*(2.*PI)*(2.*PI)*(2.*PI))) << "\n";
  //         myfile3 << real(cfield2[indexing(N_ks_x)][indexing(N_ks_y)][indexing(N_ks_z)]-cfield_ini2[indexing(N_ks_x)][indexing(N_ks_y)][indexing(N_ks_z)]) << "\n";
  //         myfile4 << imag(cfield2[indexing(N_ks_x)][indexing(N_ks_y)][indexing(N_ks_z)]-cfield_ini2[indexing(N_ks_x)][indexing(N_ks_y)][indexing(N_ks_z)]) << "\n";
  //      }
  //     }
  // }


  std::complex<double> Int;

  // int sum_index = 13265430;
  int sum_index = 2860534;
  // int sum_index = 99068;
  int i;

  double h = 0.6777;
  double upperlimit = (3./5.)*470.*PI/(4000./h);
  double lowerlimit = 38.*PI/(4000./h);

  double M_kl;
  double M_ks;
  double M_ks_counterpart;
  double dot_ls;
  double dot_ls_counterpart;
  double F2;
  double F2_counterpart;
  double f;
  double g;

  double xi_magg;
  double xj_magg;
  double xk_magg;

  double eps = 3./7.;

  // std::ofstream myfile1 ("lin/delta_real.txt");
  // std::ofstream myfile2 ("lin/delta_imag.txt");
  // std::ofstream myfile3 ("lin/Integral_real.txt");
  // std::ofstream myfile4 ("lin/Integral_imag.txt");

  for(N_kl_x=-6;N_kl_x<=6;N_kl_x=N_kl_x+2){
      for(N_kl_y=-6;N_kl_y<=6;N_kl_y=N_kl_y+2){
        for(N_kl_z=-6;N_kl_z<=6;N_kl_z=N_kl_z+2){
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
            // Int = Int + g*(cfield_ini2[indexing(int(xi_magg))][indexing(int(xj_magg))][indexing(int(xk_magg))])*(cfield_ini2[indexing(-int(xi_magg)+N_kl_x)][indexing(-int(xj_magg)+N_kl_y)][indexing(-int(xk_magg)+N_kl_z)]);

            // Int = Int + g*(cfield2[indexing(int(xi_magg))][indexing(int(xj_magg))][indexing(int(xk_magg))]-cfield_ini2[indexing(int(xi_magg))][indexing(int(xj_magg))][indexing(int(xk_magg))])*(cfield2[indexing(-int(xi_magg)+N_kl_x)][indexing(-int(xj_magg)+N_kl_y)][indexing(-int(xk_magg)+N_kl_z)]-cfield_ini2[indexing(-int(xi_magg)+N_kl_x)][indexing(-int(xj_magg)+N_kl_y)][indexing(-int(xk_magg)+N_kl_z)]);
            f=F2*matter_power_lin(k_h,p_h,M_ks*delta_k,true)+F2_counterpart*matter_power_lin(k_h,p_h,M_ks_counterpart*delta_k,true);
            Int = Int + f*g;
          }

          std::cout << 1./(Int*4.*PI*((upperlimit/h)*(upperlimit/h)*(upperlimit/h)-(lowerlimit/h)*(lowerlimit/h)*(lowerlimit/h))/(3.*sum_index*(2.*PI)*(2.*PI)*(2.*PI))) << " " << N_kl_x << " " << N_kl_y << " " << N_kl_z << " " << cfield_ini2[indexing(N_kl_x)][indexing(N_kl_y)][indexing(N_kl_z)] << "\n";



        }
      }
  }



  // cfield2.clear();
  // cfield2.shrink_to_fit();
  cfield_ini2.clear();
  cfield_ini2.shrink_to_fit();
  
  return 0.;
}
