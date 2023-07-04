#include <Rcpp.h>
//#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector fExpCpp(NumericVector x, NumericVector param, int t_diff, double max_diam){
  int n = x.size();
  NumericVector q = (max_diam - x) * (1-exp(-t_diff*(param[0]*log(x)+param[1])));
  NumericVector z(n);
  double mean_z = 0;
  int nq = 0;
  
  for(int i=0;i<n;i++){
    if(q[i]>0) {
      z[i] = log(q[i]);
      mean_z = mean_z + log(q[i]);
      nq++;
    } else {
      z[i] = 0;
    }
  }
  
  if(nq != 0){
    
    mean_z = mean_z / nq;
    
    for(int i=0;i<n;i++){
      if(z[i] == 0){
        z[i] = mean_z;
      }
    }
    
  }
  return z;
  
}

// [[Rcpp::export]]
double quadTrapezCpp_1(NumericVector q, double h, int nx){
  double int_q = sum(q) - 0.5 * (q[0]+q[nx-1]);
  return int_q * h;
}

// [[Rcpp::export]]
NumericVector quadTrapezCpp_2(NumericMatrix q, double h, int nx){
  int ncol = q.ncol();
  NumericVector int_q(ncol);
  for(int j=0; j<ncol;j++){
    int_q[j] = sum(q(_,j));
  }
  int_q = int_q - 0.5 * (q(0,_)+q(nx-1,_));
  return h*int_q;
}

// [[Rcpp::export]]
NumericVector IPMprobSurvivalCpp(NumericVector x, NumericVector param){
  return 1/(1+exp(-(param[0]*log(x)+param[1])));
}

//[[Rcpp::export]]
NumericMatrix my_dlnorm( NumericMatrix x, NumericVector means, NumericVector sds){
    int nrow = x.nrow(), ncol = x.ncol() ;
    NumericMatrix res(nrow,ncol) ;
    for( int i=0; i<nrow; i++){
      res(i,_) = dlnorm(x(i,_),means[i],sds[i]);
    }
    return res ;
}

// [[Rcpp::export]]
NumericMatrix IPMpdfTreeGrowthCpp(NumericVector x, NumericVector y, NumericVector param, int t_diff, double max_diam, int nx, NumericMatrix y_minus_x){
  NumericVector f_gr = fExpCpp(x,param,t_diff,max_diam);
  NumericVector variance = exp(param[2] + param[3]*log(x));
  NumericMatrix z = my_dlnorm(y_minus_x,f_gr,sqrt(variance));
  return z;
}

// [[Rcpp::export]]
NumericVector IPMadultSurvivalTimesGrowthCpp(NumericVector ntrees, NumericVector x, NumericVector y, NumericVector param_survival, NumericVector param_growth, int t_diff, double max_diam, double h, int nx, NumericMatrix y_minus_x){
  NumericVector q = ntrees * IPMprobSurvivalCpp(x,param_survival);
  NumericMatrix ipm_growth = IPMpdfTreeGrowthCpp(x,y,param_growth,t_diff,max_diam,nx,y_minus_x);
  NumericMatrix res(ntrees.size(),ntrees.size());
  
  for(int j=0;j<ntrees.size();j++){
    res(_,j) = q * ipm_growth(_,j);
  }
  
  NumericVector res2 = quadTrapezCpp_2(res,h,nx);
  
  return res2;
}

// [[Rcpp::export]]
NumericVector IPMIngrowthIdentityCpp(NumericVector y, NumericVector param){
  return (param[0]*exp(-y*param[0])*param[1]);
}

// [[Rcpp::export]]
double IPMSaplingsCpp(double nsapl, NumericVector param){
  return (param[0]*param[1]);
}

// [[Rcpp::export]]
NumericVector CoefTimesVarCpp(Rcpp::List param, NumericVector z, NumericVector s = 0, NumericVector spba = 0, std::string type = "survival") {
  
  int n = Rcpp::as<NumericVector>(param[0]).size();
  NumericVector res(n);
  
  for(int i = 0;i<n;i++){
    
    if (type == "ingrowth"){
      
      res[i] = as<NumericVector>(param["intercept"])[i] + 
      as<NumericVector>(param["rain"])[i]*z[0] + 
      as<NumericVector>(param["temper"])[i]*z[1] +
      as<NumericVector>(param["basal.area"])[i]*z[2] +
      as<NumericVector>(param["anom.rain"])[i]*z[3] +
      as<NumericVector>(param["anom.temper"])[i]*z[4] +
      as<NumericVector>(param["rain.anom.temper"])[i]*z[0]*z[4];
      
      if(s.size() > 1){
        res[i] += as<NumericVector>(param["saplings"])[i]*s[i];
      }
      
      res[i] = res[i] * (10000/(M_PI*25));
      
    }else if (type == "growth" | type == "survival"){
      
      res[i] = as<NumericVector>(param["intercept"])[i] + 
      as<NumericVector>(param["rain"])[i]*z[0] + 
      as<NumericVector>(param["temper"])[i]*z[1] +
      as<NumericVector>(param["basal.area"])[i]*z[2] +
      as<NumericVector>(param["anom.rain"])[i]*z[3] +
      as<NumericVector>(param["anom.temper"])[i]*z[4] +
      as<NumericVector>(param["rain.anom.temper"])[i]*z[0]*z[4];
      
    }else if (type == "recruitment"){
      
      res[i] = as<NumericVector>(param["intercept"])[i] + 
      as<NumericVector>(param["rain"])[i]*z[0] + 
      as<NumericVector>(param["temper"])[i]*z[1] +
      as<NumericVector>(param["basal.area"])[i]*z[2] +
      as<NumericVector>(param["anom.rain"])[i]*z[3] +
      as<NumericVector>(param["anom.temper"])[i]*z[4] +
      as<NumericVector>(param["rain.anom.temper"])[i]*z[0]*z[4] + 
      as<NumericVector>(param["basal.area.sp"])[i]*spba[i];
      
    }
  }
  return res;
}
