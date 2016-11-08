#include <Rcpp.h>
#include <string>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector extract_vector(NumericMatrix mat,int m){
  NumericVector mat_row(mat.ncol());
  for(int i=0;i<mat.ncol();i++){
    mat_row[i] = mat(m,i);
  }
  return mat_row;
}
// [[Rcpp::export]]
double SumC(NumericVector vector){
  int n = vector.size();
  double total = 0;
  for(int i=0;i<n;i++){
    total += vector[i];
  }
  return total;
}
// [[Rcpp::export]]
double corC(NumericVector X,NumericVector Y){
  int i;
  double temp1=0,temp2=0;
  double X_mean,Y_mean;
  double fenzi=0;
  double fenmu,coeff;
  if(X.size() != Y.size()){
    printf("error dim not the same!");
    return 0;
  }
  X_mean = SumC(X)/X.size();
  Y_mean = SumC(Y)/Y.size();
  for(i=0;i<X.size();i++){
    fenzi += (X[i]-X_mean)*(Y[i]-Y_mean);
  }
  for(i=0;i<X.size();i++){
     temp1 += pow((X[i]-X_mean),2);
  }
  for(i=0;i<Y.size();i++){
     temp2 += pow((Y[i]-Y_mean),2);
  }
  fenmu = sqrt(temp1*temp2);
  coeff = fenzi / fenmu;
  return coeff;
}

// [[Rcpp::export]]
DataFrame create_cor(NumericMatrix a,NumericMatrix b,CharacterVector name1,CharacterVector name2){
  int n = a.nrow()*b.nrow();
  NumericVector cor_vector(n);
  CharacterVector namei(n),namej(n);
  for(int i=0;i<a.nrow();i++){
    for(int j=0;j<b.nrow();j++){
      namei[i*b.nrow()+j] = name1[i];
      namej[i*b.nrow()+j] = name2[j];        
      cor_vector[i*b.nrow()+j] = corC(extract_vector(a,i),extract_vector(b,j));
    }
  }
  return DataFrame::create(_["cor_name1"]=namei, _["cor_name2"]=namej, _["cor"]=cor_vector);
}


//NumericVector create_cor(NumericMatrix a,NumericMatrix b){
//DataFrame create_cor(NumericMatrix a,NumericMatrix b,CharacterVector name1,CharacterVector name2){
