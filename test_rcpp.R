#2016-8-9 for test Rcpp:
library(Rcpp)
sourceCpp('test.cpp')
#R code
fib_R_1 <- function(n){
  if(n == 1 | n == 2){
    return(1)
  }else{
    return(fib_R_1(n-1)+fib_R_1(n-2))
  }
}
#test 30
system.time(fib_cpp_1(30))
system.time(fib_R_1(30))
#test 50
system.time(fib_R_1(50))
system.time(fib_cpp_1(50))
#--------use cppFunction-------------------
#scalar input,scalar output
cppFunction('int add(int x,int y){
return(x+y);
}')
#vector input,scalar output
cppFunction('double sumC(NumericVector x){
  int n = x.size();  // same as length()
  double total = 0;
  for(int i=0;i<n;i++){
    total += x[i];
  }
  return total;
}')
#R 内置sum函数是使用底层C编写的 效率很高
system.time(
  sum(runif(10^8))
)
system.time(
  sumC(runif(10^8))
)
#vector input vector output
cppFunction('NumericVector pdistC(double x,NumericVector ys){
  int n = ys.size();
  NumericVector out(n); //make the same length of ys
  for(int i=0;i<n;i++){
    out[i] = sqrt(pow(ys[i]-x,2.0));
  }
  return out;
}')
#matrix input vector output
cppFunction('NumericVector rowSumsC(NumericMatrix x){
              int nrow = x.nrow();
              int ncol = x.ncol();
              NumericVector out(nrow);
              for(int i = 0;i < nrow;i++){
                double total = 0;
                for(int j = 0;j < ncol;j++){
                  total += x(i,j);
                }
                out[i] = total;
              }
              return out;}')
set.seed(101)
x <- matrix(sample(100),10)
rowSums(x)
rowSumsC(x)
#-------put a R function into C++ file-------


