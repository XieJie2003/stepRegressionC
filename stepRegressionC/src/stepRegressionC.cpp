#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include<algorithm>
using namespace arma;

// [[Rcpp::export]]
arma::vec betaHat(arma::vec y, arma::mat X) {
vec betaHat, XtY;
mat XtX;
// equivalent to betaHat = solve(X,y);
XtX = X.t() * X;
XtY = X.t() * y;
// XtX * betaHat = XtY
betaHat = solve(XtX, XtY);

return betaHat;
}

// [[Rcpp::export]]
arma::mat MatrixDelete(arma::mat x,int p){
 int q=x.n_cols;
 if(p==q-1){
  return x.cols(0,q-2);
 }
 else{
  for(int i=p+1;i<q;i++){
   x.col(i-1)=x.col(i);
  }
  return x.cols(0,q-2);
 }
}

// [[Rcpp::export]]
double AIC(arma::vec y, mat x){
 double aic;
 int n = x.n_rows;
 int p = x.n_cols;
 vec betahat = betaHat(y,x);
 mat yHat=x * betahat;
 double rss=0;
 for(int i=0;i<n;i++){
  rss+=(yHat[i]-y[i])*(yHat[i]-y[i]);
 }
 aic=n*log(rss/n)+2*p;
 return aic;
}

// [[Rcpp::export]]
int search(arma::vec &deletE,int start,int end){
 int i = start, j = end, pivot=deletE[start];
 while (i < j)
	{
		while (i<j && deletE[j]>pivot) 
		{
			j--;
		}
		if (i < j)
		{
		  int temp=deletE[j];
		  deletE[j]=deletE[i];
		  deletE[i]=temp;
			i++; 
		}
		while (i < j && deletE[i] <= pivot) 
		{
			i++;
		}
		if (i < j)
		{
		  int temp=deletE[j];
		  deletE[j]=deletE[i];
		  deletE[i]=temp;
			j--;
		}
	}
	return i;
} 


// [[Rcpp::export]]
void quickSort(arma::vec &r, int low, int hight)
{
	int mid;
	if (low < hight)
	{
		mid = search(r, low, hight);  
		quickSort(r, low, mid - 1); 
		quickSort(r, mid+1, hight); 
	}
}

// [[Rcpp::export]]
int Reset(int t,arma::vec deletE){
 int p=deletE.size();
 quickSort(deletE,0,p-1);
 int r=0;
 for(int i=0;i<p;i++){
  if(deletE[i]==0){
   r++;
  }
  else{
   break;
  }
 }
 vec wholE=zeros<vec>(p);
 vec deletEr=zeros<vec>(p-r);
 for(int i=0;i<p-r;i++){
  deletEr[i]=deletE[r+i];
 }
 int q=deletEr.size();
 for(int i=0;i<q;i++){
  int w=deletEr[i];
  wholE[w-1]=w;
 }
 int v=0;
 for(int i=0;i<p;i++){
  if(wholE[i]==0){
   v++;
  }
  if(v==t){
   return i+1;
  }
 }
 return 0;
}

// [[Rcpp::export]]
arma::vec bcstepC2(arma::vec y, arma::mat x,int u,arma::vec deletE){
 int p = x.n_cols;
 double aic=AIC(y,x);
 int t=-1;
 for(int i=0;i<p;i++){
  mat x_bc=MatrixDelete(x,i);
  double aic2=AIC(y,x_bc);
  if(aic>aic2){
   aic=aic2;
   t=i;
  }
 }
 if(t!=-1){
  deletE[u]=Reset(t+1,deletE);
  u++;
  return bcstepC2(y,MatrixDelete(x,t),u,deletE);
 }
 else{
  int g=0;
  for(int i=0;i<p;i++){
   if(deletE[i]!=0){
    g++;
   }
   else{
    break;
   }
  }
  return deletE.subvec(0,g-1);
 }
}

// [[Rcpp::export]]
arma::vec bcstepC(arma::vec y, arma::mat x){
 int u=0;
 int p = x.n_cols;
 vec deletE=zeros<vec>(p);
 double aic=AIC(y,x);
 int t=-1;
 for(int i=0;i<p;i++){
  mat x_bc=MatrixDelete(x,i);
  double aic2=AIC(y,x_bc);
  if(aic>aic2){
   aic=aic2;
   t=i;
  }
 }
 if(t!=-1){
  deletE[u]=t+1;
  u++;
  return bcstepC2(y,MatrixDelete(x,t),u,deletE);
 }
 else if(t==-1){
  cout<<t<<endl;
  return deletE;
 }
}