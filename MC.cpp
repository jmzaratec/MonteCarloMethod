#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <fstream>
#include <ctime>

using namespace std;

class BlackScholes {

public:
	double normal_distribution(double x){
	// declare the number of steps
	int n;
	// declare Pi, h, minus_infty,normal_dist
	double Pi, h, minus_infty,normal_dist;
	// assign constants
	n=10000;
	Pi = 3.141592653589793;
	minus_infty = -10.;
	h = (x - minus_infty)/2./n;
	// summation and the current value of x
	double sum,x_i;
	// add the two ends y_0 and y_2n
	sum = exp(-0.5*minus_infty*minus_infty) + 
		exp(-0.5*x*x);
	// find x_i at position 1
	x_i = minus_infty + h;
	// sum = sum + 4y_1
	sum = sum + 4.*exp(-0.5*x_i*x_i);
	// now loop over i=1 to n-1 to add on all other terms
	for(int i=1;i<=n-1;i=i+1)
	{
	  // add on terms y_2, y_4, y_6 ,..., y_{n-2}
	  x_i = minus_infty + 2.*i*h;
	  sum = sum + 2.*exp(-0.5*x_i*x_i);
	  // add on terms y_3, y_5, ..., y_{n-1}
	  x_i = x_i + h;	
	  sum = sum + 4.*exp(-0.5*x_i*x_i);
	}
return 1/sqrt(2.*Pi)*h/3.*sum; }
// A function to compute the price of a European Call Option using BlackSholes formula
double call(double S0,double T, double t, double r, double q, double sigma, double X){
	double d1,d2, call;
	call=0;
	if(S0==0.){	call= max(S0-X,0.);
	return call;}
	else if (t>=T){call= max(S0-X,0.);
	return call;}
	else {
	//Black Scholes Formula
		d1 = (log(S0/X) +(r-q)*(T-t) + (sigma*sigma*0.5)*(T-t))/( sigma*sqrt(T-t));
		d2 = d1-(sigma*sqrt(T-t));
	//cout<<d1<<" "<<d2<<"\n";
		call  = S0*exp((-q)*(T-t))*normal_distribution(d1)-X*exp((-r)*(T-t))*normal_distribution(d2);
	//To avoid error in case t=>T
	return call;}
}

double put(double S0,double T, double t, double r, double q, double sigma, double X){
	double d1,d2, put;
	put=0;
	if(S0==0.){
	put= max(-S0+X,0.);
	return put;	}
	else if (t>=T){put= max(-S0+X,0.);
	return put;
	}
	else {
	//Black Scholes Formula
		d1 = (log(S0/X) +(r-q)*(T-t) + (sigma*sigma*0.5)*(T-t))/( sigma*sqrt(T-t));
		d2 = d1-(sigma*sqrt(T-t));
	//cout<<d1<<" "<<d2<<"\n";
		put  = -S0*exp((-q)*(T-t))*normal_distribution(-d1)+X*exp((-r)*(T-t))*normal_distribution(-d2);
	//To avoid error in case t=>T
return put;	}
}

double binCall(double S0,double T, double t, double r, double q, double sigma, double X){
	double d1,d2, bincall;

	if (t>=T && S0>X) {bincall= 1;
		return bincall;	}
	else if (t>=T && S0<X) {bincall= 0;
		return bincall;	}
	else if (S0==0.&& S0>X) {bincall= 1;
		return bincall;	}
	else if (S0==0.&& X==0) {bincall= 0;
		return bincall;	}
	else{
	//Black Scholes Formula
		d1 = (log(S0/X) +(r-q)*(T-t) + (sigma*sigma*0.5)*(T-t))/( sigma*sqrt(T-t));
		d2 = d1-(sigma*sqrt(T-t));
	//cout<<d1<<" "<<normal_distribution(d2)<<"\n";
		bincall  = 1*exp(-(r)*(T-t))*normal_distribution(d2);
	//To avoid error in case t=>T
	//if(S0==0.){	bincall= max(1.,0.);}
return bincall;}
}

double binPut(double S0,double T, double t, double r, double q, double sigma, double X){
	double d1,d2, binput;

	if (t>=T && S0>X) {binput= 0;
		return binput;	}
	else if (t>=T && S0<X) {binput= 1;
		return binput;	}
	else if (S0==0.&& S0>X) {binput= 0;
		return binput;	}
	else if (S0==0.&& X==0) {binput= 1;
		return binput;	}
	else{
	//Black Scholes Formula
		d1 = (log(S0/X) +(r-q)*(T-t) + (sigma*sigma*0.5)*(T-t))/( sigma*sqrt(T-t));
		d2 = d1-(sigma*sqrt(T-t));
	//cout<<d1<<" "<<normal_distribution(d2)<<"\n";
		binput  = 1*exp(-(r)*(T-t))*normal_distribution(-d2);
	//To avoid error in case t=>T
	//if(S0==0.){	bincall= max(1.,0.);}
return binput;}
}
// end of class
};
//discount factor for continuous fixed payments
double disF(double maturity, double interestRate){
	double N = 100000;
	double payment=1.;
	double dt= maturity/N;
	double r= interestRate;
	double PV=0;



	for (int i=1;i<=N;i++){	

	PV= PV+ (payment/N)*exp(-r*dt*i);
	//cout<<i<<"PV: "<<PV<<"\n";

}





	return PV;
}
class DBonds{
public:
	double NB, F, NS, S;
	double maturity,sigma,interestRate, S0, X;
// auxiliry function
double h(double N, double x){
	double h;
	
	h= 0.5+sqrt((0.25-(0.25*exp(-1*pow(x/(N+(0.3333333333333)),2)*(N+(0.16666666666666666))))));
	

	return h;
}
double f(double S,double X, double F)
{
	
	if(S<X)   { return F;}
	else {return 0.;}
}
double qstar(double maturity, double sigma, double r, double S0, double X, int N, double NB, double F){
double d1,d2;
double qstar;
 
d1 = (log(S0/X) +(r)*(maturity) + (sigma*sigma*0.5)*(maturity))/( sigma*sqrt(maturity));
d2 = (log(S0/X) +(r)*(maturity) - (sigma*sigma*0.5)*(maturity))/( sigma*sqrt(maturity));

qstar = h(N,d1);


return qstar;


}
double qnormal(double maturity, double sigma, double r, double S0, double X, int N, double NB, double F){
double d1,d2;
double qnormal;
 
d1 = (log(S0/X) +(r)*(maturity) + (sigma*sigma*0.5)*(maturity))/( sigma*sqrt(maturity));
d2 = (log(S0/X) +(r)*(maturity) - (sigma*sigma*0.5)*(maturity))/( sigma*sqrt(maturity));

qnormal=h(N,d2);

return qnormal;
}
double f_Put_heston(  double S, double X, double dt, double sigma){
	//double a;
	//a= max(S-X,0.);

	  // if nearest to node use different formula
	//which of my node is closest to strike price
	//cout<<"Value: "<<log(S)-log(X)<<" Com "<<sigma*sqrt(dt);
	if(log(S)-log(X)>(sigma*sqrt(dt))){   return 0;}
	else if(log(S)-log(X)<(-sigma*sqrt(dt))) {return  (X-(S*((exp(sigma*sqrt(dt))-exp(-sigma*sqrt(dt)))/(2*sigma*sqrt(dt)))));}
	else{return (((X*(sigma*sqrt(dt)-log((S/X))))/(2*sigma*sqrt(dt)))+(S/(2*sigma*sqrt(dt)))*(exp(-sigma*sqrt(dt))-(X/S)));}


}
//no sirve
double f_Put_hestonConvertible(  double S, double NB, double dt, double sigma, double F, double Z, double Cp){
	//double a;
	//a= max(S-X,0.);

	  // if nearest to node use different formula
	//which of my node is closest to strike price
	//cout<<"Value: "<<log(S)-log(X)<<" Com "<<sigma*sqrt(dt);
	  //min(tree[N][j]/NB,max(F,Z*tree[N][j]));
	
	  if(log(S)-log(F/Z)>(sigma*sqrt(dt))){   return  (-F/Z+(S*((exp(sigma*sqrt(dt))-exp(-sigma*sqrt(dt)))/(2*sigma*sqrt(dt))))) ;}
	  else if(log(S)-log(F/Z)<(sigma*sqrt(dt))&&log(S)-log(F/Z)>(-sigma*sqrt(dt))) {return  (((-F/Z*(sigma*sqrt(dt)+log((S/NB*F))))/(2*sigma*sqrt(dt)))-(S/(2*sigma*sqrt(dt)))*(exp(-sigma*sqrt(dt))+(NB*F/S)));}
	  else if(log(S)-log(F/Z)<(-sigma*sqrt(dt)) && log(S)-log(NB*F)>(sigma*sqrt(dt))) {return F ;}

	
	else if(log(S)-log(NB*F)>(sigma*sqrt(dt))&& log(S)-log(F/Z)<(-sigma*sqrt(dt))) {   return F;}
	else if(log(S)-log(NB*F)<(sigma*sqrt(dt))&& log(S)-log(F/Z)>(-sigma*sqrt(dt))) {return (((NB*F*(sigma*sqrt(dt)-log((S/NB*F))))/(2*sigma*sqrt(dt)))+(S/(2*sigma*sqrt(dt)))*(exp(-sigma*sqrt(dt))-(NB*F/S)));}
	else {return  (NB*F-(S*((exp(sigma*sqrt(dt))-exp(-sigma*sqrt(dt)))/(2*sigma*sqrt(dt)))));}
	



}

// Extrafunction
double binomialPut(double maturity, double sigma, double interestRate, double S0, double X, int N){

  vector<vector<double> > tree(N+1,vector<double>(N+1));
  // local variables
  double u,d,q,dt;
  // initialise variables
  dt=maturity/N;
  u = exp(sigma*sqrt(dt));
  d = exp(-sigma*sqrt(dt));
  q = (exp(interestRate*dt)-d)/(u-d);
  // calculate tree
  for(int i=0;i<tree.size();i++)
  {
	  for(int j=0;j<tree[i].size();j++)
	  {
		  tree[i][j] = S0*pow(u,j)*pow(d,i-j);// value of stock at i,j
		//  cout << "S_{"<<i<<"," << j << "} =";
	//	  cout << tree[i][j] << endl;
	  }
  }
  // create storage for the value tree
  vector<vector<double> > valueTree(N+1,vector<double>(N+1));
  // now set option value at maturity
  // set V[N][j]= payoff condition
  for (int j=0; j<valueTree[N].size();j++){
	  valueTree[N][j] = max(-tree[N][j]+X,0.);
	  // USE THE ADJUSTED PRICE
		//o.5*pow((Sj+0.5*dS-X) ,2)

	  //check values
	  //cout << "V_{"<<N<<"," << j << "} =";
	 // cout<<valueTree[N][j]<<endl;
	   }

  // use backward induction to fill in values back to (0,0)

  for(int i=N-1;i>=0;i--){

	  for(int j=0; j<=i;j++){
		   //calculate value of tree at this node
		  valueTree[i][j]= exp(-interestRate*dt)*(q*valueTree[i+1][j+1] + (1-q)*valueTree[i+1][j]);
		  //check values
	 // cout << "V_{"<<i<<"," << j << "} =";
	 // cout<<valueTree[i][j]<<endl;
	  }


  }
	
return valueTree[0][0];
}
double amePut(double maturity, double sigma, double interestRate, double S0, double X, int N, double div2){

  vector<vector<double> > tree(N+1,vector<double>(N+1));
  // local variables
  double u,d,q,dt;
  // initialise variables
  dt=maturity/(N); //dt=maturity/(N-1); in order to match results on internet
  u = exp(sigma*sqrt(dt));
  d = exp(-sigma*sqrt(dt));
  q = (exp((interestRate-div2)*dt)-d)/(u-d);
  // calculate tree
  for(int i=0;i<tree.size();i++)
  {
	  for(int j=0;j<tree[i].size();j++)
	  {
		  tree[i][j] = S0*pow(u,j)*pow(d,i-j);// value of stock at i,j
		  //tree[i][j] = max(tree[i][j],Z*tree[i][j]);// For dBonds
		 
		//  cout << "S_{"<<i<<"," << j << "} =";
		 // cout << tree[i][j] << endl;
	  }
  }
  // create storage for the value tree
  vector<vector<double> > valueTree(N+1,vector<double>(N+1));
  // now set option value at maturity
 
  // use backward induction to fill in values back to (0,0)

  for(int i=N-1;i>=0;i--){

	  for(int j=0; j<=i;j++){
		   //calculate value of tree at this node
		  double tem1=max(X-tree[i][j],0.);
			//cout<<i<<" "<<j<<"    "<<tem1<<"       "<<tree[i][j]<<"\n";
		  valueTree[i][j]= max(exp(-interestRate*dt)*(q*valueTree[i+1][j+1] + (1-q)*valueTree[i+1][j]),tem1);
		  
		  
		  //check values
	// cout << "V_{"<<i<<"," << j << "} =";
	 // cout<<valueTree[i][j]<<endl;
	  }
  }

   
	
return valueTree[0][0];
}
// DBonds
double DB_CRR(double maturity, double sigma, double interestRate, double S0, double X, int N, double NB, double F){

  vector<vector<double> > tree(N+1,vector<double>(N+1));
  // local variables
  double u,d,q,dt;
  // initialise variables
  dt=maturity/N;
  u = exp(sigma*sqrt(dt));
  d = exp(-sigma*sqrt(dt));
  q = (exp(interestRate*dt)-d)/(u-d);
  // calculate tree
  for(int i=0;i<tree.size();i++)
  {
	  for(int j=0;j<tree[i].size();j++)
	  {
		  tree[i][j] = S0*pow(u,j)*pow(d,i-j);// value of stock at i,j
		//  cout << "S_{"<<i<<"," << j << "} =";
	//	  cout << tree[i][j] << endl;
	  }
  }
  // create storage for the value tree
  vector<vector<double> > valueTree(N+1,vector<double>(N+1));
  // now set option value at maturity
  // set V[N][j]= payoff condition
  for (int j=0; j<valueTree[N].size();j++){
	  valueTree[N][j] = min(tree[N][j]/NB,F);
	   //valueTree[N][j] = max(-tree[N][j]+X,0.);
	  // USE THE ADJUSTED PRICE
		//o.5*pow((Sj+0.5*dS-X) ,2)

	  //check values
	  //cout << "V_{"<<N<<"," << j << "} =";
	 // cout<<valueTree[N][j]<<endl;
	   }

  // use backward induction to fill in values back to (0,0)

  for(int i=N-1;i>=0;i--){

	  for(int j=0; j<=i;j++){
		   //calculate value of tree at this node
		  valueTree[i][j]= exp(-interestRate*dt)*(q*valueTree[i+1][j+1] + (1-q)*valueTree[i+1][j]);
		  //check values
	 // cout << "V_{"<<i<<"," << j << "} =";
	 // cout<<valueTree[i][j]<<endl;
	  }


  }
	
return valueTree[0][0];
}
double DB_CRRdrift(double maturity, double sigma, double interestRate, double S0, double X, int N, double NB, double F){

  vector<vector<double> > tree(N+1,vector<double>(N+1));
  // local variables
  double u,d,q,dt;
   double drift = (log(X)- log(S0))/maturity;
  // initialise variables
  dt=maturity/N;
  u = exp(drift*dt + sigma*sqrt(dt));
  d = exp(drift*dt-sigma*sqrt(dt));
  q = (exp(interestRate*dt)-d)/(u-d);
  


  // calculate tree
  for(int i=0;i<tree.size();i++)
  {
	  for(int j=0;j<tree[i].size();j++)
	  {
		  tree[i][j] = S0*pow(u,j)*pow(d,i-j);// value of stock at i,j
		//  cout << "S_{"<<i<<"," << j << "} =";
	//	  cout << tree[i][j] << endl;
	  }
  }
  // create storage for the value tree
  vector<vector<double> > valueTree(N+1,vector<double>(N+1));
  // now set option value at maturity
  // set V[N][j]= payoff condition
  for (int j=0; j<valueTree[N].size();j++){
	  valueTree[N][j] = min(tree[N][j]/NB,F);
	   //valueTree[N][j] = max(-tree[N][j]+X,0.);
	  // USE THE ADJUSTED PRICE
		//o.5*pow((Sj+0.5*dS-X) ,2)

	  //check values
	  //cout << "V_{"<<N<<"," << j << "} =";
	 // cout<<valueTree[N][j]<<endl;
	   }

  // use backward induction to fill in values back to (0,0)

  for(int i=N-1;i>=0;i--){

	  for(int j=0; j<=i;j++){
		   //calculate value of tree at this node
		  valueTree[i][j]= exp(-interestRate*dt)*(q*valueTree[i+1][j+1] + (1-q)*valueTree[i+1][j]);
		  //check values
	 // cout << "V_{"<<i<<"," << j << "} =";
	 // cout<<valueTree[i][j]<<endl;
	  }


  }
	
return valueTree[0][0];
}
double DB_Leisen(double maturity, double sigma, double interestRate, double S0, double X, int N, double NB, double F){

  vector<vector<double> > tree(N+1,vector<double>(N+1));
  // local variables
  double u,d,q,dt;
    // initialise variables
  dt=maturity/N;
    u = exp(interestRate*dt)*(qstar(maturity,  sigma,  interestRate,  S0,  X,  N,  NB,  F)/qnormal( maturity,  sigma,  interestRate,  S0,  X,  N,  NB,  F)); //u = exp(sigma*sqrt(dt));
	d = (exp(interestRate*dt)-qnormal( maturity,  sigma,  interestRate,  S0,  X,  N,  NB,  F)*u)/(1-qnormal( maturity,  sigma,  interestRate,  S0,  X,  N,  NB,  F));//  d = exp(-sigma*sqrt(dt));
  q = (exp(interestRate*dt)-d)/(u-d);

  // calculate tree
  for(int i=0;i<tree.size();i++)
  {
	  for(int j=0;j<tree[i].size();j++)
	  {
		  tree[i][j] = S0*pow(u,j)*pow(d,i-j);// value of stock at i,j
		//  cout << "S_{"<<i<<"," << j << "} =";
	//	  cout << tree[i][j] << endl;
	  }
  }
  // create storage for the value tree
  vector<vector<double> > valueTree(N+1,vector<double>(N+1));
  // now set option value at maturity
  // set V[N][j]= payoff condition
  for (int j=0; j<valueTree[N].size();j++){
	  valueTree[N][j] = min(tree[N][j]/NB,F);
	   //valueTree[N][j] = max(-tree[N][j]+X,0.);
	  // USE THE ADJUSTED PRICE
		//o.5*pow((Sj+0.5*dS-X) ,2)

	  //check values
	  //cout << "V_{"<<N<<"," << j << "} =";
	 // cout<<valueTree[N][j]<<endl;
	   }

  // use backward induction to fill in values back to (0,0)

  for(int i=N-1;i>=0;i--){

	  for(int j=0; j<=i;j++){
		   //calculate value of tree at this node
		  valueTree[i][j]= exp(-interestRate*dt)*(q*valueTree[i+1][j+1] + (1-q)*valueTree[i+1][j]);
		  //check values
	 // cout << "V_{"<<i<<"," << j << "} =";
	 // cout<<valueTree[i][j]<<endl;
	  }


  }
	
return valueTree[0][0];
}
double DB_Put_Heston(double maturity, double sigma, double interestRate, double S0, double X, int N, double NB, double F){

  vector<vector<double> > tree(N+1,vector<double>(N+1));
  // local variables
  double u,d,q,dt;
  // initialise variables
  dt=maturity/N;
  u = exp(sigma*sqrt(dt));
  d = exp(-sigma*sqrt(dt));
  q = (exp(interestRate*dt)-d)/(u-d);
  // calculate tree
  for(int i=0;i<tree.size();i++)
  {
	  for(int j=0;j<tree[i].size();j++)
	  {
		  tree[i][j] = S0*pow(u,j)*pow(d,i-j);// value of stock at i,j
		//  cout << "S_{"<<i<<"," << j << "} =";
	//	  cout << tree[i][j] << endl;
	  }
  }
  // create storage for the value tree
  vector<vector<double> > valueTree(N+1,vector<double>(N+1));
  // now set option value at maturity
  // set V[N][j]= payoff condition
  for (int j=0; j<valueTree[N].size();j++){
	  valueTree[N][j] = f_Put_heston(tree[N][j],X,dt,sigma);//min(tree[N][j]/NB,F);
	   //valueTree[N][j] = max(-tree[N][j]+X,0.);
	  // USE THE ADJUSTED PRICE
		//o.5*pow((Sj+0.5*dS-X) ,2)

	  //check values
	  //cout << "V_{"<<N<<"," << j << "} =";
	 // cout<<valueTree[N][j]<<endl;
	   }

  // use backward induction to fill in values back to (0,0)

  for(int i=N-1;i>=0;i--){

	  for(int j=0; j<=i;j++){
		   //calculate value of tree at this node
		  valueTree[i][j]= exp(-interestRate*dt)*(q*valueTree[i+1][j+1] + (1-q)*valueTree[i+1][j]);
		  //check values
	 // cout << "V_{"<<i<<"," << j << "} =";
	 // cout<<valueTree[i][j]<<endl;
	  }


  }
	
return valueTree[0][0];
}
double DB_Tian(double maturity, double sigma, double interestRate, double S0, double X, int N, double NB, double F){

  vector<vector<double> > tree(N+1,vector<double>(N+1));
  // local variables
  double u,d,q,dt;
  // initialise variables
   dt=maturity/N;
   double M = exp(interestRate*dt);
   double V= exp(sigma*sigma*dt);
   u = 0.5*M*V*((V+1)+sqrt((V*V+2*V-3)));
  d =  0.5*M*V*((V+1)-sqrt((V*V+2*V-3)));
  q = (M-d)/(u-d);
 
  // calculate tree
  for(int i=0;i<tree.size();i++)
  {
	  for(int j=0;j<tree[i].size();j++)
	  {
		  tree[i][j] = S0*pow(u,j)*pow(d,i-j);// value of stock at i,j
		//  cout << "S_{"<<i<<"," << j << "} =";
	//	  cout << tree[i][j] << endl;
	  }
  }
  // create storage for the value tree
  vector<vector<double> > valueTree(N+1,vector<double>(N+1));
  // now set option value at maturity
  // set V[N][j]= payoff condition
  for (int j=0; j<valueTree[N].size();j++){
	  valueTree[N][j] = min(tree[N][j]/NB,F);
	   //valueTree[N][j] = max(-tree[N][j]+X,0.);
	  // USE THE ADJUSTED PRICE
		//o.5*pow((Sj+0.5*dS-X) ,2)

	  //check values
	  //cout << "V_{"<<N<<"," << j << "} =";
	 // cout<<valueTree[N][j]<<endl;
	   }

  // use backward induction to fill in values back to (0,0)

  for(int i=N-1;i>=0;i--){

	  for(int j=0; j<=i;j++){
		   //calculate value of tree at this node
		  valueTree[i][j]= exp(-interestRate*dt)*(q*valueTree[i+1][j+1] + (1-q)*valueTree[i+1][j]);
		  //check values
	 // cout << "V_{"<<i<<"," << j << "} =";
	 // cout<<valueTree[i][j]<<endl;
	  }


  }
	
return valueTree[0][0];
}
//CDN
double CDN_CRR(double maturity, double sigma, double interestRate, double S0, double X, int N, double NB, double F, double div){
	double e;
  vector<vector<double> > tree(N+1,vector<double>(N+1));
  // local variables
  double u,d,q,dt;
  // initialise variables
  dt=maturity/N;
  u = exp(sigma*sqrt(dt));
  d = exp(-sigma*sqrt(dt));
  q = (exp(interestRate*dt)-d)/(u-d);
  // calculate tree
  for(int i=0;i<tree.size();i++)
  {
	  for(int j=0;j<tree[i].size();j++)
	  {
		  tree[i][j] = S0*pow(u,j)*pow(d,i-j);// value of stock at i,j
		//  cout << "S_{"<<i<<"," << j << "} =";
	//	  cout << tree[i][j] << endl;
	  }
  }
  // create storage for the value tree
  vector<vector<double> > valueTree(N+1,vector<double>(N+1));
  // now set option value at maturity
  // set V[N][j]= payoff condition
  for (int j=0; j<valueTree[N].size();j++){
	  //valueTree[N][j] = min(tree[N][j]/NB,F);
	   //valueTree[N][j] = max(-tree[N][j]+X,0.);
	   //valueTree[N][j] = f(tree[N][j],X,F)-(0.05*(F/maturity)*dt);//
	     //valueTree[N][j] = f(tree[N][j]-(div*dt),X,F)-div*dt;//adding q
		valueTree[N][j] = f(tree[N][j],X,F);//adding q //this is 
		 // USE THE ADJUSTED PRICE
		//o.5*pow((Sj+0.5*dS-X) ,2)

	  //check values
	  //cout << "V_{"<<N<<"," << j << "} =";
	 // cout<<valueTree[N][j]<<endl;
	   }

  // use backward induction to fill in values back to (0,0)

  for(int i=N-1;i>=0;i--){

	  for(int j=0; j<=i;j++){
		   //calculate value of tree at this node
		  valueTree[i][j]= exp(-interestRate*dt)*(q*valueTree[i+1][j+1] + (1-q)*valueTree[i+1][j])-div*dt; // discounting the continous payment.
		  //check values
	 // cout << "V_{"<<i<<"," << j << "} =";
	 // cout<<valueTree[i][j]<<endl;
	  }


  }
	
return valueTree[0][0];//-div*maturity;
}
double CDN_CRRdrift(double maturity, double sigma, double interestRate, double S0, double X, int N, double NB, double F, double div){
	double e;
  vector<vector<double> > tree(N+1,vector<double>(N+1));
  // local variables
  double u,d,q,dt;
  double drift = (log(X)- log(S0))/maturity;
  // initialise variables
  dt=maturity/N;
  u = exp(drift*dt + sigma*sqrt(dt));
  d = exp(drift*dt-sigma*sqrt(dt));
  q = (exp(interestRate*dt)-d)/(u-d);
  

  // calculate tree
  for(int i=0;i<tree.size();i++)
  {
	  for(int j=0;j<tree[i].size();j++)
	  {
		  tree[i][j] = S0*pow(u,j)*pow(d,i-j);// value of stock at i,j
		//  cout << "S_{"<<i<<"," << j << "} =";
	//	  cout << tree[i][j] << endl;
	  }
  }
  // create storage for the value tree
  vector<vector<double> > valueTree(N+1,vector<double>(N+1));
  // now set option value at maturity
  // set V[N][j]= payoff condition
  for (int j=0; j<valueTree[N].size();j++){
	  //valueTree[N][j] = min(tree[N][j]/NB,F);
	   //valueTree[N][j] = max(-tree[N][j]+X,0.);
	   //valueTree[N][j] = f(tree[N][j],X,F)-(0.05*(F/maturity)*dt);//
	     //valueTree[N][j] = f(tree[N][j]-(div*dt),X,F)-div*dt;//adding q
		valueTree[N][j] = f(tree[N][j],X,F);//adding q //this is 
		 // USE THE ADJUSTED PRICE
		//o.5*pow((Sj+0.5*dS-X) ,2)

	  //check values
	  //cout << "V_{"<<N<<"," << j << "} =";
	 // cout<<valueTree[N][j]<<endl;
	   }


  // use backward induction to fill in values back to (0,0)

  for(int i=N-1;i>=0;i--){

	  for(int j=0; j<=i;j++){
		   //calculate value of tree at this node
		  valueTree[i][j]= exp(-interestRate*dt)*(q*valueTree[i+1][j+1] + (1-q)*valueTree[i+1][j])-div*dt; // discounting the continous payment.
		  //check values
	 // cout << "V_{"<<i<<"," << j << "} =";
	 // cout<<valueTree[i][j]<<endl;
	  }


  }
	
return valueTree[0][0];//-div*maturity;
}
double CDN_Leisen(double maturity, double sigma, double interestRate, double S0, double X, int N, double NB, double F, double div){
	double e;
  vector<vector<double> > tree(N+1,vector<double>(N+1));
 
  // local variables
  double u,d,q,dt;

  // initialise variables
  dt=maturity/N;
    u = exp(interestRate*dt)*(qstar(maturity,  sigma,  interestRate,  S0,  X,  N,  NB,  F)/qnormal( maturity,  sigma,  interestRate,  S0,  X,  N,  NB,  F)); //u = exp(sigma*sqrt(dt));
	d = (exp(interestRate*dt)-qnormal( maturity,  sigma,  interestRate,  S0,  X,  N,  NB,  F)*u)/(1-qnormal( maturity,  sigma,  interestRate,  S0,  X,  N,  NB,  F));//  d = exp(-sigma*sqrt(dt));
  q = (exp(interestRate*dt)-d)/(u-d);
  
  // calculate tree
  for(int i=0;i<tree.size();i++)
  {
	 
	  for(int j=0;j<tree[i].size();j++)
	  {
		//  Star[i][j]=0.;//setting value 0;
		  tree[i][j] = S0*pow(u,j)*pow(d,i-j);// value of stock at i,j
		// cout << "S_{"<<i<<"," << j << "} =";
		  //cout << tree[i][j] << endl;
	  }
  }


  // create storage for the value tree
  vector<vector<double> > valueTree(N+1,vector<double>(N+1));
  // now set option value at maturity
  // set V[N][j]= payoff condition
  for (int j=0; j<valueTree[N].size();j++){
	  //valueTree[N][j] = min(tree[N][j]/NB,F);
	   //valueTree[N][j] = max(-tree[N][j]+X,0.);
	  // valueTree[N][j] = f(tree[N][j],X,F)-(div*dt);
	   valueTree[N][j] = f(tree[N][j],X,F);
	   
	  // USE THE ADJUSTED PRICE


	  //check values
	 // cout << "V_{"<<N<<"," << j << "} =";
	  //cout<<valueTree[N][j]<<endl;
	   }

  

  // use backward induction to fill in values back to (0,0)
  //double check=0;
  for(int i=N-1;i>=0;i--){

	  for(int j=0; j<=i;j++){
		   //calculate value of tree at this node
		
		  valueTree[i][j]= exp(-interestRate*dt)*(q*valueTree[i+1][j+1] + (1-q)*valueTree[i+1][j])-div*dt; // discounting the continous payment.
		  //check values
	  // cout << "V_{"<<i<<"," << j << "} =";
	 // cout<<valueTree[i][j]<<endl;
	  }
 //double check = check + div*dt;
 
  }
//cout<<valueTree[0][0]<<"   "<<check<<"  "<<0.05*F<<"   ";	

return valueTree[0][0];//-div*maturity;
}
double CDN_hestonCHECK(double maturity, double sigma, double interestRate, double S0, double X, int N, double NB, double F, double div){
	double e;
  vector<vector<double> > tree(N+1,vector<double>(N+1));
 
  // local variables
  double u,d,q,dt;
  
   

  // initialise variables
  dt=maturity/N;
  u = exp(sigma*sqrt(dt));
  d = exp(-sigma*sqrt(dt));
  q = (exp(interestRate*dt)-d)/(u-d);
  // calculate tree
  for(int i=0;i<tree.size();i++)
  {
	 
	  for(int j=0;j<tree[i].size();j++)
	  {
		//  Star[i][j]=0.;//setting value 0;
		  tree[i][j] = S0*pow(u,j)*pow(d,i-j);// value of stock at i,j
		// cout << "S_{"<<i<<"," << j << "} =";
		  //cout << tree[i][j] << endl;
	  }
  }
	double Smax= tree[N][N];
	double Smin= tree[N][0];
  double dS=(Smax - Smin)/(N+1);  
 // cout<<"max: "<<Smax<<"min: "<<Smin;

  // create storage for the value tree
  vector<vector<double> > valueTree(N+1,vector<double>(N+1));
  // now set option value at maturity
  // set V[N][j]= payoff condition
  for (int j=0; j<valueTree[N].size();j++){
	  //valueTree[N][j] = min(tree[N][j]/NB,F);
	   //valueTree[N][j] = max(-tree[N][j]+X,0.);
	  // valueTree[N][j] = f(tree[N][j],X,F)-(div*dt);
	   valueTree[N][j] = f(tree[N][j],X,F);
	    
	   
	  // USE THE ADJUSTED PRICE
		//o.5*pow((Sj+0.5*dS-X) ,2)

	  //check values
	 // cout << "V_{"<<N<<"," << j << "} =";
	  //cout<<valueTree[N][j]<<endl;
	   }

   for  (int j=N; j>=0;j--){
	
	if(valueTree[N][j]==0. && valueTree[N][j-1]==F){
	 // cout<<"if";
			valueTree[N][j] =  F*0.5;
		}
	   
	  // USE THE ADJUSTED PRICE
		//o.5*pow((Sj+0.5*dS-X) ,2)

	  //check values
	//  cout << "NewV_{"<<N<<"," << j << "} =";
	 // cout<<valueTree[N][j]<<endl;
	   }

  // use backward induction to fill in values back to (0,0)
  //double check=0;
  for(int i=N-1;i>=0;i--){

	  for(int j=0; j<=i;j++){
		   //calculate value of tree at this node
		
		  valueTree[i][j]= exp(-interestRate*dt)*(q*valueTree[i+1][j+1] + (1-q)*valueTree[i+1][j])-div*dt; // discounting the continous payment.
		  //check values
	  // cout << "V_{"<<i<<"," << j << "} =";
	 // cout<<valueTree[i][j]<<endl;
	  }
 //double check = check + div*dt;
 
  }
//cout<<valueTree[0][0]<<"   "<<check<<"  "<<0.05*F<<"   ";	

return valueTree[0][0];//-div*maturity;
}
double CDN_hestonCHECK2(double maturity, double sigma, double interestRate, double S0, double X, int N, double NB, double F, double div){
	double e;
  vector<vector<double> > tree(N+1,vector<double>(N+1));
 
  // local variables
  double u,d,q,dt;
  
   

  // initialise variables
  dt=maturity/N;
  u = exp(sigma*sqrt(dt));
  d = exp(-sigma*sqrt(dt));
  q = (exp(interestRate*dt)-d)/(u-d);
  // calculate tree
  for(int i=0;i<tree.size();i++)
  {
	 
	  for(int j=0;j<tree[i].size();j++)
	  {
		//  Star[i][j]=0.;//setting value 0;
		  tree[i][j] = S0*pow(u,j)*pow(d,i-j);// value of stock at i,j
		// cout << "S_{"<<i<<"," << j << "} =";
		  //cout << tree[i][j] << endl;
	  }
  }
	double Smax= tree[N][N];
	double Smin= tree[N][0];
  double dS=(Smax - Smin)/(N+1);  
 // cout<<"max: "<<Smax<<"min: "<<Smin;

  // create storage for the value tree
  vector<vector<double> > valueTree(N+1,vector<double>(N+1));
  // now set option value at maturity
  // set V[N][j]= payoff condition
  for (int j=0; j<valueTree[N].size();j++){
	  //valueTree[N][j] = min(tree[N][j]/NB,F);
	   //valueTree[N][j] = max(-tree[N][j]+X,0.);
	  // valueTree[N][j] = f(tree[N][j],X,F)-(div*dt);
	   valueTree[N][j] = f(tree[N][j],X,F);
	    
	   
	  // USE THE ADJUSTED PRICE
		//o.5*pow((Sj+0.5*dS-X) ,2)

	  //check values
	 // cout << "V_{"<<N<<"," << j << "} =";
	  //cout<<valueTree[N][j]<<endl;
	   }

   for  (int j=N; j>=0;j--){
	
	if(valueTree[N][j]==0. && valueTree[N][j-1]==F){
	 // cout<<"if";
			valueTree[N][j] =  F*((X-tree[N][j-1])/(tree[N][j]-tree[N][j-1]));
		}
	   
	  // USE THE ADJUSTED PRICE
		//o.5*pow((Sj+0.5*dS-X) ,2)

	  //check values
	//  cout << "NewV_{"<<N<<"," << j << "} =";
	 // cout<<valueTree[N][j]<<endl;
	   }

  // use backward induction to fill in values back to (0,0)
  //double check=0;
  for(int i=N-1;i>=0;i--){

	  for(int j=0; j<=i;j++){
		   //calculate value of tree at this node
		
		  valueTree[i][j]= (exp(-interestRate*dt)*(q*valueTree[i+1][j+1] + (1-q)*valueTree[i+1][j]))-div*dt; // discounting the continous payment.
		  //check values
	  // cout << "V_{"<<i<<"," << j << "} =";
	 // cout<<valueTree[i][j]<<endl;
	  }
 //double check = check + div*dt;
 
  }
//cout<<valueTree[0][0]<<"   "<<check<<"  "<<0.05*F<<"   ";	

return valueTree[0][0];//-div*maturity;
}
double CDN_hestonCHECK3(double maturity, double sigma, double interestRate, double S0, double X, int N, double NB, double F, double div){
	double e;
  vector<vector<double> > tree(N+1,vector<double>(N+1));
 
  // local variables
  double u,d,q,dt;
  
   

  // initialise variables
  dt=maturity/N;
  u = exp(sigma*sqrt(dt));
  d = exp(-sigma*sqrt(dt));
  q = (exp(interestRate*dt)-d)/(u-d);
  // calculate tree
  for(int i=0;i<tree.size();i++)
  {
	 
	  for(int j=0;j<tree[i].size();j++)
	  {
		//  Star[i][j]=0.;//setting value 0;
		  tree[i][j] = S0*pow(u,j)*pow(d,i-j);// value of stock at i,j
		// cout << "S_{"<<i<<"," << j << "} =";
		  //cout << tree[i][j] << endl;
	  }
  }
	double Smax= tree[N][N];
	double Smin= tree[N][0];
  double dS=(Smax - Smin)/(N+1);  
 // cout<<"max: "<<Smax<<"min: "<<Smin;

  // create storage for the value tree
  vector<vector<double> > valueTree(N+1,vector<double>(N+1));
  // now set option value at maturity
  // set V[N][j]= payoff condition
  for (int j=0; j<valueTree[N].size();j++){
	  //valueTree[N][j] = min(tree[N][j]/NB,F);
	   //valueTree[N][j] = max(-tree[N][j]+X,0.);
	  // valueTree[N][j] = f(tree[N][j],X,F)-(div*dt);
	  // valueTree[N][j] = f(tree[N][j],X,F);
	   
	   
	   if(log(tree[N][j])-log(X)>(sigma*sqrt(dt))){   valueTree[N][j]= 0.;}
	else if(log(tree[N][j])-log(X)<(-sigma*sqrt(dt))) {  valueTree[N][j]=F;}
	else{ valueTree[N][j]= F*(0.5 +((X-tree[N][j])/(tree[N][j+1]-tree[N][j])));}

	
	   }



  // use backward induction to fill in values back to (0,0)
  //double check=0;
  for(int i=N-1;i>=0;i--){

	  for(int j=0; j<=i;j++){
		   //calculate value of tree at this node
		
		  valueTree[i][j]= (exp(-interestRate*dt)*(q*valueTree[i+1][j+1] + (1-q)*valueTree[i+1][j]))-div*dt; // discounting the continous payment.
		  //check values
	  // cout << "V_{"<<i<<"," << j << "} =";
	 // cout<<valueTree[i][j]<<endl;
	  }
 //double check = check + div*dt;
 
  }
//cout<<valueTree[0][0]<<"   "<<check<<"  "<<0.05*F<<"   ";	

return valueTree[0][0];//-div*maturity;
}
double CDN_Tian(double maturity, double sigma, double interestRate, double S0, double X, int N, double NB, double F, double div){
	double e;
  vector<vector<double> > tree(N+1,vector<double>(N+1));
  // local variables
  double u,d,q,dt;
  // initialise variables
  dt=maturity/N;
  double M = exp(interestRate*dt);
   double V= exp(sigma*sigma*dt);
   u = 0.5*M*V*((V+1)+sqrt((V*V+2*V-3)));
  d =  0.5*M*V*((V+1)-sqrt((V*V+2*V-3)));
  q = (M-d)/(u-d);
 
  // calculate tree
  for(int i=0;i<tree.size();i++)
  {
	  for(int j=0;j<tree[i].size();j++)
	  {
		  tree[i][j] = S0*pow(u,j)*pow(d,i-j);// value of stock at i,j
		//  cout << "S_{"<<i<<"," << j << "} =";
	//	  cout << tree[i][j] << endl;
	  }
  }
  // create storage for the value tree
  vector<vector<double> > valueTree(N+1,vector<double>(N+1));
  // now set option value at maturity
  // set V[N][j]= payoff condition
  for (int j=0; j<valueTree[N].size();j++){
	  //valueTree[N][j] = min(tree[N][j]/NB,F);
	   //valueTree[N][j] = max(-tree[N][j]+X,0.);
	   //valueTree[N][j] = f(tree[N][j],X,F)-(0.05*(F/maturity)*dt);//
	     //valueTree[N][j] = f(tree[N][j]-(div*dt),X,F)-div*dt;//adding q
		valueTree[N][j] = f(tree[N][j],X,F);//adding q //this is 
		 // USE THE ADJUSTED PRICE
		//o.5*pow((Sj+0.5*dS-X) ,2)

	  //check values
	  //cout << "V_{"<<N<<"," << j << "} =";
	 // cout<<valueTree[N][j]<<endl;
	   }

  // use backward induction to fill in values back to (0,0)

  for(int i=N-1;i>=0;i--){

	  for(int j=0; j<=i;j++){
		   //calculate value of tree at this node
		  valueTree[i][j]= exp(-interestRate*dt)*(q*valueTree[i+1][j+1] + (1-q)*valueTree[i+1][j])-div*dt; // discounting the continous payment.
		  //check values
	 // cout << "V_{"<<i<<"," << j << "} =";
	 // cout<<valueTree[i][j]<<endl;
	  }


  }
	
return valueTree[0][0];//-div*maturity;
}
//ConB
double ConB_CRR(double maturity, double sigma, double interestRate, double S0, double X, int N, double div2, double Z, double NB, double F, double Cp, double t1, double t2){

  vector<vector<double> > tree(N+1,vector<double>(N+1));
  // local variables
  double u,d,q,dt;
  // initialise variables



  dt=maturity/N;
  u = exp(sigma*sqrt(dt));
  d = exp(-sigma*sqrt(dt));
  q = (exp((interestRate-div2)*dt)-d)/(u-d);
  // calculate tree
 
  for(int i=0;i<tree.size();i++)
  {
	  for(int j=0;j<tree[i].size();j++)
	  {
		  tree[i][j] = S0*pow(u,j)*pow(d,i-j);// value of stock at i,j
		
		 
		 // cout << "S_{"<<i<<"," << j << "} =";
		  //cout << tree[i][j] << "  "<<dt*(i)<<endl;
	  }
  }
  // create storage for the value tree
  vector<vector<double> > valueTree(N+1,vector<double>(N+1));
  // now set option value at maturity
  // set V[N][j]= payoff condition
  for (int j=0; j<valueTree[N].size();j++){
	 // valueTree[N][j] = max(-tree[N][j]+X,0.); // price of a put
	
		 valueTree[N][j] = min(tree[N][j]/NB,max(F,Z*tree[N][j]));
	  // USE THE ADJUSTED PRICE
		//o.5*pow((Sj+0.5*dS-X) ,2)

	  //check values
	  //cout << "V_{"<<N<<"," << j << "} =";
	 // cout<<valueTree[N][j]<<endl;
	   }

  // use backward induction to fill in values back to (0,0)

  for(int i=N-1;i>=0;i--){

	  for(int j=0; j<=i;j++){
		   //calculate value of tree at this node
	

		   valueTree[i][j]= max(exp(-interestRate*dt)*(q*valueTree[i+1][j+1] + (1-q)*valueTree[i+1][j]),Z*tree[i][j]);
			   //valueTree[i][j]=max(valueTree[i][j],Z*tree[i][j]);
			
		   if(dt*i<= t2 && dt*i>=t1){
		valueTree[i][j] = min( valueTree[i][j],Cp);
		   }
		  
		  //check values
	  //cout << "V_{"<<i<<"," << j << "} =";
	  //cout<<valueTree[i][j]<<endl;
	  }
  }

   
	
return valueTree[0][0];
}
double ConB_CRRdrift(double maturity, double sigma, double interestRate, double S0, double X, int N, double div2, double Z, double NB, double F, double Cp, double t1, double t2){

  vector<vector<double> > tree(N+1,vector<double>(N+1));
  // local variables
  double u,d,q,dt;
  double drift = (log(X)- log(S0))/maturity;
  //double drift1 = (log(F/Z)- log(S0))/maturity;
  // initialise variables
  dt=maturity/N;
  u = exp(drift*dt + sigma*sqrt(dt));
  d = exp(drift*dt-sigma*sqrt(dt));
  q = (exp((interestRate-div2)*dt)-d)/(u-d);
  
  
  // calculate tree
 
  for(int i=0;i<tree.size();i++)
  {
	  for(int j=0;j<tree[i].size();j++)
	  {
		  tree[i][j] = S0*pow(u,j)*pow(d,i-j);// value of stock at i,j
		
		 
		 // cout << "S_{"<<i<<"," << j << "} =";
		  //cout << tree[i][j] << "  "<<dt*(i)<<endl;
	  }
  }
  // create storage for the value tree
  vector<vector<double> > valueTree(N+1,vector<double>(N+1));
  // now set option value at maturity
  // set V[N][j]= payoff condition
  for (int j=0; j<valueTree[N].size();j++){
	 // valueTree[N][j] = max(-tree[N][j]+X,0.); // price of a put
	
		 valueTree[N][j] = min(tree[N][j]/NB,max(F,Z*tree[N][j]));
	  // USE THE ADJUSTED PRICE
		//o.5*pow((Sj+0.5*dS-X) ,2)

	  //check values
	  //cout << "V_{"<<N<<"," << j << "} =";
	 // cout<<valueTree[N][j]<<endl;
	   }

  // use backward induction to fill in values back to (0,0)

  for(int i=N-1;i>=0;i--){

	  for(int j=0; j<=i;j++){
		   //calculate value of tree at this node
	

		   valueTree[i][j]= max(exp(-interestRate*dt)*(q*valueTree[i+1][j+1] + (1-q)*valueTree[i+1][j]),Z*tree[i][j]);
			   //valueTree[i][j]=max(valueTree[i][j],Z*tree[i][j]);
			
		   if(dt*i<= t2 && dt*i>=t1){
		valueTree[i][j] = min( valueTree[i][j],Cp);
		   }
		  
		  //check values
	  //cout << "V_{"<<i<<"," << j << "} =";
	  //cout<<valueTree[i][j]<<endl;
	  }
  }

   
	
return valueTree[0][0];
}
double ConB_Leisen(double maturity, double sigma, double interestRate, double S0, double X, int N, double div2, double Z, double NB, double F, double Cp, double t1, double t2){

  vector<vector<double> > tree(N+1,vector<double>(N+1));
  // local variables
  double u,d,q,dt;
  // initialise variables
  dt=maturity/N;
    u = exp((interestRate-div2)*dt)*(qstar(maturity,  sigma,  interestRate,  S0,  X,  N,  NB,  F)/qnormal( maturity,  sigma,  interestRate,  S0,  X,  N,  NB,  F)); //u = exp(sigma*sqrt(dt));
	d = (exp((interestRate-div2)*dt)-qnormal( maturity,  sigma,  interestRate,  S0,  X,  N,  NB,  F)*u)/(1-qnormal( maturity,  sigma,  interestRate,  S0,  X,  N,  NB,  F));//  d = exp(-sigma*sqrt(dt));
    q = (exp((interestRate-div2)*dt)-d)/(u-d);
  // calculate tree
 
  for(int i=0;i<tree.size();i++)
  {
	  for(int j=0;j<tree[i].size();j++)
	  {
		  tree[i][j] = S0*pow(u,j)*pow(d,i-j);// value of stock at i,j
		
		 
		 // cout << "S_{"<<i<<"," << j << "} =";
		  //cout << tree[i][j] << "  "<<dt*(i)<<endl;
	  }
  }
  // create storage for the value tree
  vector<vector<double> > valueTree(N+1,vector<double>(N+1));
  // now set option value at maturity
  // set V[N][j]= payoff condition
  for (int j=0; j<valueTree[N].size();j++){
	 // valueTree[N][j] = max(-tree[N][j]+X,0.); // price of a put
	
		 valueTree[N][j] = min(tree[N][j]/NB,max(F,Z*tree[N][j]));
	  // USE THE ADJUSTED PRICE
		//o.5*pow((Sj+0.5*dS-X) ,2)

	  //check values
	  //cout << "V_{"<<N<<"," << j << "} =";
	 // cout<<valueTree[N][j]<<endl;
	   }

  // use backward induction to fill in values back to (0,0)

  for(int i=N-1;i>=0;i--){

	  for(int j=0; j<=i;j++){
		   //calculate value of tree at this node
	

		   valueTree[i][j]= max(exp(-interestRate*dt)*(q*valueTree[i+1][j+1] + (1-q)*valueTree[i+1][j]),Z*tree[i][j]);
			   //valueTree[i][j]=max(valueTree[i][j],Z*tree[i][j]);
			
		   if(dt*i<= t2 && dt*i>=t1){
		valueTree[i][j] = min( valueTree[i][j],Cp);
		   }
		  
		  //check values
	  //cout << "V_{"<<i<<"," << j << "} =";
	  //cout<<valueTree[i][j]<<endl;
	  }
  }

   
	
return valueTree[0][0];
}
double ConB_Tian(double maturity, double sigma, double interestRate, double S0, double X, int N, double div2, double Z, double NB, double F, double Cp, double t1, double t2){

  vector<vector<double> > tree(N+1,vector<double>(N+1));
  // local variables
  double u,d,q,dt;
  // initialise variables
    dt=maturity/N;
double M = exp((interestRate-div2)*dt);
   double V= exp(sigma*sigma*dt);
   u = 0.5*M*V*((V+1)+sqrt((V*V+2*V-3)));
  d =  0.5*M*V*((V+1)-sqrt((V*V+2*V-3)));
  q = (exp((interestRate-div2)*dt)-d)/(u-d);


  // calculate tree
 
  for(int i=0;i<tree.size();i++)
  {
	  for(int j=0;j<tree[i].size();j++)
	  {
		  tree[i][j] = S0*pow(u,j)*pow(d,i-j);// value of stock at i,j
		
		 
		 // cout << "S_{"<<i<<"," << j << "} =";
		  //cout << tree[i][j] << "  "<<dt*(i)<<endl;
	  }
  }
  // create storage for the value tree
  vector<vector<double> > valueTree(N+1,vector<double>(N+1));
  // now set option value at maturity
  // set V[N][j]= payoff condition
  for (int j=0; j<valueTree[N].size();j++){
	 // valueTree[N][j] = max(-tree[N][j]+X,0.); // price of a put
	
		 valueTree[N][j] = min(tree[N][j]/NB,max(F,Z*tree[N][j]));
	  // USE THE ADJUSTED PRICE
		//o.5*pow((Sj+0.5*dS-X) ,2)

	  //check values
	  //cout << "V_{"<<N<<"," << j << "} =";
	 // cout<<valueTree[N][j]<<endl;
	   }

  // use backward induction to fill in values back to (0,0)

  for(int i=N-1;i>=0;i--){

	  for(int j=0; j<=i;j++){
		   //calculate value of tree at this node
	

		   valueTree[i][j]= max(exp(-interestRate*dt)*(q*valueTree[i+1][j+1] + (1-q)*valueTree[i+1][j]),Z*tree[i][j]);
			   //valueTree[i][j]=max(valueTree[i][j],Z*tree[i][j]);
			
		   if(dt*i<= t2 && dt*i>=t1){
		valueTree[i][j] = min( valueTree[i][j],Cp);
		   }
		  
		  //check values
	  //cout << "V_{"<<i<<"," << j << "} =";
	  //cout<<valueTree[i][j]<<endl;
	  }
  }
  
	
return valueTree[0][0];
}
double ConB_Heston(double maturity, double sigma, double interestRate, double S0, double X, int N, double div2, double Z, double NB, double F, double Cp, double t1, double t2){

  vector<vector<double> > tree(N+1,vector<double>(N+1));
  // local variables
  double u,d,q,dt;
  // initialise variables



  dt=maturity/N;
  u = exp(sigma*sqrt(dt));
  d = exp(-sigma*sqrt(dt));
  q = (exp((interestRate-div2)*dt)-d)/(u-d);
  // calculate tree
 
  for(int i=0;i<tree.size();i++)
  {
	  for(int j=0;j<tree[i].size();j++)
	  {
		  tree[i][j] = S0*pow(u,j)*pow(d,i-j);// value of stock at i,j
		
		 
		 // cout << "S_{"<<i<<"," << j << "} =";
		  //cout << tree[i][j] << "  "<<dt*(i)<<endl;
	  }
  }
  // create storage for the value tree
  vector<vector<double> > valueTree(N+1,vector<double>(N+1));
  // now set option value at maturity
  // set V[N][j]= payoff condition
  for (int j=0; j<valueTree[N].size();j++){
	 // valueTree[N][j] = max(-tree[N][j]+X,0.); // price of a put
	
		// valueTree[N][j] = min(tree[N][j]/NB,max(F,Z*tree[N][j]));
	 // valueTree[N][j] = f_Put_hestonConvertible(tree[N][j],NB,dt,sigma,F,Z,Cp);
		 
	 if(log(tree[N][j])-log(F/Z)<(sigma*sqrt(dt))&& log(tree[N][j])-log(F/Z)>(-sigma*sqrt(dt))){  
		 valueTree[N][j] = 	((((F/Z)*(sigma*sqrt(dt)-log((tree[N][j]/(F/Z)))))/(2*sigma*sqrt(dt)))-(S/(2*sigma*sqrt(dt)))*(exp(-sigma*sqrt(dt))+((F/Z)/tree[N][j])));}
	 else if(log(tree[N][j])-log(X)<(sigma*sqrt(dt))&& log(tree[N][j])-log(X)>(-sigma*sqrt(dt))){
		 valueTree[N][j] = (((X*(sigma*sqrt(dt)-log((tree[N][j]/X))))/(2*sigma*sqrt(dt)))-(S/(2*sigma*sqrt(dt)))*(exp(-sigma*sqrt(dt))+(X/tree[N][j])));}
	 else{
		 valueTree[N][j] = min(tree[N][j]/NB,max(F,Z*tree[N][j]));
	 }
	 
	
	   }

  // use backward induction to fill in values back to (0,0)

  for(int i=N-1;i>=0;i--){

	  for(int j=0; j<=i;j++){
		   //calculate value of tree at this node
	

		   valueTree[i][j]= max(exp(-interestRate*dt)*(q*valueTree[i+1][j+1] + (1-q)*valueTree[i+1][j]),Z*tree[i][j]);
			   //valueTree[i][j]=max(valueTree[i][j],Z*tree[i][j]);
			
		   if(dt*i<= t2 && dt*i>=t1){
		valueTree[i][j] = min( valueTree[i][j],Cp);
		   }
		  
		  //check values
	  //cout << "V_{"<<i<<"," << j << "} =";
	  //cout<<valueTree[i][j]<<endl;
	  }
  }

   
	
return valueTree[0][0];
}

};

double section1(){
cout. precision(15);


	// parameters
	BlackScholes o1;
	DBonds o2;
	
	double NB= 1092, F= 100, NS= 4575, S= 9.9;
	double maturity=2.4,sigma=0.37,interestRate=0.057,S0= NS*S + NB*F, X= NB*F; // if S0 = X the shape is symetric!! if we mode we dont't get symetric
	/*double maturity=.5,sigma=0.2,interestRate=0.1,S0=42, X=40., div2=0.; //To test the results
	int N=5000;
	cout<<o2.amePut(maturity,sigma,interestRate,S0,X,N,div2);*/
		

///1.2
	
 std::ofstream output;
	 output.open("C:/Users/jm_zarate_c/Documents/Visual Studio 2010/Projects/ScientificComputing/11_BinomialProject/output.csv");
	 output. precision(15);
	 //C:\Users\jm_zarate_c\Documents\Visual Studio 2010\Projects\ScientificComputing\9_class
	if(!output.is_open())
 {
     std::cout << " File not opened \n"; // To validate if the file was opened.
    
	 // stop the program here
    throw;
  }
//Writting headings 
	output<<"N"<<","<<"Analytical"<<","
					<<"DB_CRR"<<","
					<<"DB_CRRdrift"<<","
					<<"DB_Leisen"<<","
					<<"DB_Heston"<<","
					<<"DB_Tian"<<"\n";



	double olda=0., diff=0.;
	double oldclock=0., clocka=0., absdiffclock;
	double  a=0., b=0., c=0, d=0.,e=0., f=0.;
		
	//cout << "Time Elapsed(Clock ticks: " << clock() << ") Seconds:  " << double(clock()/CLOCKS_PER_SEC) << endl;
	

	
	double time2=0, oldtime=0;
	
	for( int N=100; N<1000;N=N++){ // number of steps
 //   cout<<"the value with: "<<N<< "is "<< binomialCall( maturity,  sigma,  interestRate,  S0, X,  N)<<"\n";

		clock_t startTime = clock();

	//oldtime = time2;
		//a =F*exp(-maturity*interestRate)- o1.put(S0,maturity,0,interestRate,0,sigma,X)/NB;
		
		//b = o2.DB_CRR( maturity,  sigma,  interestRate,  S0, X,  N, NB, F);		
		c = o2.DB_CRRdrift( maturity,  sigma,  interestRate,  S0, X,  N, NB, F);		
		//d = o2.DB_Leisen( maturity,  sigma,  interestRate,  S0, X,  N, NB, F);		
		//e = F*exp(-maturity*interestRate)-o2.DB_Put_Heston( maturity,  sigma,  interestRate,  S0, X,  N, NB, F)/NB;		
		//f = o2.DB_Tian( maturity,  sigma,  interestRate,  S0, X,  N, NB, F);		
		//clocka= double(clock()/CLOCKS_PER_SEC);
		
		//cout<<o1.put(S0,maturity,0,interestRate,0,sigma,X)<<"\n";
		output<<N<<","<< a<<"," 
					<<b<<","
					<<c<<","
					<<d<<","
					<<e<<","
					<<f<<",";
		
		
	
		float secsElapsed = (float)(clock() - startTime)/CLOCKS_PER_SEC;
		output << " Seconds:,  " << secsElapsed << endl;
		
		
		cout<<N<<" \n";
	//cout << " Acum Seconds:,  " << time2 <<"   individual:, "<<time2-oldtime<< endl;
	//cout << "Time Elapsed(Clock ticks: " << clock() << ") Seconds:  " << double(clock()/CLOCKS_PER_SEC) << endl;
//cout<<N<<"\n";
	 //cout<<binomialPutNormal( maturity,  sigma,  interestRate,  S0, X,  N);


	}

   std::cout << " File write successful 1 \n";
	// CLOSE file
  output.close();
  

  ///////////////////////////////////////////////////////////////////////

  /*
  
 std::ofstream output2;
	 output2.open("C:/Users/jm_zarate_c/Documents/Visual Studio 2010/Projects/ScientificComputing/11_BinomialProject/section1.1.csv");
	 
	 //C:\Users\jm_zarate_c\Documents\Visual Studio 2010\Projects\ScientificComputing\9_class
	if(!output2.is_open())
 {
     std::cout << " File not opened \n"; // To validate if the file was opened.
    
	 // stop the program here
    throw;
  }
//Writting headings 
	output2<<"B(V t)"<<","<<"Expected Payoff"<<","<<"Analytical"<<"\n";


	for( double S01=50000; S01<250000;S01=S01+5000){ // number of steps
 int N=3.;
		maturity=0.00001;
		double a = o2.DB_CRR( maturity,  sigma,  interestRate,  S01, X,  N, NB, F);
		
		
		output2<<S01<<","<< a<<","
					<<F*exp(-maturity*interestRate)- o1.put(S01,maturity,0,interestRate,0,sigma,X)/NB<<"\n";


	}


 	  std::cout << " File write successful \n";
	// CLOSE file
  output2.close();

  

  */


 return 0;




}

double section2(){
	

cout. precision(15);

BlackScholes o1;
	DBonds o2;

	double NB= 1092, F= 100, NS= 4575, S= 9.9;
	double maturity=2.4,sigma=0.37,interestRate=0.057,S0= NS*S + NB*F, X= NB*F; // if S0 = X the shape is symetric!! if we mode we dont't get symetric
	double div= 0.05*(F/maturity);
//	cout<<o2.f_Put_heston(100,120,maturity/100,sigma)<<"  ";

	//double maturity=.5,sigma=0.2,interestRate=0.1,S0=42, X=42., F=1, div=0., NB=0. ; //To test the results
		//int N=4;
		//cout<<o2.binomialCDNLeisen( maturity,  sigma,  interestRate,  S0, X,  N, NB, F,div);



 std::ofstream output2;
	 output2.open("C:/Users/jm_zarate_c/Documents/Visual Studio 2010/Projects/ScientificComputing/11_BinomialProject/output6.csv");
	 output2. precision(15);
	 //C:\Users\jm_zarate_c\Documents\Visual Studio 2010\Projects\ScientificComputing\9_class
	if(!output2.is_open())
 {
     std::cout << " File not opened \n"; // To validate if the file was opened.
    
	 // stop the program here
    throw;
  }
//Writting headings 

	output2<<"N"<<","<<"Analitical"<<","
			<<"CDN_CRR"<<","	
			<<"CDN_drift"<<","	
			<<"CDN_Leisen"<<","
			<<"CDN_Heston"<<","
			<<"CDN_Heston2"<<","
			<<"CDN_Tian"<<"\n";

double a=0., b=0.,c=0.,d=0., e=0., f=0., g=0.;
		a=F*o1.binPut(S0,maturity,0,interestRate,0,sigma,X)-div*maturity*disF(maturity,interestRate);//*exp(-interestRate*maturity);
	for( int N=9500; N<10001;N=N++){ // number of steps
 //   cout<<"the value with: "<<N<< "is "<< binomialCall( maturity,  sigma,  interestRate,  S0, X,  N)<<"\n";
			cout<<N<<"\n";	
		b = o2.CDN_CRR( maturity,  sigma,  interestRate,  S0, X,  N, NB, F,div);
		c = o2.CDN_CRRdrift( maturity,  sigma,  interestRate,  S0, X,  N, NB, F,div);
		d=o2.CDN_Leisen( maturity,  sigma,  interestRate,  S0, X,  N, NB, F,div);
		//e= o2.CDN_hestonCHECK( maturity,  sigma,  interestRate,  S0, X,  N, NB, F,div);
		f= o2.CDN_hestonCHECK3( maturity,  sigma,  interestRate,  S0, X,  N, NB, F,div);
		g=o2.CDN_Tian( maturity,  sigma,  interestRate,  S0, X,  N, NB, F,div);

		
		output2<<N<<","<< a<<","
					<<b<<","
					<<c<<","
					<<d<<","
					<<e<<","		
					<<f<<","
					<<g<<"\n";	

 //cout<<N<<"\n";
	 //cout<<binomialPutNormal( maturity,  sigma,  interestRate,  S0, X,  N);


	}
	

	  std::cout << " File write successful 2 \n";
	// CLOSE file
  output2.close();
//  */


	return 0;
}
double section3(){
	
BlackScholes o1;
	DBonds o2;

	double NB= 1092, F= 100, NS= 4575, S= 9.9;
	double maturity=2.4,sigma=0.37,interestRate=0.057,S0= NS*S + NB*F, X= NB*F; // if S0 = X the shape is symetric!! if we mode we dont't get symetric
	double div= 0.05*(F/maturity);


	//Parameters for convertible bond
	double Sc=21;// Converstion share price
	double Qration = F/Sc; // Convertion ratio
	double div2=0.056; // continuous dividend proportional to the rm value V 
	double deltaNS= Qration*NB;
	double Z= Qration/(deltaNS+NS); //dilution factor
	double Cp= 107; //If the bond may be bought back at the price CP over some time period
	double t1=0.55;
	double t2=1.;

//	cout<<o2.f_Put_heston(100,120,maturity/100,sigma)<<"  ";

	//double maturity=.5,sigma=0.2,interestRate=0.1,S0=42, X=42., F=1, div=0., NB=0. ; //To test the results
		//int N=4;
		//cout<<o2.binomialCDNLeisen( maturity,  sigma,  interestRate,  S0, X,  N, NB, F,div);
	cout<<S0<<"  ";
	

 std::ofstream output3;
	 output3.open("C:/Users/jm_zarate_c/Documents/Visual Studio 2010/Projects/ScientificComputing/11_BinomialProject/output5.csv");
	  output3. precision(15);
	 //C:\Users\jm_zarate_c\Documents\Visual Studio 2010\Projects\ScientificComputing\9_class
	if(!output3.is_open())
 {
     std::cout << " File not opened \n"; // To validate if the file was opened.
    
	 // stop the program here
    throw;
  }
//Writting headings 

	output3<<"N"<<","
			<<"CB_CRR"<<","	
			<<"CB_CRRdrift"<<","	
			<<"CB_Leisen"<<","	
			<<"CB_Tian"<<","
			<<"Heston"<<"\n";

double a=0., b=0.,c=0.,d=0., e=0., f=0.;

	for( int N=100; N<1001;N=N++){ // number of steps
 //   cout<<"the value with: "<<N<< "is "<< binomialCall( maturity,  sigma,  interestRate,  S0, X,  N)<<"\n";
		
	
		//b = o2.ConB_CRR(maturity, sigma, interestRate,S0,X,N,div2,Z,NB,F,Cp,t1,t2);
		c =  o2.ConB_CRRdrift(maturity, sigma, interestRate,S0,X,N,div2,Z,NB,F,Cp,t1,t2);
		//d= o2.ConB_Leisen(maturity, sigma, interestRate,S0,X,N,div2,Z,NB,F,Cp,t1,t2);
		//e= o2.ConB_Tian(maturity, sigma, interestRate,S0,X,N,div2,Z,NB,F,Cp,t1,t2);
		//f=o2.ConB_Heston(maturity, sigma, interestRate,S0,X,N,div2,Z,NB,F,Cp,t1,t2);
		cout<<N<<"\n";
		output3<<N<<","
					<<b<<","
					<<c<<","
					<<d<<","
					<<e<<","
					<<f<<"\n";	

 //cout<<N<<"\n";
	 //cout<<binomialPutNormal( maturity,  sigma,  interestRate,  S0, X,  N);


	}
	

	  std::cout << " File write successful 3 \n";
	// CLOSE file
  output3.close();
//  */


	



}

int main(){


//cout. precision(15);

	section3();

	system("pause");


 return 0;	


}
