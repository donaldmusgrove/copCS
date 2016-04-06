#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


// Calculate log-likelihood
// ======================================================================
double llk(mat Q, mat T){
  double trQT, logL, valQ, signQ;
  log_det(valQ,signQ,Q);
   
  trQT = trace(Q*T);  
  logL = 0.5 * ( valQ - trQT );
 
  return( logL );
}

// ======================================================================
// IPSP2 helper functions 
// ----------------------------------------------------------------------
// - Extract unique vertices from a clique subgraph
// ======================================================================
uvec Wvert(List Wi){
  int nWi = Wi.size();
  uvec V;
   
  for(int j=0; j<nWi; j++){
    uvec Vcat = Wi[j];
    V         = join_cols(V, Vcat);  
  }

  return(unique(V));
}


// Extract unique vertices from a (complement) clique subgraph
// ======================================================================
uvec WvertC(List Wi, int C){
  int nWi = Wi.size();
  uvec V;
   
  for(int j=0; j<nWi; j++ ){
    if(j != C){
      uvec Vcat = Wi[j];
      V         = join_cols(V, Vcat);  
    }
  }

  return(unique(V));
}


// Set difference
// - Return vector Ea with elements of Ec removed
// ======================================================================
uvec setdifference(uvec Ea, uvec Ec){
  int na = Ea.size();
   
  for(int j=na; j --> 0;){
    bool matchind = any(Ec == Ea(j));
    if(matchind == TRUE){
      Ea.shed_row(j);
    }
  }

  return(Ea);
}

vec setdifferencev(vec Ea, vec Ec){
  int na = Ea.size();
   
  for(int j=na; j --> 0;){
    bool matchind = any(Ec == Ea(j));
    if(matchind == TRUE){
      Ea.shed_row(j);
    }
  }

  return(Ea);
}

// Sequence from 0:(n-1)
// ======================================================================
vec sequencen(int n){
  vec seqn(n);

  for(int i=0; i<n; i++){
    seqn(i) = i;
  }

  return(seqn);
}


// Generate list of vertices in each clique partition
// - Function depends on clique partition skeleton (W)
// ======================================================================
List clique_part(List W){
  int nW = W.size();
  List Wi;
  List U(nW);
  int nWi;
  vec Ui;
   
  for(int i=0; i<nW; i++){
    Wi  = W[i];
    nWi = Wi.size();
      
    for(int j=0; j<nWi; j++){
      vec Wij = Wi[j];
      Ui      = join_cols(Ui, Wij); 
    }
    U[i] = unique(Ui);   
    Ui.clear();
  }   

  return (U);
}

// Generate complement list of vertices in each clique partition
// - Function depends on clique partition skeleton (W)
// ======================================================================
List clique_comppart(List U, int n){
  int nU = U.size();
  vec V  = sequencen(n);
  List Ubar(nU);
  vec Ubari;
   
  for(int i=0; i<nU; i++ ){
    vec Ui    = U[i];
    vec Ubari = setdifferencev(V, Ui);
    Ubar[i]   = Ubari;
  }
  
  return (Ubar);
}




// IPSP2 Algorithm of Xu et al (2014)
// ======================================================================
// [[Rcpp::export]]
SEXP IPS2(double rho, int n, mat T, List W, int maxit, double eps){
  List U    = clique_part(W);
  List Ubar = clique_comppart(U, n);

  int itcount = 0;
  double prevlogL, logL;
  bool converged;
  vec logLout(maxit);
  uvec V;
  int nW = W.size();
  List Wi;
  int nWi, nEc;   
  mat QM, QMbar, QMMbar, Psi, LambdaTemp, Sc, Scinv, Qca, Qaa, Qtemp, QT;
  vec trQT;
  
  mat Lambda(n,n);
  mat Q(n,n); 
    Q.eye();
    
  while (itcount<maxit){
    for(int i=0; i<nW; i++){
      Wi        = W[i];
      V         = Wvert(Wi);
      nWi       = Wi.size();

      uvec Uvert     = U[i];
      uvec Ubarvert  = Ubar[i];
     
      QM           = Q( Uvert, Uvert );
      QMbar        = Q( Ubarvert, Ubarvert );
      QMMbar       = Q( Uvert, Ubarvert );     
      Psi          = QMMbar * inv_sympd(QMbar) * trans(QMMbar);
      LambdaTemp   = QM - Psi;
      Lambda(V,V)  = LambdaTemp;
     
      for(int j=0; j<nWi; j++){
        uvec Ec  = Wi[j];
        uvec Ea  = WvertC(Wi, j);
        Ea       = setdifference(Ea, Ec);
        
        nEc      = Ec.size();
        arma::mat Idn(nEc,nEc); 
        Idn.eye();
        arma::mat Onesn(nEc,nEc);
        Onesn.ones();
         
        Sc     = (1-rho) * Idn + rho * Onesn;
        Scinv  = inv(Sc);
            
        Qca   = Lambda(Ec, Ea);
        Qaa   = Lambda(Ea, Ea);
       
        Lambda(Ec,Ec) = Scinv + Qca * inv(Qaa) * trans(Qca);
      }
     
      Qtemp  = Lambda(V,V) + Psi;
      Q(V,V) = Qtemp; 
    }
        
    logL             = llk(Q,T);
    logLout(itcount) = logL;  
    itcount++;

    if(itcount > 1){
      if( logL - prevlogL < eps ){
        converged = TRUE;
        break;
      } else {
        if(itcount == maxit){
          converged = FALSE;
          break;
        }
      }         
    }
    
    prevlogL = logL;
  }
   
   logLout = logLout.subvec(0,itcount-1);

   return Rcpp::List::create(Rcpp::Named("Q")         = Q,
                             Rcpp::Named("count")     = itcount,
                             Rcpp::Named("converged") = converged,
                             Rcpp::Named("logL")      = logLout);
}