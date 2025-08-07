#include "TMatrixD.h"
#include "TString.h"

#include <iostream>

void matrixWqq() {
  /*
  double h3all = 1;
  double h3gencsall = 0.289841;
  double h3genudall = 0.377007;
  double h3genxall = 0.333153;
  double h3tagcsall = 0.117846;
  double h3tagudall = 0.180028;
  double h3tagxall = 0.702127;
  double h3gencstagcs = 0.0820791;
  double h3gencstagud = 0.0243256;
  double h3gencstagx = 0.183436;
  double h3genudtagcs = 0.00609975;
  double h3genudtagud = 0.091527;
  double h3genudtagx = 0.27938;
  double h3genxtagcs = 0.0296669;
  double h3genxtagud = 0.0641751;
  double h3genxtagx = 0.23931;
  double h3all_data = 1;
  double h3tagcsall_data = 0.0952018;
  double h3tagudall_data = 0.148134;
  double h3tagxall_data = 0.756665;
*/

double h3all = 1;
double h3gencsall = 0.290962;
double h3genudall = 0.377095;
double h3genxall = 0.331943;
double h3tagcsall = 0.118699;
double h3tagudall = 0.179498;
double h3tagxall = 0.701803;
double h3gencstagcs = 0.0827145;
double h3gencstagud = 0.0243024;
double h3gencstagx = 0.183945;
double h3genudtagcs = 0.00612232;
double h3genudtagud = 0.0914332;
double h3genudtagx = 0.279539;
double h3genxtagcs = 0.0298623;
double h3genxtagud = 0.0637625;
double h3genxtagx = 0.238318;
double h3all_data = 1;
double h3tagcsall_data = 0.0952018;
double h3tagudall_data = 0.148134;
double h3tagxall_data = 0.756665;
double h3all_scaled = 1;
  // Matrix of (true x gen) efficiencies [c, x, u] x [cs, xx, ud]
  double A[3][3] =
    {{h3gencstagcs/h3gencsall, h3genxtagcs/h3genxall, h3genudtagcs/h3genudall},
     { h3gencstagx/h3gencsall,  h3genxtagx/h3genxall,  h3genudtagx/h3genudall},
     {h3gencstagud/h3gencsall, h3genxtagud/h3genxall, h3genudtagud/h3genudall}
    };

  // Tagging efficiency scale factors
  //double kcc(0.9), kuu(0.9), kxx(1.135); // diagonal
  //double kcu(0.9), kuc(0.9); // corners
  double kcc(0.84), kuu(0.85), kxx(1.09); // diagonal
  double kcu(0.84), kuc(0.84); // corners
  // Calculating neighbouring element to keep probability sum unity
  double kxc = (A[1][0]+A[0][0]*(1-kcc)+A[2][0]*(1-kcu))/A[1][0];
  double kxu = (A[1][2]+A[2][2]*(1-kuu)+A[0][2]*(1-kuc))/A[1][2];
  double fcx = (A[0][1])/(A[0][1]+A[2][1]);
  double kcx = (A[0][1]+A[1][1]*fcx*(1-kxx))/A[0][1];
  double kux = (A[2][1]+A[1][1]*(1-fcx)*(1-kxx))/A[2][1];
  // Final matrix
  double K[3][3] =
    {{kcc, kcx, kcu},
     {kxc, kxx, kxu},
     {kuc, kux, kuu}
    };
  
  // Vector of gen fractions
  double g[3] =
    {h3gencsall/h3all, h3genxall/h3all, h3genudall/h3all};

  // Vector of tag fractions
  double m[3] =
    {h3tagcsall/h3all, h3tagxall/h3all, h3tagudall/h3all};

  // Vector of data fractions
  double d[3] =
    {h3tagcsall_data/h3all_data, h3tagxall_data/h3all_data, h3tagudall_data/h3all_data};

  TMatrixD MA(3,3), MB(3,3), MC(3,3), MD(3,3);
  TMatrixD vg(3,1), vgm(3,1), vtm(3,1), vm(3,1);
  TMatrixD vgd(3,1), vtd(3,1), vd(3,1);
  TMatrixD vgr(3,1), vtr(3,1);
  for (int i = 0; i != 3; ++i) {
    for (int j = 0; j != 3; ++j) {
      MA[i][j] = A[i][j];
      MB[i][j] = A[i][j]; // A inverted
      MC[i][j] = A[i][j]*K[i][j];
      MD[i][j] = A[i][j]*K[i][j]; // C inverted
    } // for j
  } // for i
  for (int i = 0; i != 3; ++i) {
    vgm[i][0] = g[i];
    vg[i][0]  = g[i];
    //vtm[i][0] = m[i];
    vm[i][0]  = m[i];
    vtd[i][0] = d[i];
    vd[i][0]  = d[i];
    //vtr[i][0] = d[i];
  }
  
  // Check MC tag fractions with forward multiplication
  vtm.Mult(MA,vg);

  // Solve data gen fractions with matrix inversion
  MB.Invert();
  vgd.Mult(MB,vd);
  MD.Invert();
  vgr.Mult(MD,vd);

  // Check data tag fractions with forward multiplication
  vtd.Mult(MA,vgd);
  vtr.Mult(MC,vgr);

  // Copy matrices back to arrays
  double tm[3], td[3], tr[3], gm[3], gd[3], gr[3];
  for (int i = 0; i != 3; ++i) {
    gm[i] = g[i];
    gd[i] = vgd[i][0];
    gr[i] = vgr[i][0];
    tm[i] = vtm[i][0];
    td[i] = vtd[i][0];
    tr[i] = vtr[i][0];
  }
  double B[3][3], C[3][3];
    for (int i = 0; i != 3; ++i) {
    for (int j = 0; j != 3; ++j) {
      B[i][j] = MB[i][j];
      C[i][j] = MC[i][j];
    } // for j
  } // for i


  // Print fractions
  cout << "Forward pass for MC gen fractions (ref. data):" << endl;
  cout << Form("[ d_c]  [ m_c]   [  c ]   [c_cs c_xx c_ud]   [ cs ]\n");
  cout << Form("[ d_x], [ m_x] = [  x ] = [x_cs x_xx x_ud] * [ xx ]\n");
  cout << Form("[ d_u]  [ m_u]   [  u ]   [u_cs u_xx u_ud]   [ ud ]\n");
  cout << endl;
  cout << Form("[%1.2f]  [%1.2f]   [%1.2f]   [%1.2f %1.2f %1.2f]   [%1.2f]\n",
	       d[0],     m[0],     tm[0],    A[0][0],A[0][1],A[0][2],gm[0]);
  cout << Form("[%1.2f], [%1.2f] = [%1.2f] = [%1.2f %1.2f %1.2f] * [%1.2f]\n",
	       d[1],     m[1],     tm[1],    A[1][0],A[1][1],A[1][2],gm[1]);
  cout << Form("[%1.2f]  [%1.2f]   [%1.2f]   [%1.2f %1.2f %1.2f]   [%1.2f]\n",
	       d[2],     m[2],     tm[2],    A[2][0],A[2][1],A[2][2],gm[2]);
  cout << endl;
  cout << "Forward pass for data gen fractions solved by inversion (ref. MC):" << endl;
  cout << Form("[ m_c]  [ d_c]   [  c ]   [c_cs c_xx c_ud]   [ cs ]\n");
  cout << Form("[ m_x], [ d_x] = [  x ] = [x_cs x_xx x_ud] * [ xx ]\n");
  cout << Form("[ m_u]  [ d_u]   [  u ]   [u_cs u_xx u_ud]   [ ud ]\n");
  cout << endl;
  cout << Form("[%1.2f]  [%1.2f]   [%1.2f]   [%1.2f %1.2f %1.2f]   [%1.2f]\n",
	       m[0],     d[0],     td[0],    A[0][0],A[0][1],A[0][2],gd[0]);
  cout << Form("[%1.2f], [%1.2f] = [%1.2f] = [%1.2f %1.2f %1.2f] * [%1.2f]\n",
	       m[1],     d[1],     td[1],    A[1][0],A[1][1],A[1][2],gd[1]);
  cout << Form("[%1.2f]  [%1.2f]   [%1.2f]   [%1.2f %1.2f %1.2f]   [%1.2f]\n",
	       m[2],     d[2],     td[2],    A[2][0],A[2][1],A[2][2],gd[2]);
  cout << endl;
  cout << "Scaled efficiency matrix for data starting from diagonal + corners:" << endl;
  cout << Form("[c_cs c_xx c_ud]   [ kcc  kcx  kcu]\n");
  cout << Form("[x_cs x_xx x_ud] = [ kxc  kxx  kxu] (x) C_MC\n");
  cout << Form("[u_cs u_xx u_ud]   [ kuc  kux  kuu]\n");
  cout << endl;
  cout << Form("                   [%1.2f %1.2f %1.2f]  data   [%1.2f]        [%1.2f]\n",
	       K[0][0],K[0][1],K[0][2], gr[0]/gm[0], d[0]/m[0]);
  cout << Form("                   [%1.2f %1.2f %1.2f], ---- = [%1.2f] (gen), [%1.2f] (tag)\n",
	       K[1][0],K[1][1],K[1][2], gr[1]/gm[1], d[1]/m[1]);
  cout << Form("                   [%1.2f %1.2f %1.2f]   mc    [%1.2f]        [%1.2f]\n",
	       K[2][0],K[2][1],K[2][2], gr[2]/gm[2], d[2]/m[2]);
  cout << endl;
  cout << "Forward pass for data solved by inversion of scaled efficiency matrix:" << endl;
  cout << Form("[%1.2f]  [%1.2f]   [%1.2f]   [%1.2f %1.2f %1.2f]   [%1.2f]\n",
	       m[0],     d[0],     tr[0],    C[0][0],C[0][1],C[0][2],gr[0]);
  cout << Form("[%1.2f], [%1.2f] = [%1.2f] = [%1.2f %1.2f %1.2f] * [%1.2f]\n",
	       m[1],     d[1],     tr[1],    C[1][0],C[1][1],C[1][2],gr[1]);
  cout << Form("[%1.2f]  [%1.2f]   [%1.2f]   [%1.2f %1.2f %1.2f]   [%1.2f]\n",
	       m[2],     d[2],     tr[2],    C[2][0],C[2][1],C[2][2],gr[2]);


  // Interpret x as a linear kombination of cs-like (bX) and ud like (xG)
  // x = f*c + (1-f)*u <=> x-u = (c-u)*f <=> = f = (x-u)/(c-u)
  double fm[3];
  fm[0] = (A[0][1]-A[0][2])/(A[0][0]-A[0][1]);
  fm[1] = (A[1][1]-A[1][2])/(A[1][0]-A[1][1]);
  fm[2] = (A[2][1]-A[2][2])/(A[2][0]-A[2][1]);
  double fd[3];
  fd[0] = (C[0][1]-C[0][2])/(C[0][0]-C[0][1]);
  fd[1] = (C[1][1]-C[1][2])/(C[1][0]-C[1][1]);
  fd[2] = (C[2][1]-C[2][2])/(C[2][0]-C[2][1]);
  cout << endl;
  cout << "Gen-X decomposed to cs-like (bX) vs ud-like (Xg) using tagging efficiency:" << endl;
  cout << Form("           [%1.2f]       [%1.2f]\n",fm[0],fd[0]);
  cout << Form("f(cs/bX) = [%1.2f] (MC), [%1.2f] (data)\n",fm[1],fm[1]);
  cout << Form("           [%1.2f]       [%1.2f]\n",fm[2],fm[2]);
}
