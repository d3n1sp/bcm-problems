#include "stdafx.h"
#include "cgrid_el.h"


///////////////////////////////////////////////////////////////////////////////
//...auxilliary function for grid nodes of the straight line (parametric line);
void CGrid_el::grid_box1D(int N1, int shift, int id_link, double fy, double fz)
{
  int i, l; double d; 
  delete_struct(geom);
  if (N1 > 0 && bufferization(N = N1+1, NULL_STATE) 
				 && (geom = (int *)new_struct((N1+4)*sizeof(int))) != NULL) {
/////////////////////////////////////////////////////////////////
//...filling points -- vertical of homotating boxes and geometry;
      for (Y[0]   = fy, Z[0] = fz, d = 1./N1, geom[l = 0] = 1, i = 0; i < N1; i++) {
           X[i+1] = X[i]+d;
           Y[i+1] = Y[i];
           Z[i+1] = Z[i];
      }
		for (geom[++l] = GL_LINE_STRIP, geom[++l] = N1+1, i = 0; i <= N1; i++) geom[++l] = i+shift;
		hit[0] = -1; hit[N1] = -2;	if (id_link) hit[0] = hit[N1] = -1;
      buf_count = 0; return;
  }
  zero_grid();
}

///////////////////////////////////////////////////////////////////////////////
//...auxilliary function for grid nodes of rectangular cell (parametric quads);
void CGrid_el::grid_box2D(int N1, int N2, int shift, int id_link, double fz)
{
  int    i, j, l, m;
  double d, g; 
  delete_struct(geom);
  if (N1 > 0 && N2 > 0 && bufferization(N = (N1+1)*(N2+1), NULL_STATE) 
				 && (geom = (int *)new_struct((1+2*N2*(N1+2))*sizeof(int))) != NULL) {
/////////////////////////////////////////////////////////////////
//...filling points -- vertical of homotating boxes and geometry;
      for (Z[0]   = fz, d = 1./N1, geom[l = 0] = N2, m = N2*(N1+1), i = 0; i < N1; i++) {
           X[i+1] = X[i]+d;
           Y[i+1] = Y[i];
           Z[i+1] = Z[i];
           hit[i] = hit[i+1] = -1; hit[m+i] = hit[m+i+1] = -3;
           if (id_link) hit[i] = hit[i+1] = hit[m+i] = hit[m+i+1] = -1;
      }
      for (g = 1./N2,  j = 0; j < geom[0]; j++) {
           geom[++l] = GL_QUAD_STRIP;
           geom[++l] = 2*(N1+1);
           for (m = j*(N1+1), i = 0; i <= N1; i++) {
                geom[++l]  = m+i+N1+1+shift;
                geom[++l]  = m+i+shift;
                X[m+i+N1+1] = X[m+i];
                Y[m+i+N1+1] = Y[m+i]+g;
                Z[m+i+N1+1] = Z[m+i];
           }
           hit[m+N1] = -2; hit[m+N1+1] = -4;
           if (id_link) hit[m+N1] = hit[m+N1+1] = -1;
      }
      buf_count = 0; return;
  }
  zero_grid();
}

/////////////////////////////////////////////////////////////////////////////////
//...auxilliary function for grid nodes of rectangular box (parametric unit box);
void CGrid_el::grid_box3D(int N1, int N2, int N3, int shift, int id_link)
{
  int    i, j, k, m = 0, n, NN = (N1+1)*(N2+1);
  double d, g, h;
  delete_struct(geom);
  if (N1 > 0 && N2 > 0 && bufferization(N = (N1+1)*(N2+1)*(N3+1), NULL_STATE)
				 && (geom = (int *)new_struct((1+2*N2*(N1+2))*sizeof(int))) != NULL) {
/////////////////////////////////////////////////////////////////
//...filling points -- vertical of homotating boxes and geometry;
      for (Z[0]   = 0., d = 1./N1, i = 0; i < N1; i++) {
           X[i+1] = X[i]+d;
           Y[i+1] = Y[i];
           Z[i+1] = Z[i];
           hit[i] = hit[i+1] = -1; hit[m+i] = hit[m+i+1] = -3;
           if (id_link) hit[i] = hit[i+1] = hit[m+i] = hit[m+i+1] = -1;
      }
      for (g = 1./N2,  j = 0; j < N2; j++) {
           for (m = j*(N1+1), i = 0; i <= N1; i++) {//...где geom ???
                X[m+i+N1+1] = X[m+i];
                Y[m+i+N1+1] = Y[m+i]+g;
                Z[m+i+N1+1] = Z[m+i];
           }
           hit[m+N1] = -2; hit[m+N1+1] = -4;
           if (id_link) hit[m+N1] = hit[m+N1+1] = -1;
      }
      for (h = 1./N3, k = 0; k < N3; k++) {
           for (n = k*NN, i = 0; i < NN; i++) {
                X[n+i+NN] = X[n+i];
                Y[n+i+NN] = Y[n+i];
                Z[n+i+NN] = Z[n+i]+h;
           }
           hit[m+N1] = -2; hit[m+N1+1] = -4;
           if (id_link) hit[m+N1] = hit[m+N1+1] = -1;
      }
      buf_count = 0; return;
  }
  zero_grid();
}

/////////////////////////////////////////////////
//...auxilliary function for triangle grid nodes;
void CGrid_el::grid_tria_3D(int & N_elem, double * P, int shift)
{
  double d, px1, px2, py1, py2, pz1, pz2;
  int    i, j, k, l, m = (int)(sqrt((double)N_elem)-1.);
  if (N_elem > 0 && P && bufferization(N = m*(m+5)/2+3, NULL_STATE) &&
	  (geom   = (int *)new_struct((1+(m+1)*(m+5))*sizeof(int))) != NULL) {
      N_elem = (m+1)*(m+1);

/////////////////////////////////////////////////////////
//...shifting vector homothetic triangles in subdivision;
      px1 = (P[6]-P[0])*(d = 1./(1.+m)); px2 = (P[3]-P[0])*d;
      py1 = (P[7]-P[1])* d;              py2 = (P[4]-P[1])*d;
      pz1 = (P[8]-P[2])* d;              pz2 = (P[5]-P[2])*d;

//////////////////////////////////////////////////////
//...filling points -- apices of homothetic triangles;
      for (X[0]  = P[0], Y[0]  = P[1],  Z[0]  = P[2],
           nX[0] = P[9], nY[0] = P[10], nZ[0] = P[11],
           hit[0]    = -1, geom[l = 0] =m+1, i = 1; i <= geom[0]; i++) {
           geom[++l] = GL_TRIANGLE_STRIP;
           geom[++l] = 2*i+1;
           for (j = i*(i+1)/2, k = 0; k < i; k++) {
                X[j+k]    = X[j+k-i]+px1;
                Y[j+k]    = Y[j+k-i]+py1;
                Z[j+k]    = Z[j+k-i]+pz1;
                nX[j+k]   = P[9];
                nY[j+k]   = P[10];
                nZ[j+k]   = P[11];
                hit[j+k]  = -1;
                geom[++l] = j+k+shift;
                geom[++l] = j+k-i+shift;
           }
           X[j+i]    = X[j-1]+px2;
           Y[j+i]    = Y[j-1]+py2;
           Z[j+i]    = Z[j-1]+pz2;
           nX[j+i]   = P[9];
           nY[j+i]   = P[10];
           nZ[j+i]   = P[11];
           hit[j+i]  = -1;
           geom[++l] = j+i+shift;
      }
      buf_count = 0; return;
  }
  zero_grid();
}

////////////////////////////////////////////////////////
//...auxilliary function for grid nodes surface of hexa;
void CGrid_el::grid_hexa_refine(int N_elem, double * P)
{
/////////////
//...face Z;
	grid_quad_3D_refine(N_elem, N_elem, P+12);	swap(P[3],  P[9]);  swap(P[4],  P[10]);  swap(P[5],  P[11]);
/////////////
//...face Z';
	grid_quad_3D_refine(N_elem, N_elem, P);		swap(P[15], P[21]); swap(P[16], P[22]);  swap(P[17], P[23]);
																swap(P[3],  P[12]); swap(P[4],  P[13]);  swap(P[5],  P[14]);
																swap(P[6],  P[21]); swap(P[7],  P[22]);  swap(P[8],  P[23]);
/////////////
//...face Y';
	grid_quad_3D_refine(N_elem, N_elem, P+12);	swap(P[3],  P[9]);  swap(P[4],  P[10]);  swap(P[5],  P[11]);
/////////////
//...face Y;
	grid_quad_3D_refine(N_elem, N_elem, P);		swap(P[3],  P[9]);  swap(P[4],  P[10]);  swap(P[5],  P[11]);
																swap(P[6],  P[15]); swap(P[7],  P[16]);  swap(P[8],  P[17]);
																swap(P[9],  P[12]); swap(P[10], P[13]);  swap(P[11], P[14]);
/////////////
//...face X';
	grid_quad_3D_refine(N_elem, N_elem, P);		swap(P[15], P[21]); swap(P[16], P[22]);  swap(P[17], P[23]);
/////////////
//...face X;
	grid_quad_3D_refine(N_elem, N_elem, P+12);	swap(P[6],  P[21]); swap(P[7],  P[22]);  swap(P[8],  P[23]);
																swap(P[6],  P[15]); swap(P[7],  P[16]);  swap(P[8],  P[17]);
																swap(P[3],  P[12]); swap(P[4],  P[13]);  swap(P[5],  P[14]);
}

/////////////////////////////////////////////////////////
//...auxilliary function for grid nodes surface of tetra;
void CGrid_el::grid_tetra_refine(int N_elem, double * P)
{
///////////////
//...lateral 1;
													swap(P[0], P[3]); swap(P[1], P[4]);  swap(P[2], P[5]); 
	grid_tria_3D_refine(N_elem, P);		swap(P[0], P[3]); swap(P[1], P[4]);  swap(P[2], P[5]); 
///////////////
//...lateral 2;
	grid_tria_3D_refine(N_elem, P+3);	swap(P[6], P[9]); swap(P[7], P[10]); swap(P[8], P[11]); 
///////////////
//...lateral 4;
	grid_tria_3D_refine(N_elem, P);		swap(P[6], P[9]); swap(P[7], P[10]); swap(P[8], P[11]); 
///////////////
//...lateral 3;
													swap(P[3], P[9]); swap(P[4], P[10]); swap(P[5], P[11]); 
	grid_tria_3D_refine(N_elem, P);		swap(P[3], P[9]); swap(P[4], P[10]); swap(P[5], P[11]); 
}

/////////////////////////////////////////////////////
//...auxilliary function for basic sphere grid nodes;
int  CGrid_el::grid_sph_2D(int N_elem, CMap * mp)
{
  if (N_elem > 0 && size_of_map(mp) > 6) {
      double r0 = fabs(mp[7]), v, h;
      int    i;
      if (bufferization(N = N_elem, NULL_STATE)) {
          for (h = M_PI*2.0/N_elem, v = 0., i = 0; i < N; i++, v += h) {
               X[i]   = r0*cos(v)+mp[1];
               Y[i]   = r0*sin(v)+mp[2];
               Z[i]   = mp[3];
               nX[i]  = -sin(v);
               nY[i]  =  cos(v);
               nZ[i]  = 0.;
               hit[i] = -1;
          }
          buf_count = 0; return(1);
      }
  }
  zero_grid(); return(0);
}

/////////////////////////////////////////////////////
//...auxilliary function for basic sphere grid nodes;
int  CGrid_el::grid_sph_3D(int N_elem, CMap * mp)
{
  if (N_elem > 0 && mp && ID_MAP(2, SPHERE_GENUS) == mp[0]) {
      double r0 = fabs(mp[7]), f0, f, HY, v, u, h;
      int    NX, NY, k, l, i = 0;

///////////////////////////////
//...defining number of points;
      f  = sqrt(4.*M_PI/(double)N_elem);
      NY = f > 0. ? (int)ceil(M_PI/f)+1 : 4;
      f0 = f > EE ?  M_PI/f : 4.;
      HY = M_PI/NY;

      for (v  = 0., N = k = 0; k <= NY; k++, v += HY)
      if ((NX = max(1, (int)(f0*sin(v)))) > 1) N += 2*NX;
      else                                     N += 2*NX-1;

//////////////////////////////
//...filling values of points;
      if (bufferization(N, NULL_STATE)) {
          for (v = 0, k = 0; k <= NY; k++, v += HY) {
               h = M_PI/(NX = max(1, (int)(f0*sin(v)))); u = 0.;
               if (NX > 1) {
                   X[i]     = r0*cos(u)*sin(v)+mp[1];
                   Y[i]     = r0*sin(u)*sin(v)+mp[2];
                   Z[i]     = r0*cos(v)+mp[3];
                   nX[i]    = cos(u)*sin(v);
                   nY[i]    = sin(u)*sin(v);
                   nZ[i]    = cos(v);
                   hit[i++] = -1;
               }
               for (u += h, l = 1-NX; l < NX; l++, u += h) {
                   X[i]     = r0*cos(u)*sin(v)+mp[1];
                   Y[i]     = r0*sin(u)*sin(v)+mp[2];
                   Z[i]     = r0*cos(v)+mp[3];
                   nX[i]    = cos(u)*sin(v);
                   nY[i]    = sin(u)*sin(v);
                   nZ[i]    = cos(v);
                   hit[i++] = -1;
               }
          }
          buf_count = 0; return(1);
      }
  }
  zero_grid(); return(0);
}

///////////////////////////////////////////////////////
//...auxilliary function for basic cylinder grid nodes;
int  CGrid_el::grid_cyl_3D(int N_elem, CMap * mp, int id_flag)
{
  if (N_elem > 0 && mp && ID_MAP(2, CYL_GENUS) == mp[0]) {
      double r0 = fabs(mp[7]), f0, f, HY, HZ, v, u, h;
      int    NX, NY, NZ, k, l, m, i = 0;

///////////////////////////////
//...defining number of points;
      f  = sqrt(M_PI/(double)N_elem);
      NY = f > 0. ? (int)ceil(M_PI_2/f)+1 : 4;
      f0 = f > EE ?  M_PI/f : 4.;
      HY = M_PI_2/NY;

      for (v  = 0., k = 0; k <= NY; k++, v += HY)
      if ((NX = max(1, (int)(f0*sin(v)))) > 1) N +=  4*NX;
      else                                     N += (2*NX-1)*2;

      NZ = r0 > EE ? (int)(fabs(mp[9]-mp[8])/(M_PI*r0)*NX) : 0;
      HZ = (mp[8]-mp[9])/(NZ+1);
      N += 2*NX*NZ;

//////////////////////////////
//...filling values of points;
      if (bufferization(N, NULL_STATE)) {

///////////////////////////////////////////
//...filling points for ending semispheres;
          for (v = 0, k = 0; k < NY; k++, v += HY) {
               h = M_PI/(NX = max(1, (int)(f0*sin(v)))); u = 0.;
               if (NX > 1) {
                   X[i]     = r0*cos(u)*sin(v)+mp[1];
                   Y[i]     = r0*sin(u)*sin(v)+mp[2];
                   Z[i]     = mp[3]+mp[9]+(id_flag ? 0. : r0*cos(v));
                   nX[i]    = id_flag ? 0. : cos(u)*sin(v);
                   nY[i]    = id_flag ? 0. : sin(u)*sin(v);
                   nZ[i]    = id_flag ? 1. : cos(v);
                   hit[i++] = -1;

                   X[i]     = X[i-1];
                   Y[i]     = Y[i-1];
                   Z[i]     = mp[3]*2.+mp[8]+mp[9]-Z[i-1];
                   nX[i]    = id_flag ?  0. : cos(u)*sin(v);
                   nY[i]    = id_flag ?  0. : sin(u)*sin(v);
                   nZ[i]    = id_flag ? -1. : -cos(v);
                   hit[i++] = -1;
               }
               for (u += h, l = 1-NX; l < NX; l++, u += h) {
                   X[i]     = r0*cos(u)*sin(v)+mp[1];
                   Y[i]     = r0*sin(u)*sin(v)+mp[2];
                   Z[i]     = mp[3]+mp[9]+(id_flag ? 0. : r0*cos(v));
                   nX[i]    = id_flag ? 0. : cos(u)*sin(v);
                   nY[i]    = id_flag ? 0. : sin(u)*sin(v);
                   nZ[i]    = id_flag ? 1. : cos(v);
                   hit[i++] = -1;

                   X[i]     = X[i-1];
                   Y[i]     = Y[i-1];
                   Z[i]     = mp[3]*2.+mp[8]+mp[9]-Z[i-1];
                   nX[i]    = id_flag ?  0. : cos(u)*sin(v);
                   nY[i]    = id_flag ?  0. : sin(u)*sin(v);
                   nZ[i]    = id_flag ? -1. : -cos(v);
                   hit[i++] = -1;
               }
          }

////////////////////////////////////////////////////
//...filling points for lateral surface of cylinder;
          h = M_PI/(NX = max(1, (int)(f0*sin(v)))); u = 0.;
          if (NX > 1) {
              X[i]     = r0*cos(u)*sin(v)+mp[1];
              Y[i]     = r0*sin(u)*sin(v)+mp[2];
              Z[i]     = mp[3]+mp[9]+(id_flag ? 0. : r0*cos(v));
              nX[i]    = cos(u);
              nY[i]    = sin(u);
              hit[i++] = -1;
              for (m = 0; m <= NZ; m++) {
                   X[i]     = X[i-1];
                   Y[i]     = Y[i-1];
                   Z[i]     = Z[i-1]+HZ;
                   nX[i]    = cos(u);
                   nY[i]    = sin(u);
                   hit[i++] = -1;
              }
          }
          for (u += h, l = 1-NX; l < NX; l++, u += h) {
              X[i]     = r0*cos(u)*sin(v)+mp[1];
              Y[i]     = r0*sin(u)*sin(v)+mp[2];
              Z[i]     = mp[3]+mp[9]+(id_flag ? 0. : r0*cos(v));
              nX[i]    = cos(u);
              nY[i]    = sin(u);
              hit[i++] = -1;
              for (m = 0; m <= NZ; m++) {
                   X[i]     = X[i-1];
                   Y[i]     = Y[i-1];
                   Z[i]     = Z[i-1]+HZ;
                   nX[i]    = cos(u);
                   nY[i]    = sin(u);
                   hit[i++] = -1;
              }
          }
          buf_count = 0; return(1);
      }
  }
  zero_grid(); return(0);
}

////////////////////////////////////////////////////////////////////////////////
//...auxilliary function for basic box or sphere grid nodes (in some direction);
int  CGrid_el::grid_box_3D(int N_elem, CMap * mp, int dir)
{
  int k;
  if (N_elem > 0 && mp && (mp[k = size_of_map(mp)] == (CMap)  BASIS_BOX ||
                           mp[k]                   == (CMap)SPECIAL_BOX ||
                           mp[k]                   == (CMap)NULL_CELL)) {
      double f, X1, X2, Y1, Y2, Z1, ZZ, dX, dY, dZ, P[3], pp[3], HX, HY, HZ;
      int    i, j, m, n, l, NX, NY, NZ;

//////////////////////////////////////
//...defining parameters of direction;
      n  = (i = (m = dir)/2+1)%3+1; 
      j = n%3+1;

//////////////////////////////////
//...operating with spherical box;
      if (mp[k] == (CMap)NULL_CELL) {
         if  (grid_sph_3D(2*N_elem, mp))
         for (k = 0; k < N; k++) {
             P[0] = X[k]-mp[1];
             P[1] = Y[k]-mp[2];
             P[2] = Z[k]-mp[3];
             if (P[i-1]*((dir % 2)*2.-1.) < EE_ker) hit[k] = 0;
         }
         else {
            zero_grid(); return(0);
         }
         return(1);
      }

////////////////////////////////////////
//...defining coordinates of box center;
      X1 = mp[k+n+3]-mp[k+n];
      X2 = mp[k+n+3]+mp[k+n];
      Y1 = mp[k+j+3]-mp[k+j];
      Y2 = mp[k+j+3]+mp[k+j]; HX = mp[k+n]; HY = mp[k+j]; HZ = mp[k+i];
      Z1 = mp[k+i+3];
      ZZ = mp[k+i+3]+mp[k+i]*(2.*(m%2)-1.); if (ZZ < Z1) Z1 = ZZ;

///////////////////////////////
//...defining number of points;
      f  = sqrt(4.*(HX*HY+2.*HY*HZ+2.*HX*HZ)/N_elem);
      NX = 1+(int)ceil((dX = 2.*HX)/f);     dX /= NX;
      NY = 1+(int)ceil((dY = 2.*HY)/f);     dY /= NY;
      NZ = 1+(int)ceil((dZ =    HZ)/f);     dZ /= NZ;

//////////////////////////////
//...filling values of points;
      if (bufferization(N = (NX+1)*(NY+1)+2*(NY+1)*(NZ+1)+2*(NX+1)*(NZ+1), NULL_STATE)) {

//////////////////////////////////
//...filling upper part of surface;
          pp[0] = pp[1] = pp[2] = 0.;
          for (P[i-1] = ZZ, pp[i-1] = (dir % 2)*2.-1., k = l = 0; l <= NY; l++)
          for (P[j-1] = Y1+l*dY,                           m = 0; m <= NX; m++) {
               P[n-1] = X1+m*dX;
               X[k]     = P[0];
               Y[k]     = P[1];
               Z[k]     = P[2];
               nX[k]    = pp[0];
               nY[k]    = pp[1];
               nZ[k]    = pp[2];
               hit[k++] = -1;
          }

/////////////////////////////
//...filling lateral surface;
          pp[0] = pp[1] = pp[2] = 0.;
          for (P[j-1] = Y1, pp[j-1] = -1., l = 0; l <= NZ; l++)
          for (P[i-1] = Z1+l*dZ,           m = 0; m <= NX; m++) {
               P[n-1] = X1+m*dX;
               X[k]     = P[0];
               Y[k]     = P[1];
               Z[k]     = P[2];
               nX[k]    = pp[0];
               nY[k]    = pp[1];
               nZ[k]    = pp[2];
               hit[k++] = -1;
          }
          pp[0] = pp[1] = pp[2] = 0.;
          for (P[j-1] = Y2, pp[j-1] = 1., l = 0; l <= NZ; l++)
          for (P[i-1] = Z1+l*dZ,          m = 0; m <= NX; m++) {
               P[n-1] = X1+m*dX;
               X[k]     = P[0];
               Y[k]     = P[1];
               Z[k]     = P[2];
               nX[k]    = pp[0];
               nY[k]    = pp[1];
               nZ[k]    = pp[2];
               hit[k++] = -1;
          }
          pp[0] = pp[1] = pp[2] = 0.;
          for (P[n-1] = X1, pp[n-1] = -1., l = 0; l <= NZ; l++)
          for (P[i-1] = Z1+l*dZ,           m = 0; m <= NY; m++) {
               P[j-1] = Y1+m*dY;
               X[k]     = P[0];
               Y[k]     = P[1];
               Z[k]     = P[2];
               nX[k]    = pp[0];
               nY[k]    = pp[1];
               nZ[k]    = pp[2];
               hit[k++] = -1;
          }
          pp[0] = pp[1] = pp[2] = 0.;
          for (P[n-1] = X2, pp[n-1] = 1., l = 0; l <= NZ; l++)
          for (P[i-1] = Z1+l*dZ,          m = 0; m <= NY; m++) {
               P[j-1] = Y1+m*dY;
               X[k]     = P[0];
               Y[k]     = P[1];
               Z[k]     = P[2];
               nX[k]    = pp[0];
               nY[k]    = pp[1];
               nZ[k]    = pp[2];
               hit[k++] = -1;
          }
          buf_count = 0; return(1);
      }
  }
  zero_grid(); return(0);
}

//////////////////////////////////////////////////////////////////
//...auxilliary function for grid nodes of line between two boxes;
int CGrid_el::grid_box_box_2D(int N_elem, CMap * mp, CMap * mp_ext)
{
  int k, l, i, m, n, NX;
  if (N_elem > 0 && mp     && mp    [k = size_of_map(mp)]     == (CMap)SPECIAL_SHEET &&
                    mp_ext && mp_ext[l = size_of_map(mp_ext)] == (CMap)SPECIAL_SHEET) {
      double f, X1, X2, X3, X4, YY, P[2], dX;
      int    dir = -1;

//////////////////////////////////
//...direction to the nearest box;
      if ((m = dir) < 0) {
          if (fabs(fabs(f = mp[k+3]-mp_ext[l+3])-(mp[k+1]+mp_ext[l+1])) < EE_ker) m = f > 0. ? 0 : 1; else
          if (fabs(fabs(f = mp[k+4]-mp_ext[l+4])-(mp[k+2]+mp_ext[l+2])) < EE_ker) m = f > 0. ? 2 : 3; dir = m;
      }
      if (dir < 0.) goto err;

///////////////////////////////////////////////////////////////
//...definition parameters of interlace line between two boxes;
      n  = (i = m/2+1)%2+1;
      X1 = mp[k+n+2]-mp[k+n]; X3 = mp_ext[l+n+2]-mp_ext[l+n]; X1 = max(X1, X3);
      X2 = mp[k+n+2]+mp[k+n]; X4 = mp_ext[l+n+2]+mp_ext[l+n]; X2 = min(X2, X4);
      YY = mp[k+i+2]+mp[k+i]*(2.*(m%2)-1.);

///////////////////////////////
//...defining number of points;
      f  = 4.*(mp[k+1]+mp[k+2])/N_elem;
      NX = 1+(int)ceil((dX = X2-X1)/f); dX /= NX;

//////////////////////////////
//...filling values of points;
      if (bufferization(N = NX+1), NULL_STATE) {
          for (P[i-1] = YY,  k = l = 0; l <= NX; l++) {
               P[n-1] = X1+l*dX;
               X[k]     = P[0];
               Y[k]     = P[1];
               nX[k]    = dir/2 == 0 ? (2.*(dir%2)-1.) : 0.;
               nY[k]    = dir/2 == 1 ? (2.*(dir%2)-1.) : 0.;
               hit[k++] = -1;
          }
          buf_count = 0; return(1);
      }
  }
err:
  zero_grid(); return(0);
}

////////////////////////////////////////////////////////////////////////////////////
//...auxilliary function for grid nodes of interlace region between two boxes in 3D;
void CGrid_el::grid_box_box_3D(int N_elem, CMap * mp, CMap * mp_ext, int id_flag, int & dir)
{
  int k, l, i, j, m, n, NX, NY, NZ;
  if (N_elem > 0 && mp     && (mp    [k = size_of_map(mp)]     == (CMap)SPECIAL_BOX || mp    [k] == (CMap)BASIS_BOX) &&
                    mp_ext && (mp_ext[l = size_of_map(mp_ext)] == (CMap)SPECIAL_BOX || mp_ext[l] == (CMap)BASIS_BOX)) {
      double f, X1, X2, X3, X4, Y1, Y2, Y3, Y4, Z1, Z2, ZZ, P[3], dX, dY, dZ;

//////////////////////////////////////////////////////////////////
//...®ЇаҐ¤Ґ«пҐ¬ ­ Їа ў«Ґ­ЁҐ а бЇ®«®¦Ґ­Ёп б®бҐ¤­ҐЈ® Їап¬®гЈ®«м­ЁЄ ;
      if (fabs(fabs(f = mp[k+4]-mp_ext[l+4])-(mp[k+1]+mp_ext[l+1])) < EE_ker) m = f > 0. ? 0 : 1; else
      if (fabs(fabs(f = mp[k+5]-mp_ext[l+5])-(mp[k+2]+mp_ext[l+2])) < EE_ker) m = f > 0. ? 2 : 3; else
      if (fabs(fabs(f = mp[k+6]-mp_ext[l+6])-(mp[k+3]+mp_ext[l+3])) < EE_ker) m = f > 0. ? 4 : 5; dir = m;

///////////////////////////////////////////////////////////////////////////////////////
//...®ЇаҐ¤Ґ«пҐ¬ Ї а ¬Ґвал ­ Їа ў«Ґ­Ёп Ё Є®®а¤Ё­ вл жҐ­ва  Ї®¤­Ё¬ Ґ¬®Ј® Ї а ««Ґ«ҐЇЁЇҐ¤ ;
      n  = (i = m/2+1)%3+1; j = n%3+1;
      X1 = mp[k+n+3]-mp[k+n]; X3 = mp_ext[l+n+3]-mp_ext[l+n]; X1 = max(X1, X3);
      X2 = mp[k+n+3]+mp[k+n]; X4 = mp_ext[l+n+3]+mp_ext[l+n]; X2 = min(X2, X4);
      Y1 = mp[k+j+3]-mp[k+j]; Y3 = mp_ext[l+j+3]-mp_ext[l+j]; Y1 = max(Y1, Y3);
      Y2 = mp[k+j+3]+mp[k+j]; Y4 = mp_ext[l+j+3]+mp_ext[l+j]; Y2 = min(Y2, Y4);
      Z1 = mp[k+i+3]+mp[k+i]*(2.*(m%2)-1.); ZZ = Z1;
      Z2 = mp_ext[l+i+3]; if (Z2 < Z1) swap(Z1, Z2);

//////////////////////////////////////
//...®ЇаҐ¤Ґ«пҐ¬ зЁб«® в®зҐЄ а §ЎЁҐ­Ёп;
      f  = id_flag ? sqrt(fabs((X2-X1)*(Y2-Y1))/N_elem) :
                     sqrt(8.*(mp[k+1]*mp[k+2]+mp[k+2]*mp[k+3]+mp[k+1]*mp[k+3])/N_elem);
      NX = 1+(int)ceil((dX = X2-X1)/f); dX /= NX;
      NY = 1+(int)ceil((dY = Y2-Y1)/f); dY /= NY;
      NZ = 1+(int)ceil((dZ = Z2-Z1)/f); dZ /= NZ;

//////////////////////////////
//...filling values of points;
      if (bufferization(N = (NX+1)*(NY+1)+(id_flag ? 0 : 2*(NX+1)*(NZ+1)+2*(NY+1)*(NZ+1)), NULL_STATE)) {
          for (P[i-1] = id_flag ? ZZ : mp_ext[l+i+3], k = l = 0; l <= NY; l++)
          for (P[j-1] = Y1+l*dY,                          m = 0; m <= NX; m++) {
               P[n-1] = X1+m*dX;
               X[k]     = P[0];
               Y[k]     = P[1];
               Z[k]     = P[2];
               hit[k++] = -1;
          }

//////////////////////////////////////////////
//...§ Ї®«­пҐ¬ Ў®Є®ўЁ­ЄЁ Ї®¤­пв®© Ї®ўҐае­®бвЁ;
          if (! id_flag) {
              for (P[j-1] = Y1,        l = 0; l <= NZ; l++)
              for (P[i-1] = Z1+l*dZ,   m = 0; m <= NX; m++) {
                   P[n-1] = X1+m*dX;
                   X[k]     = P[0];
                   Y[k]     = P[1];
                   Z[k]     = P[2];
                   hit[k++] = -1;
              }
              for (P[j-1] = Y2,        l = 0; l <= NZ; l++)
              for (P[i-1] = Z1+l*dZ,   m = 0; m <= NX; m++) {
                   P[n-1] = X1+m*dX;
                   X[k]     = P[0];
                   Y[k]     = P[1];
                   Z[k]     = P[2];
                   hit[k++] = -1;
              }
              for (P[n-1] = X1,        l = 0; l <= NZ; l++)
              for (P[i-1] = Z1+l*dZ,   m = 0; m <= NY; m++) {
                   P[j-1] = Y1+m*dY;
                   X[k]     = P[0];
                   Y[k]     = P[1];
                   Z[k]     = P[2];
                   hit[k++] = -1;
              }
              for (P[n-1] = X2,        l = 0; l <= NZ; l++)
              for (P[i-1] = Z1+l*dZ,   m = 0; m <= NY; m++) {
                   P[j-1] = Y1+m*dY;
                   X[k]     = P[0];
                   Y[k]     = P[1];
                   Z[k]     = P[2];
                   hit[k++] = -1;
              }
          }
          buf_count = 0; return;
      }
  }
  zero_grid();
}

////////////////////////////////////////////////////////////////////////////////
//...auxilliary function for grid nodes of interlace region between two circles;
void CGrid_el::grid_sph_sph_2D(int N_elem, CMap * mp, CMap * mp_ext)
{
  if (N_elem > 0 && mp     && ID_MAP(2, SPHERE_GENUS) == mp[0]     &&
                    mp_ext && ID_MAP(2, SPHERE_GENUS) == mp_ext[0]) {
      double r0 = fabs(mp[7]), r1 = fabs(mp_ext[7]), R0 = abs_point(mp+1, mp_ext+1),
              x = R0 > EE_dop ? .5*(R0+(r0*r0-r1*r1)/R0) : r0, f0, f, HX, HY,
              P[3], P0[3] = {mp[1], mp[2], 0.}, f1;
      int     NX, i;
      if ((f0 = r0*r0-x*x) > EE) {
///////////////////////////////
//...defining number of points;
           f0 = sqrt(f0);
           f1 = arg2(comp(x, f0));
           f  = (HX = 2.*f1)/N_elem;
           NX = f > 0. ? (int)ceil(HX/f)+1 : 4; HX /= NX-1;

///////////////////////////////////////////////////////
//...filling values of points (between seating points);
           if (bufferization(N = NX-1), NULL_STATE) {
               for (i = 0; i < N; i++) {
                   nX[i]  = cos(f = HX*(i+.5)-f1);
                   nY[i]  = sin(f);
                   X[i]   = r0*nX[i];
                   Y[i]   = r0*nY[i];
                   hit[i] = -1;
               }

///////////////////////////////////////////////////
//...shift and rotate points into interlace region;
               f = arg0(comp(HX = mp_ext[1]-mp[1], HY = mp_ext[2]-mp[2]));
               for (HX = cos(f), HY = sin(f), i = 0; i < N; i++) {
                    P[0] = X[i];
                    P[1] = Y[i]; P[2] = 0.; point_iso(P, P0, HX, HY, 1., 0., 1., 0.);
                    X[i] = P[0];
                    Y[i] = P[1];
               }
               buf_count = 0; return;
           }
      }
  }
  zero_grid();
}

//////////////////////////////////////////////////////////////////////////////////////
//...auxilliary function for grid nodes of interlace region between two spheres in 3D;
void CGrid_el::grid_sph_sph_3D(int N_elem, CMap * mp, CMap * mp_ext)
{
  if (N_elem > 0 && mp     && ID_MAP(2, SPHERE_GENUS) == mp[0]     &&
                    mp_ext && ID_MAP(2, SPHERE_GENUS) == mp_ext[0]) {
      double r0 = fabs(mp[7]), r1 = fabs(mp_ext[7]), R0 = abs_point(mp+1, mp_ext+1),
              x = R0 > EE_dop ? .5*(R0+(r0*r0-r1*r1)/R0) : r0, f0, f, d, HX, HY, HZ,
              P[3], P0[3] = {mp[1], mp[2], mp[3]}, v, u, h;
      int    NX, NY, k, l, i = 0;
      if ((f0 = r0*r0-x*x) > EE) {
///////////////////////////////
//...defining number of points;
           f0 = sqrt(f0);
           HY = arg2(comp(x, f0));
           f  = sqrt(M_PI*(1.-cos(HY))/N_elem);
           NY = f > 0. ? (int)ceil(HY/f)+1 : 4; HY /= NY;
           f0 = f > EE ?  M_PI/f : 4.;

           for (v  = 0., k = 0; k <= NY; k++, v += HY)
           if ((NX = max(1, (int)(f0*sin(v)))) > 1) N += 2*NX;
           else                                     N += 2*NX-1;

//////////////////////////////
//...filling values of points;
           if (bufferization(N, NULL_STATE)) {
              for (v = 0, k = 0; k <= NY; k++, v += HY) {
                   h = M_PI/(NX = max(1, (int)(f0*sin(v)))); u = 0.;
                   if (NX > 1) {
                       X[i]     = r0*cos(u)*sin(v);
                       Y[i]     = r0*sin(u)*sin(v);
                       Z[i]     = r0*cos(v);
                       hit[i++] = -1;
                   }
                   for (u += h, l = 1-NX; l < NX; l++, u += h) {
                       X[i]     = r0*cos(u)*sin(v);
                       Y[i]     = r0*sin(u)*sin(v);
                       Z[i]     = r0*cos(v);
                       hit[i++] = -1;
                   }
              }

///////////////////////////////////////////////////
//...shift and rotate points into interlace region;
              f  = arg0(comp(HX = mp_ext[1]-mp[1], HY = mp_ext[2]-mp[2]));
              d  = arg0(comp(HZ = mp_ext[3]-mp[3], sqrt(HX*HX+HY*HY)));
              for (HX = cos(f), HY = sin(f),
                   HZ = cos(d), f0 = sin(d), i = 0; i < N; i++) {
                   P[0] = X[i];
                   P[1] = Y[i];
                   P[2] = Z[i]; point_iso(P, P0, HX, HY, HZ, f0, 1., 0.);
                   X[i] = P[0];
                   Y[i] = P[1];
                   Z[i] = P[2];
              }
              buf_count = 0; return;
           }
      }
  }
  zero_grid();
}

//////////////////////////////////////////////////////////////////////////////////////
//...auxilliary function for grid nodes of interlace region between two spheres in 3D;
int CGrid_el::grid_sph_sph_3D_diameter(int N_elem, CMap * mp, CMap * mp_ext, double theta0)
{
  if (N_elem > 0 && mp     && ID_MAP(2, SPHERE_GENUS) == mp[0]     &&
                    mp_ext && ID_MAP(2, SPHERE_GENUS) == mp_ext[0]) {

      double r0 = fabs(mp[7]), r1 = fabs(mp_ext[7]), R0 = abs_point(mp+1, mp_ext+1),
              x = R0 > EE_dop ? .5*(R0+(r0*r0-r1*r1)/R0) : r0, f0 = r0*r0-x*x, f, d, HX, HY, HZ,
              P[3], P0[3] = {mp[1], mp[2], mp[3]}, v, u, h, ff, theta = f0 > EE ? arg2(comp(x, f0 = sqrt(f0))): M_PI;
      int    NX, NY, k, l, i = 0;

      if ((theta *= 180./M_PI) > theta0) {
///////////////////////////////
//...defining number of points;
           HY = arg2(comp(R0, r1)); r0 = sqrt(R0*R0+r1*r1);
           f  = sqrt(M_PI*(1.-cos(HY))/N_elem);
           NY = f > 0. ? (int)ceil(HY/f)+1 : 4; HY /= NY;
           f0 = f > EE ?  M_PI/f : 4.;

           for (v  = 0., k = 0; k <= NY; k++, v += HY)
           if ((NX = max(1, (int)(f0*sin(v)))) > 1) N += 2*NX;
           else                                     N += 2*NX-1;

//////////////////////////////
//...filling values of points;
           if (bufferization(N, NULL_STATE)) {
              for (v = 0, k = 0; k <= NY; k++, v += HY) {
                   h = M_PI/(NX = max(1, (int)(f0*sin(v)))); u = 0.;
                   ff = R0/(r0*cos(v));
                   if (NX > 1) {
                       X[i]     = r0*cos(u)*sin(v)*ff;
                       Y[i]     = r0*sin(u)*sin(v)*ff;
                       Z[i]     = R0;
                       hit[i++] = -1;
                   }
                   for (u += h, l = 1-NX; l < NX; l++, u += h) {
                       X[i]     = r0*cos(u)*sin(v)*ff;
                       Y[i]     = r0*sin(u)*sin(v)*ff;
                       Z[i]     = R0;
                       hit[i++] = -1;
                   }
              }

///////////////////////////////////////////////////
//...shift and rotate points into interlace region;
              f  = arg0(comp(HX = mp_ext[1]-mp[1], HY = mp_ext[2]-mp[2]));
              d  = arg0(comp(HZ = mp_ext[3]-mp[3], sqrt(HX*HX+HY*HY)));
              for (HX = cos(f), HY = sin(f),
                   HZ = cos(d), f0 = sin(d), i = 0; i < N; i++) {
                   P[0] = X[i];
                   P[1] = Y[i];
                   P[2] = Z[i]; point_iso(P, P0, HX, HY, HZ, f0, 1., 0.);
                   X[i] = P[0];
                   Y[i] = P[1];
                   Z[i] = P[2];
              }
              buf_count = 0; return(1);
           }
      }
  }
  return(0);
//zero_grid();
}

///////////////////////////////////////////////////////////////////////////////////////////////
//...вспомогательная функция для построения 2D сетки в четырехугольнике (возможно вырожденном);
int CGrid_el::grid_quad(complex * w, double *& X, double *& Y, double h1, double h2,
                           int *& geom, int & N, int shift, int Max_N1, int Max_N2)
{
	double f1 = abs(w[2]-w[0]), f2 = abs(w[3]-w[1]), f3, f4;
	int    N1, N2, m1, m2, j, l; if (f1 < EE_ker) f1 = 0.; if (f2 < EE_ker) f2 = 0.;
	complex z, z1, z2;

////////////////////////////////////////////////////////////
//...определяем количество узлов сетки для четыpехугольника;
	f3 = 0.5*(f1+f2);
	f4 = 0.5*(abs(w[1]-w[0])+abs(w[3]-w[2]));
	N1 = h1 > 0. ? 2+(int)ceil(f3/h1) : 1+abs(Max_N1);
	N2 = h2 > 0. ? 2+(int)ceil(f4/h2) : 1+abs(Max_N2);
	if (f1 == 0. && f2 == 0.) {
		X = Y = NULL; geom = NULL;
		N = 0; return(1);
	}

///////////////////////////////////////////////////////////////////////
//...распределяем массивы точек для всего контура и описание геометрии;
	m1 = f1 == 0. ? N1-1 : 0; m2 = f2 == 0. ? N1-1 : 0; N = N1*N2-m1-m2;
	X    = (double *)new_struct(N*sizeof(double));
	Y    = (double *)new_struct(N*sizeof(double));
	geom = (int *)new_struct((1+2*(N2-1)*(N1+1)-m1-m2)*sizeof(int));
	if (! X || ! Y || ! geom) {
		delete_struct(X);
		delete_struct(Y);
		delete_struct(geom); N = 0; return(0);
	}

///////////////////////////////////////////////////////////
//...инициализиpуем узлы сетки на границе области и внутри;
	for (j = 0; j < N1; j++) {
		z1 = w[0]+j*(w[2]-w[0])/(N1-1.);
		z2 = w[1]+j*(w[3]-w[1])/(N1-1.);
		for (l = 0; l < N2; l++) {
				z = z1+l*(z2-z1)/(N2-1.);
				if (l == 0 && f1 == 0.) {
				  X[0] = real(z); Y[0] = imag(z);
				} 
				else
				if (l == N2-1 && f2 == 0.) {
				  X[N-1] = real(z); Y[N-1] = imag(z);
				} 
				else {
				  X[l*N1+j-m1] = real(z); Y[l*N1+j-m1] = imag(z);
				}
		}
	}

/////////////////////////////////////////////////////////////////////////////
//...заполняем массив, определяющий геометрию области и выходим из программы;
	set_geom(geom, N1, N2, m1, N, shift, f1, f2); 
	return(1);
}

///////////////////////////////////////////////////////////////////////
//...вспомогательная функция для построения 2D сетки в "кривом колене";
int CGrid_el::grid_knee(complex * w, double *& X, double *& Y, double h1, double h2,
                           int *& geom, int & N, int shift, int Max_N1, int Max_N2)
{
	double r1 = abs(w[2]-w[0]), r2 = abs(w[5]-w[3]), f1, f2, f3, f4;
	int    N1, N2, m1, m2, j, l; if (r1 < EE_ker) r1 = 0.; if (r2 < EE_ker) r2 = 0.;
	complex z, z1, z2;

////////////////////////////////////////////////////////////
//...определяем количество узлов сетки для "кривого колена";
	f1 = arc_length(w[0], w[1], w[2], 1);
	f2 = arc_length(w[3], w[4], w[5], 0);
	f3 = 0.5*(f1+f2);
	f4 = 0.5*(abs(w[5]-w[1])+abs(w[4]-w[2]));
	N1 = h1 > 0. ? 2+(int)ceil(f3/h1) : 1+abs(Max_N1);
	N2 = h2 > 0. ? 2+(int)ceil(f4/h2) : 1+abs(Max_N2);
	if (r1 == 0. && r2 == 0.) {
		X = Y = NULL; geom = NULL; N = 0; return(1);
	}
	if (r1 == 0.) f1 = 0.; else f1 /= r1;
	if (r2 == 0.) f2 = 0.; else f2 /= r2;

///////////////////////////////////////////////////////////////////////
//...распределяем массивы точек для всего контура и описание геометрии;
	m1 = r1 == 0. ? N1-1 : 0; m2 = r2 == 0. ? N1-1 : 0; N = N1*N2-m1-m2;
	X    = (double *)new_struct(N*sizeof(double));
	Y    = (double *)new_struct(N*sizeof(double));
	geom = (int *)new_struct((1+2*(N2-1)*(N1+1)-m1-m2)*sizeof(int));
	if (! X || ! Y || ! geom) {
		delete_struct(X);
		delete_struct(Y);
		delete_struct(geom); N = 0; return(0);
	}

///////////////////////////////////////////////////////////
//...инициализиpуем узлы сетки на границе области и внутри;
	for (j = 0; j < N1; j++) {
		if (w[0] == w[1]) z1 = w[0]+j*(w[2]-w[0])/(N1-1.);
		else              z1 = w[0]+(w[1]-w[0])*polar(1., -j*f1/(N1-1.));
		if (w[3] == w[4]) z2 = w[5]+j*(w[3]-w[5])/(N1-1.);
		else              z2 = w[3]+(w[5]-w[3])*polar(1., -j*f2/(N1-1.));
		for (l = 0; l < N2; l++) {
				z = z1+l*(z2-z1)/(N2-1.);
				if (l == 0 && r1 == 0.) {
				  X[0] = real(z); Y[0] = imag(z);
				} 
				else
				if (l == N2-1 && r2 == 0.) {
				  X[N-1] = real(z); Y[N-1] = imag(z);
				} 
				else {
				  X[l*N1+j-m1] = real(z); Y[l*N1+j-m1] = imag(z);
				}
		}
	}

/////////////////////////////////////////////////////////////////////////////
//...заполняем массив, определяющий геометрию области и выходим из программы;
	set_geom(geom, N1, N2, m1, N, shift, r1, r2); 
	return(1);
}

/////////////////////////////////////////////////////////////////////////
//...рабочая функция, заполняющая массив, определяющий геометрию области;
void CGrid_el::set_geom(int *& geom, int &N1, int &N2, int &m1, int &N, int &shift, double &f1, double &f2)
{
  int k = 0, j, l;

  geom[k++] = N2-1;
  for (l = 1; l < N2; l++) {
      if (l == 1 && f1 == 0.) {
         geom[k++] = GL_TRIANGLE_FAN;
         geom[k++] = N1+1; geom[k++] = shift;
         for (j = 0; j < N1; j++) geom[k++] = l*N1+j-m1+shift;
      } 
		else
      if (l == N2-1 && f2 == 0.) {
         geom[k++] = GL_TRIANGLE_FAN;
         geom[k++] = N1+1; geom[k++] = N-1+shift;
         for (j = 0; j < N1; j++) geom[k++] = (l-1)*N1+(N1-j-1)-m1+shift;
      } 
		else {
         geom[k++] = GL_QUAD_STRIP;
         geom[k++] = 2*N1;
         for (j = 0; j < N1; j++) {
             geom[k++] = (l-1)*N1-m1+j+shift;
             geom[k++] = l*N1+j-m1+shift;
         }
      }
  }
  return;
}

///////////////////////////////////////
//...test for internal point of circle;
int CGrid_el::in_circle(CMap * mp, double X, double Y, double eps)
{
	if (mp && ID_MAP(1, SPHERE_GENUS) == mp[0] && fabs(mp[7]) > fabs(eps)) {
		double A = fabs(mp[7]), A_inv = 1./A, 
				CZ = cos(mp[4]), SZ = sin(mp[4]), CY = cos(mp[5]), SY = sin(mp[5]),
				CX = cos(mp[6]), SX = sin(mp[6]), P[3];
			P[0] = X-mp[1];
			P[1] = Y-mp[2];
			P[2] = 0.;
			point_iso<double>(P, NULL, CX, -SX, CY, -SY, CZ, -SZ);
			P[0] *= A_inv;
			P[1] *= A_inv;
			if (P[0]*P[0]+P[1]*P[1] < 1.+eps) return(1);
	}
	return(0);
}

////////////////////////////////////////
//...test for internal point of ellipse;
int CGrid_el::in_ellipse(CMap * mp, double X, double Y, double eps)
{
	if (mp && ID_MAP(1, CYL_GENUS) == mp[0] && fabs(mp[7]) > fabs(eps) && fabs(mp[8]) > fabs(eps)) {
		double A = fabs(mp[7]), B = fabs(mp[8]), A_inv = 1./A, B_inv = 1./B, 
				CZ = cos(mp[4]), SZ = sin(mp[4]), CY = cos(mp[5]), SY = sin(mp[5]),
				CX = cos(mp[6]), SX = sin(mp[6]), P[3];
			P[0] = X-mp[1];
			P[1] = Y-mp[2];
			P[2] = 0.;
			point_iso<double>(P, NULL, CX, -SX, CY, -SY, CZ, -SZ);
			P[0] *= A_inv;
			P[1] *= B_inv;
			if (P[0]*P[0]+P[1]*P[1] < 1.+eps) return(1);
	}
	return(0);
}

///////////////////////////////////////
//...test for internal point of sphere;
int CGrid_el::in_sphere(CMap * mp, double X, double Y, double Z, double eps)
{
	if (mp && ID_MAP(2, SPHERE_GENUS) == mp[0] && fabs(mp[7]) > fabs(eps)) {
		double A = fabs(mp[7]), A_inv = 1./A, 
				CZ = cos(mp[4]), SZ = sin(mp[4]), CY = cos(mp[5]), SY = sin(mp[5]),
				CX = cos(mp[6]), SX = sin(mp[6]), dd, P[3];
			P[0] = X-mp[1];
			P[1] = Y-mp[2];
			P[2] = Z-mp[3];
			point_iso<double>(P, NULL, CX, -SX, CY, -SY, CZ, -SZ);
			P[0] *= A_inv;
			P[1] *= A_inv;
			P[2] *= A_inv;

			dd = P[0]*P[0]+P[1]*P[1]+P[2]*P[2];
			if (mp[7] > 0. && dd < 1.+eps || 
				 mp[7] < 0. && dd > 1.-eps) return(1);
	}
	return(0);
}

/////////////////////////////////////////
//...test for internal point of spheroid;
int CGrid_el::in_spheroid(CMap * mp, double X, double Y, double Z, double eps)
{
	if (mp && ID_MAP(2, SPHEROID_GENUS) == mp[0] && fabs(mp[7]) > fabs(eps) && fabs(mp[8]) > fabs(eps)) {
		double A = fabs(mp[7]), B = fabs(mp[8]), A_inv = 1./A, B_inv = 1./B, 
				CZ = cos(mp[4]), SZ = sin(mp[4]), CY = cos(mp[5]), SY = sin(mp[5]),
				CX = cos(mp[6]), SX = sin(mp[6]), P[3];
			P[0] = X-mp[1];
			P[1] = Y-mp[2];
			P[2] = Z-mp[3];
			point_iso<double>(P, NULL, CX, -SX, CY, -SY, CZ, -SZ);
			P[0] *= B_inv;
			P[1] *= B_inv;
			P[2] *= A_inv;
			if (P[0]*P[0]+P[1]*P[1]+P[2]*P[2] < 1.+eps) return(1);
	}
	return(0);
}

//////////////////////////////////////
//...test for internal point of sheet;
int CGrid_el::in_sheet(CMap * mp, double X, double Y, double eps)
{
	int m;
	if (mp && ID_MAP(2, NULL_GENUS) == mp[0] && mp[m = size_of_map(2, NULL_GENUS)] == (CMap)SHEET_CELL) {
		double CZ = cos(mp[4]), SZ = sin(mp[4]), CY = cos(mp[5]), SY = sin(mp[5]),
				 CX = cos(mp[6]), SX = sin(mp[6]), P[3];
		P[0] = X-mp[1];
		P[1] = Y-mp[2];
		P[2] = 0.;
		point_iso<double>(P, NULL, CX, -SX, CY, -SY, CZ, -SZ);
		if (-eps < P[0]+mp[m+1]*.5 && P[0]-mp[m+1]*.5 < eps && 
			 -eps < P[1]+mp[m+2]*.5 && P[1]-mp[m+2]*.5 < eps) return(1);
	}
	return(0);
}

//////////////////////////////////////
//...test for internal point of box;
int CGrid_el::in_box(CMap * mp, double X, double Y, double Z, double eps)
{
	int m;
	if (mp && ID_MAP(3, NULL_GENUS) == mp[0] && mp[m = size_of_map(2, NULL_GENUS)] == (CMap)BOX_CELL) {
		double CZ = cos(mp[4]), SZ = sin(mp[4]), CY = cos(mp[5]), SY = sin(mp[5]),
				 CX = cos(mp[6]), SX = sin(mp[6]), P[3];
		P[0] = X-mp[1];
		P[1] = Y-mp[2];
		P[1] = Z-mp[3];
		point_iso<double>(P, NULL, CX, -SX, CY, -SY, CZ, -SZ);
		if (-eps < P[0]+mp[m+1]*.5 && P[0]-mp[m+1]*.5 < eps && 
			 -eps < P[1]+mp[m+2]*.5 && P[1]-mp[m+2]*.5 < eps && 
			 -eps < P[2]+mp[m+3]*.5 && P[2]-mp[m+3]*.5 < eps) return(1);
	}
	return(0);
}

