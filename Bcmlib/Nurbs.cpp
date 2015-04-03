#include "stdafx.h"
#include "nurbs.h"

////////////////////////////////////////////////////////////////////
//...функция, вычисляющая значение сплайновых функций в одной точке;
void nurbs_spl(double t, double * t_a, int N_t, int q, int & i, double * N)
{
	memset(N, 0x0000, q*sizeof(double));
	N[0] = 1.;
	i    = 0;
/////////////////////////
//...определяем интервал;
	while (i < N_t && t >= t_a[i]) i++; i = min(i, N_t-q); if (! i) return;

//////////////////////////////////////////////////////////
//...вычисляем примыкающие к интервалу сплайновые функции;
#ifdef ANOTHER_ALGORITHM
	double temp; int k, j;
	for (k = 1; k < q; k++) {
		 N[0] = (temp = N[0])*(t_a[i]-t)/(t_a[i]-t_a[i-k]);
		 for (j = 1; j < k; j++)
		 N[j] =  temp *   (t-t_a[i+j-k-1])/(t_a[i+j-1]-t_a[i+j-k-1])+
				  (temp = N[j])*(t_a[i+j]-t)/(t_a[i+j]-t_a[i+j-k]);
		 N[k] =  temp * (t-t_a[i-1])/(t_a[i+k-1]-t_a[i-1]);
	}
#else
	double temp1, temp2;
	for (int k = 1; k < q; k++) {
		N[0] = (temp1 = N[0])*(t_a[i]-t)/(t_a[i]-t_a[i-k]);
		for (int j = 1; j < k; j++) {
			N[j] =  temp1*(t-t_a[i+j-k-1])/(t_a[i+j-1]-t_a[i+j-k-1])+
					 (temp2 = N[j])*(t_a[i+j]-t)/(t_a[i+j]-t_a[i+j-k]);
			temp1 = temp2;
		}
		N[k] = temp1*(t-t_a[i-1])/(t_a[i+k-1]-t_a[i-1]);
	}
#endif
	return;
}

/////////////////////////////////////////////////////////////////////////////////////
//...функция, вычисляющая значение сплайновых функций и их производных в одной точке;
void nurbs_spl(double t, double * t_a, int N_t, int q, int & i, double * N, double * D)
{
	memset(N, 0x0000, q*sizeof(double));
	memset(D, 0x0000, q*sizeof(double));
	N[0] = 1.;
	i    = 0;
/////////////////////////
//...определяем интервал;
	while (i < N_t && t >= t_a[i]) i++; i = min(i, N_t-q); if (! i) return;

///////////////////////////////////////////////////////////////////////////
//...вычисляем примыкающие к интервалу сплайновые функции и их производные;
	double temp1, temp2, d1, d2, f, g;
	for (int k = 1; k < q; k++) {
		 D[0] = ((d1   = D[0])*(t_a[i]-t)-N[0])*(f = 1./(t_a[i]-t_a[i-k]));
		 N[0] = (temp1 = N[0])*(t_a[i]-t)*f;
		 for (int j  = 1; j < k; j++) {
			 D[j] =  (d1 * (t-t_a[i+j-k-1]) + temp1)*(f = 1./(t_a[i+j-1]-t_a[i+j-k-1]))+
					  ((d2 = D[j])*(t_a[i+j]-t)- N[j])*(g = 1./(t_a[i+j]-t_a[i+j-k]));
			 N[j] =  temp1*(t-t_a[i+j-k-1])*f+(temp2 = N[j])*(t_a[i+j]-t)*g;
			d1 = d2;	temp1 = temp2;
		 }
		 D[k] =  (d1 * (t-t_a[i-1])+temp1)*(f = 1./(t_a[i+k-1]-t_a[i-1]));
		 N[k] =  temp1*(t-t_a[i-1])*f;
	}
	return;
}


////////////////////////////////////////////////////////////////////////////////////////////
//...функция, вычисляющая значение сплайновых функций и их вторых производных в одной точке;
void nurbs_spl(double t, double * t_a, int N_t, int q, int & i, double * N, double * D, double * K)
{
	memset(N, 0x0000, q*sizeof(double));
	memset(D, 0x0000, q*sizeof(double));
	memset(K, 0x0000, q*sizeof(double));
	N[0] = 1.;
	i    = 0;
/////////////////////////
//...определяем интервал;
	while (i < N_t && t >= t_a[i]) i++; i = min(i, N_t-q); if (! i) return;

///////////////////////////////////////////////////////////////////////////
//...вычисляем примыкающие к интервалу сплайновые функции и их производные;
	double temp1, temp2, d1, d2, h1, h2, f, g;
	for (int k = 1; k < q; k++) {
		 K[0] = ((h1   = K[0])*(t_a[i]-t)-D[0]*2.)*(f = 1./(t_a[i]-t_a[i-k]));
		 D[0] = ((d1   = D[0])*(t_a[i]-t)-N[0])*f;
		 N[0] = (temp1 = N[0])*(t_a[i]-t)*f;
		 for (int j  = 1; j < k; j++) {
			 K[j] = (h1*(t-t_a[i+j-k-1])+d1*2.)*(f = 1./(t_a[i+j-1]-t_a[i+j-k-1]))+
					 ((h2 = K[j])*(t_a[i+j]-t)-D[j]*2.)*(g = 1./(t_a[i+j]-t_a[i+j-k]));
			 D[j] = (d1*(t-t_a[i+j-k-1])+temp1)*f+((d2 = D[j])*(t_a[i+j]-t)-N[j])*g;
			 N[j] = temp1*(t-t_a[i+j-k-1])*f+(temp2 = N[j])*(t_a[i+j]-t)*g;
			h1 = h2; d1 = d2; temp1 = temp2;
		 }
		 K[k] = (h1*(t-t_a[i-1])+d1*2)*(f = 1./(t_a[i+k-1]-t_a[i-1]));
		 D[k] = (d1*(t-t_a[i-1])+temp1)*f;
		 N[k] = temp1*(t-t_a[i-1])*f;
	}
	return;
}

////////////////////////////////////////////////////////////////////////////////
//...функция, вычисляющая значение сплайна в 3D по заданным сплайновым функциям;
void nurbs_P(double * N, double * p, int q, int i, int dim, double * P)
{
	memset(P, 0x0000, dim*sizeof(double));
	for (int k = 0; k < q;   k++)
	for (int l = 0; l < dim; l++) P[l] += N[k]*p[(i+k)*dim+l];
	return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//...функция, вычисляющая значение сплайна и его производной в 3D по заданным сплайновым функциям;
void nurbs_P(double * N, double * D, double * p, int q, int i, int dim, double * P, double * Q)
{
	memset(P, 0x0000, dim*sizeof(double));
	memset(Q, 0x0000, dim*sizeof(double));
	for (int k = 0; k < q;   k++)
	for (int l = 0; l < dim; l++) {
		P[l] += N[k]*p[(i+k)*dim+l];
		Q[l] += D[k]*p[(i+k)*dim+l];
	}
	return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
//...функция, вычисляющая значение сплайна и его второй производной в 3D по заданным сплайновым функциям;
void nurbs_P(double * N, double * D, double * K, double * p, int q, int i, int dim, double * P, double * Q, double * R)
{
	memset(P, 0x0000, dim*sizeof(double));
	memset(Q, 0x0000, dim*sizeof(double));
	memset(R, 0x0000, dim*sizeof(double));
	for (int k = 0; k < q;   k++)
	for (int l = 0; l < dim; l++) {
		P[l] += N[k]*p[(i+k)*dim+l];
		Q[l] += D[k]*p[(i+k)*dim+l];
		R[l] += K[k]*p[(i+k)*dim+l];
	}
	return;
}

///////////////////////////////////////////////////////////////////////////////////////////
//...функция, вычисляющая значение двумерного сплайна в 3D по заданным сплайновым функциям;
void nurbs_P(double * N_u, double * N_v, double * p, int q_u, int q_v, int icp, int i, int j,
                                                     int dim, double * P)
{
	memset(P, 0x0000, dim*sizeof(double));
	for (int k = 0; k < q_u; k++)
	for (int m = 0; m < q_v; m++)
	for (int l = 0; l < dim; l++) P[l] += N_u[k]*N_v[m]*p[((j+m)*icp+i+k)*dim+l];
	return;
}

/////////////////////////////////////////////////////////////////////////////
//...функция, вычисляющая значение двумерного сплайна и его производных в 3D;
void nurbs_P(double * N_u, double * D_u, double * N_v, double * D_v, double * p, int q_u, int q_v,
                           int icp, int i, int j, int dim, double * P, double * Q, double * R)
{
	memset(P, 0x0000, dim*sizeof(double));
	memset(Q, 0x0000, dim*sizeof(double));
	memset(R, 0x0000, dim*sizeof(double));
	for (int k = 0; k < q_u; k++)
	for (int m = 0; m < q_v; m++)
	for (int l = 0; l < dim; l++) {
		P[l] += N_u[k]*N_v[m]*p[((j+m)*icp+i+k)*dim+l];
		Q[l] += D_u[k]*N_v[m]*p[((j+m)*icp+i+k)*dim+l];
		R[l] += N_u[k]*D_v[m]*p[((j+m)*icp+i+k)*dim+l];
	}
	return;
}

////////////////////////////////////////////////////////////////////////////////////
//...функция, вычисляющая значение двумерного сплайна и его вторых производных в 3D;
void nurbs_P(double * N_u, double * D_u, double * K_u, double * N_v, double * D_v, double * K_v, double * p, int q_u, int q_v,
                           int icp, int i, int j, int dim, double * P, double * Q, double * R, double * S, double * T, double * W)
{
	memset(P, 0x0000, dim*sizeof(double));
	memset(Q, 0x0000, dim*sizeof(double));
	memset(R, 0x0000, dim*sizeof(double));
	memset(S, 0x0000, dim*sizeof(double));
	memset(T, 0x0000, dim*sizeof(double));
	memset(W, 0x0000, dim*sizeof(double));
	for (int k = 0; k < q_u; k++)
	for (int m = 0; m < q_v; m++)
	for (int l = 0; l < dim; l++) {
		P[l] += N_u[k]*N_v[m]*p[((j+m)*icp+i+k)*dim+l];
		Q[l] += D_u[k]*N_v[m]*p[((j+m)*icp+i+k)*dim+l];
		R[l] += N_u[k]*D_v[m]*p[((j+m)*icp+i+k)*dim+l];
		S[l] += K_u[k]*N_v[m]*p[((j+m)*icp+i+k)*dim+l];
		T[l] += D_u[k]*D_v[m]*p[((j+m)*icp+i+k)*dim+l];
		W[l] += N_u[k]*K_v[m]*p[((j+m)*icp+i+k)*dim+l];
	}
	return;
}

////////////////////////////////////////////////////////////////////////
//...функция, вычисляющая значение числителя рационального сплайна в 3D;
void nurbs_wP(double * N, double * p, int q, int i, int dim, double * P)
{
	memset(P, 0x0000, (dim-1)*sizeof(double));
	for (int k = 0; k < q;     k++)
	for (int l = 0; l < dim-1; l++) P[l] += N[k]*p[(i+k)*dim+l]*p[(i+k+1)*dim-1];
	return;
}

//////////////////////////////////////////////////////////////////////////////////////////
//...функция, вычисляющая значение числителя рационального сплайна в 3D и его производной;
void nurbs_wP(double * N, double * D, double * p, int q, int i, int dim, double * P, double * Q)
{
	memset(P, 0x0000, (dim-1)*sizeof(double));
	memset(Q, 0x0000, (dim-1)*sizeof(double));
	for (int k = 0; k < q;     k++)
	for (int l = 0; l < dim-1; l++) {
		P[l] += N[k]*p[(i+k)*dim+l]*p[(i+k+1)*dim-1];
		Q[l] += D[k]*p[(i+k)*dim+l]*p[(i+k+1)*dim-1];
	}
	return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
//...функция, вычисляющая значение числителя рационального сплайна в 3D и его второй производной;
void nurbs_wP(double * N, double * D, double * K, double * p, int q, int i, int dim, double * P, double * Q, double * R)
{
	memset(P, 0x0000, (dim-1)*sizeof(double));
	memset(Q, 0x0000, (dim-1)*sizeof(double));
	memset(R, 0x0000, (dim-1)*sizeof(double));
	for (int k = 0; k < q;     k++)
	for (int l = 0; l < dim-1; l++) {
		P[l] += N[k]*p[(i+k)*dim+l]*p[(i+k+1)*dim-1];
		Q[l] += D[k]*p[(i+k)*dim+l]*p[(i+k+1)*dim-1];
		R[l] += D[k]*p[(i+k)*dim+l]*p[(i+k+1)*dim-1];
	}
	return;
}

///////////////////////////////////////////////////////////////////////////////////
//...функция, вычисляющая значение числителя двумерного рационального сплайна в 3D;
void nurbs_wP(double * N_u, double * N_v, double * p, int q_u, int q_v, int icp, int i, int j,
                                                      int dim, double * P)
{
	memset(P, 0x0000, (dim-1)*sizeof(double));
	for (int k = 0; k < q_u;   k++)
	for (int m = 0; m < q_v;   m++)
	for (int l = 0; l < dim-1; l++)
		P[l] += N_u[k]*N_v[m]*p[((j+m)*icp+i+k)*dim+l]*p[((j+m)*icp+i+k+1)*dim-1];
	return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
//...функция, вычисляющая значение числителя двумерного рационального сплайна в 3D и производных;
void nurbs_wP(double * N_u, double * D_u, double * N_v, double * D_v, double * p, int q_u, int q_v,
                            int icp, int i, int j, int dim, double * P, double * Q, double * R)
{
	memset(P, 0x0000, (dim-1)*sizeof(double));
	memset(Q, 0x0000, (dim-1)*sizeof(double));
	memset(R, 0x0000, (dim-1)*sizeof(double));
	for (int k = 0; k < q_u;   k++)
	for (int m = 0; m < q_v;   m++)
	for (int l = 0; l < dim-1; l++) {
		P[l] += N_u[k]*N_v[m]*p[((j+m)*icp+i+k)*dim+l]*p[((j+m)*icp+i+k+1)*dim-1];
		Q[l] += D_u[k]*N_v[m]*p[((j+m)*icp+i+k)*dim+l]*p[((j+m)*icp+i+k+1)*dim-1];
		R[l] += N_u[k]*D_v[m]*p[((j+m)*icp+i+k)*dim+l]*p[((j+m)*icp+i+k+1)*dim-1];
	}
	return;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//...функция, вычисляющая значение числителя двумерного рационального сплайна в 3D и вторых производных;
void nurbs_wP(double * N_u, double * D_u, double * K_u, double * N_v, double * D_v, double * K_v, double * p, int q_u, int q_v,
                            int icp, int i, int j, int dim, double * P, double * Q, double * R, double * S, double * T, double * W)
{
	memset(P, 0x0000, (dim-1)*sizeof(double));
	memset(Q, 0x0000, (dim-1)*sizeof(double));
	memset(R, 0x0000, (dim-1)*sizeof(double));
	memset(S, 0x0000, (dim-1)*sizeof(double));
	memset(T, 0x0000, (dim-1)*sizeof(double));
	memset(W, 0x0000, (dim-1)*sizeof(double));
	for (int k = 0; k < q_u;   k++)
	for (int m = 0; m < q_v;   m++)
	for (int l = 0; l < dim-1; l++) {
		P[l] += N_u[k]*N_v[m]*p[((j+m)*icp+i+k)*dim+l]*p[((j+m)*icp+i+k+1)*dim-1];
		Q[l] += D_u[k]*N_v[m]*p[((j+m)*icp+i+k)*dim+l]*p[((j+m)*icp+i+k+1)*dim-1];
		R[l] += N_u[k]*D_v[m]*p[((j+m)*icp+i+k)*dim+l]*p[((j+m)*icp+i+k+1)*dim-1];
		S[l] += K_u[k]*N_v[m]*p[((j+m)*icp+i+k)*dim+l]*p[((j+m)*icp+i+k+1)*dim-1];
		T[l] += D_u[k]*D_v[m]*p[((j+m)*icp+i+k)*dim+l]*p[((j+m)*icp+i+k+1)*dim-1];
		W[l] += N_u[k]*K_v[m]*p[((j+m)*icp+i+k)*dim+l]*p[((j+m)*icp+i+k+1)*dim-1];
	}
	return;
}

//////////////////////////////////////////////////////////////////////////
//...функция, вычисляющая значение знаменателя рационального сплайна в 3D;
void nurbs_w(double * N, double * p, int q, int i, int dim, double & f)
{
	f = 0.;
	for (int k = 0; k < q; k++) f += N[k]*p[(i+k+1)*dim-1];
	return;
}

////////////////////////////////////////////////////////////////////////////////////////////
//...функция, вычисляющая значение знаменателя рационального сплайна в 3D и его производной;
void nurbs_w(double * N, double * D, double * p, int q, int i, int dim, double & f, double & d)
{
	f = d = 0.;
	for (int k = 0; k < q; k++) {
		f += N[k]*p[(i+k+1)*dim-1];
		d += D[k]*p[(i+k+1)*dim-1];
	}
	return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//...функция, вычисляющая значение знаменателя рационального сплайна в 3D и его второй производной;
void nurbs_w(double * N, double * D, double * K, double * p, int q, int i, int dim, double & f, double & d, double & h)
{
	f = d = h = 0.;
	for (int k = 0; k < q; k++) {
		f += N[k]*p[(i+k+1)*dim-1];
		d += D[k]*p[(i+k+1)*dim-1];
		h += K[k]*p[(i+k+1)*dim-1];
	}
	return;
}

/////////////////////////////////////////////////////////////////////////////////////
//...функция, вычисляющая значение знаменателя двумерного рационального сплайна в 3D;
void nurbs_w(double * N_u, double * N_v, double * p, int q_u, int q_v, int icp, int i, int j,
                                                     int dim, double & f)
{
	f = 0.;
	for (int k = 0; k < q_u;   k++)
	for (int m = 0; m < q_v;   m++) f += N_u[k]*N_v[m]*p[((j+m)*icp+i+k+1)*dim-1];
	return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//...функция, вычисляющая значение знаменателя двумерного рационального сплайна в 3D и производных;
void nurbs_w(double * N_u, double * D_u, double * N_v, double * D_v, double * p, int q_u, int q_v,
                           int icp, int i, int j, int dim, double & f, double & d, double & r)
{
	f = d = r = 0.;
	for (int k = 0; k < q_u;   k++)
	for (int m = 0; m < q_v;   m++) {
		f += N_u[k]*N_v[m]*p[((j+m)*icp+i+k+1)*dim-1];
		d += D_u[k]*N_v[m]*p[((j+m)*icp+i+k+1)*dim-1];
		r += N_u[k]*D_v[m]*p[((j+m)*icp+i+k+1)*dim-1];
	}
	return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//...функция, вычисляющая значение знаменателя двумерного рационального сплайна в 3D и вторых производных;
void nurbs_w(double * N_u, double * D_u, double * K_u, double * N_v, double * D_v, double * K_v, double * p, int q_u, int q_v,
                           int icp, int i, int j, int dim, double & f, double & d, double & r, double & h, double & q, double & w)
{
	f = d = r = h = q = w = 0.;
	for (int k = 0; k < q_u;   k++)
	for (int m = 0; m < q_v;   m++) {
		f += N_u[k]*N_v[m]*p[((j+m)*icp+i+k+1)*dim-1];
		d += D_u[k]*N_v[m]*p[((j+m)*icp+i+k+1)*dim-1];
		r += N_u[k]*D_v[m]*p[((j+m)*icp+i+k+1)*dim-1];
		h += K_u[k]*N_v[m]*p[((j+m)*icp+i+k+1)*dim-1];
		q += D_u[k]*D_v[m]*p[((j+m)*icp+i+k+1)*dim-1];
		w += N_u[k]*K_v[m]*p[((j+m)*icp+i+k+1)*dim-1];
	}
	return;
}

//////////////////////////////////////////////////////////////
//...функция, вычисляющая значение рационального сплайна в 2D;
void nurbs2D(double t, double * t_a, double * p, int N_t, int q, int dim, double * P)
{
	double N[MAX_SPLINE], f;
	if (q <= MAX_SPLINE) {
		int i;
		nurbs_spl(t, t_a, N_t, q, i, N); if (i < q) return;
		if (dim == 3) {
			nurbs_wP(N, p, q, i-q, dim, P);
			nurbs_w (N, p, q, i-q, dim, f);
			P[0] *= (f = 1./f);
			P[1] *= f;
			P[2] *= f;
		}
		else nurbs_P(N, p, q, i-q, dim, P);
	}
	return;
}

/////////////////////////////////////////////////
//...функция, вычисляющая сплайновую кривую в 2D;
void nurbs2D(double t, CMap * mp, double * P, int id_map)
{
	if (mp   && ID_MAP(1, NURBS_GENUS)  == mp[0] &&
		mp[size_of_map(1, NURBS_GENUS)] == (CMap)B1SPLINE_CELL && P) {
		int m1 = (int)(mp[8]-mp[9]), m2 = (int)mp[9], m3 = (int)mp[7]/2;
		if (! id_map) nurbs2D(t, mp+11, mp+m1+m2+13, m1+m2, m1, m3, P); 
		else {
			memset(P, 0x0000, (m3-1)*sizeof(double));
			if (mp[m1+m2+12] > mp[m1+m2+11]+EE_ker)
			nurbs2D(mp[m1+m2+11]+t*(mp[m1+m2+12]-mp[m1+m2+11]), 
								mp+11, mp+m1+m2+13, m1+m2, m1, m3, P);
		}
	}
	return;
}

//////////////////////////////////////////////////////////////
//...функция, вычисляющая значение рационального сплайна в 3D;
void nurbs3D(double t, double * t_a, double * p, int N_t, int q, int dim, double * P)
{
	double N[MAX_SPLINE], f;
	if (q <= MAX_SPLINE) {
		int i;
		nurbs_spl(t, t_a, N_t, q, i, N); if (i < q) return;
		if (dim == 4) {
			nurbs_wP(N, p, q, i-q, dim, P);
			nurbs_w (N, p, q, i-q, dim, f);
			P[0] *= (f = 1./f);
			P[1] *= f;
			P[2] *= f;
		}
		else nurbs_P(N, p, q, i-q, dim, P);
	}
	return;
}

/////////////////////////////////////////////////////////////////////////////////////////////
//...функция, вычисляющая значение рационального сплайна в 3D и его производной по параметру;
void nurbs3D(double t, double * t_a, double * p, int N_t, int q, int dim, double * P, double * Q)
{
	double N[MAX_SPLINE], D[MAX_SPLINE], f, d;
	if (q <= MAX_SPLINE) {
		int i;
		nurbs_spl(t, t_a, N_t, q, i, N, D); if (i < q) return;
		if (dim == 4) {
			nurbs_wP(N, D, p, q, i-q, dim, P, Q);
			nurbs_w (N, D, p, q, i-q, dim, f, d);
			P[0] *= (f = 1./f); Q[0] *= f; Q[0] -= P[0]*(d *= f);
			P[1] *= f;          Q[1] *= f; Q[1] -= P[1]*d;
			P[2] *= f;          Q[2] *= f; Q[2] -= P[2]*d;
		}
		else nurbs_P(N, D, p, q, i-q, dim, P, Q);
	}
	return;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
//...функция, вычисляющая значение рационального сплайна в 3D и его второй производной по параметру;
void nurbs3D(double t, double * t_a, double * p, int N_t, int q, int dim, double * P, double * Q, double * R)
{
	double N[MAX_SPLINE], D[MAX_SPLINE], K[MAX_SPLINE], f, d, h;
	if (q <= MAX_SPLINE) {
		int i;
		nurbs_spl(t, t_a, N_t, q, i, N, D, K); if (i < q) return;
		if (dim == 4) {
			nurbs_wP(N, D, K, p, q, i-q, dim, P, Q, R);
			nurbs_w (N, D, K, p, q, i-q, dim, f, d, h);
			P[0] *= (f = 1./f); Q[0] *= f; Q[0] -= P[0]*(d *= f); R[0] *= f; R[0] -= P[0]*(h *= f)+2.*Q[0]*d;
			P[1] *= f;          Q[1] *= f; Q[1] -= P[1]*d;        R[1] *= f; R[1] -= P[1]*h+2.*Q[1]*d;
			P[2] *= f;          Q[2] *= f; Q[2] -= P[2]*d;        R[2] *= f; R[2] -= P[2]*h+2.*Q[2]*d;
		}
		else nurbs_P(N, D, K, p, q, i-q, dim, P, Q, R);
	}
	return;
}

/////////////////////////////////////////////////
//...функция, вычисляющая сплайновую кривую в 3D;
void nurbs3D(double t, CMap * mp, double * P, int id_map)
{
	if (mp   && ID_MAP(1, NURBS_GENUS)  == mp[0] &&
		mp[size_of_map(1, NURBS_GENUS)] == (CMap)B1SPLINE_CELL && P) {
		int m1 = (int)(mp[8]-mp[9]), m2 = (int)mp[9], m3 = (int)mp[7]/2;
		if (! id_map) nurbs3D(t, mp+11, mp+m1+m2+13, m1+m2, m1, m3, P); 
		else {
			memset(P, 0x0000, (m3-1)*sizeof(double));
			if (mp[m1+m2+12] > mp[m1+m2+11]+EE_ker)
			nurbs3D(mp[m1+m2+11]+t*(mp[m1+m2+12]-mp[m1+m2+11]), 
								mp+11, mp+m1+m2+13, m1+m2, m1, m3, P);
		}
	}
	return;
}

///////////////////////////////////////////////////////////////
//...функция, вычисляющая сплайновую кривую в 3D и касательную;
void nurbs3D(double t, CMap * mp, double * P, double * ort, int id_map)
{
	if (mp   && ID_MAP(1, NURBS_GENUS)  == mp[0] &&
		mp[size_of_map(1, NURBS_GENUS)] == (CMap)B1SPLINE_CELL && P && ort) {
		int m1 = (int)(mp[8]-mp[9]), m2 = (int)mp[9], m3 = (int)mp[7]/2;
		if (! id_map) nurbs3D(t, mp+11, mp+m1+m2+13, m1+m2, m1, m3, P, ort); 
		else {
			memset(P,   0x0000, (m3-1)*sizeof(double));
			memset(ort, 0x0000, (m3-1)*sizeof(double)); ort[2] = 1.;
			if (mp[m1+m2+12] > mp[m1+m2+11]+EE_ker)
			nurbs3D(mp[m1+m2+11]+t*(mp[m1+m2+12]-mp[m1+m2+11]), 
								mp+11, mp+m1+m2+13, m1+m2, m1, m3, P, ort);
		}
		double f = sqrt(ort[0]*ort[0]+ort[1]*ort[1]+ort[2]*ort[2]);
		ort[0] *= (f = f > EE ? 1./f : 1.);
		ort[1] *= f;
		ort[2] *= f;
	}
	return;
}

/////////////////////////////////////////////////////////////////////////////////////////////
//...функция, вычисляющая сплайновую кривую в 3D, касательную и обратный радиус кривизны R_t[0];
void nurbs3D(double t, CMap * mp, double * P, double * ort, double * R_t, int id_map)
{
	if (mp   && ID_MAP(1, NURBS_GENUS)  == mp[0] &&
		mp[size_of_map(1, NURBS_GENUS)] == (CMap)B1SPLINE_CELL && P && ort) {
		int m1 = (int)(mp[8]-mp[9]), m2 = (int)mp[9], m3 = (int)mp[7]/2;
		if (! id_map) nurbs3D(t, mp+11, mp+m1+m2+13, m1+m2, m1, m3, P, ort, R_t); 
		else {
			memset(P,   0x0000, (m3-1)*sizeof(double));
			memset(ort, 0x0000, (m3-1)*sizeof(double)); ort[2] = 1.;
			memset(R_t, 0x0000, (m3-1)*sizeof(double));
			if (mp[m1+m2+12] > mp[m1+m2+11]+EE_ker)
			nurbs3D(mp[m1+m2+11]+t*(mp[m1+m2+12]-mp[m1+m2+11]), 
								mp+11, mp+m1+m2+13, m1+m2, m1, m3, P, ort, R_t);
		}
		double f = sqrt(ort[0]*ort[0]+ort[1]*ort[1]+ort[2]*ort[2]), 
				 d = sqrt(R_t[0]*R_t[0]+R_t[1]*R_t[1]+R_t[2]*R_t[2]), 
				 h =		 ort[0]*R_t[0]+ort[1]*R_t[1]+ort[2]*R_t[2];
		if (f > EE) R_t[0] = sqrt(sqr(f)*sqr(d)-h*h)/(f*f*f); else R_t[0] = 0.;
		ort[0] *= (f = f > EE ? 1./f : 1.);
		ort[1] *= f;
		ort[2] *= f;
	}
	return;
}

///////////////////////////////////////////////////////////
//...функция, вычисляющая значение двумерного сплайна в 3D;
void nurbs3D(double u, double * t_u, int n_u, int q_u,
             double v, double * t_v, int n_v, int q_v, double * p, int icp, int dim, double * P)
{
	double N_u[MAX_SPLINE], N_v[MAX_SPLINE], f;
	if (q_u <= MAX_SPLINE && q_v <= MAX_SPLINE) {
		int i, j;
		nurbs_spl(u, t_u, n_u, q_u, i, N_u); if (i < q_u) return;
		nurbs_spl(v, t_v, n_v, q_v, j, N_v); if (j < q_v) return;
		if (dim == 4) {
			nurbs_wP(N_u, N_v, p, q_u, q_v, icp, i-q_u, j-q_v, dim, P);
			nurbs_w (N_u, N_v, p, q_u, q_v, icp, i-q_u, j-q_v, dim, f);
			P[0] *= (f = 1./f);
			P[1] *= f;
			P[2] *= f;
		}
		else nurbs_P(N_u, N_v, p, q_u, q_v, icp, i-q_u, j-q_v, dim, P);
	}
	return;
}
	
///////////////////////////////////////////////////////////////////////////////////////////
//...функция, вычисляющая значение двумерного сплайна в 3D и его производные по параметрам;
void nurbs3D(double u, double * t_u, int n_u, int q_u,
             double v, double * t_v, int n_v, int q_v, double * p, int icp, int dim,
             double * P, double * Q, double * R)
{
	double N_u[MAX_SPLINE], D_u[MAX_SPLINE], N_v[MAX_SPLINE], D_v[MAX_SPLINE], f, d, r;
	if (q_u <= MAX_SPLINE && q_v <= MAX_SPLINE) {
		int i, j;
		nurbs_spl(u, t_u, n_u, q_u, i, N_u, D_u); if (i < q_u) return;
		nurbs_spl(v, t_v, n_v, q_v, j, N_v, D_v); if (j < q_v) return;
		if (dim == 4) {
			nurbs_wP(N_u, D_u, N_v, D_v, p, q_u, q_v, icp, i-q_u, j-q_v, dim, P, Q, R);
			nurbs_w (N_u, D_u, N_v, D_v, p, q_u, q_v, icp, i-q_u, j-q_v, dim, f, d, r);
			P[0] *= (f = 1./f); Q[0] *= f; Q[0] -= P[0]*(d *= f); R[0] *= f; R[0] -= P[0]*(r *= f);
			P[1] *= f;          Q[1] *= f; Q[1] -= P[1]*d;        R[1] *= f; R[1] -= P[1]*r;
			P[2] *= f;          Q[2] *= f; Q[2] -= P[2]*d;        R[2] *= f; R[2] -= P[2]*r;
		}
		else nurbs_P(N_u, D_u, N_v, D_v, p, q_u, q_v, icp, i-q_u, j-q_v, dim, P, Q, R);
	}
	return;
}

//////////////////////////////////////////////////////////////////////////////////////////////////
//...функция, вычисляющая значение двумерного сплайна в 3D и его вторые производные по параметрам;
void nurbs3D(double u, double * t_u, int n_u, int q_u,
             double v, double * t_v, int n_v, int q_v, double * p, int icp, int dim,
             double * P, double * Q, double * R, double * S, double * T, double * W)
{
	double N_u[MAX_SPLINE], D_u[MAX_SPLINE], K_u[MAX_SPLINE], N_v[MAX_SPLINE], D_v[MAX_SPLINE], K_v[MAX_SPLINE], 
		f, d, r, h, q, w;
	if (q_u <= MAX_SPLINE && q_v <= MAX_SPLINE) {
		int i, j;
		nurbs_spl(u, t_u, n_u, q_u, i, N_u, D_u, K_u); if (i < q_u) return;
		nurbs_spl(v, t_v, n_v, q_v, j, N_v, D_v, K_v); if (j < q_v) return;
		if (dim == 4) {
			nurbs_wP(N_u, D_u, K_u, N_v, D_v, K_v, p, q_u, q_v, icp, i-q_u, j-q_v, dim, P, Q, R, S, T, W);
			nurbs_w (N_u, D_u, K_u, N_v, D_v, K_v, p, q_u, q_v, icp, i-q_u, j-q_v, dim, f, d, r, h, q, w);
			P[0] *= (f = 1./f); Q[0] *= f; Q[0] -= P[0]*(d *= f); S[0] *= f; S[0] -= P[0]*(h *= f)+2.*Q[0]*d;
			P[1] *= f;          Q[1] *= f; Q[1] -= P[1]*d;        S[1] *= f; S[1] -= P[1]*h+2.*Q[1]*d;
			P[2] *= f;          Q[2] *= f; Q[2] -= P[2]*d;        S[2] *= f; S[2] -= P[2]*h+2.*Q[2]*d;
			R[0] *= f; R[0] -= P[0]*(r *= f); W[0] *= f; W[0] -= P[0]*(w *= f)+2.*R[0]*r;
			R[1] *= f; R[1] -= P[1]*r;        W[1] *= f; W[1] -= P[1]*w+2.*R[1]*r;
			R[2] *= f; R[2] -= P[2]*r;        W[2] *= f; W[2] -= P[2]*w+2.*R[2]*r;
			T[0] *= f; T[0] -= P[0]*(q *= f)+Q[0]*r+R[0]*d;
			T[1] *= f; T[1] -= P[1]*q+Q[1]*r+R[1]*d;
			T[2] *= f; T[2] -= P[2]*q+Q[2]*r+R[2]*d;
		}
		else nurbs_P(N_u, D_u, K_u, N_v, D_v, K_v, p, q_u, q_v, icp, i-q_u, j-q_v, dim, P, Q, R, S, T, W);
	}
	return;
}

//////////////////////////////////////////////////////
//...функция, вычисляющая сплайновую поверхность в 3D;
void nurbs3D(double u, double v, CMap * mp, double * P, int id_map)
{
	if (mp   && ID_MAP(2, NURBS_GENUS)  == mp[0] &&
		mp[size_of_map(2, NURBS_GENUS)] == (CMap)B2SPLINE_CELL && P) {
		int m1 = (int)(mp[8]-mp[10]), m2 = (int)mp[8],
			 m3 = (int)(mp[9]-mp[11]), m4 = (int)mp[9], m5 = (int)mp[7]/4;
		if (! id_map) nurbs3D(u, mp+13, m2, m1, v, mp+m2+15, m4, m3, mp+m2+m4+17, m2-m1, m5, P);
		else {
			memset(P,   0x0000, (m5-1)*sizeof(double));
			if (mp[m2+14] > mp[m2+13]+EE_ker && mp[m2+m4+16] > mp[m2+m4+15]+EE_ker)
			nurbs3D(mp[m2+13]+u*(mp[m2+14]-mp[m2+13]), mp+13, m2, m1, 
					mp[m2+m4+15]+v*(mp[m2+m4+16]-mp[m2+m4+15]),
					mp+m2+15, m4, m3, mp+m2+m4+17, m2-m1, m5, P);
		}
	}
	return;
}

/////////////////////////////////////////////////////////////////////
//...функция, вычисляющая сплайновую поверхность в 3D, а также репер;
void nurbs3D(double u, double v, CMap * mp, double * P, double * ort, int id_map)
{
	if (mp  && ID_MAP(2, NURBS_GENUS) == mp[0] &&
		mp[size_of_map(2, NURBS_GENUS)] == (CMap)B2SPLINE_CELL && P && ort) {
		int m1 = (int)(mp[8]-mp[10]), m2 = (int)mp[8],
			 m3 = (int)(mp[9]-mp[11]), m4 = (int)mp[9], m5 = (int)mp[7]/4;
		if (! id_map) nurbs3D(u, mp+13, m2, m1, v, mp+m2+15, m4, m3, mp+m2+m4+17, m2-m1, m5, P, ort+6, ort+3);
		else {
			memset(P,   0x0000, (m5-1)*sizeof(double));
			memset(ort, 0x0000,      9*sizeof(double)); ort[8] = ort[4] = 1.;
			if (mp[m2+14] > mp[m2+13]+EE_ker && mp[m2+m4+16] > mp[m2+m4+15]+EE_ker) {
				double f, d;
				nurbs3D(mp[m2+13]+u*(f = mp[m2+14]-mp[m2+13]), mp+13, m2, m1,
						mp[m2+m4+15]+v*(d = mp[m2+m4+16]-mp[m2+m4+15]),
						mp+m2+15, m4, m3, mp+m2+m4+17, m2-m1, m5, P, ort+6, ort+3);
				ort[6] *= f; ort[7] *= f; ort[8] *= f;
				ort[3] *= d; ort[4] *= d; ort[5] *= d;
			}
		}
		double f =  sqrt(ort[6]*ort[6]+ort[7]*ort[7]+ort[8]*ort[8]),
				 d =  sqrt(ort[3]*ort[3]+ort[4]*ort[4]+ort[5]*ort[5]);
		ort[6] *= (f = f > EE ? 1./f : 1.); ort[3] *= (d = d > EE ? 1./d : 1.);
		ort[7] *=  f;                       ort[4] *=  d;
		ort[8] *=  f;                       ort[5] *=  d;
		ort[0] = ort[5]*ort[7]-ort[4]*ort[8];
		ort[1] = ort[3]*ort[8]-ort[5]*ort[6];
		ort[2] = ort[4]*ort[6]-ort[3]*ort[7];
	}
	return;
}

////////////////////////////////////////////////////////////////////////////////////////////////
//...функция, вычисляющая сплайновую поверхность в 3D, а также репер и обратный радиус кривизны;
void nurbs3D(double u, double v, CMap * mp, double * P, double * ort, double * R_t, int id_map)
{
	if (mp  && ID_MAP(2, NURBS_GENUS) == mp[0] &&
		mp[size_of_map(2, NURBS_GENUS)] == (CMap)B2SPLINE_CELL && P && ort && R_t) {
		int m1 = (int)(mp[8]-mp[10]), m2 = (int)mp[8],
			 m3 = (int)(mp[9]-mp[11]), m4 = (int)mp[9], m5 = (int)mp[7]/4;
		if (! id_map) nurbs3D(u, mp+13, m2, m1, v, mp+m2+15, m4, m3, mp+m2+m4+17, m2-m1, m5, P, ort+6, ort+3, R_t, R_t+3, R_t+6);
		else {
			memset(P,   0x0000, (m5-1)*sizeof(double));
			memset(ort, 0x0000,      9*sizeof(double)); ort[8] = ort[4] = 1.;
			memset(R_t, 0x0000,      9*sizeof(double));
			if (mp[m2+14] > mp[m2+13]+EE_ker && mp[m2+m4+16] > mp[m2+m4+15]+EE_ker) {
				double f, d, h; 
				nurbs3D(mp[m2+13]+u*(f = mp[m2+14]-mp[m2+13]), mp+13, m2, m1,
						mp[m2+m4+15]+v*(d = mp[m2+m4+16]-mp[m2+m4+15]),
						mp+m2+15, m4, m3, mp+m2+m4+17, m2-m1, m5, P, ort+6, ort+3, R_t, R_t+3, R_t+6);
				ort[6] *= f; ort[7] *= f; ort[8] *= f;
				ort[3] *= d; ort[4] *= d; ort[5] *= d; h = f*d; f *= f; d *= d;
				R_t[0] *= f; R_t[1] *= f; R_t[2] *= f;
				R_t[3] *= h; R_t[4] *= h; R_t[5] *= h;
				R_t[6] *= d; R_t[7] *= d; R_t[8] *= d;
			}
		}
		double f = sqrt(ort[6]*ort[6]+ort[7]*ort[7]+ort[8]*ort[8]),
				 d = sqrt(ort[3]*ort[3]+ort[4]*ort[4]+ort[5]*ort[5]),
				 h = sqrt(R_t[0]*R_t[0]+R_t[1]*R_t[1]+R_t[2]*R_t[2]), 
				 q = sqrt(R_t[6]*R_t[6]+R_t[7]*R_t[7]+R_t[8]*R_t[8]), 
				 r =		 ort[6]*R_t[0]+ort[7]*R_t[1]+ort[8]*R_t[2],
				 w =		 ort[3]*R_t[6]+ort[4]*R_t[7]+ort[5]*R_t[8];
		if (f > EE) R_t[0] = sqrt(sqr(f)*sqr(h)-r*r)/(f*f*f); else R_t[0] = 0.;
		if (d > EE) R_t[1] = sqrt(sqr(d)*sqr(q)-w*w)/(d*d*d); else R_t[1] = 0.;
		ort[6] *= (f = f > EE ? 1./f : 1.); ort[3] *= (d = d > EE ? 1./d : 1.);
		ort[7] *=  f;                       ort[4] *=  d; R_t[2] = f;
		ort[8] *=  f;                       ort[5] *=  d; R_t[3] = d;
		ort[0] = ort[5]*ort[7]-ort[4]*ort[8];
		ort[1] = ort[3]*ort[8]-ort[5]*ort[6];
		ort[2] = ort[4]*ort[6]-ort[3]*ort[7];
	}
	return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//...функция, вычисляющая значения обобщенного сплайна (включая трактовки для простейших элементов)
void Nurbs3D(double u, double v, CMap * mp, double * P, int id_map)
{
	if (mp) {
		switch ((int)mp[0]) {
		case    (int)ID_MAP(2, SPHERE_GENUS):
					P[0] = P[1] = 0.;
					P[2] = fabs(mp[7]);
					point_iso  (P, NULL, u, M_PI_2-v);
		break;
		case    (int)ID_MAP(2,  NURBS_GENUS): nurbs3D(u, v, mp, P, id_map); break;
		}
	}
	return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//...вычисление обобщенного сплайна и его репера (включая трактовки для простейших элементов)
void Nurbs3D(double u, double v, CMap * mp, double * P, double * ort, int id_map)
{
	double CZ, SZ, CY, SY;
	if (mp) {
		switch ((int)mp[0]) {
		case    (int)ID_MAP(2, SPHERE_GENUS):
					P[0]   = P[1] = 0.;
					P[2]   = fabs(mp[7]);
					ort[0] = ort[1] = ort[3] = ort[5] = ort[7] = ort[8] = 0.;
					ort[2] = ort[4] = ort[6] = 1.;
					point_iso<double>(P, NULL, CZ = cos(u), SZ = sin(u),
										   CY = sin(v), SY = cos(v), 1., 0.);
					point_iso<double>(ort,   NULL, CZ, SZ, CY, SY, 1., 0.);
					point_iso<double>(ort+3, NULL, CZ, SZ, CY, SY, 1., 0.);
					point_iso<double>(ort+6, NULL, CZ, SZ, CY, SY, 1., 0.);
		break;
		case    (int)ID_MAP(2,  NURBS_GENUS): nurbs3D(u, v, mp, P, ort, id_map); break;
		}
	}
	return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//...вспомогательная процедура - извлечение параметрического описания из ячейки (если оно есть)
CMap * get_pm_2D(Cells * ce, int Max_N)
{
	CMap * pm = ce && ce->pm && ce->graph && Max_N > 0 ?
		  get_map(1, NURBS_GENUS, ce->graph[1]*(Max_N+1)*2) : NULL;
	if (pm) {
		double P[3], f = 1./Max_N, t;
		int    i, j, k, l, m = size_of_map(1, NURBS_GENUS);
		pm[i  = m] = (CMap)B1POLY_CELL;
		pm[7] = (CMap)4;
		pm[8] = (CMap)(ce->graph[1]*(Max_N+1)+2);
		pm[9] = (CMap)(ce->graph[1]*(Max_N+1));

//////////////////////////////////////////////////////////////////////////////////////////
//...заполняем параметрическое описание с учетом типа кривой, входящей в описание контура;
		for (l = 0; l <  ce->graph[1]; l++)
		if  (ID_MAP(1, NURBS_GENUS) == ce->pm[l][0])
		switch ((int)ce->pm[l][m]) {
			  case (int)B1SPLINE_CELL :
			  if (fabs(ce->pm[l][(int)ce->pm[l][8]+12]-ce->pm[l][(int)ce->pm[l][8]+11]) > EE_ker)
			  for (t = 0., k = 0; k <= Max_N;   k++) {
					 nurbs2D(t, ce->pm[l], P, 1); t += f;
					 pm[++i] = P[0];
					 pm[++i] = P[1];
			  }
			  else {
					 pm[8] -= (CMap)(Max_N+1);
					 pm[9] -= (CMap)(Max_N+1);
			  }    break;
			  case (int)B1POLY_CELL : {
					  add_new_maps(pm, (int)pm[9]*2+m, (k = (j = (int)ce->pm[l][9])-(Max_N+1))*2);
					  pm[8] += (CMap)k;
					  pm[9] += (CMap)k;
					  for (k = m;  j > 0; j--) {
							 pm[++i] = ce->pm[l][k++];
							 pm[++i] = ce->pm[l][k++];
					  }
			  }     break;
		}

/////////////////////////////////////////////////////////////////
//...для надежности соединяем начало и конец многоугольной линии;
		pm[i-1] = pm[k = m+1];
		pm[i]   = pm[k+1];
	}
	return pm;
}

/////////////////////////////////////////////////////////////////////////
//...преобразование параметрического описания к области параметров сферы;
void pm_sph_convert(CMap * pm, CMap * mp, CMap * sph)
{
  int k, i;
  if (pm  &&      ID_MAP(1,  NURBS_GENUS)  == pm[0]              &&
      pm[i = size_of_map(1,  NURBS_GENUS)] == (CMap)B1POLY_CELL   &&
      mp  &&      ID_MAP(2,  NURBS_GENUS)  == mp[0]              &&
      mp[    size_of_map(2,  NURBS_GENUS)] == (CMap)B2SPLINE_CELL &&
      sph &&      ID_MAP(2, SPHERE_GENUS)  == sph[0]) {
      double P[3] = { 0., 0., 0.}, u, v, f = 2.*M_PI;

////////////////////////////////////////////////////
//...пересчитываем точки в область параметров сферы;
      for (k = (int)pm[9]; k > 0; k--) {
           u = pm[++i];
           v = pm[++i];

           nurbs3D(u, v, mp, P, 0);
           make_local(P, sph);

           u = pm[i-1] = arg0(comp(P[0], P[1]));
           v = pm[i]   = M_PI_2-arg0(comp(P[2], sqrt(P[0]*P[0]+P[1]*P[1])));
      }

#ifdef ___PARAMETRIZATION_CONSTRUCTION___
//////////////////////////////////////////////////////////////////
//...корректируем точки, приводя их к одному листу многозначности;
      u0 = pm[++(i = size_of_map(1,  NURBS_GENUS))];
      v0 = pm[++i];
      for (k = (int)pm[9]-1; k > 0; k--) {
           u = pm[++i];
           v = pm[++i];

/////////////////////////////////////////////////
//...нормируем точки и запоминаем новое значение;
           while (u-u0 >  f-EE_dop) u -= f;
           while (u-u0 < -f+EE_dop) u += f;
           while (v-v0 >  M_PI-EE_dop) v -= M_PI;
           while (v-v0 < -M_PI+EE_dop) v += M_PI; 

           pm[i-1] = u0 = u; 
           pm[i]   = v0 = v;
      }
#endif
  }
  return;
}

///////////////////////////////////////////////////////////////
//...функция движения по трассе на сплайновой поверхности в 3D;
void Trassa3D(double & u, double & v, CMap * mp, double * P, double * ort, int MAX_NEWTON)
{
	if (mp  && ID_MAP(2, NURBS_GENUS) == mp[0] &&
		mp[size_of_map(2, NURBS_GENUS)] == (CMap)B2SPLINE_CELL && P && ort) {
		int m1 = (int)(mp[8]-mp[10]), m2 = (int)mp[8],
			 m3 = (int)(mp[9]-mp[11]), m4 = (int)mp[9], m5 = (int)mp[7]/4;
		double R_t[9], P_temp[3], f, d, h, q, r, w, alpha, betta, J;
		if (mp[m2+14] > mp[m2+13]+EE_ker && mp[m2+m4+16] > mp[m2+m4+15]+EE_ker) {
/////////////////////////////////
//...сажаем точку на поверхность;
			for (int k = 0; k < MAX_NEWTON; k++) {
				nurbs3D(mp[m2+13]+u*(f = mp[m2+14]-mp[m2+13]), mp+13, m2, m1,
						mp[m2+m4+15]+v*(d = mp[m2+m4+16]-mp[m2+m4+15]),
						mp+m2+15, m4, m3, mp+m2+m4+17, m2-m1, m5, P_temp, ort+3, ort+6);
				ort[3] *= f; ort[4] *= f; ort[5] *= f;
				ort[6] *= d; ort[7] *= d; ort[8] *= d;
				f = ort[3]*ort[3]+ort[4]*ort[4]+ort[5]*ort[5];
				d = ort[6]*ort[6]+ort[7]*ort[7]+ort[8]*ort[8];
				//if (f > EE && d > EE) {
				//	u += ((P[0]-P_temp[0])*ort[3]+(P[1]-P_temp[1])*ort[4]+(P[2]-P_temp[2])*ort[5])/f; u = max(min(u, 1.), 0.);
				//	v += ((P[0]-P_temp[0])*ort[6]+(P[1]-P_temp[1])*ort[7]+(P[2]-P_temp[2])*ort[8])/d; v = max(min(v, 1.), 0.);
				//}
////////////////////////////////////
//...точный алгоритм метода Ньютона;
				J = f*d-sqr(h = ort[3]*ort[6]+ort[4]*ort[7]+ort[5]*ort[8]); f /= J; d /= J; h /= J;
				alpha = (P[0]-P_temp[0])*ort[3]+(P[1]-P_temp[1])*ort[4]+(P[2]-P_temp[2])*ort[5];
				betta = (P[0]-P_temp[0])*ort[6]+(P[1]-P_temp[1])*ort[7]+(P[2]-P_temp[2])*ort[8];
				u += alpha*d-betta*h; u = max(min(u, 1.), 0.);
				v += betta*f-alpha*h; v = max(min(v, 1.), 0.);
			}
///////////////////////////////////////
//...вычисляем скорректированную точку;
			nurbs3D(mp[m2+13]+u*(f = mp[m2+14]-mp[m2+13]), mp+13, m2, m1,
					mp[m2+m4+15]+v*(d = mp[m2+m4+16]-mp[m2+m4+15]),
					mp+m2+15, m4, m3, mp+m2+m4+17, m2-m1, m5, P, ort+3, ort+6, R_t, R_t+3, R_t+6);
			ort[3] *= f; ort[4] *= f; ort[5] *= f;
			ort[6] *= d; ort[7] *= d; ort[8] *= d; h = f*d; f *= f; d *= d;
			R_t[0] *= f; R_t[1] *= f; R_t[2] *= f;
			R_t[3] *= h; R_t[4] *= h; R_t[5] *= h;
			R_t[6] *= d; R_t[7] *= d; R_t[8] *= d;
///////////////////////
//...вычисляем нормаль;
			if ((f = ort[3]*ort[3]+ort[4]*ort[4]+ort[5]*ort[5]) > EE) h = 1./sqrt(f); else h = 1.;
			if ((d = ort[6]*ort[6]+ort[7]*ort[7]+ort[8]*ort[8]) > EE) h = h /sqrt(d);
			q = (ort[4]*ort[8]-ort[5]*ort[7])*h;
			r = (ort[5]*ort[6]-ort[3]*ort[8])*h;
			w = (ort[3]*ort[7]-ort[4]*ort[6])*h;
			if ((h = q*q+r*r+w*w) > EE) {
				q *= (h = 1./sqrt(h)); r *= h; w *= h;
			}
			else {
				r = 1.; q = w = 0.;
			}
//////////////////////////////////////////////////////////////////////////////////////////
//...вычисляем производные первого и второго порядка вдоль скорректированного направления;
			J = f*d-sqr(h = ort[3]*ort[6]+ort[4]*ort[7]+ort[5]*ort[8]); f /= J; d /= J; h /= J;
			alpha  = ort[3]*ort[0]+ort[4]*ort[1]+ort[5]*ort[2]; J = alpha;
			betta  = ort[6]*ort[0]+ort[7]*ort[1]+ort[8]*ort[2];
			alpha  = J*d-betta*h;
			betta  = betta*f-J*h;
			if ((J = sqrt(sqr(alpha)+sqr(betta))) > EE) {
				alpha *= (J = 1./J); betta *= J;
			}
			else {
				alpha = 1.; betta = 0.;
			}
			ort[6] = alpha*ort[3]+betta*ort[6]; ort[3] = -q;
			ort[7] = alpha*ort[4]+betta*ort[7]; ort[4] = -r;
			ort[8] = alpha*ort[5]+betta*ort[8]; ort[5] = -w;
			R_t[0] = sqr(alpha)*R_t[0]+sqr(betta)*R_t[6]+2.*alpha*betta*R_t[3];
			R_t[1] = sqr(alpha)*R_t[1]+sqr(betta)*R_t[7]+2.*alpha*betta*R_t[4];
			R_t[2] = sqr(alpha)*R_t[2]+sqr(betta)*R_t[8]+2.*alpha*betta*R_t[5];
/////////////////////////////////////////////////////////////////////////////////////////////////////
//...вычисляем касательную вдоль траектории и нормаль кривизны (с нормой обратного радиуса кривизны);
			f = ort[6]*ort[6]+ort[7]*ort[7]+ort[8]*ort[8];
			if (f > EE) {
				ort[0] = ort[6]*(f = 1./sqrt(f)); ort[1] = ort[7]*f; ort[2] = ort[8]*f;
			} 
			else {
				ort[0] = (f = 1.); ort[1] = ort[2] = 0.;
			}
			R_t[3] = (R_t[1]*ort[2]-R_t[2]*ort[1])*sqr(f);
			R_t[4] = (R_t[2]*ort[0]-R_t[0]*ort[2])*sqr(f);
			R_t[5] = (R_t[0]*ort[1]-R_t[1]*ort[0])*sqr(f);
			ort[6] = R_t[4]*ort[2]-R_t[5]*ort[1];
			ort[7] = R_t[5]*ort[0]-R_t[3]*ort[2];
			ort[8] = R_t[3]*ort[1]-R_t[4]*ort[0];
		}
	}
	return;
}


/////////////////////////////////////////////
//...функция движения по скелету трассы в 3D;
void Trassa3D_skel(double & u, CMap * mp, double * P, double * ort, int MAX_NEWTON)
{
	if (mp  && ID_MAP(1, NURBS_GENUS) == mp[0] &&
		mp[size_of_map(1, NURBS_GENUS)] == (CMap)B1SPLINE_CELL && P && ort) {
		int m1 = (int)(mp[8]-mp[9]), m2 = (int)mp[9], m3 = (int)mp[7]/2;
		double P_temp[3], f, h, alpha;
		if (mp[m1+m2+12] > mp[m1+m2+11]+EE_ker) {
////////////////////////////
//...сажаем точку на кривую;
			for (int k = 0; k < MAX_NEWTON; k++) {
				nurbs3D(mp[m1+m2+11]+u*(f = mp[m1+m2+12]-mp[m1+m2+11]), 
									mp+11, mp+m1+m2+13, m1+m2, m1, m3, P_temp, ort+3);
				ort[3] *= f; ort[4] *= f; ort[5] *= f;
				alpha = (P[0]-P_temp[0])*ort[3]+(P[1]-P_temp[1])*ort[4]+(P[2]-P_temp[2])*ort[5];
				u += alpha/(ort[3]*ort[3]+ort[4]*ort[4]+ort[5]*ort[5]); u = max(min(u, 1.), 0.);
			}
///////////////////////////////////////
//...вычисляем скорректированную точку;
			nurbs3D(mp[m1+m2+11]+u*(f = mp[m1+m2+12]-mp[m1+m2+11]), 
								mp+11, mp+m1+m2+13, m1+m2, m1, m3, P, ort+3, ort+6);
			ort[3] *= f; ort[4] *= f; ort[5] *= f; h = f*f;
			ort[6] *= h; ort[7] *= h; ort[8] *= h;
/////////////////////////////////////////////////////////////////////////////////////////////////////
//...вычисляем касательную вдоль траектории и нормаль кривизны (с нормой обратного радиуса кривизны);
			f = ort[3]*ort[3]+ort[4]*ort[4]+ort[5]*ort[5];
			if (f > EE) {
				ort[0] = ort[3]*(f = 1./sqrt(f)); ort[1] = ort[4]*f; ort[2] = ort[5]*f;
			} 
			else {
				ort[0] = (f = 1.); ort[1] = ort[2] = 0.;
			}
			ort[3] = (ort[7]*ort[2]-ort[8]*ort[1])*sqr(f);
			ort[4] = (ort[8]*ort[0]-ort[6]*ort[2])*sqr(f);
			ort[5] = (ort[6]*ort[1]-ort[7]*ort[0])*sqr(f);
			ort[6] = ort[4]*ort[2]-ort[5]*ort[1];
			ort[7] = ort[5]*ort[0]-ort[3]*ort[2];
			ort[8] = ort[3]*ort[1]-ort[4]*ort[0];
		}
/////////////////////////////////////////
//...вычисляем вспомогательную конормаль;
			if ((f = ort[6]*ort[6]+ort[7]*ort[7]+ort[8]*ort[8]) > EE) h = 1./sqrt(f); else h = 1.;
			ort[3] = (ort[1]*ort[8]-ort[2]*ort[7])*h;
			ort[4] = (ort[2]*ort[6]-ort[0]*ort[8])*h;
			ort[5] = (ort[0]*ort[7]-ort[1]*ort[6])*h;
	}
	return;
}

//////////////////////////////////////////////////////////////////
//...функция, вычисляющая параметры заданной точки на поверхности;
double nurbs3D_u_par(double v, CMap * mp, double * P, double eps, int M_iter, int N_ini, double fnorm)
{
	double u_par = -1.;
	if (mp  && ID_MAP(2, NURBS_GENUS) == mp[0] &&
		mp[size_of_map(2, NURBS_GENUS)] == (CMap)B2SPLINE_CELL && P) {
		int m1 = (int)(mp[8]-mp[10]), m2 = (int)mp[8],
			 m3 = (int)(mp[9]-mp[11]), m4 = (int)mp[9], m5 = (int)mp[7]/4, k;
		double P_temp[6], ort[6], f, d, q, r, u, u0, h = 1./N_ini, w = h/N_ini;
		if (mp[m2+14] > mp[m2+13]+EE_ker && mp[m2+m4+16] > mp[m2+m4+15]+EE_ker) {
////////////////////////////////////////////////////
//...ищем ближайшую точку на параметрической кривой;
			u = 0.;
			nurbs3D(mp[m2+13]+u*(mp[m2+14]-mp[m2+13]), mp+13, m2, m1,
					mp[m2+m4+15]+v*(mp[m2+m4+16]-mp[m2+m4+15]),
					mp+m2+15, m4, m3, mp+m2+m4+17, m2-m1, m5, P_temp); 
			d = abs_point(P, P_temp); set_point(P_temp+3, P_temp); r = 0.;
			for (k = 0; k < N_ini; k++) {
				nurbs3D(mp[m2+13]+(u += h)*(mp[m2+14]-mp[m2+13]), mp+13, m2, m1,
						mp[m2+m4+15]+v*(mp[m2+m4+16]-mp[m2+m4+15]),
						mp+m2+15, m4, m3, mp+m2+m4+17, m2-m1, m5, P_temp); 
				q = fnorm*abs_point(P_temp, P_temp+3);	 
				if (! k && d < q) {
					r = q; u_par = u-h; 
				}
				if ((f = abs_point(P, P_temp)) < d && f < q) {
					d = f; r = q; u_par = u; 
				}
				set_point(P_temp+3, P_temp);
			}
			if (d > r) return(u_par = -1);
/////////////////////////////////////////////
//...ищем окрестность для алгоритма бисекции;
			u0 = u_par; do {
				if ((u0 -= w) < 0.) u0 = 0.; 
				nurbs3D(mp[m2+13]+u0*(f = mp[m2+14]-mp[m2+13]), mp+13, m2, m1,
						mp[m2+m4+15]+v*(d = mp[m2+m4+16]-mp[m2+m4+15]),
						mp+m2+15, m4, m3, mp+m2+m4+17, m2-m1, m5, P_temp, ort, ort+3);
				ort[0] *= f; ort[1] *= f; ort[2] *= f;
				ort[3] *= d; ort[4] *= d; ort[5] *= d;
				d = (P[0]-P_temp[0])*ort[0]+(P[1]-P_temp[1])*ort[1]+(P[2]-P_temp[2])*ort[2];
			}
			while (u0 > 0. && d < 0.); 
			u = u_par; do {
				if ((u += w) > 1.) u = 1.; 
				nurbs3D(mp[m2+13]+u*(f = mp[m2+14]-mp[m2+13]), mp+13, m2, m1,
						mp[m2+m4+15]+v*(d = mp[m2+m4+16]-mp[m2+m4+15]),
						mp+m2+15, m4, m3, mp+m2+m4+17, m2-m1, m5, P_temp, ort, ort+3);
				ort[0] *= f; ort[1] *= f; ort[2] *= f;
				ort[3] *= d; ort[4] *= d; ort[5] *= d;
				d = (P[0]-P_temp[0])*ort[0]+(P[1]-P_temp[1])*ort[1]+(P[2]-P_temp[2])*ort[2];
			}
			while (u < 1. && d > 0.); 
///////////////////////
//...алгоритм бисекции;
			k = 0; do {
				u_par = (u+u0)*.5; 
				nurbs3D(mp[m2+13]+u_par*(f = mp[m2+14]-mp[m2+13]), mp+13, m2, m1,
						mp[m2+m4+15]+v*(d = mp[m2+m4+16]-mp[m2+m4+15]),
						mp+m2+15, m4, m3, mp+m2+m4+17, m2-m1, m5, P_temp, ort, ort+3);
				ort[0] *= f; ort[1] *= f; ort[2] *= f;
				ort[3] *= d; ort[4] *= d; ort[5] *= d;
				d = (P[0]-P_temp[0])*ort[0]+(P[1]-P_temp[1])*ort[1]+(P[2]-P_temp[2])*ort[2];
				if (d < 0.) u  = u_par; else 
				if (d > 0.) u0 = u_par; else u = u0 = u_par;
			}
			while (++k < M_iter && fabs(u-u0) >= EE); 
			if ((d = abs_point(P, P_temp)) < eps) 
				  u_par = (u+u0)*.5;
			else u_par = -1.;
		}
	}
	return u_par;
}

double nurbs3D_v_par(double u, CMap * mp, double * P, double eps, int M_iter, int N_ini, double fnorm)
{
	double v_par = -1.;
	if (mp  && ID_MAP(2, NURBS_GENUS) == mp[0] &&
		mp[size_of_map(2, NURBS_GENUS)] == (CMap)B2SPLINE_CELL && P) {
		int m1 = (int)(mp[8]-mp[10]), m2 = (int)mp[8],
			 m3 = (int)(mp[9]-mp[11]), m4 = (int)mp[9], m5 = (int)mp[7]/4, k;
		double P_temp[6], ort[6], f, d, q, r, v, v0, h = 1./N_ini, w = h/N_ini;
		if (mp[m2+14] > mp[m2+13]+EE_ker && mp[m2+m4+16] > mp[m2+m4+15]+EE_ker) {
////////////////////////////////////////////////////
//...ищем ближайшую точку на параметрической кривой;
			v = 0.;
			nurbs3D(mp[m2+13]+u*(mp[m2+14]-mp[m2+13]), mp+13, m2, m1,
					mp[m2+m4+15]+v*(mp[m2+m4+16]-mp[m2+m4+15]),
					mp+m2+15, m4, m3, mp+m2+m4+17, m2-m1, m5, P_temp); 
			d = abs_point(P, P_temp); set_point(P_temp+3, P_temp); r = 0.;
			for (k = 0; k < N_ini; k++) {
				nurbs3D(mp[m2+13]+u*(mp[m2+14]-mp[m2+13]), mp+13, m2, m1,
						mp[m2+m4+15]+(v += h)*(mp[m2+m4+16]-mp[m2+m4+15]),
						mp+m2+15, m4, m3, mp+m2+m4+17, m2-m1, m5, P_temp); 
				q = fnorm*abs_point(P_temp, P_temp+3);
				if (! k && d < q) {
					r = q; v_par = v-h; 
				}
				if ((f = abs_point(P, P_temp)) < d && f < q) {
					d = f; r = q; v_par = v; 
				}
				set_point(P_temp+3, P_temp); 
			}
			if (d > r) return(v_par = -1);
/////////////////////////////////////////////
//...ищем окрестность для алгоритма бисекции;
			v0 = v_par; do {
				if ((v0 -= w) < 0.) v0 = 0.; 
				nurbs3D(mp[m2+13]+u*(f = mp[m2+14]-mp[m2+13]), mp+13, m2, m1,
						mp[m2+m4+15]+v0*(d = mp[m2+m4+16]-mp[m2+m4+15]),
						mp+m2+15, m4, m3, mp+m2+m4+17, m2-m1, m5, P_temp, ort+3, ort);
				ort[3] *= f; ort[4] *= f; ort[5] *= f;
				ort[0] *= d; ort[1] *= d; ort[2] *= d;
				d = (P[0]-P_temp[0])*ort[0]+(P[1]-P_temp[1])*ort[1]+(P[2]-P_temp[2])*ort[2];
			}
			while (v0 > 0. && d < 0.); 
			v = v_par; do {
				if ((v += w) > 1.) v = 1.; 
				nurbs3D(mp[m2+13]+u*(f = mp[m2+14]-mp[m2+13]), mp+13, m2, m1,
						mp[m2+m4+15]+v*(d = mp[m2+m4+16]-mp[m2+m4+15]),
						mp+m2+15, m4, m3, mp+m2+m4+17, m2-m1, m5, P_temp, ort+3, ort);
				ort[3] *= f; ort[4] *= f; ort[5] *= f;
				ort[0] *= d; ort[1] *= d; ort[2] *= d;
				d = (P[0]-P_temp[0])*ort[0]+(P[1]-P_temp[1])*ort[1]+(P[2]-P_temp[2])*ort[2];
			}
			while (v < 1. && d > 0.); 
///////////////////////
//...алгоритм бисекции;
			k = 0; do {
				v_par = (v+v0)*.5; 
				nurbs3D(mp[m2+13]+u*(f = mp[m2+14]-mp[m2+13]), mp+13, m2, m1,
						mp[m2+m4+15]+v_par*(d = mp[m2+m4+16]-mp[m2+m4+15]),
						mp+m2+15, m4, m3, mp+m2+m4+17, m2-m1, m5, P_temp, ort+3, ort);
				ort[3] *= f; ort[4] *= f; ort[5] *= f;
				ort[0] *= d; ort[1] *= d; ort[2] *= d;
				d = (P[0]-P_temp[0])*ort[0]+(P[1]-P_temp[1])*ort[1]+(P[2]-P_temp[2])*ort[2];
				if (d < 0.) v  = v_par; else 
				if (d > 0.) v0 = v_par; else v = v0 = v_par;
			}
			while (++k < M_iter && fabs(v-v0) >= EE); 
			if ((d = abs_point(P, P_temp)) < eps) 
				  v_par = (v+v0)*.5;
			else v_par = -1.;
		}
	}
	return v_par;
}

/////////////////////////////////////////////////////
//...зачитывание обобщенного сплайна из формата IGES;
CMap * Nurbs3D_IGES(char * ch_IGES, int IGES_cod, unsigned long id_num)
{
	const int STR_SIZE = 80;
	CMap * mp = NULL;
	char * id_IGES = read_struct_ascii(ch_IGES), temp = 0, temp_dop = 0;
	if (id_IGES) {
		char one_line[STR_SIZE+1], * pchar;	int cod = 0;
		unsigned long ppos_cur = 0, upper_limit, count_cod = 0, beg_data = 0, weight_data = 0;
		user_Count (id_IGES, 0, upper_limit, '\x0');

////////////////////////////////////////////////////////////////////////////////
//...пропускаем стартовый и глобальных параметров разделы (раздел S и раздел G);
		do {
			ONE_LINE(id_IGES, pchar, ppos_cur, upper_limit, one_line, STR_SIZE);
		}
		while (ppos_cur < upper_limit && one_line[72] != 'D');

//////////////////////////////////////////////////////////////////////////////
//...ищем вхождение нужного элемента в разделе вхождений элементов (раздел D);
		if (ppos_cur < upper_limit && one_line[72] == 'D') {
			swap(one_line[8], temp); cod = atoi(one_line);
			swap(one_line[8], temp);
			if (cod == IGES_cod) {
				count_cod++;
				if (count_cod != id_num)
				ONE_LINE(id_IGES, pchar, ppos_cur, upper_limit, one_line, STR_SIZE);
			}
		}
		while (ppos_cur < upper_limit && one_line[72] == 'D' && count_cod != id_num) {
			ONE_LINE(id_IGES, pchar, ppos_cur, upper_limit, one_line, STR_SIZE);
			swap(one_line[8], temp); cod = atoi(one_line);
			swap(one_line[8], temp);
			if (cod == IGES_cod) {
				count_cod++;
				if (count_cod != id_num)
				ONE_LINE(id_IGES, pchar, ppos_cur, upper_limit, one_line, STR_SIZE);
			}
		}

///////////////////////////////////////////////////////////////////////////
//...определяем начало и конец зачитываемого элемента в разделе параметров;
		if (one_line[72] == 'D' && count_cod == id_num) {
			swap(one_line[16], temp); beg_data = atoi(one_line+8); 
			swap(one_line[16], temp);
			ONE_LINE(id_IGES, pchar, ppos_cur, upper_limit, one_line, STR_SIZE);
			swap(one_line[32], temp); weight_data = atoi(one_line+24); 
			swap(one_line[32], temp);
		}

/////////////////////////////////////////////////////
//...пропускаем остаток раздела вхождений (раздел D);
		do {
			ONE_LINE(id_IGES, pchar, ppos_cur, upper_limit, one_line, STR_SIZE);
		}
		while (ppos_cur < upper_limit && one_line[72] == 'D');

////////////////////////////////////////////////////////////////////////////////////////
//...настраиваемся на начало вхождения нужного элемента в разделе параметров (раздел P);
		while (ppos_cur < upper_limit && one_line[72] == 'P' && --beg_data != 0)
			ONE_LINE(id_IGES, pchar, ppos_cur, upper_limit, one_line, STR_SIZE);

////////////////////////////////////////////////////////////////////
//...зачитываем параметры элемента из раздела параметров (раздел P);
		if (one_line[72] == 'P' && count_cod == id_num && beg_data == 0) {
			int K1, K2 = 0, M1, M2 = 0, i, j;
			swap(one_line[64], temp); 
			cod = atoi(one_line); pchar = strstr(one_line, ",");
			K1 = atoi(pchar+1); pchar = strstr(pchar+1, ",");
			M1 = atoi(pchar+1); pchar = strstr(pchar+1, ",");
			if (cod == 128) {
				K2 = M1;
				M1 = atoi(pchar+1); pchar = strstr(pchar+1, ",");
				M2 = atoi(pchar+1); pchar = strstr(pchar+1, ",");
			}
			pchar = strstr(pchar+1, ",");
			pchar = strstr(pchar+1, ",");
			pchar = strstr(pchar+1, ",");
			pchar = strstr(pchar+1, ",");
			if (cod == 128) pchar = strstr(pchar+1, ",");
		
///////////////////////////////////////////////////////
//...распределяем соответствующую геометрическую карту;
			if (cod == 126) {
				mp = get_map(1, NURBS_GENUS);
				mp[7] = 6;
				mp[8] = K1+M1+2;
				mp[9] = K1+1;
				add_new_cell_maps(mp, B1SPLINE_CELL);
			}
			else if (cod == 128) {
				mp = get_map(2, NURBS_GENUS);
				mp[7] = 12;
				mp[8] = K1+M1+2;
				mp[9] = K2+M2+2;
				mp[10] = K1+1;
				mp[11] = K2+1;
				add_new_cell_maps(mp, B2SPLINE_CELL);
			}
		
//////////////////////////////////////////////////////////////
//...зачитываем все оставшиеся параметры геометрической карты;
			if (mp && pchar) {
				i = size_of_map(mp);
				for (j = K1+M1+1; pchar && j >= 0; j--) {
					mp[++i] = user_strtod(pchar+1);
					pchar = strstr(pchar+1, ",");
					if (! pchar || (! strstr(pchar+1, ",") && ! strstr(pchar+1, ";"))) {	
						ONE_LINE(id_IGES, pchar, ppos_cur, upper_limit, one_line, STR_SIZE);
						pchar = one_line-1;
					}
				}
				i += 2;
				if (cod == 128) {
					for (j = K2+M2+1; pchar && j >= 0; j--) {
						mp[++i] = user_strtod(pchar+1);
						pchar = strstr(pchar+1, ",");
						if (! pchar || (! strstr(pchar+1, ",") && ! strstr(pchar+1, ";"))) {	
							ONE_LINE(id_IGES, pchar, ppos_cur, upper_limit, one_line, STR_SIZE);
							pchar = one_line-1;
						}
					}
					i += 2;
				}

//////////////////////////////////////////////////////////////////////////
//...пропускаем веса рационального сплайна и зачитываем контрольные точки;
				if (cod == 126) {
					for (j = (K1+1); pchar && j > 0; j--) {
						pchar = strstr(pchar+1, ",");
						if (! pchar || (! strstr(pchar+1, ",") && ! strstr(pchar+1, ";"))) {	
							ONE_LINE(id_IGES, pchar, ppos_cur, upper_limit, one_line, STR_SIZE);
							pchar = one_line-1;
						}
					}
					for (j = 3*(K1+1); pchar && j > 0; j--) {
						mp[++i] = user_strtod(pchar+1);
						pchar = strstr(pchar+1, ",");
						if (! pchar || (! strstr(pchar+1, ",") && ! strstr(pchar+1, ";"))) {	
							ONE_LINE(id_IGES, pchar, ppos_cur, upper_limit, one_line, STR_SIZE);
							pchar = one_line-1;
						}
					}
					i = size_of_map(mp)+K1+M1+2;
					for (j = 2; pchar && j > 0; j--) { //...зачитываем интервал параметризации;
						mp[++i] = user_strtod(pchar+1);
						pchar = strstr(pchar+1, ",");
						if (! pchar || (! strstr(pchar+1, ",") && ! strstr(pchar+1, ";"))) {	
							ONE_LINE(id_IGES, pchar, ppos_cur, upper_limit, one_line, STR_SIZE);
							pchar = one_line-1;
						}
					}
				}
				if (cod == 128) {
					for (j = (K1+1)*(K2+1); pchar && j > 0; j--) {
						pchar = strstr(pchar+1, ",");
						if (! pchar || (! strstr(pchar+1, ",") && ! strstr(pchar+1, ";"))) {	
							ONE_LINE(id_IGES, pchar, ppos_cur, upper_limit, one_line, STR_SIZE);
							pchar = one_line-1;
						}
					}
					for (j = 3*(K1+1)*(K2+1); pchar && j > 0; j--) {
						mp[++i] = user_strtod(pchar+1);
						pchar = strstr(pchar+1, ",");
						if (! pchar || (! strstr(pchar+1, ",") && ! strstr(pchar+1, ";"))) {	
							ONE_LINE(id_IGES, pchar, ppos_cur, upper_limit, one_line, STR_SIZE);
							pchar = one_line-1;
						}
					}
					i = size_of_map(mp)+K1+M1+2;
					for (j = 2; pchar && j > 0; j--) { //...зачитываем интервал параметризации;
						mp[++i] = user_strtod(pchar+1);
						pchar = strstr(pchar+1, ",");
						if (! pchar || (! strstr(pchar+1, ",") && ! strstr(pchar+1, ";"))) {	
							ONE_LINE(id_IGES, pchar, ppos_cur, upper_limit, one_line, STR_SIZE);
							pchar = one_line-1;
						}
					}
					i += K2+M2+2;
					for (j = 2; pchar && j > 0; j--) { //...зачитываем второй интервал параметризации;
						mp[++i] = user_strtod(pchar+1);
						pchar = strstr(pchar+1, ",");
						if (! pchar || (! strstr(pchar+1, ",") && ! strstr(pchar+1, ";"))) {	
							ONE_LINE(id_IGES, pchar, ppos_cur, upper_limit, one_line, STR_SIZE);
							pchar = one_line-1;
						}
					}
				}
			}
		}
	}
	delete_struct(id_IGES);
	return mp;
}

////////////////////////////////////////////////////////////
//...запись на диск и чтение в двоичном формате всей трассы;
void Trassa_write(char * ch_TRASSA, CMap ** trassa, int max_trassa)
{
	FILE * device = fopen(ch_TRASSA, "w+b");
	if (device) {
		fwrite(& max_trassa, sizeof(int), 1, device);
		for (int j = 0; j < max_trassa; j++) {
			int m = size_of_map(trassa[j])+size_of_dop(trassa[j])+1;
			fwrite(& m, sizeof(int), 1, device);
			fwrite(trassa[j], sizeof(CMap), m, device);
		}
		fclose (device);
	}
}

int Trassa_read(char * ch_TRASSA, CMap ** trassa, int MAX_TRASSA, int * max_trassa)
{
	FILE * device = fopen(ch_TRASSA, "r+b"); int m, max_common = 0, k = 0;
	if (device) {
      size_t res = fread(& k, sizeof(int), 1, device); if (k > MAX_TRASSA) k = MAX_TRASSA;
		if (max_trassa) max_trassa[0] = k; max_common += 1;
		for (int j = 0; j < k; j++) {
			res = fread(& m, sizeof(int), 1, device);
			trassa[j] = (CMap *)new_struct(m*sizeof(CMap));
			res = fread(trassa[j], sizeof(CMap), m, device);
			max_common += m+1;
		}
		fclose (device);
	}
	return(max_common);
}

int Trassa_convert(char * ch_TRASSA, char * ch_TRASSA_OUT, CMap ** trassa, int MAX_TRASSA, int cod)
{
	int k;
	for (k = 0; k < MAX_TRASSA; k++) {
		trassa[k] = Nurbs3D_IGES(ch_TRASSA, cod, k+1);
		if (! trassa[k]) break;
	}
	Trassa_write(ch_TRASSA_OUT, trassa, k);
	return(k);
}

/////////////////////////////////////////////////////////////////////////
//...переформатриование трассы вместе с нормировкой и линкованием патчей;
void Trassa_link(char * ch_TRASSA, char * ch_TRASSA_OUT, CMap ** mp_trassa, int MAX_TRASSA, double f)
{
	int i, j, k, l;
/////////////////////////////////////////////
//...образуем поверхностную структуру трассы;
	Trassa_read(ch_TRASSA, mp_trassa, MAX_TRASSA, & k);
	CCells * trassa = new CCells(k); trassa->graph[0] = k;
	for (k = 0; k < trassa->graph[0]; k++) {
		  trassa->ce[k] = new CCells;
		  trassa->ce[k]->cells_new(0, 2);
		  trassa->ce[k]->mp = mp_trassa[k];
		  trassa->ce[k]->ce = trassa->ce;
		  trassa->ce[k]->graph[0] = trassa->graph[0];
	}

///////////////////////////////////////////////
//...нормировка трассы заданным коэффициентом;
	for (j = 0; j < k; j++) 
	if ( ID_MAP(2, NURBS_GENUS) == trassa->ce[j]->mp[0] && trassa->ce[j]->mp[size_of_map(2, NURBS_GENUS)] == (CMap)B2SPLINE_CELL) {
		int m1 = (int)trassa->ce[j]->mp[10], m2 = (int)trassa->ce[j]->mp[8],  
			 m3 = (int)trassa->ce[j]->mp[11], m4 = (int)trassa->ce[j]->mp[9], m5 = (int)trassa->ce[j]->mp[7]/4;
		CMap * p = trassa->ce[j]->mp+m2+m4+17;
		for (l = m1*m3; l > 0; l--, p = (m5 == 4 ? p+1 : p))  
		for (i = 0;     i < 3; i++, p++) p[0] *= f;
	}

////////////////////////
//... линкование трассы;
	int M_iter = 5, N_ini = 40, M_ini = 2;	
	double  u_par, v_par, P[3], eps = 1e-1;
	for (j = 0; j < k; j++) {
		for (l = 1; l < M_ini; l++) {
			nurbs3D(u_par = 0., v_par = (double)l/M_ini, trassa->ce[j]->mp, P, OK_STATE); 
			for (i = 0; i < k; i++) if (i != j) {
				if ((v_par = nurbs3D_v_par(u_par = 0., trassa->ce[i]->mp, P, eps, M_iter, N_ini)) >= 0.) {
					if (! topo_id(trassa->ce[j]->graph, i)) topo_insert(trassa->ce[j]->graph, i); 
					if (! topo_id(trassa->ce[i]->graph, j)) topo_insert(trassa->ce[i]->graph, j); 
					break;
				}
				if ((v_par = nurbs3D_v_par(u_par = 1., trassa->ce[i]->mp, P, eps, M_iter, N_ini)) >= 0.) {
					if (! topo_id(trassa->ce[j]->graph, i)) topo_insert(trassa->ce[j]->graph, i); 
					if (! topo_id(trassa->ce[i]->graph, j+k)) topo_insert(trassa->ce[i]->graph, j+k); 
					break;
				}
				if ((u_par = nurbs3D_u_par(v_par = 0., trassa->ce[i]->mp, P, eps, M_iter, N_ini)) >= 0.) {
					if (! topo_id(trassa->ce[j]->graph, i)) topo_insert(trassa->ce[j]->graph, i); 
					if (! topo_id(trassa->ce[i]->graph, j+k*2)) topo_insert(trassa->ce[i]->graph, j+k*2); 
					break;
				}
				if ((u_par = nurbs3D_u_par(v_par = 1., trassa->ce[i]->mp, P, eps, M_iter, N_ini)) >= 0.) {
					if (! topo_id(trassa->ce[j]->graph, i)) topo_insert(trassa->ce[j]->graph, i); 
					if (! topo_id(trassa->ce[i]->graph, j+k*3)) topo_insert(trassa->ce[i]->graph, j+k*3); 
					break;
				}
			}
		}
		for (l = 1; l < M_ini; l++) {
			nurbs3D(u_par = 1., v_par = (double)l/M_ini, trassa->ce[j]->mp, P, OK_STATE); 
			for (i = 0; i < k; i++) if (i != j) {
				if (i == 31 || i == 18)
					 i = i;
				if ((v_par = nurbs3D_v_par(u_par = 0., trassa->ce[i]->mp, P, eps, M_iter, N_ini)) >= 0.) {
					if (! topo_id(trassa->ce[j]->graph, i+k)) topo_insert(trassa->ce[j]->graph, i+k); 
					if (! topo_id(trassa->ce[i]->graph, j)) topo_insert(trassa->ce[i]->graph, j); 
					break;
				}
				if ((v_par = nurbs3D_v_par(u_par = 1., trassa->ce[i]->mp, P, eps, M_iter, N_ini)) >= 0.) {
					if (! topo_id(trassa->ce[j]->graph, i+k)) topo_insert(trassa->ce[j]->graph, i+k); 
					if (! topo_id(trassa->ce[i]->graph, j+k)) topo_insert(trassa->ce[i]->graph, j+k); 
					break;
				}
				if ((u_par = nurbs3D_u_par(v_par = 0., trassa->ce[i]->mp, P, eps, M_iter, N_ini)) >= 0.) {
					if (! topo_id(trassa->ce[j]->graph, i+k)) topo_insert(trassa->ce[j]->graph, i+k); 
					if (! topo_id(trassa->ce[i]->graph, j+k*2)) topo_insert(trassa->ce[i]->graph, j+k*2); 
					break;
				}
				if ((u_par = nurbs3D_u_par(v_par = 1., trassa->ce[i]->mp, P, eps, M_iter, N_ini)) >= 0.) {
					if (! topo_id(trassa->ce[j]->graph, i+k)) topo_insert(trassa->ce[j]->graph, i+k); 
					if (! topo_id(trassa->ce[i]->graph, j+k*3)) topo_insert(trassa->ce[i]->graph, j+k*3); 
					break;
				}
			}
		}
		for (l = 1; l < M_ini; l++) {
			nurbs3D(u_par = (double)l/M_ini, v_par = 0., trassa->ce[j]->mp, P, OK_STATE); 
			for (i = 0; i < k; i++) if (i != j) {
				if (i == 31 || i == 18)
					 i = i;
				if ((v_par = nurbs3D_v_par(u_par = 0., trassa->ce[i]->mp, P, eps, M_iter, N_ini)) >= 0.) {
					if (! topo_id(trassa->ce[j]->graph, i+k*2)) topo_insert(trassa->ce[j]->graph, i+k*2); 
					if (! topo_id(trassa->ce[i]->graph, j)) topo_insert(trassa->ce[i]->graph, j); 
					break;
				}
				if ((v_par = nurbs3D_v_par(u_par = 1., trassa->ce[i]->mp, P, eps, M_iter, N_ini)) >= 0.) {
					if (! topo_id(trassa->ce[j]->graph, i+k*2)) topo_insert(trassa->ce[j]->graph, i+k*2); 
					if (! topo_id(trassa->ce[i]->graph, j+k)) topo_insert(trassa->ce[i]->graph, j+k); 
					break;
				}
				if ((u_par = nurbs3D_u_par(v_par = 0., trassa->ce[i]->mp, P, eps, M_iter, N_ini)) >= 0.) {
					if (! topo_id(trassa->ce[j]->graph, i+k*2)) topo_insert(trassa->ce[j]->graph, i+k*2); 
					if (! topo_id(trassa->ce[i]->graph, j+k*2)) topo_insert(trassa->ce[i]->graph, j+k*2); 
					break;
				}
				if ((u_par = nurbs3D_u_par(v_par = 1., trassa->ce[i]->mp, P, eps, M_iter, N_ini)) >= 0.) {
					if (! topo_id(trassa->ce[j]->graph, i+k*2)) topo_insert(trassa->ce[j]->graph, i+k*2); 
					if (! topo_id(trassa->ce[i]->graph, j+k*3)) topo_insert(trassa->ce[i]->graph, j+k*3); 
					break;
				}
			}
		}
		for (l = 1; l < M_ini; l++) {
			nurbs3D(u_par = (double)l/M_ini, v_par = 1., trassa->ce[j]->mp, P, OK_STATE); 
			for (i = 0; i < k; i++) if (i != j) {
				if (i == 31 || i == 18)
					 i = i;
				if ((v_par = nurbs3D_v_par(u_par = 0., trassa->ce[i]->mp, P, eps, M_iter, N_ini)) >= 0.) {
					if (! topo_id(trassa->ce[j]->graph, i+k*3)) topo_insert(trassa->ce[j]->graph, i+k*3); 
					if (! topo_id(trassa->ce[i]->graph, j)) topo_insert(trassa->ce[i]->graph, j); 
					break;
				}
				if ((v_par = nurbs3D_v_par(u_par = 1., trassa->ce[i]->mp, P, eps, M_iter, N_ini)) >= 0.) {
					if (! topo_id(trassa->ce[j]->graph, i+k*3)) topo_insert(trassa->ce[j]->graph, i+k*3); 
					if (! topo_id(trassa->ce[i]->graph, j+k)) topo_insert(trassa->ce[i]->graph, j+k); 
					break;
				}
				if ((u_par = nurbs3D_u_par(v_par = 0., trassa->ce[i]->mp, P, eps, M_iter, N_ini)) >= 0.) {
					if (! topo_id(trassa->ce[j]->graph, i+k*3)) topo_insert(trassa->ce[j]->graph, i+k*3); 
					if (! topo_id(trassa->ce[i]->graph, j+k*2)) topo_insert(trassa->ce[i]->graph, j+k*2); 
					break;
				}
				if ((u_par = nurbs3D_u_par(v_par = 1., trassa->ce[i]->mp, P, eps, M_iter, N_ini)) >= 0.) {
					if (! topo_id(trassa->ce[j]->graph, i+k*3)) topo_insert(trassa->ce[j]->graph, i+k*3); 
					if (! topo_id(trassa->ce[i]->graph, j+k*3)) topo_insert(trassa->ce[i]->graph, j+k*3); 
					break;
				}
			}
		}
	}

//////////////////////////////////////////////////////////////////////////////////
//... запись отлинкованной трассы на диск в двоичном формате и выход из программы;
	FILE * device = fopen(ch_TRASSA_OUT, "w+b");
	if (device) {
		fwrite(& k, sizeof(int), 1, device);
		for (int j = 0; j < k; j++) {
			int m = trassa->ce[j]->graph[1]+2;
			fwrite(& m, sizeof(int), 1, device);
			fwrite(trassa->ce[j]->graph, sizeof(Topo), m, device);

			m = size_of_map(trassa->ce[j]->mp)+size_of_dop(trassa->ce[j]->mp)+1;
			fwrite(& m, sizeof(int), 1, device);
			fwrite(trassa->ce[j]->mp, sizeof(CMap), m, device);
		}
		fclose(device);
	}
	trassa->zero_cells(); delete trassa;
	return;
}

/////////////////////////////////////////////////////////////////////////////////
//...переформатриование скелета трассы вместе с нормировкой и линкованием кривых;
void Trassa_link_skel(char * ch_TRASSA, char * ch_TRASSA_OUT, CMap ** mp_trassa, int MAX_TRASSA, double f)
{
	int i, j, k, l;
/////////////////////////////////////////////
//...образуем поверхностную структуру трассы;
	Trassa_read(ch_TRASSA, mp_trassa, MAX_TRASSA, & k);
	CCells * trassa = new CCells(k); trassa->graph[0] = k;
	for (k = 0; k < trassa->graph[0]; k++) {
		  trassa->ce[k] = new CCells;
		  trassa->ce[k]->cells_new(0, 2);
		  trassa->ce[k]->mp = mp_trassa[k];
		  trassa->ce[k]->ce = trassa->ce;
		  trassa->ce[k]->graph[0] = trassa->graph[0];
	}

///////////////////////////////////////////////
//...нормировка трассы заданным коэффициентом;
	for (j = 0; j < k; j++) 
	if ( ID_MAP(1, NURBS_GENUS) == trassa->ce[j]->mp[0] && trassa->ce[j]->mp[size_of_map(1, NURBS_GENUS)] == (CMap)B1SPLINE_CELL) {
		int m1 = (int)trassa->ce[j]->mp[9], m2 = (int)trassa->ce[j]->mp[8], m3 = (int)trassa->ce[j]->mp[7]/2;
		CMap * p = trassa->ce[j]->mp+m2+13;
		for (l = m1; l > 0; l--, p = (m3 == 4 ? p+1 : p))  
		for (i = 0;  i < 3; i++, p++) p[0] *= f;
	}

////////////////////////
//... линкование трассы;
	double u_par, P[3], P_temp[3], eps = 1e-1;
	for (j = 0; j < k; j++) {
		nurbs3D(u_par = 0., trassa->ce[j]->mp, P, OK_STATE); 
		for (i = 0; i < k; i++) if (i != j) {
			nurbs3D(u_par = 0., trassa->ce[i]->mp, P_temp, OK_STATE);
			if (abs_point(P, P_temp) < eps) {
				if (! topo_id(trassa->ce[j]->graph, i)) topo_insert(trassa->ce[j]->graph, i); 
				if (! topo_id(trassa->ce[i]->graph, j)) topo_insert(trassa->ce[i]->graph, j); 
				break;
			}
			nurbs3D(u_par = 1., trassa->ce[i]->mp, P_temp, OK_STATE);
			if (abs_point(P, P_temp) < eps) {
				if (! topo_id(trassa->ce[j]->graph, i)) topo_insert(trassa->ce[j]->graph, i); 
				if (! topo_id(trassa->ce[i]->graph, j+k)) topo_insert(trassa->ce[i]->graph, j+k); 
				break;
			}
		}
		nurbs3D(u_par = 1., trassa->ce[j]->mp, P, OK_STATE); 
		for (i = 0; i < k; i++) if (i != j) {
			nurbs3D(u_par = 0., trassa->ce[i]->mp, P_temp, OK_STATE);
			if (abs_point(P, P_temp) < eps) {
				if (! topo_id(trassa->ce[j]->graph, i+k)) topo_insert(trassa->ce[j]->graph, i+k); 
				if (! topo_id(trassa->ce[i]->graph, j)) topo_insert(trassa->ce[i]->graph, j); 
				break;
			}
			nurbs3D(u_par = 1., trassa->ce[i]->mp, P_temp, OK_STATE);
			if (abs_point(P, P_temp) < eps) {
				if (! topo_id(trassa->ce[j]->graph, i+k)) topo_insert(trassa->ce[j]->graph, i+k); 
				if (! topo_id(trassa->ce[i]->graph, j+k)) topo_insert(trassa->ce[i]->graph, j+k); 
				break;
			}
		}
	}

//////////////////////////////////////////////////////////////////////////////////
//... запись отлинкованной трассы на диск в двоичном формате и выход из программы;
	FILE * device = fopen(ch_TRASSA_OUT, "w+b");
	if (device) {
		fwrite(& k, sizeof(int), 1, device);
		for (int j = 0; j < k; j++) {
			int m = trassa->ce[j]->graph[1]+2;
			fwrite(& m, sizeof(int), 1, device);
			fwrite(trassa->ce[j]->graph, sizeof(Topo), m, device);

			m = size_of_map(trassa->ce[j]->mp)+size_of_dop(trassa->ce[j]->mp)+1;
			fwrite(& m, sizeof(int), 1, device);
			fwrite(trassa->ce[j]->mp, sizeof(CMap), m, device);
		}
		fclose(device);
	}
	trassa->zero_cells(); delete trassa;
	return;
}

/////////////////////////////////////////////////////////////////////////
//...чтение отлинкованной трассы (или скелета трассы) в двоичном формате;
CCells * Trassa_read(char * ch_TRASSA)
{
	FILE * device = fopen(ch_TRASSA, "r+b"); int m, k = 0;
	CCells * trassa = NULL;
	if (device) {
      size_t res = fread(& k, sizeof(int), 1, device);
		trassa = new CCells(k); trassa->graph[0] = k;
		for (k = 0; k < trassa->graph[0]; k++) {
			res = fread(& m, sizeof(int), 1, device);
			trassa->ce[k] = new CCells;
			trassa->ce[k]->cells_new(0, m);
			trassa->ce[k]->ce = trassa->ce;
			res = fread(trassa->ce[k]->graph, sizeof(Topo), m, device);

			res = fread(& m, sizeof(int), 1, device);
			trassa->ce[k]->mp = (CMap *)new_struct(m*sizeof(CMap));
			res = fread(trassa->ce[k]->mp, sizeof(CMap), m, device);
		}
		fclose (device);
	}
	return(trassa);
}

//////////////////////////////////////
//...запись траектории в формате IGES;
void track_igs(char * ch_TRACK, CGrid * nd_track)
{
	FILE * device = fopen(ch_TRACK, "w+b"); int i, j, k, l;
	if (device) {
		char buff[100], probel[100];
///////////////////////////////
//...запись начальных разделов;
		fprintf(device, "3D InterOp IGES (Version 20 0 3), Spatial Corp. Copyright (c) 1999-2007.S      1\n");
		fprintf(device, "1H,,1H;,6HNoname,9HTrack.igs,13HSpatial Corp.,20H3D InterOp ACIS/IGES,  G      1\n");
		fprintf(device, "32,38, 6,308,15,6HNoname,1.000,2,2HMM,1,1.000,15H20121029.143951,       G      2\n");
		fprintf(device, "1.0e-006,0.00,6HNoname,6HNoname,11,0,15H20121029.143951;                G      3\n");
//////////////////////////////
//...запись раздела вхождений;
		for (k = 1; k < nd_track->N; k++) { 
			fprintf(device, "     110 %7i       0       0       0       0       0       000010001D%7i\n", 2*k-1, 2*k-1);
			fprintf(device, "     110       0       0       2       0       0       0               0D%7i\n", 2*k);
		}
		sprintf(buff, "102,%i,", nd_track->N-1); l = (int)strlen(buff);
		for (i = k = 0; k < nd_track->N-1; k++) {
			if (k < nd_track->N-2) sprintf(buff, "%i,", 2*k+1); 
			else						  sprintf(buff, "%i;", 2*k+1); 
			if (l+(j = (int)strlen(buff)) <= 64) l += j; 
			else {
				i++; l = j;
			}
		}
		fprintf(device, "     102 %7i       0       0       0       0       0       000000001D%7i\n", 2*k+1, 2*k+1);
		fprintf(device, "     102       0       0 %7i       0       0       0               0D%7i\n", i+1, 2*k+2);

///////////////////////////////
//...запись раздела параметров;
		for (k = 1; k < nd_track->N; k++) { 
			fprintf(device, "110,%19g,%19g,%19g, %7iP%7i\n", nd_track->X[k-1], nd_track->Y[k-1], nd_track->Z[k-1], 2*k-1, 2*k-1);
			fprintf(device, "%19g,%19g,%19g;     %7iP%7i\n", nd_track->X[k], nd_track->Y[k], nd_track->Z[k],2*k-1, 2*k);
		}
		sprintf(buff, "102,%i,", nd_track->N-1); l = (int)strlen(buff);
		fprintf(device, "%s", buff);
		for (i = k = 0; k < nd_track->N-1; k++) {
			if (k < nd_track->N-2) sprintf(buff, "%i,", 2*k+1); 
			else						  sprintf(buff, "%i;", 2*k+1); 
			if (l+(j = (int)strlen(buff)) <= 64) {
				l += j; fprintf(device, "%s", buff);
			}
			else {
				sprintf(probel, "                                                                "); probel[64-l] = 0;
				fprintf(device, "%s%8iP%7i\n", probel, 2*nd_track->N-1, 2*nd_track->N-1+i); i++;
				fprintf(device, "%s", buff); l = j;
			}
		}
		sprintf(probel, "                                                                "); probel[64-l] = 0;
		fprintf(device, "%s%8iP%7i\n", probel, 2*nd_track->N-1, 2*nd_track->N-1+i);
/////////////////////////////////
//...запись завершающего раздела;
		fprintf(device, "S      1G      3D%7iP%7i                                        T      1\n", 2*nd_track->N, 2*nd_track->N-1+i);
		fclose(device);
	}
	return;
}

/////////////////////////////////////////////////
//...запись совокупности маркеров в формате IGES;
void marker_igs(char * ch_MARKER, CGrid * nd_marker)
{
	FILE * device = fopen(ch_MARKER, "w+b"); int i, j, k, l;
	if (device) {
		char buff[100], probel[100];
///////////////////////////////
//...запись начальных разделов;
		fprintf(device, "3D InterOp IGES (Version 20 0 3), Spatial Corp. Copyright (c) 1999-2007.S      1\n");
		fprintf(device, "1H,,1H;,6HNoname,9HTrack.igs,13HSpatial Corp.,20H3D InterOp ACIS/IGES,  G      1\n");
		fprintf(device, "32,38, 6,308,15,6HNoname,1.000,2,2HMM,1,1.000,15H20121029.143951,       G      2\n");
		fprintf(device, "1.0e-006,0.00,6HNoname,6HNoname,11,0,15H20121029.143951;                G      3\n");
//////////////////////////////
//...запись раздела вхождений;
		for (k = 1; k < nd_marker->N; k++) { 
			fprintf(device, "     110 %7i       0       0       0       0       0       000010001D%7i\n", 2*k-1, 2*k-1);
			fprintf(device, "     110       0       0       2       0       0       0               0D%7i\n", 2*k);
		}
		sprintf(buff, "102,%i,", nd_marker->N-1); l = (int)strlen(buff);
		for (i = k = 0; k < nd_marker->N-1; k++) {
			if (k < nd_marker->N-2) sprintf(buff, "%i,", 2*k+1); 
			else							sprintf(buff, "%i;", 2*k+1); 
			if (l+(j = (int)strlen(buff)) <= 64) l += j; 
			else {
				i++; l = j;
			}
		}
		fprintf(device, "     102 %7i       0       0       0       0       0       000000001D%7i\n", 2*k+1, 2*k+1);
		fprintf(device, "     102       0       0 %7i       0       0       0               0D%7i\n", i+1, 2*k+2);

///////////////////////////////
//...запись раздела параметров;
		for (k = 1; k < nd_marker->N; k++) { 
			fprintf(device, "110,%19g,%19g,%19g, %7iP%7i\n", nd_marker->X[k-1], nd_marker->Y[k-1], nd_marker->Z[k-1], 2*k-1, 2*k-1);
			if (nd_marker->hit[k-1] == nd_marker->hit[k])
				fprintf(device, "%19g,%19g,%19g;     %7iP%7i\n", nd_marker->X[k], nd_marker->Y[k], nd_marker->Z[k],2*k-1, 2*k); else
				fprintf(device, "%19g,%19g,%19g;     %7iP%7i\n", nd_marker->X[k-1], nd_marker->Y[k-1], nd_marker->Z[k-1],2*k-1, 2*k);
		}
		sprintf(buff, "102,%i,", nd_marker->N-1); l = (int)strlen(buff);
		fprintf(device, "%s", buff);
		for (i = k = 0; k < nd_marker->N-1; k++) {
			if (k < nd_marker->N-2) sprintf(buff, "%i,", 2*k+1); 
			else							sprintf(buff, "%i;", 2*k+1); 
			if (l+(j = (int)strlen(buff)) <= 64) {
				l += j; fprintf(device, "%s", buff);
			}
			else {
				sprintf(probel, "                                                                "); probel[64-l] = 0;
				fprintf(device, "%s%8iP%7i\n", probel, 2*nd_marker->N-1, 2*nd_marker->N-1+i); i++;
				fprintf(device, "%s", buff); l = j;
			}
		}
		sprintf(probel, "                                                                "); probel[64-l] = 0;
		fprintf(device, "%s%8iP%7i\n", probel, 2*nd_marker->N-1, 2*nd_marker->N-1+i);
/////////////////////////////////
//...запись завершающего раздела;
		fprintf(device, "S      1G      3D%7iP%7i                                        T      1\n", 2*nd_marker->N, 2*nd_marker->N-1+i);
		fclose(device);
	}
	return;
}
