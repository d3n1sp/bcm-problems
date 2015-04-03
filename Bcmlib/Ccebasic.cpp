#include "stdafx.h"
#include "ccebasic.h"

//////////////////////////////////////////////////////////
//...прямое формирование плоской фасеты в 3D пространстве;
void CCeBasic::get_facet_directly(double * P, int N, int id_plane, int id_fast, int id_dbl_point, double eps)
{ //...id_fast == OK_STATE - плоский элемент с нормалью и точками (без прямых линий); 
	if (P && N >= 3) {
		int   i, k, l, m = id_fast == OK_STATE ? N : N*2;
		cells_new(m+1, N+2, (l = size_of_map(2, NULL_GENUS))+1+size_of_dop(FACET_CELL));

////////////////////////
//...вычисление нормали;
      double p0[3], fN, fi, theta, CZ = 0., SZ = 0., CY = 0., SY = 0.,
				ort[3] = { (P[4]-P[1])*(P[N*3-1]-P[2])-(P[5]-P[2])*(P[N*3-2]-P[1]),
							  (P[5]-P[2])*(P[N*3-3]-P[0])-(P[3]-P[0])*(P[N*3-1]-P[2]),
							  (P[3]-P[0])*(P[N*3-2]-P[1])-(P[4]-P[1])*(P[N*3-3]-P[0])}; 

///////////////////////////////////////////////////////////////////////////////////////
//...прямая запись нормали или углов Эйлера (сферических углов) в геометрическую карту;
      if (id_fast == NULL_STATE) {
			CZ = cos(fi = arg0(comp(ort[0], ort[1])));
			SZ = sin(fi);
			CY = cos(theta = arg0(comp(ort[2], sqrt(ort[0]*ort[0]+ort[1]*ort[1]))));
			SY = sin(theta);
		}
		else {
			if ((fi = sqrt(sqr(ort[0])+sqr(ort[1])+sqr(ort[2]))) > EE) {
				ort[0] *= (fi = 1./fi);
				ort[1] *=  fi;
				ort[2] *=  fi;
			}
			else {
				ort[0] = 0.;
				ort[1] = 0.;
				ort[2] = 1.;
			}
			if ((SY = sqrt(ort[0]*ort[0]+ort[1]*ort[1])) > EE) {
				CZ = ort[0]*(fi = 1./SY);
				SZ = ort[1]* fi;
			}
			else {
				CZ = 1.;
				SZ = 0.;
			}
			CY = ort[2];
      }

/////////////////////////////
//...проверка плоской фигуры;
		if (id_plane == OK_STATE)
		for (i = k = 0; k < N; k++) {
			p0[0] = P[i++]-P[0];
			p0[1] = P[i++]-P[1];
			p0[2] = P[i++]-P[2];
			point_iso<double>(p0, NULL, 1., 0., CY, -SY, CZ, -SZ);
			if  (fabs(p0[2]) >= eps) goto err;
		}

/////////////////////////////////
//...запись геометрической карты;
		if (mp) {
			mp[0] = ID_MAP(2, NULL_GENUS);
			mp[1] = (CMap)P[0];
			mp[2] = (CMap)P[1];
			mp[3] = (CMap)P[2];
         mp[l] = (CMap)FACET_CELL;
			if (id_fast == NULL_STATE) {
				mp[4] = (CMap)fi;
				mp[5] = (CMap)theta;
				mp[6] = (CMap)0.;
			}
			else {
				mp[4] = (CMap)ort[0];
				mp[5] = (CMap)ort[1];
				mp[6] = (CMap)ort[2];
				for (p0[0] = p0[1] = p0[2] = 0, fN = 1./N, i = k = 0; k < N; k++) { //...центроид;
					p0[0] += P[i++]*fN;
					p0[1] += P[i++]*fN;
					p0[2] += P[i++]*fN;
				}
				mp[1] = (CMap)p0[0];
				mp[2] = (CMap)p0[1];
				mp[3] = (CMap)p0[2];
			}
		}

//////////////////////////////////////////////////
//...создаем топологию пространственного полигона;
		if (graph) {
			for ( l = size_of_map(0, NULL_GENUS), k = N-1; ce && k >= 0; k--) {
				  ce[k] = new CCells(-1); if (! ce[k]) goto err;
				  ce[k]->cells_new(0, 2, l+1);
				  ce[k]->mp[l]		= (CMap)NULL_CELL;
				  ce[k]->ce			= ce;
				  ce[k]->graph[0] = m;
				  graph[k+2]		= m-N+k;
			}
			if (id_fast != OK_STATE)
			for ( l = size_of_map(1, NULL_GENUS), k = N; ce && k < m; k++) {
				  ce[k] = new CCells(-1); if (! ce[k]) goto err;
				  ce[k]->cells_new(0, 4, l+1);
				  ce[k]->mp[l]    = (CMap)NULL_CELL;
				  ce[k]->ce       = ce;
				  ce[k]->graph[0] = m;
				  ce[k]->graph[1] = 2;
				  ce[k]->graph[2] = k-N;
				  ce[k]->graph[3] = (k-N+1)%N;
			}
			graph[0] = m;
			graph[1] = N;
		}

//////////////////////////////////////////////////
//...создаем геометрию пространственного полигона;
      for (i = k = 0; k < N; k++) {
           ce[k]->mp[0] = ID_MAP(0, NULL_GENUS);
           ce[k]->mp[1] = P[i++];
           ce[k]->mp[2] = P[i++];
           ce[k]->mp[3] = P[i++];
      }
		if (id_fast != OK_STATE)
      for (; k < m; k++) {
           ce[k]->mp[0] = ID_MAP(1, NULL_GENUS);
           ce[k]->line_correct(id_fast);
      }

/////////////////////////////////
//...удаляем повторяющиеся точки;
		if (id_dbl_point == OK_STATE)
      for (k = N; k < m; k++) ce[k]->cell_correct_double_point(eps);
      return;
  }
err:
  zero_cells();
}

//////////////////////////////
//...forming plane tria-facet;
void CCeBasic::get_tria_facet(double * P1, double * P2, double * P3, int id_fast)
{
	double P[9] = { P1[0], P1[1], P1[2], 
						 P2[0], P2[1], P2[2], 
						 P3[0], P3[1], P3[2]	}, A, B, E, R;
	get_facet_directly(P, 3, NULL_STATE, id_fast, OK_STATE);
	A = abs_point(P,   P+3);
	B = abs_point(P+3, P+6);
	E = abs_point(P,   P+6); R = (A+B+E)*.5;
	mp[size_of_map(mp)+1] = (CMap)sqrt(fabs(R*(R-A)*(R-B)*(R-E)));
}

//////////////////////////////
//...forming plane quad-facet;
void CCeBasic::get_quad_facet(double * P1, double * P2, double * P3, double * P4, int id_fast)
{
	double P[12] = {P1[0], P1[1], P1[2], 
						 P2[0], P2[1], P2[2], 
						 P3[0], P3[1], P3[2], 
						 P4[0], P4[1], P4[2]}, A, B, C, D, E, Q, R;
	get_facet_directly(P, 4, NULL_STATE, id_fast, OK_STATE);
	A = abs_point(P,   P+3);
	B = abs_point(P+3, P+6);
	E = abs_point(P,   P+6); R = (A+B+E)*.5;
	C = abs_point(P+6, P+9);
	D = abs_point(P+9, P);   Q = (C+D+E)*.5;
	mp[size_of_map(mp)+1] = (CMap)sqrt(fabs(R*(R-A)*(R-B)*(R-E)))+sqrt(fabs(Q*(Q-C)*(Q-D)*(Q-E)));
}

/////////////////////////////////////////////////////
//...setting additional parameters for a plane facet;
int CCeBasic::SetFacetParam(int k, double square, double * pp, int id_property)
{
	int l, j;
	if (ce && graph && 0 <= k && k < graph[0] && ce[k]) {
		CMap * mp = ce[k]->mp;
      if (mp && ID_MAP(2, NULL_GENUS)  == mp[0] &&
          mp[  l  = size_of_map(2, NULL_GENUS)] == (CMap)FACET_CELL) {
          mp[++l] = (CMap)square;

          double p1, p2, p3, p4;
          if (id_property == DEFAULT_BND) p1 = p2 = p3 = p4 = 0.; else
          if (id_property == RIGID_BND)   p1 = p2 = p3 = p4 = 0.; else
          if (id_property == PN_BND && pp) {
					p1 = 0.;
					p2 = MIN_HIT;
					p3 = pp[2];
					p4 = pp[3];
          }
          else
          if (id_property == VEL_BND && pp) {
					p1 = p2 = 0.;
					p3 = (pp[0]*mp[4]+pp[1]*mp[5]+pp[2]*mp[6]);
					p4 = (pp[3]*mp[4]+pp[4]*mp[5]+pp[5]*mp[6]);
          }
          else
          if (id_property == VN_BND && pp) {
					p1 = p2 = 0.;
					p3 = pp[2];
					p4 = pp[3];
          }
          else
          if (id_property == VNTABLE_BND && pp) {
					p1 = p2 = p3 = 0.;
					p4 = MIN_HIT;
          }
          else
          if (id_property == ABSORB_BND && pp) {
					double f, w1 = (f = sqrt(fabs(1.-pp[0])))*cos(pp[1]*M_PI),
								 w2 =                          f*sin(pp[1]*M_PI);
					p2 = (-2.*w2/(f = (1.+w1)*(1.+w1)+w2*w2));
					p1 = ((1.-w1*w1-w2*w2)/f);//...because of external normal vector on boundary;
					p3 = pp[2];
					p4 = pp[3];
          }
          else
          if (id_property == ABSORBTABLE_BND) {
					double f, w1 = (f = sqrt(fabs(1.-pp[0])))*cos(pp[1]*M_PI),
								 w2 =                          f*sin(pp[1]*M_PI);
					p2 = (-2.*w2/(f = (1.+w1)*(1.+w1)+w2*w2));
					p1 = ((1.-w1*w1-w2*w2)/f);//...because of external normal vector on boundary;
					p3 = pp[2];
					p4 = MAX_HIT;
          }
          else
          if (id_property == EXTERN_BND) {
					p1 = MIN_HIT;
					p2 = p3 = p4 = 0.;
          }
          else
          if (id_property == DISPLN_BND) {
					p1 = pp[0];
					p2 = p3 = 0.;
					p4 = MIN_HIT;
          }
          else
          if (id_property == FORCEN_BND) {
					p1 = pp[0];
					p2 = p3 = 0.;
					p4 = 1.;
          }
          else
          if (id_property == ADHESION_BND) {
					p1 = p2 = p3 = 0.;
					p4 = 0.;
          }
            else
          if (id_property == PRESSURE_BND) {
					p1 = pp[0];
					p2 = p3 = 0.;
					p4 = MAX_HIT;
          }
          else
          if (id_property == VELOCITY_BND) {
					p1 = pp[0];
					p2 = pp[1];
					p3 = pp[2];
					p4 = 1.;
          }
          else
            if (id_property == VELOCITY_N_BND) {
					p1 = pp[0];
					p2 = p3 = 0.;
					p4 = 2.;
          }
          else
          if (id_property == RIGID_DISPL3D_BND) {
					p1 = pp[0];
					p2 = pp[1];
					p3 = pp[2];
					p4 = MIN_HIT;
          }
          else
          if (id_property == DISPL3D_BND) {
					p1 = pp[0];
					p2 = pp[1];
					p3 = pp[2];
					p4 = MAX_HIT;
          }
          else
          if (id_property == FORCE3D_BND) {
					p1 = pp[0];
					p2 = pp[1];
					p3 = pp[2];
					p4 = MARKER1_BND-SPECIAL_BND;
          }
          else
          if (id_property == MOMENT3D_BND) {
					p1 = pp[0];
					p2 = pp[1];
					p3 = pp[2];
					p4 = MOMENTS_BND-SPECIAL_BND;
          }
          else
          if (id_property == FNORMS_BND) {
					p1 = pp[0];
					p2 = pp[1];
					p3 = pp[2];
					p4 = NORMS_BND-SPECIAL_BND;
          }
          else
          if (id_property == FSKEWS_BND) {
					p1 = pp[0];
					p2 = p3 = 0.;
					p4 = SKEWS_BND-SPECIAL_BND;
          }
          else
          if (id_property == UPTAKES_BND) {
					p1 = pp[0];
					p2 = p3 = 0.;
					p4 = UPTAKE_BND-SPECIAL_BND;
          }
          else
          if (id_property == SPECIAL_BND) {
					p1 = p2 = p3 = 0.;
					p4 = NUMS_BND;
          }
          else {
				 for (j = MARKER1_BND; j < NUMS_BND; j++) {
					  if (id_property == j) {
							p1 = p2 = p3 = 0.;
							p4 = j-SPECIAL_BND+1.;
							break;
					  }
				 }
				 if (j == NUMS_BND) p1 = p2 = p3 = p4 = 0.;
			 }
			mp[++l] = (CMap)p1;
			mp[++l] = (CMap)p2;
			mp[++l] = (CMap)p3;
			mp[++l] = (CMap)p4;

			return(1);
		}
	}
	return(0);
}

///////////////////////////////////////////////////////
//...setting additional parameters for a boundary node;
int CCeBasic::SetNodeParam(double * pp, int id_property)
{
	if (mp && ID_MAP(0, NULL_GENUS) == mp[0] && add_new_cell_maps(mp, BOUNDARY_NODE)) {
		double p1, p2, p3, p4;
		if (id_property == DEFAULT_BND) p1 = p2 = p3 = p4 = 0.; else
		if (id_property == RIGID_BND) { //...жесткое закрепление перемещений и углов;
			p1 = p2 = p3 = 0.;
			p4 = MIN_HIT;
		}
		else
		if (id_property == RIGID_DISPL3D_BND) { //...жесткий изгиб (закрепление перемещений и углов);
			p1 = pp[0];
			p2 = pp[1];
			p3 = pp[2];
			p4 = MIN_HIT;
		}
		else
		if (id_property == DISPL3D_BND) { //...изгиб перемещением со свободным моментом;
			p1 = pp[0];
			p2 = pp[1];
			p3 = pp[2];
			p4 = MAX_HIT;
		}
		else
		if (id_property == FORCE3D_BND) { //...изгиб заданной силой cо свободным моментом;
			p1 = pp[0];
			p2 = pp[1];
			p3 = pp[2];
			p4 = 1.;
		}
		else
		if (id_property == MOMENT3D_BND) { //...изгиб заданным моментом со свободной силой;
			p1 = pp[0];
			p2 = pp[1];
			p3 = pp[2];
			p4 = 20.;
		}
		else {
			int  j;
			for (j = MARKER1_BND; j < NUMS_BND; j++) { //...расстановка маркеров точек приложения внещней силы;
				if (id_property == j) {
					p1 = p2 = p3 = 0.;
					p4 = j-SPECIAL_BND+1.;
					break;
				}
			}
			if (j == NUMS_BND) p1 = p2 = p3 = p4 = 0.;
		}
		int   l = size_of_map(0, NULL_GENUS);
		mp[++l] = (CMap)p1;
		mp[++l] = (CMap)p2;
		mp[++l] = (CMap)p3;
		mp[++l] = (CMap)p4;

		return(1);
	}
	return(0);
}

/////////////////////////////////////////////////////////////////
//...additional parameters for a plane restriction boundary node;
void CCeBasic::SetNodeXParam(double X0, double * pp, int id_property) 
{
	for (int k = 0; graph && k < graph[0]; k++)
	for (int num = 0; num < ce[k]->graph[1]; num++)
		if ( fabs(ce[ce[k]->graph[2+num]]->mp[1]-X0) < EE)
	((CCeBasic *)ce[ce[k]->graph[2+num]])->SetNodeParam(pp, id_property);
}

void CCeBasic::SetNodeYParam(double Y0, double * pp, int id_property) 
{
	for (int k = 0; graph && k < graph[0]; k++)
	for (int num = 0; num < ce[k]->graph[1]; num++)
		if ( fabs(ce[ce[k]->graph[2+num]]->mp[2]-Y0) < EE)
	((CCeBasic *)ce[ce[k]->graph[2+num]])->SetNodeParam(pp, id_property);
}

void CCeBasic::SetNodeZParam(double Z0, double * pp, int id_property) 
{
	for (int k = 0; graph && k < graph[0]; k++)
	for (int num = 0; num < ce[k]->graph[1]; num++)
		if ( fabs(ce[ce[k]->graph[2+num]]->mp[3]-Z0) < EE)
	((CCeBasic *)ce[ce[k]->graph[2+num]])->SetNodeParam(pp, id_property);
}

////////////////////////////////////////////////////////////////
//...additional parameters for a edge restriction boundary node;
void CCeBasic::SetNodeXYParam(double X0, double Y0, double * pp, int id_property) 
{
	for (int k = 0; graph && k < graph[0]; k++)
	for (int num = 0; num < ce[k]->graph[1]; num++)
		if ( fabs(ce[ce[k]->graph[2+num]]->mp[1]-X0) < EE &&
			  fabs(ce[ce[k]->graph[2+num]]->mp[2]-Y0) < EE)
	((CCeBasic *)ce[ce[k]->graph[2+num]])->SetNodeParam(pp, id_property);
}

void CCeBasic::SetNodeXZParam(double X0, double Z0, double * pp, int id_property) 
{
	for (int k = 0; graph && k < graph[0]; k++)
	for (int num = 0; num < ce[k]->graph[1]; num++)
		if ( fabs(ce[ce[k]->graph[2+num]]->mp[1]-X0) < EE &&
			  fabs(ce[ce[k]->graph[2+num]]->mp[3]-Z0) < EE)
	((CCeBasic *)ce[ce[k]->graph[2+num]])->SetNodeParam(pp, id_property);
}

void CCeBasic::SetNodeYZParam(double Y0, double Z0, double * pp, int id_property) 
{
	for (int k = 0; graph && k < graph[0]; k++)
	for (int num = 0; num < ce[k]->graph[1]; num++)
		if ( fabs(ce[ce[k]->graph[2+num]]->mp[2]-Y0) < EE &&
			  fabs(ce[ce[k]->graph[2+num]]->mp[3]-Z0) < EE)
	((CCeBasic *)ce[ce[k]->graph[2+num]])->SetNodeParam(pp, id_property);
}

//////////////////////////////////////////////////////////////////
//...additional parameters for a vertex restriction boundary node;
void CCeBasic::SetNodeXYZParam(double X0, double Y0, double Z0, double * pp, int id_property) 
{
	for (int k = 0; graph && k < graph[0]; k++)
	for (int num = 0; num < ce[k]->graph[1]; num++)
		if ( fabs(ce[ce[k]->graph[2+num]]->mp[1]-X0) < EE &&
			  fabs(ce[ce[k]->graph[2+num]]->mp[2]-Y0) < EE &&
			  fabs(ce[ce[k]->graph[2+num]]->mp[3]-Z0) < EE)
	((CCeBasic *)ce[ce[k]->graph[2+num]])->SetNodeParam(pp, id_property);
}

/////////////////////////////////////////////////////////////////////////
//...forming plane facet composition (element) on the base OpenGL format;
void CCeBasic::get_nd_bar_directly(CGrid * nd, int k)
{
	if (! nd->geom || ! nd->geom_ptr || ! nd->X || ! nd->Y || ! nd->Z || nd->geom[0] < k || k < 0) return;
	int      i, m1, m2, m3, m4, m5, m6, m7, m8, l;
	double   P[12], A, B, C, D, E, R, Q, pp[4] = {0., 0., 0., 0.}; 
	CCeBasic * ce;
  
/////////////////////////////////////////
//...forming surface geometry of element;
	cells_new (1, 2, 0);
	switch (nd->geom[i = nd->geom_ptr[k]]) {
		case GL_LINE_STRIP: if (nd->geom[i+1] >= 4) { //...composite line;
			for (l = nd->geom[i+1]-2; l >= 2; l--) {
				m1 = nd->geom[i+l+2];
				m2 = nd->geom[i+l+3];

				P[0] = nd->X[m1]; P[1] = nd->Y[m1]; P[2] = nd->Z[m1];
				P[3] = nd->X[m2]; P[4] = nd->Y[m2]; P[5] = nd->Z[m2];
				ce = new CCeBasic;
				ce->get_line(P, P+3, NULL_STATE, OK_STATE);
				i = bar_add(ce);
			}
		}  break;
		case GL_TRIANGLES: if (nd->geom[i+1] == 5) {
			m1   = nd->geom[i+4];
			m2   = nd->geom[i+5];
			m3   = nd->geom[i+6];
////////////////
//...tria facet;
			P[0] = nd->X[m1]; P[1]  = nd->Y[m1]; P[2]  = nd->Z[m1];
			P[3] = nd->X[m2]; P[4]  = nd->Y[m2]; P[5]  = nd->Z[m2];
			P[6] = nd->X[m3]; P[7]  = nd->Y[m3]; P[8]  = nd->Z[m3];
			ce = new CCeBasic;
			ce->get_facet_directly(P, 3, NULL_STATE, ZERO_STATE, NULL_STATE);
			i = bar_add(ce);
			A = abs_point(P,   P+3);
			B = abs_point(P+3, P+6);
			E = abs_point(P,   P+6); R = (A+B+E)*.5;
			SetFacetParam(i-1, sqrt(fabs(R*(R-A)*(R-B)*(R-E))), pp);
		}  break;
		case GL_QUADS: if (nd->geom[i+1] == 6) {
			m1   = nd->geom[i+4];
			m2   = nd->geom[i+5];
			m3   = nd->geom[i+6];
			m4   = nd->geom[i+7];
////////////////
//...quad facet;
			P[0] = nd->X[m1]; P[1]  = nd->Y[m1]; P[2]  = nd->Z[m1];
			P[3] = nd->X[m2]; P[4]  = nd->Y[m2]; P[5]  = nd->Z[m2];
			P[6] = nd->X[m3]; P[7]  = nd->Y[m3]; P[8]  = nd->Z[m3];
			P[9] = nd->X[m4]; P[10] = nd->Y[m4]; P[11] = nd->Z[m4];
			ce = new CCeBasic;
			ce->get_facet_directly(P, 4, NULL_STATE, ZERO_STATE, NULL_STATE);
			i = bar_add(ce);
			A = abs_point(P,   P+3);
			B = abs_point(P+3, P+6);
			E = abs_point(P,   P+6); R = (A+B+E)*.5;
			C = abs_point(P+6, P+9);
			D = abs_point(P+9, P);   Q = (C+D+E)*.5;
			SetFacetParam(i-1, sqrt(fabs(R*(R-A)*(R-B)*(R-E)))+sqrt(fabs(Q*(Q-C)*(Q-D)*(Q-E))), pp);
		}  break;
		case GL_BOXS: if (nd->geom[i+1] == 10) {
			m1   = nd->geom[i+4];
			m2   = nd->geom[i+5];
			m3   = nd->geom[i+6];
			m4   = nd->geom[i+7];
			m5   = nd->geom[i+8];
			m6   = nd->geom[i+9];
			m7   = nd->geom[i+10];
			m8   = nd->geom[i+11];
/////////////
//...face X';
			P[0] = nd->X[m1]; P[1]  = nd->Y[m1]; P[2]  = nd->Z[m1];
			P[3] = nd->X[m5]; P[4]  = nd->Y[m5]; P[5]  = nd->Z[m5];
			P[6] = nd->X[m8]; P[7]  = nd->Y[m8]; P[8]  = nd->Z[m8];
			P[9] = nd->X[m4]; P[10] = nd->Y[m4]; P[11] = nd->Z[m4];
			ce = new CCeBasic;
			ce->get_facet_directly(P, 4, NULL_STATE, OK_STATE, NULL_STATE);
			i = bar_add(ce);
			A = abs_point(P,   P+3);
			B = abs_point(P+3, P+6);
			E = abs_point(P,   P+6); R = (A+B+E)*.5;
			C = abs_point(P+6, P+9);
			D = abs_point(P+9, P);   Q = (C+D+E)*.5;
			SetFacetParam(i-1, sqrt(fabs(R*(R-A)*(R-B)*(R-E)))+sqrt(fabs(Q*(Q-C)*(Q-D)*(Q-E))), pp);
/////////////
//...face X;
			P[0] = nd->X[m2]; P[1]  = nd->Y[m2]; P[2]  = nd->Z[m2];
			P[3] = nd->X[m3]; P[4]  = nd->Y[m3]; P[5]  = nd->Z[m3];
			P[6] = nd->X[m7]; P[7]  = nd->Y[m7]; P[8]  = nd->Z[m7];
			P[9] = nd->X[m6]; P[10] = nd->Y[m6]; P[11] = nd->Z[m6];
			ce = new CCeBasic;
			ce->get_facet_directly(P, 4, NULL_STATE, OK_STATE, NULL_STATE);
			i = bar_add(ce);
			A = abs_point(P,   P+3);
			B = abs_point(P+3, P+6);
			E = abs_point(P,   P+6); R = (A+B+E)*.5;
			C = abs_point(P+6, P+9);
			D = abs_point(P+9, P);   Q = (C+D+E)*.5;
			SetFacetParam(i-1, sqrt(fabs(R*(R-A)*(R-B)*(R-E)))+sqrt(fabs(Q*(Q-C)*(Q-D)*(Q-E))), pp);
/////////////
//...face Y';
			P[0] = nd->X[m1]; P[1]  = nd->Y[m1]; P[2]  = nd->Z[m1];
			P[3] = nd->X[m2]; P[4]  = nd->Y[m2]; P[5]  = nd->Z[m2];
			P[6] = nd->X[m6]; P[7]  = nd->Y[m6]; P[8]  = nd->Z[m6];
			P[9] = nd->X[m5]; P[10] = nd->Y[m5]; P[11] = nd->Z[m5];
			ce = new CCeBasic;
			ce->get_facet_directly(P, 4, NULL_STATE, OK_STATE, NULL_STATE);
			i = bar_add(ce);
			A = abs_point(P,   P+3);
			B = abs_point(P+3, P+6);
			E = abs_point(P,   P+6); R = (A+B+E)*.5;
			C = abs_point(P+6, P+9);
			D = abs_point(P+9, P);   Q = (C+D+E)*.5;
			SetFacetParam(i-1, sqrt(fabs(R*(R-A)*(R-B)*(R-E)))+sqrt(fabs(Q*(Q-C)*(Q-D)*(Q-E))), pp);
/////////////
//...face Y;
			P[0] = nd->X[m4]; P[1]  = nd->Y[m4]; P[2]  = nd->Z[m4];
			P[3] = nd->X[m8]; P[4]  = nd->Y[m8]; P[5]  = nd->Z[m8];
			P[6] = nd->X[m7]; P[7]  = nd->Y[m7]; P[8]  = nd->Z[m7];
			P[9] = nd->X[m3]; P[10] = nd->Y[m3]; P[11] = nd->Z[m3];
			ce = new CCeBasic;
			ce->get_facet_directly(P, 4, NULL_STATE, OK_STATE, NULL_STATE);
			i = bar_add(ce);
			A = abs_point(P,   P+3);
			B = abs_point(P+3, P+6);
			E = abs_point(P,   P+6); R = (A+B+E)*.5;
			C = abs_point(P+6, P+9);
			D = abs_point(P+9, P);   Q = (C+D+E)*.5;
			SetFacetParam(i-1, sqrt(fabs(R*(R-A)*(R-B)*(R-E)))+sqrt(fabs(Q*(Q-C)*(Q-D)*(Q-E))), pp);
/////////////
//...face Z';
			P[0] = nd->X[m1]; P[1]  = nd->Y[m1]; P[2]  = nd->Z[m1];
			P[3] = nd->X[m4]; P[4]  = nd->Y[m4]; P[5]  = nd->Z[m4];
			P[6] = nd->X[m3]; P[7]  = nd->Y[m3]; P[8]  = nd->Z[m3];
			P[9] = nd->X[m2]; P[10] = nd->Y[m2]; P[11] = nd->Z[m2];
			ce = new CCeBasic;
			ce->get_facet_directly(P, 4, NULL_STATE, OK_STATE, NULL_STATE);
			i = bar_add(ce);
			A = abs_point(P,   P+3);
			B = abs_point(P+3, P+6);
			E = abs_point(P,   P+6); R = (A+B+E)*.5;
			C = abs_point(P+6, P+9);
			D = abs_point(P+9, P);   Q = (C+D+E)*.5;
			SetFacetParam(i-1, sqrt(fabs(R*(R-A)*(R-B)*(R-E)))+sqrt(fabs(Q*(Q-C)*(Q-D)*(Q-E))), pp);
/////////////
//...face Z;
			P[0] = nd->X[m5]; P[1]  = nd->Y[m5]; P[2]  = nd->Z[m5];
			P[3] = nd->X[m6]; P[4]  = nd->Y[m6]; P[5]  = nd->Z[m6];
			P[6] = nd->X[m7]; P[7]  = nd->Y[m7]; P[8]  = nd->Z[m7];
			P[9] = nd->X[m8]; P[10] = nd->Y[m8]; P[11] = nd->Z[m8];
			ce = new CCeBasic;
			ce->get_facet_directly(P, 4, NULL_STATE, OK_STATE, NULL_STATE);
			i = bar_add(ce);
			A = abs_point(P,   P+3);
			B = abs_point(P+3, P+6);
			E = abs_point(P,   P+6); R = (A+B+E)*.5;
			C = abs_point(P+6, P+9);
			D = abs_point(P+9, P);   Q = (C+D+E)*.5;
			SetFacetParam(i-1, sqrt(fabs(R*(R-A)*(R-B)*(R-E)))+sqrt(fabs(Q*(Q-C)*(Q-D)*(Q-E))), pp);
		}  break;
		case GL_PENTA: if (nd->geom[i+1] == 8) {
			m1   = nd->geom[i+4];
			m2   = nd->geom[i+5];
			m3   = nd->geom[i+6];
			m4   = nd->geom[i+7];
			m5   = nd->geom[i+8];
			m6   = nd->geom[i+9];
////////////////////
//...face lateral 1;
			P[0] = nd->X[m1]; P[1]  = nd->Y[m1]; P[2]  = nd->Z[m1];
			P[3] = nd->X[m2]; P[4]  = nd->Y[m2]; P[5]  = nd->Z[m2];
			P[6] = nd->X[m5]; P[7]  = nd->Y[m5]; P[8]  = nd->Z[m5];
			P[9] = nd->X[m4]; P[10] = nd->Y[m4]; P[11] = nd->Z[m4];
			ce = new CCeBasic;
			ce->get_facet_directly(P, 4, NULL_STATE, OK_STATE, NULL_STATE);
			i = bar_add(ce);
			A = abs_point(P,   P+3);
			B = abs_point(P+3, P+6);
			E = abs_point(P,   P+6); R = (A+B+E)*.5;
			C = abs_point(P+6, P+9);
			D = abs_point(P+9, P);   Q = (C+D+E)*.5;
			SetFacetParam(i-1, sqrt(fabs(R*(R-A)*(R-B)*(R-E)))+sqrt(fabs(Q*(Q-C)*(Q-D)*(Q-E))), pp);
////////////////////
//...face lateral 2;
			P[0] = nd->X[m2]; P[1]  = nd->Y[m2]; P[2]  = nd->Z[m2];
			P[3] = nd->X[m3]; P[4]  = nd->Y[m3]; P[5]  = nd->Z[m3];
			P[6] = nd->X[m6]; P[7]  = nd->Y[m6]; P[8]  = nd->Z[m6];
			P[9] = nd->X[m5]; P[10] = nd->Y[m5]; P[11] = nd->Z[m5];
			ce = new CCeBasic;
			ce->get_facet_directly(P, 4, NULL_STATE, OK_STATE, NULL_STATE);
			i = bar_add(ce);
			A = abs_point(P,   P+3);
			B = abs_point(P+3, P+6);
			E = abs_point(P,   P+6); R = (A+B+E)*.5;
			C = abs_point(P+6, P+9);
			D = abs_point(P+9, P);   Q = (C+D+E)*.5;
			SetFacetParam(i-1, sqrt(fabs(R*(R-A)*(R-B)*(R-E)))+sqrt(fabs(Q*(Q-C)*(Q-D)*(Q-E))), pp);
////////////////////
//...face lateral 3;
			P[0] = nd->X[m3]; P[1]  = nd->Y[m3]; P[2]  = nd->Z[m3];
			P[3] = nd->X[m1]; P[4]  = nd->Y[m1]; P[5]  = nd->Z[m1];
			P[6] = nd->X[m4]; P[7]  = nd->Y[m4]; P[8]  = nd->Z[m4];
			P[9] = nd->X[m6]; P[10] = nd->Y[m6]; P[11] = nd->Z[m6];
			ce = new CCeBasic;
			ce->get_facet_directly(P, 4, NULL_STATE, OK_STATE, NULL_STATE);
			i = bar_add(ce);
			A = abs_point(P,   P+3);
			B = abs_point(P+3, P+6);
			E = abs_point(P,   P+6); R = (A+B+E)*.5;
			C = abs_point(P+6, P+9);
			D = abs_point(P+9, P);   Q = (C+D+E)*.5;
			SetFacetParam(i-1, sqrt(fabs(R*(R-A)*(R-B)*(R-E)))+sqrt(fabs(Q*(Q-C)*(Q-D)*(Q-E))), pp);
/////////////
//...face Z';
			P[0] = nd->X[m1]; P[1] = nd->Y[m1]; P[2] = nd->Z[m1];
			P[3] = nd->X[m3]; P[4] = nd->Y[m3]; P[5] = nd->Z[m3];
			P[6] = nd->X[m2]; P[7] = nd->Y[m2]; P[8] = nd->Z[m2];
			ce = new CCeBasic;
			ce->get_facet_directly(P, 3, NULL_STATE, OK_STATE, NULL_STATE);
			i = bar_add(ce);
			A = abs_point(P,   P+3);
			B = abs_point(P+3, P+6);
			E = abs_point(P,   P+6); R = (A+B+E)*.5;
			SetFacetParam(i-1, sqrt(fabs(R*(R-A)*(R-B)*(R-E))), pp);
/////////////
//...face Z;
			P[0] = nd->X[m4]; P[1] = nd->Y[m4]; P[2] = nd->Z[m4];
			P[3] = nd->X[m5]; P[4] = nd->Y[m5]; P[5] = nd->Z[m5];
			P[6] = nd->X[m6]; P[7] = nd->Y[m6]; P[8] = nd->Z[m6];
			ce = new CCeBasic;
			ce->get_facet_directly(P, 3, NULL_STATE, OK_STATE, NULL_STATE);
			i = bar_add(ce);
			A = abs_point(P,   P+3);
			B = abs_point(P+3, P+6);
			E = abs_point(P,   P+6); R = (A+B+E)*.5;
			SetFacetParam(i-1, sqrt(fabs(R*(R-A)*(R-B)*(R-E))), pp);
		}  break;
		case GL_TETRA: if (nd->geom[i+1] == 6 || nd->geom[i+1] == 12) {
			m1   = nd->geom[i+4];
			m2   = nd->geom[i+5];
			m3   = nd->geom[i+6];
			m4   = nd->geom[i+7];
///////////////
//...lateral 1;
			P[0] = nd->X[m1]; P[1] = nd->Y[m1]; P[2] = nd->Z[m1];
			P[3] = nd->X[m3]; P[4] = nd->Y[m3]; P[5] = nd->Z[m3];
			P[6] = nd->X[m2]; P[7] = nd->Y[m2]; P[8] = nd->Z[m2];
			ce = new CCeBasic;
			ce->get_facet_directly(P, 3, NULL_STATE, OK_STATE, NULL_STATE);
			i = bar_add(ce);
			A = abs_point(P,   P+3);
			B = abs_point(P+3, P+6);
			E = abs_point(P,   P+6); R = (A+B+E)*.5;
			SetFacetParam(i-1, sqrt(fabs(R*(R-A)*(R-B)*(R-E))), pp);
///////////////
//...lateral 2;
			P[0] = nd->X[m2]; P[1] = nd->Y[m2]; P[2] = nd->Z[m2];
			P[3] = nd->X[m3]; P[4] = nd->Y[m3]; P[5] = nd->Z[m3];
			P[6] = nd->X[m4]; P[7] = nd->Y[m4]; P[8] = nd->Z[m4];
			ce = new CCeBasic;
			ce->get_facet_directly(P, 3, NULL_STATE, OK_STATE, NULL_STATE);
			i = bar_add(ce);
			A = abs_point(P,   P+3);
			B = abs_point(P+3, P+6);
			E = abs_point(P,   P+6); R = (A+B+E)*.5;
			SetFacetParam(i-1, sqrt(fabs(R*(R-A)*(R-B)*(R-E))), pp);
///////////////
//...lateral 3;
			P[0] = nd->X[m3]; P[1] = nd->Y[m3]; P[2] = nd->Z[m3];
			P[3] = nd->X[m1]; P[4] = nd->Y[m1]; P[5] = nd->Z[m1];
			P[6] = nd->X[m4]; P[7] = nd->Y[m4]; P[8] = nd->Z[m4];
			ce = new CCeBasic;
			ce->get_facet_directly(P, 3, NULL_STATE, OK_STATE, NULL_STATE);
			i = bar_add(ce);
			A = abs_point(P,   P+3);
			B = abs_point(P+3, P+6);
			E = abs_point(P,   P+6); R = (A+B+E)*.5;
			SetFacetParam(i-1, sqrt(fabs(R*(R-A)*(R-B)*(R-E))), pp);
///////////////
//...lateral 4;
			P[0] = nd->X[m1]; P[1] = nd->Y[m1]; P[2] = nd->Z[m1];
			P[3] = nd->X[m2]; P[4] = nd->Y[m2]; P[5] = nd->Z[m2];
			P[6] = nd->X[m4]; P[7] = nd->Y[m4]; P[8] = nd->Z[m4];
			ce = new CCeBasic;
			ce->get_facet_directly(P, 3, NULL_STATE, OK_STATE, NULL_STATE);
			i = bar_add(ce);
			A = abs_point(P,   P+3);
			B = abs_point(P+3, P+6);
			E = abs_point(P,   P+6); R = (A+B+E)*.5;
			SetFacetParam(i-1, sqrt(fabs(R*(R-A)*(R-B)*(R-E))), pp);
		}  break;
		case GL_PYRAMID: if (nd->geom[i+1] == 7) {
			m1   = nd->geom[i+4];
			m2   = nd->geom[i+5];
			m3   = nd->geom[i+6];
			m4   = nd->geom[i+7];
			m5   = nd->geom[i+8];
///////////////
//...lateral 1;
			P[0] = nd->X[m1]; P[1] = nd->Y[m1]; P[2] = nd->Z[m1];
			P[3] = nd->X[m2]; P[4] = nd->Y[m2]; P[5] = nd->Z[m2];
			P[6] = nd->X[m5]; P[7] = nd->Y[m5]; P[8] = nd->Z[m5];
			ce = new CCeBasic;
			ce->get_facet_directly(P, 3, NULL_STATE, OK_STATE, NULL_STATE);
			i = bar_add(ce);
			A = abs_point(P,   P+3);
			B = abs_point(P+3, P+6);
			E = abs_point(P,   P+6); R = (A+B+E)*.5;
			SetFacetParam(i-1, sqrt(fabs(R*(R-A)*(R-B)*(R-E))), pp);
///////////////
//...lateral 2;
			P[0] = nd->X[m2]; P[1] = nd->Y[m2]; P[2] = nd->Z[m2];
			P[3] = nd->X[m3]; P[4] = nd->Y[m3]; P[5] = nd->Z[m3];
			P[6] = nd->X[m5]; P[7] = nd->Y[m5]; P[8] = nd->Z[m5];
			ce = new CCeBasic;
			ce->get_facet_directly(P, 3, NULL_STATE, OK_STATE, NULL_STATE);
			i = bar_add(ce);
			A = abs_point(P,   P+3);
			B = abs_point(P+3, P+6);
			E = abs_point(P,   P+6); R = (A+B+E)*.5;
			SetFacetParam(i-1, sqrt(fabs(R*(R-A)*(R-B)*(R-E))), pp);
///////////////
//...lateral 3;
			P[0] = nd->X[m3]; P[1] = nd->Y[m3]; P[2] = nd->Z[m3];
			P[3] = nd->X[m4]; P[4] = nd->Y[m4]; P[5] = nd->Z[m4];
			P[6] = nd->X[m5]; P[7] = nd->Y[m5]; P[8] = nd->Z[m5];
			ce = new CCeBasic;
			ce->get_facet_directly(P, 3, NULL_STATE, OK_STATE, NULL_STATE);
			i = bar_add(ce);
			A = abs_point(P,   P+3);
			B = abs_point(P+3, P+6);
			E = abs_point(P,   P+6); R = (A+B+E)*.5;
			SetFacetParam(i-1, sqrt(fabs(R*(R-A)*(R-B)*(R-E))), pp);
///////////////
//...lateral 4;
			P[0] = nd->X[m4]; P[1] = nd->Y[m4]; P[2] = nd->Z[m4];
			P[3] = nd->X[m1]; P[4] = nd->Y[m1]; P[5] = nd->Z[m1];
			P[6] = nd->X[m5]; P[7] = nd->Y[m5]; P[8] = nd->Z[m5];
			ce = new CCeBasic;
			ce->get_facet_directly(P, 3, NULL_STATE, OK_STATE, NULL_STATE);
			i = bar_add(ce);
			A = abs_point(P,   P+3);
			B = abs_point(P+3, P+6);
			E = abs_point(P,   P+6); R = (A+B+E)*.5;
			SetFacetParam(i-1, sqrt(fabs(R*(R-A)*(R-B)*(R-E))), pp);
/////////////
//...basis 5;
			P[0] = nd->X[m1]; P[1] =  nd->Y[m1]; P[2] =  nd->Z[m1];
			P[3] = nd->X[m4]; P[4] =  nd->Y[m4]; P[5] =  nd->Z[m4];
			P[6] = nd->X[m3]; P[7] =  nd->Y[m3]; P[8] =  nd->Z[m3];
			P[9] = nd->X[m2]; P[10] = nd->Y[m2]; P[11] = nd->Z[m2];
			ce = new CCeBasic;
			ce->get_facet_directly(P, 4, NULL_STATE, OK_STATE, NULL_STATE);
			i = bar_add(ce);
			A = abs_point(P,   P+3);
			B = abs_point(P+3, P+6);
			E = abs_point(P,   P+6); R = (A+B+E)*.5;
			C = abs_point(P+6, P+9);
			D = abs_point(P+9, P);   Q = (C+D+E)*.5;
			SetFacetParam(i-1, sqrt(fabs(R*(R-A)*(R-B)*(R-E)))+sqrt(fabs(Q*(Q-C)*(Q-D)*(Q-E))), pp);
		}  break;
	}
	bar_ord();
}

/////////////////////////////////////////////////////////////
//...setting boundary conditions on the conjunction elements;
void CCeBasic::set_nd_bar_condit(int k, int id_property, double * pp, int * id_pp, double S)
{
   int j, l;
   if (ce && graph)
   for (j = 0; j < graph[0]; j++)
   if  (ce[j]->mp && ID_MAP(2, NULL_GENUS)  == ce[j]->mp[0] &&
        ce[j]->mp[l = size_of_map(2, NULL_GENUS)] == (CMap)FACET_CELL &&
        ce[j]->mp[++l] == (CMap)(-k-1)) {

        SetFacetParam(j, S, pp+6*(id_property-1), id_pp[2*id_property]);
        return;
   }
}

//////////////////////////////////////////////////////////////////////////////////
//...setting boundary conditions on the surface of the plane composition geometry;
int CCeBasic::set_nd_bar_condit(int k, CGrid * nd, int k1, int k2, int id_property, double * pp, int * id_pp, double & S, int id_facet)
{
	if (! nd->geom || ! nd->geom_ptr) return(0);
	int m, m1, m2, m3, m4, m5, m6, m7, m8, i, j;
	double P[12], A, B, C, D, E, R, Q; 
  
////////////////////////////////////
//...setting new boundary condition;
	switch (nd->geom[i = nd->geom_ptr[k]]) {
		case GL_BOXS: if (nd->geom[i+1] == 10) {
			m1 = nd->geom[i+4];
			m2 = nd->geom[i+5];
			m3 = nd->geom[i+6];
			m4 = nd->geom[i+7];
			m5 = nd->geom[i+8];
			m6 = nd->geom[i+9];
			m7 = nd->geom[i+10];
			m8 = nd->geom[i+11];
/////////////
//...face X';
			if (id_facet == 0 || geom_link_id4(m1, m5, m8, m4, k1, k2)) {
				m = id_property; //...property number;
				j = (int)(S = ce[0]->mp[size_of_map(2, NULL_GENUS)+1]); //...conjugate facet;
				if (j >= 0) SetFacetParam(0, S, pp+6*(m-1), id_pp[2*m]);
				else { //...устанавливаем граничное условие на разорванную связь;
					P[0] = nd->X[m1]; P[1]  = nd->Y[m1]; P[2]  = nd->Z[m1];
					P[3] = nd->X[m5]; P[4]  = nd->Y[m5]; P[5]  = nd->Z[m5];
					P[6] = nd->X[m8]; P[7]  = nd->Y[m8]; P[8]  = nd->Z[m8];
					P[9] = nd->X[m4]; P[10] = nd->Y[m4]; P[11] = nd->Z[m4];
					A = abs_point(P,   P+3);
					B = abs_point(P+3, P+6);
					E = abs_point(P,   P+6); R = (A+B+E)*.5;
					C = abs_point(P+6, P+9);
					D = abs_point(P+9, P);   Q = (C+D+E)*.5;
					SetFacetParam(0, S = sqrt(fabs(R*(R-A)*(R-B)*(R-E)))+sqrt(fabs(Q*(Q-C)*(Q-D)*(Q-E))), pp+6*(m-1), id_pp[2*m]);
				}
				if (j < 0) return(j);
				else       return(0);  
			}
/////////////
//...face X;
			if (id_facet == 1 || geom_link_id4(m2, m3, m7, m6, k1, k2)) {
				m = id_property; //...property number;
				j = (int)(S = ce[1]->mp[size_of_map(2, NULL_GENUS)+1]); //...conjugate facet;
				if (j >= 0) SetFacetParam(1, S, pp+6*(m-1), id_pp[2*m]);
				else { //...устанавливаем граничное условие на разорванную связь;
					P[0] = nd->X[m2]; P[1]  = nd->Y[m2]; P[2]  = nd->Z[m2];
					P[3] = nd->X[m3]; P[4]  = nd->Y[m3]; P[5]  = nd->Z[m3];
					P[6] = nd->X[m7]; P[7]  = nd->Y[m7]; P[8]  = nd->Z[m7];
					P[9] = nd->X[m6]; P[10] = nd->Y[m6]; P[11] = nd->Z[m6];
					A = abs_point(P,   P+3);
					B = abs_point(P+3, P+6);
					E = abs_point(P,   P+6); R = (A+B+E)*.5;
					C = abs_point(P+6, P+9);
					D = abs_point(P+9, P);   Q = (C+D+E)*.5;
					SetFacetParam(1, S = sqrt(fabs(R*(R-A)*(R-B)*(R-E)))+sqrt(fabs(Q*(Q-C)*(Q-D)*(Q-E))), pp+6*(m-1), id_pp[2*m]);
				}
				if (j < 0) return(j);
				else       return(0);  
			}
/////////////
//...face Y';
			if (id_facet == 2 || geom_link_id4(m1, m2, m6, m5, k1, k2)) {
				m = id_property; //...property number;
				j = (int)(S = ce[2]->mp[size_of_map(2, NULL_GENUS)+1]); //...conjugate facet;
				if (j >= 0) SetFacetParam(2, S, pp+6*(m-1), id_pp[2*m]);
				else { //...устанавливаем граничное условие на разорванную связь;
					P[0] = nd->X[m1]; P[1]  = nd->Y[m1]; P[2]  = nd->Z[m1];
					P[3] = nd->X[m2]; P[4]  = nd->Y[m2]; P[5]  = nd->Z[m2];
					P[6] = nd->X[m6]; P[7]  = nd->Y[m6]; P[8]  = nd->Z[m6];
					P[9] = nd->X[m5]; P[10] = nd->Y[m5]; P[11] = nd->Z[m5];
					A = abs_point(P,   P+3);
					B = abs_point(P+3, P+6);
					E = abs_point(P,   P+6); R = (A+B+E)*.5;
					C = abs_point(P+6, P+9);
					D = abs_point(P+9, P);   Q = (C+D+E)*.5;
					SetFacetParam(2, S = sqrt(fabs(R*(R-A)*(R-B)*(R-E)))+sqrt(fabs(Q*(Q-C)*(Q-D)*(Q-E))), pp+6*(m-1), id_pp[2*m]);
				}
				if (j < 0) return(j);
				else       return(0);  
			}
/////////////
//...face Y;
			if (id_facet == 3 || geom_link_id4(m4, m8, m7, m3, k1, k2)) {
				m = id_property; //...property number;
				j = (int)(S = ce[3]->mp[size_of_map(2, NULL_GENUS)+1]); //...conjugate facet;
				if (j >= 0) SetFacetParam(3, S, pp+6*(m-1), id_pp[2*m]);
				else { //...устанавливаем граничное условие на разорванную связь;
					P[0] = nd->X[m4]; P[1]  = nd->Y[m4]; P[2]  = nd->Z[m4];
					P[3] = nd->X[m8]; P[4]  = nd->Y[m8]; P[5]  = nd->Z[m8];
					P[6] = nd->X[m7]; P[7]  = nd->Y[m7]; P[8]  = nd->Z[m7];
					P[9] = nd->X[m3]; P[10] = nd->Y[m3]; P[11] = nd->Z[m3];
					A = abs_point(P,   P+3);
					B = abs_point(P+3, P+6);
					E = abs_point(P,   P+6); R = (A+B+E)*.5;
					C = abs_point(P+6, P+9);
					D = abs_point(P+9, P);   Q = (C+D+E)*.5;
					SetFacetParam(3, S = sqrt(fabs(R*(R-A)*(R-B)*(R-E)))+sqrt(fabs(Q*(Q-C)*(Q-D)*(Q-E))), pp+6*(m-1), id_pp[2*m]);
				}
				if (j < 0) return(j);
				else       return(0);  
			}
/////////////
//...face Z';
			if (id_facet == 4 || geom_link_id4(m1, m4, m3, m2, k1, k2)) {
				m = id_property; //...property number;
				j = (int)(S = ce[4]->mp[size_of_map(2, NULL_GENUS)+1]); //...conjugate facet;
				if (j >= 0) SetFacetParam(4, S, pp+6*(m-1), id_pp[2*m]);
				else { //...устанавливаем граничное условие на разорванную связь;
					P[0] = nd->X[m1]; P[1]  = nd->Y[m1]; P[2]  = nd->Z[m1];
					P[3] = nd->X[m4]; P[4]  = nd->Y[m4]; P[5]  = nd->Z[m4];
					P[6] = nd->X[m3]; P[7]  = nd->Y[m3]; P[8]  = nd->Z[m3];
					P[9] = nd->X[m2]; P[10] = nd->Y[m2]; P[11] = nd->Z[m2];
					A = abs_point(P,   P+3);
					B = abs_point(P+3, P+6);
					E = abs_point(P,   P+6); R = (A+B+E)*.5;
					C = abs_point(P+6, P+9);
					D = abs_point(P+9, P);   Q = (C+D+E)*.5;
					SetFacetParam(4, S = sqrt(fabs(R*(R-A)*(R-B)*(R-E)))+sqrt(fabs(Q*(Q-C)*(Q-D)*(Q-E))), pp+6*(m-1), id_pp[2*m]);
				}
				if (j < 0) return(j);
				else       return(0);  
			}
/////////////
//...face Z;
			if (id_facet == 5 || geom_link_id4(m5, m6, m7, m8, k1, k2)) {
				m = id_property; //...property number;
				j = (int)(S = ce[5]->mp[size_of_map(2, NULL_GENUS)+1]); //...conjugate facet;
				if (j >= 0) SetFacetParam(5, S, pp+6*(m-1), id_pp[2*m]);
				else { //...устанавливаем граничное условие на разорванную связь;
					P[0] = nd->X[m5]; P[1]  = nd->Y[m5]; P[2]  = nd->Z[m5];
					P[3] = nd->X[m6]; P[4]  = nd->Y[m6]; P[5]  = nd->Z[m6];
					P[6] = nd->X[m7]; P[7]  = nd->Y[m7]; P[8]  = nd->Z[m7];
					P[9] = nd->X[m8]; P[10] = nd->Y[m8]; P[11] = nd->Z[m8];
					A = abs_point(P,   P+3);
					B = abs_point(P+3, P+6);
					E = abs_point(P,   P+6); R = (A+B+E)*.5;
					C = abs_point(P+6, P+9);
					D = abs_point(P+9, P);   Q = (C+D+E)*.5;
					SetFacetParam(5, S = sqrt(fabs(R*(R-A)*(R-B)*(R-E)))+sqrt(fabs(Q*(Q-C)*(Q-D)*(Q-E))), pp+6*(m-1), id_pp[2*m]);
				}
				if (j < 0) return(j);
				else       return(0);  
			}
		}  break;
		case GL_PENTA: if (nd->geom[i+1] == 8) {
			m1 = nd->geom[i+4];
			m2 = nd->geom[i+5];
			m3 = nd->geom[i+6];
			m4 = nd->geom[i+7];
			m5 = nd->geom[i+8];
			m6 = nd->geom[i+9];
////////////////////
//...face lateral 1;
			if (id_facet == 0 || geom_link_id4(m1, m2, m5, m4, k1, k2)) {
				m = id_property; //...property number;
				j = (int)(S = ce[0]->mp[size_of_map(2, NULL_GENUS)+1]); //...conjugate facet;
				if (j >= 0) SetFacetParam(0, S, pp+6*(m-1), id_pp[2*m]);
				else { //...устанавливаем граничное условие на разорванную связь;
					P[0] = nd->X[m1]; P[1]  = nd->Y[m1]; P[2]  = nd->Z[m1];
					P[3] = nd->X[m2]; P[4]  = nd->Y[m2]; P[5]  = nd->Z[m2];
					P[6] = nd->X[m5]; P[7]  = nd->Y[m5]; P[8]  = nd->Z[m5];
					P[9] = nd->X[m4]; P[10] = nd->Y[m4]; P[11] = nd->Z[m4];
					A = abs_point(P,   P+3);
					B = abs_point(P+3, P+6);
					E = abs_point(P,   P+6); R = (A+B+E)*.5;
					C = abs_point(P+6, P+9);
					D = abs_point(P+9, P);   Q = (C+D+E)*.5;
					SetFacetParam(0, S = sqrt(fabs(R*(R-A)*(R-B)*(R-E)))+sqrt(fabs(Q*(Q-C)*(Q-D)*(Q-E))), pp+6*(m-1), id_pp[2*m]);
				}
				if (j < 0) return(j);
				else       return(0);  
			}
////////////////////
//...face lateral 2;
			if (id_facet == 1 || geom_link_id4(m2, m3, m6, m5, k1, k2)) {
				m = id_property; //...property number;
				j = (int)(S = ce[1]->mp[size_of_map(2, NULL_GENUS)+1]); //...conjugate facet;
				if (j >= 0) SetFacetParam(1, S, pp+6*(m-1), id_pp[2*m]);
				else { //...устанавливаем граничное условие на разорванную связь;
					P[0] = nd->X[m2]; P[1]  = nd->Y[m2]; P[2]  = nd->Z[m2];
					P[3] = nd->X[m3]; P[4]  = nd->Y[m3]; P[5]  = nd->Z[m3];
					P[6] = nd->X[m6]; P[7]  = nd->Y[m6]; P[8]  = nd->Z[m6];
					P[9] = nd->X[m5]; P[10] = nd->Y[m5]; P[11] = nd->Z[m5];
					A = abs_point(P,   P+3);
					B = abs_point(P+3, P+6);
					E = abs_point(P,   P+6); R = (A+B+E)*.5;
					C = abs_point(P+6, P+9);
					D = abs_point(P+9, P);   Q = (C+D+E)*.5;
					SetFacetParam(1, S = sqrt(fabs(R*(R-A)*(R-B)*(R-E)))+sqrt(fabs(Q*(Q-C)*(Q-D)*(Q-E))), pp+6*(m-1), id_pp[2*m]);
				}
				if (j < 0) return(j);
				else       return(0);  
			}
////////////////////
//...face lateral 3;
			if (id_facet == 2 || geom_link_id4(m3, m1, m4, m6, k1, k2)) {
				m = id_property; //...property number;
				j = (int)(S = ce[2]->mp[size_of_map(2, NULL_GENUS)+1]); //...conjugate facet;
				if (j >= 0) SetFacetParam(2, S, pp+6*(m-1), id_pp[2*m]);
				else { //...устанавливаем граничное условие на разорванную связь;
					P[0] = nd->X[m3]; P[1]  = nd->Y[m3]; P[2]  = nd->Z[m3];
					P[3] = nd->X[m1]; P[4]  = nd->Y[m1]; P[5]  = nd->Z[m1];
					P[6] = nd->X[m4]; P[7]  = nd->Y[m4]; P[8]  = nd->Z[m4];
					P[9] = nd->X[m6]; P[10] = nd->Y[m6]; P[11] = nd->Z[m6];
					A = abs_point(P,   P+3);
					B = abs_point(P+3, P+6);
					E = abs_point(P,   P+6); R = (A+B+E)*.5;
					C = abs_point(P+6, P+9);
					D = abs_point(P+9, P);   Q = (C+D+E)*.5;
					SetFacetParam(2, S = sqrt(fabs(R*(R-A)*(R-B)*(R-E)))+sqrt(fabs(Q*(Q-C)*(Q-D)*(Q-E))), pp+6*(m-1), id_pp[2*m]);
				}
				if (j < 0) return(j);
				else       return(0);  
			}
/////////////
//...face Z';
			if (id_facet == 3 || geom_link_id3(m1, m3, m2, k1)) {
				m = id_property; //...property number;
				j = (int)(S = ce[3]->mp[size_of_map(2, NULL_GENUS)+1]); //...conjugate facet;
				if (j >= 0) SetFacetParam(3, S, pp+6*(m-1), id_pp[2*m]);
				else { //...устанавливаем граничное условие на разорванную связь;
					P[0] = nd->X[m1]; P[1] = nd->Y[m1]; P[2] = nd->Z[m1];
					P[3] = nd->X[m3]; P[4] = nd->Y[m3]; P[5] = nd->Z[m3];
					P[6] = nd->X[m2]; P[7] = nd->Y[m2]; P[8] = nd->Z[m2];
					A = abs_point(P,   P+3);
					B = abs_point(P+3, P+6);
					E = abs_point(P,   P+6); R = (A+B+E)*.5;
					SetFacetParam(3, S = sqrt(fabs(R*(R-A)*(R-B)*(R-E))), pp+6*(m-1), id_pp[2*m]);
				}
				if (j < 0) return(j);
				else       return(0);  
			}
/////////////
//...face Z;
			if (id_facet == 4 || geom_link_id3(m4, m5, m6, k1)) {
				m = id_property; //...property number;
				j = (int)(S = ce[4]->mp[size_of_map(2, NULL_GENUS)+1]); //...conjugate facet;
				if (j >= 0) SetFacetParam(4, S, pp+6*(m-1), id_pp[2*m]);
				else { //...устанавливаем граничное условие на разорванную связь;
					P[0] = nd->X[m4]; P[1] = nd->Y[m4]; P[2] = nd->Z[m4];
					P[3] = nd->X[m5]; P[4] = nd->Y[m5]; P[5] = nd->Z[m5];
					P[6] = nd->X[m6]; P[7] = nd->Y[m6]; P[8] = nd->Z[m6];
					A = abs_point(P,   P+3);
					B = abs_point(P+3, P+6);
					E = abs_point(P,   P+6); R = (A+B+E)*.5;
					SetFacetParam(4, S = sqrt(fabs(R*(R-A)*(R-B)*(R-E))), pp+6*(m-1), id_pp[2*m]);
				}
				if (j < 0) return(j);
				else       return(0);  
			}
		}  break;
		case GL_TETRA: if (nd->geom[i+1] == 6) {
			m1 = nd->geom[i+4];
			m2 = nd->geom[i+5];
			m3 = nd->geom[i+6];
			m4 = nd->geom[i+7];
///////////////
//...lateral 1;
			if (id_facet == 0 || k2 == m4 && (k1 == m1 || k1 == m2 || k1 == m3)) {
				m = id_property; //...property number;
				j = (int)(S = ce[0]->mp[size_of_map(2, NULL_GENUS)+1]); //...conjugate facet;
				if (j >= 0) SetFacetParam(0, S, pp+6*(m-1), id_pp[2*m]);
				else { //...устанавливаем граничное условие на разорванную связь;
					P[0] = nd->X[m1]; P[1] = nd->Y[m1]; P[2] = nd->Z[m1];
					P[3] = nd->X[m3]; P[4] = nd->Y[m3]; P[5] = nd->Z[m3];
					P[6] = nd->X[m2]; P[7] = nd->Y[m2]; P[8] = nd->Z[m2];
					A = abs_point(P,   P+3);
					B = abs_point(P+3, P+6);
					E = abs_point(P,   P+6); R = (A+B+E)*.5;
					SetFacetParam(0, S = sqrt(fabs(R*(R-A)*(R-B)*(R-E))), pp+6*(m-1), id_pp[2*m]);
				}
				if (j < 0) return(j);
				else       return(0);  
			}
///////////////
//...lateral 2;
			if (id_facet == 1 || k2 == m1 && (k1 == m2 || k1 == m3 || k1 == m4)) {
				m = id_property; //...property number;
				j = (int)(S = ce[1]->mp[size_of_map(2, NULL_GENUS)+1]); //...conjugate facet;
				if (j >= 0) SetFacetParam(1, S, pp+6*(m-1), id_pp[2*m]);
				else { //...устанавливаем граничное условие на разорванную связь;
					P[0] = nd->X[m2]; P[1] = nd->Y[m2]; P[2] = nd->Z[m2];
					P[3] = nd->X[m3]; P[4] = nd->Y[m3]; P[5] = nd->Z[m3];
					P[6] = nd->X[m4]; P[7] = nd->Y[m4]; P[8] = nd->Z[m4];
					A = abs_point(P,   P+3);
					B = abs_point(P+3, P+6);
					E = abs_point(P,   P+6); R = (A+B+E)*.5;
					SetFacetParam(1, S = sqrt(fabs(R*(R-A)*(R-B)*(R-E))), pp+6*(m-1), id_pp[2*m]);
				}
				if (j < 0) return(j);
				else       return(0);  
			}
///////////////
//...lateral 3;
			if (id_facet == 2 || k2 == m2 && (k1 == m1 || k1 == m3 || k1 == m4)) {
				m = id_property; //...property number;
				j = (int)(S = ce[2]->mp[size_of_map(2, NULL_GENUS)+1]); //...conjugate facet;
				if (j >= 0) SetFacetParam(2, S, pp+6*(m-1), id_pp[2*m]);
				else { //...устанавливаем граничное условие на разорванную связь;
					P[0] = nd->X[m3]; P[1] = nd->Y[m3]; P[2] = nd->Z[m3];
					P[3] = nd->X[m1]; P[4] = nd->Y[m1]; P[5] = nd->Z[m1];
					P[6] = nd->X[m4]; P[7] = nd->Y[m4]; P[8] = nd->Z[m4];
					A = abs_point(P,   P+3);
					B = abs_point(P+3, P+6);
					E = abs_point(P,   P+6); R = (A+B+E)*.5;
					SetFacetParam(2, S = sqrt(fabs(R*(R-A)*(R-B)*(R-E))), pp+6*(m-1), id_pp[2*m]);
				}
				if (j < 0) return(j);
				else       return(0);  
			}
///////////////
//...lateral 4;
			if (id_facet == 3 || k2 == m3 && (k1 == m1 || k1 == m2 || k1 == m4)) {
				m = id_property; //...property number;
				j = (int)(S = ce[3]->mp[size_of_map(2, NULL_GENUS)+1]); //...conjugate facet;
				if (j >= 0) SetFacetParam(3, S, pp+6*(m-1), id_pp[2*m]);
				else { //...устанавливаем граничное условие на разорванную связь;
					P[0] = nd->X[m1]; P[1] = nd->Y[m1]; P[2] = nd->Z[m1];
					P[3] = nd->X[m2]; P[4] = nd->Y[m2]; P[5] = nd->Z[m2];
					P[6] = nd->X[m4]; P[7] = nd->Y[m4]; P[8] = nd->Z[m4];
					A = abs_point(P,   P+3);
					B = abs_point(P+3, P+6);
					E = abs_point(P,   P+6); R = (A+B+E)*.5;
					SetFacetParam(3, S = sqrt(fabs(R*(R-A)*(R-B)*(R-E))), pp+6*(m-1), id_pp[2*m]);
				}
				if (j < 0) return(j);
				else       return(0);  
			}
		}  break;
		case GL_PYRAMID: if (nd->geom[i+1] == 7) {
			m1 = nd->geom[i+4];
			m2 = nd->geom[i+5];
			m3 = nd->geom[i+6];
			m4 = nd->geom[i+7];
			m5 = nd->geom[i+8];
///////////////
//...lateral 1;
			if (id_facet == 0 || geom_link_id3(m1, m2, m5, k1)) { //...не проверена идентификация?????
				m = id_property; //...property number;
				j = (int)(S = ce[0]->mp[size_of_map(2, NULL_GENUS)+1]); //...conjugate facet;
				if (j >= 0) SetFacetParam(0, S, pp+6*(m-1), id_pp[2*m]);
				else { //...устанавливаем граничное условие на разорванную связь;
					P[0] = nd->X[m1]; P[1] = nd->Y[m1]; P[2] = nd->Z[m1];
					P[3] = nd->X[m2]; P[4] = nd->Y[m2]; P[5] = nd->Z[m2];
					P[6] = nd->X[m5]; P[7] = nd->Y[m5]; P[8] = nd->Z[m5];
					A = abs_point(P,   P+3);
					B = abs_point(P+3, P+6);
					E = abs_point(P,   P+6); R = (A+B+E)*.5;
					SetFacetParam(0, S = sqrt(fabs(R*(R-A)*(R-B)*(R-E))), pp+6*(m-1), id_pp[2*m]);
				}
				if (j < 0) return(j);
				else       return(0);  
			}
///////////////
//...lateral 2;
			if (id_facet == 1 || geom_link_id3(m2, m3, m5, k1)) { //...идентификация?????
				m = id_property; //...property number;
				j = (int)(S = ce[1]->mp[size_of_map(2, NULL_GENUS)+1]); //...conjugate facet;
				if (j >= 0) SetFacetParam(1, S, pp+6*(m-1), id_pp[2*m]);
				else { //...устанавливаем граничное условие на разорванную связь;
					P[0] = nd->X[m2]; P[1] = nd->Y[m2]; P[2] = nd->Z[m2];
					P[3] = nd->X[m3]; P[4] = nd->Y[m3]; P[5] = nd->Z[m3];
					P[6] = nd->X[m5]; P[7] = nd->Y[m5]; P[8] = nd->Z[m5];
					A = abs_point(P,   P+3);
					B = abs_point(P+3, P+6);
					E = abs_point(P,   P+6); R = (A+B+E)*.5;
					SetFacetParam(1, S = sqrt(fabs(R*(R-A)*(R-B)*(R-E))), pp+6*(m-1), id_pp[2*m]);
				}
				if (j < 0) return(j);
				else       return(0);  
			}
///////////////
//...lateral 3;
			if (id_facet == 2 || geom_link_id3(m3, m4, m5, k1)) { //...идентификация?????
				m = id_property; //...property number;
				j = (int)(S = ce[2]->mp[size_of_map(2, NULL_GENUS)+1]); //...conjugate facet;
				if (j >= 0) SetFacetParam(2, S, pp+6*(m-1), id_pp[2*m]);
				else { //...устанавливаем граничное условие на разорванную связь;
					P[0] = nd->X[m3]; P[1] = nd->Y[m3]; P[2] = nd->Z[m3];
					P[3] = nd->X[m4]; P[4] = nd->Y[m4]; P[5] = nd->Z[m4];
					P[6] = nd->X[m5]; P[7] = nd->Y[m5]; P[8] = nd->Z[m5];
					A = abs_point(P,   P+3);
					B = abs_point(P+3, P+6);
					E = abs_point(P,   P+6); R = (A+B+E)*.5;
					SetFacetParam(2, S = sqrt(fabs(R*(R-A)*(R-B)*(R-E))), pp+6*(m-1), id_pp[2*m]);
				}
				if (j < 0) return(j);
				else       return(0);  
			}
///////////////
//...lateral 4;
			if (id_facet == 3 || geom_link_id3(m4, m1, m5, k1)) { //...идентификация?????
				m = id_property; //...property number;
				j = (int)(S = ce[3]->mp[size_of_map(2, NULL_GENUS)+1]); //...conjugate facet;
				if (j >= 0) SetFacetParam(3, S, pp+6*(m-1), id_pp[2*m]);
				else { //...устанавливаем граничное условие на разорванную связь;
					P[0] = nd->X[m4]; P[1] = nd->Y[m4]; P[2] = nd->Z[m4];
					P[3] = nd->X[m1]; P[4] = nd->Y[m1]; P[5] = nd->Z[m1];
					P[6] = nd->X[m5]; P[7] = nd->Y[m5]; P[8] = nd->Z[m5];
					A = abs_point(P,   P+3);
					B = abs_point(P+3, P+6);
					E = abs_point(P,   P+6); R = (A+B+E)*.5;
					SetFacetParam(3, S = sqrt(fabs(R*(R-A)*(R-B)*(R-E))), pp+6*(m-1), id_pp[2*m]);
				}
				if (j < 0) return(j);
				else       return(0);  
			}
/////////////
//...basis 5;
			if (id_facet == 4 || geom_link_id4(m1, m4, m3, m2, k1, k2)) { //...идентификация?????
				m = id_property; //...property number;
				j = (int)(S = ce[4]->mp[size_of_map(2, NULL_GENUS)+1]); //...conjugate facet;
				if (j >= 0) SetFacetParam(4, S, pp+6*(m-1), id_pp[2*m]);
				else { //...устанавливаем граничное условие на разорванную связь;
					P[0] = nd->X[m1]; P[1] =  nd->Y[m1]; P[2] =  nd->Z[m1];
					P[3] = nd->X[m4]; P[4] =  nd->Y[m4]; P[5] =  nd->Z[m4];
					P[6] = nd->X[m3]; P[7] =  nd->Y[m3]; P[8] =  nd->Z[m3];
					P[9] = nd->X[m2]; P[10] = nd->Y[m2]; P[11] = nd->Z[m2];
					A = abs_point(P,   P+3);
					B = abs_point(P+3, P+6);
					E = abs_point(P,   P+6); R = (A+B+E)*.5;
					C = abs_point(P+6, P+9);
					D = abs_point(P+9, P);   Q = (C+D+E)*.5;
					SetFacetParam(4, S = sqrt(fabs(R*(R-A)*(R-B)*(R-E)))+sqrt(fabs(Q*(Q-C)*(Q-D)*(Q-E))), pp+6*(m-1), id_pp[2*m]);
				}
				if (j < 0) return(j);
				else       return(0);  
			}
		}  break;
	}
	return(0);
}

/*===================================================================*/
/*                   LIBRARY OF BASIC ELEMENTS                       */
/*===================================================================*/
//////////////////////////////////////
//...forming geometrical cell "point";
void CCeBasic::get_point(double * P)
{
	get_point(P[0], P[1], P[2]);
}
///////////////////////////////////////////////////
//...forming geometrical cell "point" (coordinate);
void CCeBasic::get_point(double X, double Y, double Z)
{
	int l;
	cells_new(1, 2, (l = size_of_map(0, NULL_GENUS))+1);
	if (mp) {
		mp[0] = ID_MAP(0, NULL_GENUS);
		mp[1] = (CMap)X;
		mp[2] = (CMap)Y;
		mp[3] = (CMap)Z;
		mp[l] = (CMap)NULL_CELL;
	}
}

///////////////////////////////////////
//...forming geometrical cell "circle";
void CCeBasic::get_circle(double R)
{
	int l;
	cells_new(1, 2, (l = size_of_map(1, SPHERE_GENUS))+1);
	if (mp) {
      mp[0] = ID_MAP(1, SPHERE_GENUS);
      mp[7] = (CMap)fabs(R);
      mp[8] = (CMap)M_PI;
      mp[l] = (CMap)NULL_CELL;
	}
}

/////////////////////////////////////
//...forming geometrical cell "disk";
void CCeBasic::get_disk(double R)
{
	CCeBasic * ce = new CCeBasic; 
	CMap     * mp = get_map(2, NULL_GENUS);

	ce->get_circle(R); 
	cells_new(1, 2, 0);

	bar_add((CCells *)ce);
	bar_span(mp);
}

/////////////////////////////////////
//...forming geometrical cell "ring";
void CCeBasic::get_ring(double rad, double R)
{
	CCeBasic * ce = new CCeBasic; 
	CMap     * mp = get_map(2, NULL_GENUS);

	ce->get_circle(rad); 
	cells_new(1, 2, 0);
	if (ce->mp) ce->mp[5] = M_PI;

	bar_add((CCells *)ce); ce = new CCeBasic; ce->get_circle(R);
	bar_add((CCells *)ce);
	bar_invers();
	bar_span(mp);
}

/////////////////////////////////////////
//...forming geometrical cell "cylinder";
void CCeBasic::get_cylinder(double R)
{
	int l;
	cells_new(1, 2, (l = size_of_map(2, CYL_GENUS))+1);
	if (mp) {
		mp[0] = ID_MAP(2, CYL_GENUS);
		mp[7] = (CMap)R; //...R < 0 - internal side of surface;
		mp[l] = (CMap)NULL_CELL;
	}
}

///////////////////////////////////////
//...forming geometrical cell "sphere";
void CCeBasic::get_sphere(double R)
{
	int l;
	cells_new(1, 2, (l = size_of_map(2, SPHERE_GENUS))+1);
	if (mp) {
		mp[0] = ID_MAP(2, SPHERE_GENUS);
		mp[7] = (CMap)R; //...R < 0 - internal side of surface;
      mp[l] = (CMap)NULL_CELL;
	}
}

/////////////////////////////////////////
//...forming geometrical cell "spheroid";
void CCeBasic::get_spheroid(double A, double B)
{
	int l;
	cells_new(1, 2, (l = size_of_map(2, SPHEROID_GENUS))+1);
	if (mp) {
		mp[0] = ID_MAP(2, SPHEROID_GENUS);
		mp[7] = (CMap)A;
		mp[8] = (CMap)B;
      mp[l] = (CMap)NULL_CELL;
	}
}
//////////////////////////////////////////
//...forming geometrical cell "ellipsoid";
void CCeBasic::get_ellipsoid(double A, double B, double C)
{
	int l;
	cells_new(1, 2, (l = size_of_map(2, ELLIPSOID_GENUS))+1);
	if (mp) {
		mp[0] = ID_MAP(2, ELLIPSOID_GENUS);
		mp[7] = (CMap)A;
		mp[8] = (CMap)B;
		mp[9] = (CMap)B;
      mp[l] = (CMap)NULL_CELL;
	}
}

/////////////////////////////////////
//...forming geometrical cell "ball";
void CCeBasic::get_ball(double R)
{
  CCeBasic * ce = new CCeBasic; 
  CMap     * mp = get_map(3, NULL_GENUS);

  ce->get_sphere(fabs(R)); 
  cells_new(1, 2, 0);

  bar_add((CCells *)ce);
  bar_span(mp);
}

/////////////////////////////////////
//...forming geometrical cell "roll";
void CCeBasic::get_roll(double R)
{
  CCeBasic * ce = new CCeBasic; 
  CMap     * mp = get_map(3, NULL_GENUS);

  ce->get_cylinder(fabs(R)); 
  cells_new(1, 2, 0);

  bar_add((CCells *)ce);
  bar_span(mp);
}

/////////////////////////////////////////////////////////
//...forming geometrical cell "space of dimension N_dim";
void CCeBasic::get_space(int N_dim)
{
  cells_new(1, 2, size_of_map(N_dim = abs(N_dim) % MAX_DIM, NULL_GENUS)+1);
  if (mp) {
      mp[0] = ID_MAP(N_dim, NULL_GENUS);
      mp[7] = (CMap)NULL_CELL;
  }
}

/////////////////////////////////////
//...forming geometrical cell "line";
void CCeBasic::get_line(double * P1, double * P2, int id_cell, int id_fast)
{
	CCeBasic * p1 = new CCeBasic,
				* p2 = new CCeBasic;
	CMap     * mp = get_map(1, NULL_GENUS);
	double    X  = P2[0]-P1[0], Y = P2[1]-P1[1],
				 Z  = P2[2]-P1[2], fi;
	p1->get_point(P1);
	p2->get_point(P2);

///////////////////////////////////////
//...setting geometrical card and cell;
	if (mp) {
		mp[1] = (CMap)(P1[0]+P2[0])*.5;
		mp[2] = (CMap)(P1[1]+P2[1])*.5;
		mp[3] = (CMap)(P1[2]+P2[2])*.5;
	}
	if (bar_add(p2, id_cell) && bar_add(p1, id_cell) && bar_span(mp)) {

//////////////////////////////////////////////
//...directly normal or spherical coordinates;
      if (id_fast == NULL_STATE) {
			this->mp[4] = (CMap)arg0(comp(X, Y));
			this->mp[5] = (CMap)arg0(comp(Z, sqrt(X*X+Y*Y)));
		}
		else {
          if ((fi = sqrt(sqr(X)+sqr(Y)+sqr(Z))) > EE) {
               this->mp[4] = X*(fi = 1./fi);
               this->mp[5] = Y* fi;
               this->mp[6] = Z* fi;
          }
          else {
					this->mp[4] = 0.;
					this->mp[5] = 0.;
					this->mp[6] = 1.;
          }
      }
	}
	delete_struct(mp);
}

//////////////////////////////////////////////
//...forming geometrical cell "arc of circle";
void CCeBasic::get_arc(double R, double f0, double f1, int id_cell)
{
  CCeBasic * p1 = new CCeBasic,
           * p2 = new CCeBasic;
  CMap     * mp = get_map(1, SPHERE_GENUS);
  double  P1[3] = {(R = fabs(R))*cos(f0), R*sin(f0), 0.},
          P2[3] = { R*cos(f1), R*sin(f1), 0.};
  p1->get_point(P1);
  p2->get_point(P2);

///////////////////////////////////////
//...setting geometrical card and cell;
  if (mp) {
      mp[5] = (CMap)(f1 < f0 ?  M_PI : 0.); //...mp[5] == M_PI - positive turning;
      mp[6] = (CMap)(mp[5]+cos(mp[5])*(f0+f1)*.5);
      mp[7] = (CMap)R;
      mp[8] = (CMap)fabs(f1-f0)*.5;         //...mp[8] - semi-length circular arc in radian;
  }
  bar_add(p2, id_cell);
  bar_add(p1, id_cell);
  bar_span(mp);
}

/////////////////////////////////////////////////////////////
//...forming geometrical cell "circle" on two complex points;
void CCeBasic::get_arc(double R, complex z1, complex z2, int topo)
{
  double  P1[3] = {real(z1),    imag(z1),    0.},
          P2[3] = {real(z2),    imag(z2),    0.},
          P0[3] = {P2[0]-P1[0], P2[1]-P1[1], 0.}, f = sqrt(P0[0]*P0[0]+P0[1]*P0[1]),
          R0    = fabs(f) < EE_ker ? 0. : fabs(R), f0 = R0 == 0. ? 0. : R*R-f*f*.25;

  if (! topo || (f0 = filtr_(f0, EE_ker)) < 0.) get_line(P1, P2);
  else {
      double beta  = arg0(comp(sqrt(f0), f*.5)),
             fi    = arg0(comp(P0[0], P0[1]))-M_PI_2,
             sg    = topo > 0 ? -1. : 1.,
             theta = topo > 0 ?  0. : M_PI; if (2 == abs(topo)) beta = M_PI-beta;
      P0[0] = P1[0]+sg*R0*cos(fi+sg*beta);
      P0[1] = P1[1]+sg*R0*sin(fi+sg*beta);

      get_arc  (R0, -beta, beta);   //...operation set   "returning" into mp[6];
      cells_iso(P0, fi, theta, 0.); //...operation reset "returning" into mp[6];
  }
}

/////////////////////////////////////////////////////
//...forming special geometrical cell "elliptic arc";
void CCeBasic::get_ellipt_arc(double A, double B, double f0, double f1, int id_cell)
{
  double C = .25*(A*A-B*B), R, Co, ff = .25*B*B; A *= .5; B *= .5;
  if (A < EE_ker || B < EE_ker || C < 0.) return; C = sqrt(C);

  CCeBasic * p1 = new CCeBasic,
           * p2 = new CCeBasic;
  CMap     * mp = get_map(1, CYL_GENUS, ELLIPT_ARC);
  double  P1[3] = {(R = ff/(A+C*(Co = cos(f0))))*Co+C, R*sin(f0), 0.},
          P2[3] = {(R = ff/(A+C*(Co = cos(f1))))*Co+C, R*sin(f1), 0.};
  int    k;
  p1->get_point(P1);
  p2->get_point(P2);

///////////////////////////////////////
//...setting geometrical card and cell;
  if (mp) {
      mp[5]   = (CMap)(f1 < f0 ?  M_PI : 0.); //...mp[5] == M_PI - positive turning;
      mp[7]   = (CMap)A;
      mp[8]   = (CMap)B;
      mp[++(k = size_of_map(mp))] = (CMap)(mp[5] == M_PI ? mp[5]-f0 : f0);
      mp[++k] = (CMap)(mp[5] == M_PI ? mp[5]-f1 : f1);
  }
  bar_add(p2, id_cell);
  bar_add(p1, id_cell);
  bar_span(mp);
}

//////////////////////////////////////
//...forming geometrical cell "sheet";
void CCeBasic::get_sheet(double A, double B)
{
	double P1[3] = {-.5*A, 0., 0.}, P2[3] = {.5*A, 0., 0.},
			P3[3] = {0., -.5*B, 0.}, P4[3] = {0., .5*B, 0.};
	int    k;
	if (A < EE_ker) get_line(P3, P4); else
	if (B < EE_ker) get_line(P1, P2); 
	else {
		CCeBasic * l1 = new CCeBasic,
					* l2 = new CCeBasic,
					* l3 = new CCeBasic,
					* l4 = new CCeBasic;
		CMap     * mp = get_map(2, NULL_GENUS, SHEET_CELL);
		l1->get_line(P1, P2);
		l2->get_line(P3, P4);
		l3->get_line(P2, P1);
		l4->get_line(P4, P3);

////////////////////////////////////////////
//...setting parameters of geometrical card;
		if (mp) {
			mp[++(k = size_of_map(mp))] = (CMap)A;
			mp[++k] = (CMap)B;
		}

//////////////////////////////////////////////////////
//...forming surface and removing degenerate elements;
		l4->cells_iso(P1); l2->cells_iso(P2);
		l1->cells_iso(P3); l3->cells_iso(P4);

		bar_add(l4); bar_add(l3);
		bar_add(l2); bar_add(l1);
		bar_generate();

//////////////////////////////////////////
//...spanning geometrical card and return;
		bar_span(mp);
	}
}

/////////////////////////////////////////////////////
//...forming geometrical cell "sheet with intrusion";
void CCeBasic::get_sheet_intrusion(double A, double B, double rad)
{
	if (rad > A*.5 || rad > B*.5) return;
	double P1[3] = {-.5*A, 0., 0.}, P2[3] = {.5*A, 0., 0.},
			 P3[3] = {0., -.5*B, 0.}, P4[3] = {0., .5*B, 0.};
	int k;
	if (A < EE_ker) get_line(P3, P4); else
	if (B < EE_ker) get_line(P1, P2); else {
		CCeBasic * l1 = new CCeBasic,
					* l2 = new CCeBasic,
					* l3 = new CCeBasic,
					* l4 = new CCeBasic,
					* cc = new CCeBasic;
		CMap * mp = get_map(2, NULL_GENUS, SHT_INTRUSION_CELL);
		l1->get_line(P1, P2);
		l2->get_line(P3, P4);
		l3->get_line(P2, P1);
		l4->get_line(P4, P3);
		cc->get_circle (rad);

////////////////////////////////////////////
//...setting parameters of geometrical card;
		if (mp) {
			mp[++(k = size_of_map(mp))] = (CMap)A;
			mp[++k] = (CMap)B;
			mp[++k] = (CMap)rad;
		}

//////////////////////////////////////////////////////
//...forming surface and removing degenerate elements;
		l4->cells_iso(P1); l2->cells_iso(P2);
		l1->cells_iso(P3); l3->cells_iso(P4);

		bar_add(cc);
		bar_add(l4); bar_add(l3); bar_add(l2); bar_add(l1);
		bar_generate();

//////////////////////////////////////////
//...spanning geometrical card and return;
		bar_span(mp);
	}
}

////////////////////////////////////
//...forming geometrical cell "box";
void CCeBasic::get_box(double A, double B, double C)
{
  CCeBasic * f1 = new CCeBasic,
           * f2 = new CCeBasic,
           * f3 = new CCeBasic,
           * f4 = new CCeBasic,
           * f5 = new CCeBasic,
           * f6 = new CCeBasic;
  CMap     * mp = get_map(3, NULL_GENUS, BOX_CELL);
  int     k;
  double  P1[3] = {-.5*A, 0., 0.},
          P2[3] = {0., -.5*B, 0.},
          P3[3] = {0., 0., -.5*C};
  f1->get_sheet(C, B);
  f2->get_sheet(C, B);
  f3->get_sheet(C, A);
  f4->get_sheet(C, A);
  f5->get_sheet(A, B);
  f6->get_sheet(A, B);

////////////////////////////////////////////
//...setting parameters of geometrical card;
     if (mp) {
         mp[++(k = size_of_map(mp))] = (CMap)A;
         mp[++k] = (CMap)B;
         mp[++k] = (CMap)C;
     }

///////////////////////////////////////////////////
//...forming surface and spanning geometrical card;
  f1->cells_iso(P1, 0.,     -M_PI_2); P1[0] = -P1[0]; f2->cells_iso(P1, 0.,     M_PI_2);
  f3->cells_iso(P2, M_PI_2, -M_PI_2); P2[1] = -P2[1]; f4->cells_iso(P2, M_PI_2, M_PI_2);
  f5->cells_iso(P3, 0.,        M_PI); P3[2] = -P3[2]; f6->cells_iso(P3, 0.,     0.);

  bar_add(f6); bar_add(f5); bar_add(f4);
  bar_add(f3); bar_add(f2); bar_add(f1);
  bar_span(mp);
}

//////////////////////////////////////////////////////
//...forming geometrical cell "sphere with intrusion";
void CCeBasic::get_sph_intrusion(double R, double L)
{
	if (fabs(R) < L || fabs(R) > M_SQRT2*L) return;
	CCeBasic * f1 = new CCeBasic,
			   * f2 = new CCeBasic,
			   * f3 = new CCeBasic,
			   * f4 = new CCeBasic,
			   * f5 = new CCeBasic,
			   * f6 = new CCeBasic;
	CMap * mp = get_map(2, SPHERE_GENUS, SPH_INTRUSION_CELL);
	int    k;
	double P1[3] = {-L, 0., 0.},
			 P2[3] = {0., -L, 0.},
			 P3[3] = {0., 0., -L}, rad = sqrt(R*R-L*L);
	f1->get_circle(rad);
	f2->get_circle(rad);
	f3->get_circle(rad);
	f4->get_circle(rad);
	f5->get_circle(rad);
	f6->get_circle(rad);

////////////////////////////////////////////
//...setting parameters of geometrical card;
	if (mp) {
		mp[7] = (CMap)R; //...R < 0 - internal side of surface;
		mp[++(k = size_of_map(mp))] = (CMap)2.*M_PI;
		mp[++k] = (CMap)0.;
		mp[++k] = (CMap)M_PI;
		mp[++k] = (CMap)L;
	}

///////////////////////////////////////////////////
//...forming surface and spanning geometrical card;
	f1->cells_iso(P1, 0.,     -M_PI_2); P1[0] = -P1[0]; f2->cells_iso(P1, 0.,     M_PI_2);
	f3->cells_iso(P2, M_PI_2, -M_PI_2); P2[1] = -P2[1]; f4->cells_iso(P2, M_PI_2, M_PI_2);
	f5->cells_iso(P3, 0.,        M_PI); P3[2] = -P3[2]; f6->cells_iso(P3, 0.,     0.);

	bar_add(f6); bar_add(f5); bar_add(f4);
	bar_add(f3); bar_add(f2); bar_add(f1);
	bar_span(mp);
}

///////////////////////////////////////////////////////////////////////////////
//...directly forming plane polygon (reverse ordering for boundary elements !);
void CCeBasic::get_polygon_directly(double * X, double * Y, int N, int id_dbl, double eps)
{
  if (! X || ! Y || N < 3) return;
  int k, l;
  cells_new(2*(N-1)+1, 2, 0); if (graph) graph[0] = 2*(N-1);

///////////////////////////////
//...creating polygon topology;
  for (k = N-1; k < 2*N-2; k++) { //...points;
       ce[k] = new CCells; if (! ce[k]) goto err;
       ce[k]->cells_new(0, 2, (l = size_of_map(0, NULL_GENUS))+1);
       ce[k]->mp[l] = (CMap)NULL_CELL;
       ce[k]->ce    = ce;
  }
  for (k = N-2; k >= 0; k--) { //...lines;
       ce[k] = new CCells; if (! ce[k]) goto err;
       ce[k]->cells_new(0, 4, (l = size_of_map(1, NULL_GENUS))+1);
       ce[k]->mp[l]    = (CMap)NULL_CELL;
       ce[k]->ce       = ce;
       if (ce[k]->graph) {
           ce[k]->graph[0] =
           ce[k]->graph[1] = 2;
           ce[k]->graph[2] = 2*N-3-k;
           ce[k]->graph[3] = k ? 2*N-2-k : N-1;
       }
  }

///////////////////////////////
//...creating polygon geometry;
  for (k = N-1; k < 2*N-2; k++) if (ce[k]->mp) {
       ce[k]->mp[0] = ID_MAP(0, NULL_GENUS);
       ce[k]->mp[1] = X[k-N+1];
       ce[k]->mp[2] = Y[k-N+1];
  }
  for (k = N-2; k >= 0; k--) if (ce[k]->mp) {
       ce[k]->mp[0] = ID_MAP(1, NULL_GENUS);
       ce[k]->line_correct();
  }
  install_struct();

//////////////////////////////////////
//...removing double numbering points;
  if (id_dbl) bar_correct_double_point(eps);
  return;
err:
  zero_cells();
}

////////////////////////////////////////////////////////////////////////////////////////
//...directly forming plane circular polygon (reverse ordering for boundary elements !);
void CCeBasic::get_arc_polygon_directly(double * X, double * Y, double * R, int * topo, int N, int id_dbl, double eps)
{
  if (! X || ! Y || ! R || ! topo || N < 2) return;
  int k, size_map = max(size_of_map(1, NULL_GENUS), size_of_map(1, SPHERE_GENUS));
  cells_new(2*(N-1)+1, 2, 0); if (graph) graph[0] = 2*(N-1);

///////////////////////////////
//...creating polygon topology;
  for (k = N-1; k < 2*N-2; k++) { //...points;
       ce[k] = new CCells; if (! ce[k]) goto err;
       ce[k]->cells_new(0, 2, size_of_map(0, NULL_GENUS)+1);
       ce[k]->ce = ce;
  }
  for (k = N-2; k >= 0; k--) { //...arcs;
       ce[k] = new CCells; if (! ce[k]) goto err;
       ce[k]->cells_new(0, 4, size_map+1);
       ce[k]->ce = ce;
       if (ce[k]->graph) {
           ce[k]->graph[0] =
           ce[k]->graph[1] = 2;
           ce[k]->graph[2] = 2*N-3-k;
           ce[k]->graph[3] = k ? 2*N-2-k : N-1;
       }
  }

///////////////////////////////
//...creating polygon geometry;
  for (k = N-1; k < 2*N-2; k++) if (ce[k]->mp) {
       ce[k]->mp[0] = ID_MAP(0, NULL_GENUS);
       ce[k]->mp[1] = (CMap)X[k-N+1];
       ce[k]->mp[2] = (CMap)Y[k-N+1];
       ce[k]->mp[4] = (CMap)NULL_CELL;
  }
  for (k = N-2; k >= 0; k--) if (ce[k]->mp)
  if  (topo[N-2-k]) {
       ce[k]->mp[0] = ID_MAP(1, SPHERE_GENUS);
       ce[k]->mp[9] = (CMap)NULL_CELL;
       ce[k]->circ_correct(R[N-2-k], topo[N-2-k]);
  }
  else {
       ce[k]->mp[0] = ID_MAP(1, NULL_GENUS);
       ce[k]->mp[7] = (CMap)NULL_CELL;
       ce[k]->line_correct();
  }
  install_struct();

//////////////////////////////////////
//...removing double numbering points;
  if (id_dbl) bar_correct_double_point(eps);
  return;
err:
  zero_cells();
}

///////////////////////////////////////////
//...directly forming space polygonal line;
void CCeBasic::get_line_strip_directly(double * X, double * Y, double * Z, int * geom, int * mask, int shift, int id_fast, int id_dbl, double eps)
{
  int k, l, m, N, i;
  if (X && Y && Z && geom && (N = geom[0]-shift) >= 2) {
      cells_new(2*N, 2, 0); if (graph) graph[0] = 2*N-1;

//////////////////////////////////
//...creating space line topology;
      for (l = size_of_map(0, NULL_GENUS), m = 2*N-1, k = N-1; ce && k < m; k++) {
           ce[k] = new CCells; if (! ce[k]) goto err;
           ce[k]->cells_new(0, 2, l+1);
           ce[k]->mp[l] = (CMap)NULL_CELL;
           ce[k]->ce    = ce;
      }
      for (l = size_of_map(1, NULL_GENUS), k = N-2; ce && k >= 0; k--) {
           ce[k] = new CCells; if (! ce[k]) goto err;
           ce[k]->cells_new(0, 4, l+1);
           ce[k]->mp[l]    = (CMap)NULL_CELL;
           ce[k]->ce       = ce;
           if (ce[k]->graph) {
               ce[k]->graph[0] =
               ce[k]->graph[1] = 2;
               ce[k]->graph[2] = k+N-1;
               ce[k]->graph[3] = k+N;
           }
      }

/////////////////////////////////////
//...creating geometry of space line;
      for (k = N-1; k < m; k++) {
           ce[k]->mp[0] = ID_MAP(0, NULL_GENUS);
           ce[k]->mp[1] = X ? X[i = mask ? mask[geom[k-N+2+shift]] : geom[k-N+2+shift]] : 0.;
           ce[k]->mp[2] = Y ? Y[i] : 0.;
           ce[k]->mp[3] = Z ? Z[i] : 0.;
      }
      for (k = N-2; k >= 0; k--) {
           ce[k]->mp[0] = ID_MAP(1, NULL_GENUS);
           ce[k]->line_correct(id_fast);
      }
		install_struct();

//////////////////////////////////////
//...removing double numbering points;
      if (id_dbl) bar_correct_double_point(eps);
      return;
  }
err:
  zero_cells();
}

////////////////////////////////////////////////////////
//...beam forming on the base of plane circular polygon;
void CCeBasic::get_beam(CCells * f0, double length)
{
	if (! f0 || f0->mp && f0->mp[0] != ID_MAP(2, NULL_GENUS) || length < EE_ker) return;

	double P[3] = {0., 0., -.5*length};
	int k, arc, prev;
	if (! f0->mp) {
		CMap * mp = get_map (2, NULL_GENUS);
		f0->bar_correct_double_point();
		f0->bar_generate();
		f0->bar_span(mp); delete_struct(mp);
	}

///////////////////////////////////////
//...корректируем прямолинейные звенья;
	for (k = 2; k <= f0->graph[1]+1; k++)
		f0->ce[f0->graph[k]]->line_correct();

////////////////////////////////
//...create first beam butt-end;
	f0->cells_iso(P); P[2] = length;
	bar_add(f0); f0 = bar_cpy(1);

///////////////////////////////////////////////////////////////////////////////////////////
//...create second beam butt-end (and reverse ordering of boundary elements of 1 butt_end);
	ce[0]->cells_invers();
	f0->cells_iso(P);

	bar_add(f0);
	bar_ord();

////////////////////////////////////////////////////////////////////////
//...create beam boundary surface (from elements -- dimension segments);
	for (k = 2, prev = ce[1]->graph[ce[1]->graph[1]+1]; k <= ce[1]->graph[1]+1; k++, prev = arc)
	if (ce[arc = ce[1]->graph[k]]->mp[0] == ID_MAP(1, NULL_GENUS)) {
			complex z1 = ce[arc]->get_arc_point(0),
					  z2 = ce[arc]->get_arc_point(1),
					  z0 = (z1+z2)*.5;
			P[0] = real(z0);
			P[1] = imag(z0);
			P[2] = ce[1]->mp[3]-.5*length;
			if ( ! ce[prev]->topo_id(ce[arc]->graph[2])) swap(z1, z2);

			f0 = new CCells; ((CCeBasic *)f0)->get_sheet(length, abs(z2-z1));
			f0->cells_iso(P, arg2(z2-z1)-M_PI_2, M_PI_2);
			bar_add(f0);
	}
	else
	if (ce[arc]->mp[0] == ID_MAP(1, SPHERE_GENUS)) {
			double  t1 =      ce[arc]->mp[4],
					  t2 =      ce[arc]->mp[5],
					  t3 =      ce[arc]->mp[6],
					  fi =      ce[arc]->mp[8]*cos(t2),
					  R  = fabs(ce[arc]->mp[7]);
			complex z0 =      ce[arc]->get_arc_center();
			P[0] = real(z0);
			P[1] = imag(z0);
			P[2] = ce[1]->mp[3]-.5*length;

			f0 = new CCells; ((CCeBasic *)f0)->get_cyl_segment(R, -fi, fi, -.5*length, .5*length);
			f0->cells_iso(P, t1, t2, t3);
			bar_add(f0);
	}
	bar_ord();
}

////////////////////////////////////////////////////////
//...beam forming on the base of plane circular polygon;
void CCeBasic::get_cyl_beam(CCells * f0, double fi)
{
	if (! f0 || f0->mp && f0->mp[0] != ID_MAP(2, NULL_GENUS) || fi < EE_ker) return;

	int  k, arc, prev;
	if (! f0->mp) {
		CMap * mp = get_map (2, NULL_GENUS);
		f0->bar_correct_double_point();
		f0->bar_generate();
		f0->bar_span(mp); delete_struct(mp);
	}

////////////////////////////////
//...create first beam butt-end;
	f0->cells_invers();
	bar_add(f0); f0 = bar_cpy(1);

///////////////////////////////////////////////////////////////////////////////////////////
//...create second beam butt-end (and reverse ordering of boundary elements of 1 butt_end);
	ce[0]->cells_invers();
	f0->cells_iso(NULL, 0., (fi = fi*M_PI/180.));
	bar_add(f0);
	bar_ord();

////////////////////////////////////////////////////////////////////////
//...create beam boundary surface (from elements -- dimension segments);
	for (k = 2, prev = ce[0]->graph[ce[0]->graph[1]+1]; k <= ce[0]->graph[1]+1; k++, prev = arc)
	if (ce[arc = ce[1]->graph[k]]->mp[0] == ID_MAP(1, NULL_GENUS)) {
			complex	z1 = ce[arc]->get_arc_point(0),
						z2 = ce[arc]->get_arc_point(1),
						z0 = (z1+z2)*.5;
			if (! ce[prev]->topo_id(ce[arc]->graph[2])) swap(z1, z2);

			f0 = new CCells; 
			if (fabs(real(z2)-real(z1)) < EE_ker) {
				((CCeBasic *)f0)->get_cyl_segment(real(z0), 0., fi, imag(z1), imag(z2));
				f0->cells_iso(NULL, M_PI_2, M_PI_2, -M_PI_2);
			}
			else
			if (fabs(imag(z2)-imag(z1)) < EE_ker) {
				double pp[] = {0., 0., imag(z0)};
				((CCeBasic *)f0)->get_ring_segment(0., fi, real(z1), real(z2));
				f0->cells_iso(pp);
				f0->cells_iso(NULL, M_PI_2, M_PI_2, -M_PI_2);
			}
			else {
//void get_cone_segment (double theta, double f0, double f1, double t0, double t1);
			}
			bar_add(f0);
	}
	else
	if (ce[arc]->mp[0] == ID_MAP(1, SPHERE_GENUS)) {
			complex	z0 =      ce[arc]->get_arc_center();
			double	t1 =      ce[arc]->mp[4],
						t2 =      ce[arc]->mp[5],
						t3 =      ce[arc]->mp[6],
						ff =      ce[arc]->mp[8],
						rr = fabs(ce[arc]->mp[7])*cos(t2), 
						tt = cos(t2) > 0. ? t1+t3+M_PI : t1-t3, pp[] = {0., 0., imag(z0)};
			f0 = new CCells; ((CCeBasic *)f0)->get_torus_segment(real(z0), rr, 0., fi, -ff+tt, ff+tt);
			f0->cells_iso(pp);
			f0->cells_iso(NULL, M_PI_2, M_PI_2, -M_PI_2);
			bar_add(f0);
	}
	bar_ord();
	cells_iso(NULL, 0., -fi*.5);
}

////////////////////////////////////////////////////////////////////////////
//...forming ring segment (0 > (f1-f0)*(t1-t0) - additional turn of normal);
void CCeBasic::get_ring_segment(double f0, double f1, double t0, double t1)
{
  if (t0 < 0. || t1 < 0.) return;
  int     k;
  double  P1[3] = {t0, 0., 0.},
          P2[3] = {t1, 0., 0.};
  CCeBasic * l1 = new CCeBasic,
           * l2 = new CCeBasic,
           * c1 = new CCeBasic,
           * c2 = new CCeBasic;
  CMap      * mp = get_map(2, NULL_GENUS, RING_SEGMENT);
  l1->get_line(P1, P2);
  l2->get_line(P2, P1);
  c1->get_arc(t1, f0, f1);
  c2->get_arc(t0, f1, f0);

////////////////////////////////////////////
//...setting parameters of geometrical card;
  if (mp) {
      mp[5] = (CMap)((f1 < f0 ? -1 : 1)*(t1 < t0 ? -1 : 1) == -1 ?  M_PI : 0.);
      mp[6] = (CMap)(mp[5]+cos(mp[5])*(f0+f1)*.5);
      mp[++(k = size_of_map(mp))] = (CMap)fabs(f1-f0);
      mp[++k] = (CMap)fabs(t0);
      mp[++k] = (CMap)fabs(t1);
  }

///////////////////////////////////////////////////
//...forming surface and spanning geometrical card;
  l1->cells_iso(NULL, f0);
  l2->cells_iso(NULL, f1);

  bar_add(c2); bar_add(l2);
  bar_add(c1); bar_add(l1); bar_generate();
  bar_span(mp);
}

//////////////////////////////////////////////////////////////////////////////
//...forming cylindrical segment (R > 0 - outer normal, R < 0 - inner normal);
void CCeBasic::get_cyl_segment(double R, double f0, double f1, double t0, double t1)
{
  int     k;
  double  P1[3] = {R = fabs(R), 0., t0},
          P2[3] = {    fabs(R), 0., t1};
  CCeBasic * l1 = new CCeBasic,
           * l2 = new CCeBasic,
           * c1 = new CCeBasic,
           * c2 = new CCeBasic;
  CMap     * mp = get_map(2, CYL_GENUS, CYL_SEGMENT);
  l1->get_line(P2, P1);
  l2->get_line(P1, P2);
  c1->get_arc(R, f0, f1);
  c2->get_arc(R, f1, f0);

////////////////////////////////////////////
//...setting parameters of geometrical card;
  if (mp) {
      mp[3] = (CMap)((t0+t1)*.5);
      mp[6] = (CMap)((f0+f1)*.5);
      mp[7] = (CMap)(R*(f1 < f0 ? -1. : 1.)*(t1 < t0 ? -1. : 1.));
      mp[++(k = size_of_map(mp))] = (CMap)fabs(f1-f0);
      mp[++k] = (CMap)(t1-t0);
  }

///////////////////////////////////////////////////
//...forming surface and spanning geometrical card;
  l1->cells_iso(NULL, f0);
  l2->cells_iso(NULL, f1); P1[0] = P2[0] = 0.;

  c1->cells_iso(P1);
  c2->cells_iso(P2);

  bar_add(c2); bar_add(l2);
  bar_add(c1); bar_add(l1); bar_generate();
  bar_span(mp);
}

//////////////////////////////////////////////////////////////
//...forming space cylindrical segment with blending on bands;
void CCeBasic::get_blend_cyl_segment(double beta, double R, double L, int * mm)
{
  double P[3] = {0., 0., L};
  int id_full = 0;
  CCeBasic * r0 = new CCeBasic,
           * r1 = new CCeBasic,
           * c2 = new CCeBasic,
           * p3 = new CCeBasic,
           * p4 = new CCeBasic;
  if (mm && mm[0]) r0->get_ring_segment( M_PI_2*beta, -M_PI_2*beta, 0., R);
  else             r0->get_sph_segment ( R, -M_PI_2*beta, M_PI_2*beta, M_PI_2, M_PI);
  if (mm && mm[1]) r1->get_ring_segment(-M_PI_2*beta,  M_PI_2*beta, 0., R);
  else             r1->get_sph_segment ( R, -M_PI_2*beta, M_PI_2*beta, 0., M_PI_2);
                   c2->get_cyl_segment ( R, -M_PI_2*beta, M_PI_2*beta, 0., L);
  if (id_full)     p3->get_sheet(L, R);
  if (id_full)     p4->get_sheet(L, R);

//////////////////////////////////////////////
//...setting element in geometrical structure;
  r1->cells_iso(P); P[1] = R*.5; P[2] = L*.5;

  p3->cells_iso(P, 0., -M_PI_2); p3->cells_iso(NULL, M_PI_2*(beta-1.)); P[1] = -P[1];
  p4->cells_iso(P, 0.,  M_PI_2); p4->cells_iso(NULL, M_PI_2*(1.-beta));

  bar_add(r0); bar_add(r1); bar_add(c2);
  if (id_full) {
      bar_add(p3);
      bar_add(p4);
  }
  bar_ord(); //bar_generate();
}

////////////////////////////////////////////////////////////////////////////
//...forming spherical segment (R > 0 - outer normal, R < 0 - inner normal);
void CCeBasic::get_sph_segment(double R, double f0, double f1, double t0, double t1)
{
  if (t0 < 0. || t1 < 0. || t0 > M_PI || t1 > M_PI) return;
  int      k;
  double   P[3] = {0., 0., (R = fabs(R))*cos(t1)};
  CCeBasic * c1 = new CCeBasic,
           * c2 = new CCeBasic,
           * c3 = new CCeBasic,
           * c4 = new CCeBasic;
  CMap     * mp = get_map(2, SPHERE_GENUS, SPH_SEGMENT);
  c1->get_arc(R, t0, t1);
  c2->get_arc(R*sin(t1), f0, f1);
  c3->get_arc(R, t1, t0);
  c4->get_arc(R*sin(t0), f1, f0); //...right away returning;

////////////////////////////////////////////
//...setting parameters of geometrical card;
  if (mp) {
      mp[6] = (CMap)((f0+f1)*.5);
      mp[7] = (CMap)(R*(f1 < f0 ? -1. : 1.)*(t1 < t0 ? -1. : 1.));
      mp[++(k = size_of_map(mp))] = (CMap)fabs(f1-f0);
      mp[++k] = (CMap)t0;
      mp[++k] = (CMap)t1;
  }

///////////////////////////////////////////////////
//...forming surface and spanning geometrical card;
  c1->cells_iso(NULL, -M_PI_2, -M_PI_2); c1->cells_iso(NULL, f0);
  c3->cells_iso(NULL, -M_PI_2, -M_PI_2); c3->cells_iso(NULL, f1);

  c2->cells_iso(P); P[2] = R*cos(t0);
  c4->cells_iso(P);

  bar_add(c4); bar_add(c3);
  bar_add(c2); bar_add(c1); bar_generate();
  bar_span(mp);
}

///////////////////////////////////////////////////////////////////////////////
//...forming cone segment (mp[7] > 0 - outer normal, mp[7] < 0 - inner normal);
void CCeBasic::get_cone_segment(double theta, double f0, double f1, double t0, double t1)
{
  if (theta < EE_ker || theta > M_PI-EE_ker) return;
  int     k;
  double  Co    = cos(theta), Si = sin(theta),
          P1[3] = {t0*Si, 0., t0*Co},
          P2[3] = {t1*Si, 0., t1*Co};
  CCeBasic * l1 = new CCeBasic,
           * l2 = new CCeBasic,
           * c1 = new CCeBasic,
           * c2 = new CCeBasic;
  CMap     * mp = get_map(2, CONE_GENUS, CONE_SEGMENT);
  l1->get_line(P2, P1);
  l2->get_line(P1, P2);
  c1->get_arc(P1[0], f0, f1);
  c2->get_arc(P2[0], f1, f0);

////////////////////////////////////////////
//...setting parameters of geometrical card;
  if (mp) {
      mp[6] = (CMap)((f0+f1)*.5);
      mp[7] = (CMap)(theta*(f1 < f0 ? -1. : 1.)*(t1 < t0 ? -1. : 1.));
      mp[++(k = size_of_map(mp))] = (CMap)fabs(f1-f0);
      mp[++k] = (CMap)t0;
      mp[++k] = (CMap)t1;
  }

///////////////////////////////////////////////////
//...forming surface and spanning geometrical card;
  l1->cells_iso(NULL, f0);
  l2->cells_iso(NULL, f1); P1[0] = P2[0] = 0.;

  c1->cells_iso(P1);
  c2->cells_iso(P2);

  bar_add(c2); bar_add(l2);
  bar_add(c1); bar_add(l1); bar_generate();
  bar_span(mp);
}

////////////////////////////////////////////////////////////////////////////////
//...forming torus segment (mp[7] > 0 - outer normal, mp[7] < 0 - inner normal);
void CCeBasic::get_torus_segment(double r0, double r1, double f0, double f1, double t0, double t1)
{
  if ((r0 = fabs(r0)) < (r1 = fabs(r1))+EE_ker) return;
  int      k;
  double   P[3] = {0., 0., -r1*sin(t0)};
  CCeBasic * c1 = new CCeBasic,
           * c2 = new CCeBasic,
           * c3 = new CCeBasic,
           * c4 = new CCeBasic;
  CMap      * mp = get_map(2, TORUS_GENUS, TORUS_SEGMENT);
  c1->get_arc(r0-r1*cos(t0), f1, f0);
  c3->get_arc(r0-r1*cos(t1), f0, f1);
  c2->get_arc(r1, -M_PI_2-t0, -M_PI_2-t1);
  c4->get_arc(r1, -M_PI_2-t1, -M_PI_2-t0);

////////////////////////////////////////////
//...setting parameters of geometrical card;
  if (mp) {
      mp[6] = (CMap)((f0+f1)*.5);
      mp[7] = (CMap)(r0*(f1 < f0 ? -1. : 1.)*(t1 < t0 ? -1. : 1.));
      mp[8] = (CMap)r1;
      mp[++(k = size_of_map(mp))] = (CMap)fabs(f1-f0);
      mp[++k] = (CMap)t0;
      mp[++k] = (CMap)t1;
  }

///////////////////////////////////////////////////
//...forming surface and spanning geometrical card;
  c1->cells_iso(P); P[2] = -r1*sin(t1);
  c3->cells_iso(P); P[2] = 0.; P[0] = r0;

  c2->cells_iso(P, -M_PI_2, -M_PI_2); c2->cells_iso(NULL, f0);
  c4->cells_iso(P, -M_PI_2, -M_PI_2); c4->cells_iso(NULL, f1);

  bar_add(c4); bar_add(c3);
  bar_add(c2); bar_add(c1); bar_generate();
  bar_span(mp);
}

///////////////////////////////////////////////
//...фоpмиpование объединенного блока для кpуга;
void CCeBasic::get_circle_profile(double r)
{
  if (r < EE_dop) return;

////////////////////////////////////////////////////////////
//...определяем положение образца в общей системе координат;
  complex z0 = comp(r),
          z1 = polar(r, 2.*M_PI/3.),
          z2 = conj(z1);
  CCeBasic * ce;

//////////////////////////////
//...стpоим объединенный блок;
  cells_new(1, 2, 0);
  ce = new CCeBasic;  ce->get_arc(r, z0, z1, 1);  bar_add(ce);
  ce = new CCeBasic;  ce->get_arc(r, z1, z2, 1);  bar_add(ce);
  ce = new CCeBasic;  ce->get_arc(r, z2, z0, 1);  bar_add(ce);
  bar_generate();
  bar_ord();
}


////////////////////////////////////////////////////
//...forming surface structure for L-type structure;
void CCeBasic::get_ugolok(double beta, double radc, double rad, double A1, double B1, double A2, double B2)
{
  if (beta <= 0 || A1 < -EE_dop || B1 <= 0.  || A2 < -EE_dop || B2 <= 0.
                                || radc < 0. || rad < 0.) return;
  if (beta < 1.) {
      radc *= -1.; rad *= -1.;
  }
  if (beta == 1.) {
      B1 = .5*(B1+B2); B2 = B1; radc = 0.; rad = 0.;
  }

//////////////////////////
//...auxilliary variables;
  double  C = cos(M_PI*beta), S = sin(M_PI*beta), T = tan(M_PI_2*(beta-1.)), t1 = 0., t2 = 0.;
  double  RR[9], X[9], Y[9], R0 = radc*T;
  int     t[9], m = beta < 1. ? 1 : -1;
  complex z, z0 = polar(1., M_PI_2*beta);
  if (fabs(S) > EE_dop) {
      t1 = (B2+B1*C)/S; t2 = (B1+B2*C)/S;
  }

////////////////////////////////////
//...setting conotur prizmatic beam;
  z = A2*z0;                       X[0] = real(z); Y[0] = imag(z); RR[0] = 0.;         t[0] =  0;
  z = R0*z0;                       X[1] = real(z); Y[1] = imag(z); RR[1] = fabs(radc); t[1] =  m;
  z = comp(X[1], -Y[1]);           X[2] = real(z); Y[2] = imag(z); RR[2] = 0;          t[2] =  0;
  z = A1*conj(z0);                 X[3] = real(z); Y[3] = imag(z); RR[3] = 0;          t[3] =  0;
  z = comp(A1, B1)*conj(z0);       X[4] = real(z); Y[4] = imag(z); RR[4] = 0;          t[4] =  0;
  z = comp(t1+rad*T, B1)*conj(z0); X[5] = real(z); Y[5] = imag(z); RR[5] = fabs(rad);  t[5] = -m;
  z = comp(t2+rad*T, -B2)*z0;      X[6] = real(z); Y[6] = imag(z); RR[6] = 0;          t[6] =  0;
  z = comp(A2, -B2)*z0;            X[7] = real(z); Y[7] = imag(z); RR[7] = 0;          t[7] =  0;
  z = comp(X[0], Y[0]);            X[8] = real(z); Y[8] = imag(z); RR[8] = 0;          t[8] =  0;

////////////////////////////////////////////////////////
//...checking correctness and forming contour structure;
  if (t1+rad*T > A1+EE_dop || t2+rad*T > A2+EE_dop
                           || R0 > A1+EE_dop || R0 > A2+EE_dop) return;

  get_arc_polygon_directly(X, Y, RR, t, 9);
}

/////////////////////////////////////////////////////
//...forming spacing cell, bounded by L-type contour;
void CCeBasic::get_ugolok_cell(double A1,  double A2,   double B1, double B2, double radc,
                               double rad, double beta, double x0, double y0, double f0)
{
  int    k;
  double P[3] = {x0, y0, 0.};
  CMap   * mp  = get_map(2, NULL_GENUS, UGOLOK_CELL);

////////////////////////////////////////////
//...setting parameters of geometrical card;
  if (mp) {
      mp[++(k = size_of_map(mp))] = (CMap)A1;
      mp[++k] = (CMap)A2;
      mp[++k] = (CMap)B1;
      mp[++k] = (CMap)B2;
      mp[++k] = (CMap)radc;
      mp[++k] = (CMap)rad;
      mp[++k] = (CMap)beta;
  }

///////////////////////////
//...forming cell "ugolok";
  get_ugolok(beta/180., radc, rad, A2, B2, A1, B1);
  bar_correct_double_point();
  bar_generate();
  bar_span(mp);

/////////////////////////////
//...setting ugolok in space;
  cells_iso(P, f0/180.*M_PI);
}
