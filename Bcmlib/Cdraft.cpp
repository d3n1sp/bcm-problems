#include "stdafx.h"
#include "cdraft.h"

#include "cgrid_el.h"
#include "ccebasic.h"


#ifdef ___CONSTRUCTION___

/////////////////////////////////////////////////////
//...existing of the boundary element in block links;
int CDraft::link_id_boundary(Block & B)
{
  if (B.link) {
     int j = B.link[0];

     while (0 < j && 0 <= B.link[j]) j--;
     return(0 < j);
  }
  return(1);
}

///////////////////////////////////////////////////////////////////
//...existing of the boundary element in the block (direct search);
int CDraft::link_id_boundary(int i)
{
  if (i < 0 || N <= i || ! B || ! B[i].link) return(1);
  int j = 1;

  while (j <= B[i].link[0] && 0 <= B[i].link[j]) j++;
  return(j <= B[i].link[0]);
}

////////////////////////////////////////////////////////////////////
//...number of degree of freedom for all neiborhoods of the element;
int CDraft::freedom_id(int element)
{
  int    m,  i; 
  for (  m = i = 0; i < N; i++)
  if (link_id(B[i], element))
         m += B[i].shape ? B[i].shape->NN : 0;
  return(m);
}

////////////////////////////////////////////////////////
//...number of degree of freedom for all special blocks;
int CDraft::freedom_id1(int N_links)
{
  int    m,  i; 
  for (  m = i = 0; i < N; i++)
  if (B[i].link && B[i].link[0] == N_links)
         m += B[i].shape ? B[i].shape->NN : 0;
  return(m);
}

///////////////////////////////////////////////
//...number of degree of freedom for the block;
int CDraft::freedom_block(int k)
{
  int    m = k < N && B[k].shape ? B[k].shape->NN : 0;
  return(m);
}

/////////////////////////////////////////////////////////
//...number of degree of freedom for all boundary blocks;
int CDraft::freedom_boundary()
{
  int    m,  i; 
  for (  m = i = 0; i < N; i++) if (link_id_boundary(B[i]))
         m += B[i].shape ? B[i].shape->NN : 0;
  return(m);
}

////////////////////////////////////////
//...presence point into block geometry;
int CDraft::map_into(int k, double X, double Y, double Z, int id_special)
{
    if (0 <= k && k < N && B && B[k].mp) return ::map_into(B[k].mp, X-B[k].mp[1], Y-B[k].mp[2], Z-B[k].mp[3], id_special);
    else                                 return(0);
}

////////////////////////////////////////////////////////////
//...presence point into block geometry (for visualization);
int CDraft::map1into(int k, double X, double Y, double Z)
{
    if (0 <= k && k < N && B && B[k].mp) return ::map1into(B[k].mp, X-B[k].mp[1], Y-B[k].mp[2], Z-B[k].mp[3]);
    else                                 return(0);
}


/////////////////////////////////////////////////////
//...block structure on the base of the plane sample;
int CDraft::ConvertToBeamStruct(int * NN, double * LL, int m_cell, int id_phase)
{
	if (! NN || ! LL || NN[0] <= 0) return(0);

	double pp_cond[] = {0., 0., 0.}, dd, leng;
	int i, j, k, l, m, n, j_shift, N_arcs, Max_phase, shift = bar ? 2 : 0, N0 = N;

//////////////////////////////////////////
//...construction of the blocks structure;
//	for (                           k = 0; k < N0; k++) B[k].bar->bar_invers();
	for (Max_phase = 0,             k = 0; k < N0; k++) if (B[k].link[NUM_PHASE] < Max_phase) Max_phase = B[k].link[NUM_PHASE];
	for (leng = 0., j_shift = 0,    n = 0; n < NN[0]; leng = LL[n], n++)
	for (dd = (LL[n]-leng)/NN[n+1], j = 0; j < NN[n+1]; j_shift++,  j++) if (j || n)
	for (pp_cond[2] = leng+dd*(j+.5),	  k = 0; k < N0;						 k++) {

/////////////////////////////////
//...образование геометрии блока;
		CCeBasic * beam = (CCeBasic *)CreateCells(BASIC_CELLS);

		beam->get_beam(m_cell ? B[k].bar->bar_cpy(m_cell) : B[k].bar->bar_cpy(), dd);
//		beam->cells_out("bar");

		beam->segms_bar();
		beam->cells_iso(pp_cond);
//		beam->cells_out("bar");

////////////////////////////////////
//...достраиваем торцевые сегменты;
		if ((B[k].type & ERR_MASK) == SUB_UGOLOK_BLOCK && beam->ce[0]->mp[m = size_of_map(beam->ce[0]->mp)] == NULL_CELL) {
			double A1, B1, X1, Y1, X2, Y2, rr, f0;
			int k_line, k_point1, k_point2, k_circ;

			k_line   = beam->ce[0]->graph[2];
			k_point1 = beam->ce[k_line]->graph[2];
			k_point2 = beam->ce[k_line]->graph[3];

			X1 = beam->ce[k_point1]->mp[1]-beam->ce[k_point2]->mp[1]; 
			Y1 = beam->ce[k_point1]->mp[2]-beam->ce[k_point2]->mp[2]; 
			B1 = sqrt(sqr(X1)+sqr(Y1));

			k_line   = beam->ce[0]->graph[3];
			k_point1 = beam->ce[k_line]->graph[2];
			k_point2 = beam->ce[k_line]->graph[3];

			X2 = beam->ce[k_point2]->mp[1]-beam->ce[k_point1]->mp[1]; 
			Y2 = beam->ce[k_point2]->mp[2]-beam->ce[k_point1]->mp[2];
			A1 = sqrt(sqr(X2)+sqr(Y2));

			k_circ = beam->ce[0]->graph[7];
			rr = fabs(beam->ce[k_circ]->mp[7]);

			if (beam->ce[beam->ce[0]->graph[2]]->topo_id(k_point2)) { iSWAP(k_point1, k_point2); f0 = -M_PI; }
			map_iso(beam->ce[0]->mp, NULL, f0 += arg0(comp(X2, Y2))+M_PI*.25);
			map_iso(beam->ce[1]->mp, NULL, f0);
			beam->ce[0]->mp[1] = beam->ce[1]->mp[1] = beam->ce[k_point2]->mp[1]-B1*(cos(f0-M_PI*.25)-sin(f0-M_PI*.25));
			beam->ce[0]->mp[2] = beam->ce[1]->mp[2] = beam->ce[k_point2]->mp[2]-B1*(cos(f0-M_PI*.25)+sin(f0-M_PI*.25));

			if (A1 > B1+rr && add_new_maps(beam->ce[0]->mp, (m = size_of_map(beam->ce[0]->mp))+1, 7)) {
				beam->ce[0]->mp[  m] = (CMap)UGOLOK_CELL;
				beam->ce[0]->mp[++m] = (CMap)(A1-B1);
				beam->ce[0]->mp[++m] = (CMap)(A1-B1);
				beam->ce[0]->mp[++m] = (CMap)B1;
				beam->ce[0]->mp[++m] = (CMap)B1;
				beam->ce[0]->mp[++m] = (CMap)rr;
				beam->ce[0]->mp[++m] = (CMap)0.;
				beam->ce[0]->mp[++m] = (CMap)270.;
			}
			if (A1 > B1+rr && add_new_maps(beam->ce[1]->mp, (m = size_of_map(beam->ce[1]->mp))+1, 7)) {
				beam->ce[1]->mp[  m] = (CMap)UGOLOK_CELL;
				beam->ce[1]->mp[++m] = (CMap)(A1-B1);
				beam->ce[1]->mp[++m] = (CMap)(A1-B1);
				beam->ce[1]->mp[++m] = (CMap)B1;
				beam->ce[1]->mp[++m] = (CMap)B1;
				beam->ce[1]->mp[++m] = (CMap)rr;
				beam->ce[1]->mp[++m] = (CMap)0.;
				beam->ce[1]->mp[++m] = (CMap)270.;
			}
		}
		if ((B[k].type & ERR_MASK) == SUB_SREZKA_BLOCK && beam->ce[0]->mp[m = size_of_map(beam->ce[0]->mp)] == NULL_CELL) {
			double A1, B1, X1, Y1, X2, Y2, rr, f0;
			int k_line, k_point1, k_point2, k_circ;

			k_line   = beam->ce[0]->graph[2];
			k_point1 = beam->ce[k_line]->graph[2];
			k_point2 = beam->ce[k_line]->graph[3];

			X1 = beam->ce[k_point1]->mp[1]-beam->ce[k_point2]->mp[1]; 
			Y1 = beam->ce[k_point1]->mp[2]-beam->ce[k_point2]->mp[2]; 
			B1 = sqrt(sqr(X1)+sqr(Y1));

			k_line   = beam->ce[0]->graph[3];
			k_point1 = beam->ce[k_line]->graph[2];
			k_point2 = beam->ce[k_line]->graph[3];

			X2 = beam->ce[k_point2]->mp[1]-beam->ce[k_point1]->mp[1]; 
			Y2 = beam->ce[k_point2]->mp[2]-beam->ce[k_point1]->mp[2];
			A1 = sqrt(sqr(X2)+sqr(Y2));

			k_circ = beam->ce[0]->graph[5];
			rr = fabs(beam->ce[k_circ]->mp[7]);

			if (beam->ce[beam->ce[0]->graph[2]]->topo_id(k_point2)) { iSWAP(k_point1, k_point2); f0 = -M_PI; }
			map_iso(beam->ce[0]->mp, NULL, f0 += arg0(comp(X2, Y2))+M_PI*.75);
			map_iso(beam->ce[1]->mp, NULL, f0);
			beam->ce[0]->mp[1] = beam->ce[1]->mp[1] = beam->ce[k_point1]->mp[1]+B1*cos(f0+M_PI*.75)-A1*sin(f0+M_PI*.75);
			beam->ce[0]->mp[2] = beam->ce[1]->mp[2] = beam->ce[k_point1]->mp[2]+B1*cos(f0+M_PI*.75)+A1*sin(f0+M_PI*.75);

			if (B[k].type & ERR_CODE) dSWAP(A1, B1);
			if (add_new_maps(beam->ce[0]->mp, (m = size_of_map(beam->ce[0]->mp))+1, 7)) {
				beam->ce[0]->mp[  m] = (CMap)UGOLOK_CELL;
				beam->ce[0]->mp[++m] = (CMap)A1;
				beam->ce[0]->mp[++m] = (CMap)B1;
				beam->ce[0]->mp[++m] = (CMap)B1;
				beam->ce[0]->mp[++m] = (CMap)A1;
				beam->ce[0]->mp[++m] = (CMap)rr;
				beam->ce[0]->mp[++m] = (CMap)0.;
				beam->ce[0]->mp[++m] = (CMap)90.;
			}
			if (add_new_maps(beam->ce[1]->mp, (m = size_of_map(beam->ce[1]->mp))+1, 7)) {
				beam->ce[1]->mp[  m] = (CMap)UGOLOK_CELL;
				beam->ce[1]->mp[++m] = (CMap)B1;
				beam->ce[1]->mp[++m] = (CMap)A1;
				beam->ce[1]->mp[++m] = (CMap)A1;
				beam->ce[1]->mp[++m] = (CMap)B1;
				beam->ce[1]->mp[++m] = (CMap)rr;
				beam->ce[1]->mp[++m] = (CMap)0.;
				beam->ce[1]->mp[++m] = (CMap)90.;
			}
		}
		if ((B[k].type & ERR_MASK) == SUB_UGOLOK2BLOCK && beam->ce[0]->mp[m = size_of_map(beam->ce[0]->mp)] == NULL_CELL) {
			double A1, B1, X1, Y1, X2, Y2, rr, f0;
			int k_line, k_point1, k_point2, k_circ;

			k_line   = beam->ce[0]->graph[2];
			k_point1 = beam->ce[k_line]->graph[2];
			k_point2 = beam->ce[k_line]->graph[3];

			X1 = beam->ce[k_point1]->mp[1]-beam->ce[k_point2]->mp[1]; 
			Y1 = beam->ce[k_point1]->mp[2]-beam->ce[k_point2]->mp[2]; 
			B1 = sqrt(sqr(X1)+sqr(Y1));

			k_line   = beam->ce[0]->graph[3];
			k_point1 = beam->ce[k_line]->graph[2];
			k_point2 = beam->ce[k_line]->graph[3];

			X2 = beam->ce[k_point2]->mp[1]-beam->ce[k_point1]->mp[1]; 
			Y2 = beam->ce[k_point2]->mp[2]-beam->ce[k_point1]->mp[2];
			A1 = sqrt(sqr(X2)+sqr(Y2));

			k_circ = beam->ce[0]->graph[6];
			rr = fabs(beam->ce[k_circ]->mp[7]);

			if (beam->ce[beam->ce[0]->graph[2]]->topo_id(k_point2)) { iSWAP(k_point1, k_point2); f0 = -M_PI; }
			map_iso(beam->ce[0]->mp, NULL, f0 += arg0(comp(X2, Y2))+M_PI*.25);
			map_iso(beam->ce[1]->mp, NULL, f0);
			beam->ce[0]->mp[1] = beam->ce[1]->mp[1] = beam->ce[k_point2]->mp[1]-B1*(cos(f0-M_PI*.25)-sin(f0-M_PI*.25));
			beam->ce[0]->mp[2] = beam->ce[1]->mp[2] = beam->ce[k_point2]->mp[2]-B1*(cos(f0-M_PI*.25)+sin(f0-M_PI*.25));

			if (A1 > B1+rr && add_new_maps(beam->ce[0]->mp, (m = size_of_map(beam->ce[0]->mp))+1, 7)) {
				beam->ce[0]->mp[  m] = (CMap)UGOLOK_CELL;
				beam->ce[0]->mp[++m] = (CMap)(A1-B1);
				beam->ce[0]->mp[++m] = (CMap)(A1-B1);
				beam->ce[0]->mp[++m] = (CMap)B1;
				beam->ce[0]->mp[++m] = (CMap)B1;
				beam->ce[0]->mp[++m] = (CMap)rr;
				beam->ce[0]->mp[++m] = (CMap)0.;
				beam->ce[0]->mp[++m] = (CMap)270.;
			}
			if (A1 > B1+rr && add_new_maps(beam->ce[1]->mp, (m = size_of_map(beam->ce[1]->mp))+1, 7)) {
				beam->ce[1]->mp[  m] = (CMap)UGOLOK_CELL;
				beam->ce[1]->mp[++m] = (CMap)(A1-B1);
				beam->ce[1]->mp[++m] = (CMap)(A1-B1);
				beam->ce[1]->mp[++m] = (CMap)B1;
				beam->ce[1]->mp[++m] = (CMap)B1;
				beam->ce[1]->mp[++m] = (CMap)rr;
				beam->ce[1]->mp[++m] = (CMap)0.;
				beam->ce[1]->mp[++m] = (CMap)270.;
			}
		}

/////////////////////////////////
//...образование блока в образце;
		add_block(B[k].type);
		B[k+j_shift*N0].bar = beam;
		set_block(B[k+j_shift*N0], SPHERE_GENUS, 1.);

//////////////////////////////////
//...линкование блочной структуры;
		set_link (B[k+j_shift*N0], max(NUM_PHASE, B[k].link[0]));
		if (j_shift > 1)
		B[k+(j_shift-1)*N0].link[2] = k+(j_shift+0)*N0; 
		B[k+(j_shift+0)*N0].link[1] = k+(j_shift-1)*N0; 
		if (m_cell) {
			for (m = min(NUM_PHASE-2, B[k].link[0]), i = 0; i < m; i++)
				if (B[k].link[i+1] >= 0) B[k+j_shift*N0].link[i+3] = B[k].link[i+1]+j_shift*N0; else
				if (B[k].link[i+1] <= SRF_STATE) {
					B[k+j_shift*N0].link[i+3] = B[k].link[i+1]-shift;	
					for (l = NUM_PHASE; l < B[k].link[0]; l++) //...коррекция многофазных сред;
						if (B[k].link[l+1] == SRF_STATE-B[k].link[i+1]) {
							B[k+j_shift*N0].link[i+3] += shift-j_shift*N0; break;
						}
				}
				else B[k+j_shift*N0].link[i+3] = B[k].link[i+1];
		}
		else {
			for (m = min(NUM_PHASE-2, B[k].bar->arcs_number()), i = 0; i < m; i++)
				if (B[k].link[i+1] >= 0) B[k+j_shift*N0].link[m-i+2] = B[k].link[i+1]+j_shift*N0; else
				if (B[k].link[i+1] <= SRF_STATE) {
					B[k+j_shift*N0].link[m-i+2] = B[k].link[i+1]-shift;
					for (l = NUM_PHASE; l < B[k].link[0]; l++) //...коррекция многофазных сред;
						if (B[k].link[l+1] == SRF_STATE-B[k].link[i+1]) {
							B[k+j_shift*N0].link[m-i+2] += shift-j_shift*N0; break;
						}
				}
				else B[k+j_shift*N0].link[m-i+2] = B[k].link[i+1];
		}
		for (i = NUM_PHASE; i < B[k+j_shift*N0].link[0]; i++) {
			if (B[k].link[i+1] >= 0)	B[k+j_shift*N0].link[i+1] = B[k].link[i+1]+j_shift*N0;
			else								B[k+j_shift*N0].link[i+1] = B[k].link[i+1];
		}
		if (B[k].link[0] < NUM_PHASE)	B[k+j_shift*N0].link[0] = 2+B[k].link[0];
		if (id_phase)						B[k+j_shift*N0].link[NUM_PHASE] = B[k].link[NUM_PHASE]+n*(Max_phase-1);
		else									B[k+j_shift*N0].link[NUM_PHASE] = B[k].link[NUM_PHASE];
	}

///////////////////////////////////////////////////////////////
//...преобразование блоков в нулевом слое и окаймляющего блока;
	for (dd = LL[0]/NN[1], pp_cond[2] = dd*.5, k = 0; k < N0; k++) {
		CCeBasic * beam = (CCeBasic *)CreateCells(BASIC_CELLS);

		beam->get_beam(m_cell ? B[k].bar->bar_cpy(m_cell) : B[k].bar->bar_cpy(), dd);
		beam->segms_bar();
		beam->cells_iso(pp_cond); N_arcs = m_cell ? B[k].link[0] : B[k].bar->arcs_number(); delete B[k].bar;

////////////////////////////////////
//...достраиваем торцевые сегменты;
		if ((B[k].type & ERR_MASK) == SUB_UGOLOK_BLOCK && beam->ce[0]->mp[m = size_of_map(beam->ce[0]->mp)] == NULL_CELL) {
			double A1, B1, X1, Y1, X2, Y2, rr, f0;
			int k_line, k_point1, k_point2, k_circ;

			k_line   = beam->ce[0]->graph[2];
			k_point1 = beam->ce[k_line]->graph[2];
			k_point2 = beam->ce[k_line]->graph[3];

			X1 = beam->ce[k_point1]->mp[1]-beam->ce[k_point2]->mp[1]; 
			Y1 = beam->ce[k_point1]->mp[2]-beam->ce[k_point2]->mp[2]; 
			B1 = sqrt(sqr(X1)+sqr(Y1));

			k_line   = beam->ce[0]->graph[3];
			k_point1 = beam->ce[k_line]->graph[2];
			k_point2 = beam->ce[k_line]->graph[3];

			X2 = beam->ce[k_point2]->mp[1]-beam->ce[k_point1]->mp[1]; 
			Y2 = beam->ce[k_point2]->mp[2]-beam->ce[k_point1]->mp[2];
			A1 = sqrt(sqr(X2)+sqr(Y2));

			k_circ = beam->ce[0]->graph[7];
			rr = fabs(beam->ce[k_circ]->mp[7]);

			if (beam->ce[beam->ce[0]->graph[2]]->topo_id(k_point2)) { iSWAP(k_point1, k_point2); f0 = -M_PI; }
			map_iso(beam->ce[0]->mp, NULL, f0 += arg0(comp(X2, Y2))+M_PI*.25);
			map_iso(beam->ce[1]->mp, NULL, f0);
			beam->ce[0]->mp[1] = beam->ce[1]->mp[1] = beam->ce[k_point2]->mp[1]-B1*(cos(f0-M_PI*.25)-sin(f0-M_PI*.25));
			beam->ce[0]->mp[2] = beam->ce[1]->mp[2] = beam->ce[k_point2]->mp[2]-B1*(cos(f0-M_PI*.25)+sin(f0-M_PI*.25));

			if (A1 > B1+rr && add_new_maps(beam->ce[0]->mp, (m = size_of_map(beam->ce[0]->mp))+1, 7)) {
				beam->ce[0]->mp[  m] = (CMap)UGOLOK_CELL;
				beam->ce[0]->mp[++m] = (CMap)(A1-B1);
				beam->ce[0]->mp[++m] = (CMap)(A1-B1);
				beam->ce[0]->mp[++m] = (CMap)B1;
				beam->ce[0]->mp[++m] = (CMap)B1;
				beam->ce[0]->mp[++m] = (CMap)rr;
				beam->ce[0]->mp[++m] = (CMap)0.;
				beam->ce[0]->mp[++m] = (CMap)270.;
			}
			if (A1 > B1+rr && add_new_maps(beam->ce[1]->mp, (m = size_of_map(beam->ce[1]->mp))+1, 7)) {
				beam->ce[1]->mp[  m] = (CMap)UGOLOK_CELL;
				beam->ce[1]->mp[++m] = (CMap)(A1-B1);
				beam->ce[1]->mp[++m] = (CMap)(A1-B1);
				beam->ce[1]->mp[++m] = (CMap)B1;
				beam->ce[1]->mp[++m] = (CMap)B1;
				beam->ce[1]->mp[++m] = (CMap)rr;
				beam->ce[1]->mp[++m] = (CMap)0.;
				beam->ce[1]->mp[++m] = (CMap)270.;
			}
		}
		if ((B[k].type & ERR_MASK) == SUB_SREZKA_BLOCK && beam->ce[0]->mp[m = size_of_map(beam->ce[0]->mp)] == NULL_CELL) {
			double A1, B1, X1, Y1, X2, Y2, rr, f0;
			int k_line, k_point1, k_point2, k_circ;

			k_line   = beam->ce[0]->graph[2];
			k_point1 = beam->ce[k_line]->graph[2];
			k_point2 = beam->ce[k_line]->graph[3];

			X1 = beam->ce[k_point1]->mp[1]-beam->ce[k_point2]->mp[1]; 
			Y1 = beam->ce[k_point1]->mp[2]-beam->ce[k_point2]->mp[2]; 
			B1 = sqrt(sqr(X1)+sqr(Y1));

			k_line   = beam->ce[0]->graph[3];
			k_point1 = beam->ce[k_line]->graph[2];
			k_point2 = beam->ce[k_line]->graph[3];

			X2 = beam->ce[k_point2]->mp[1]-beam->ce[k_point1]->mp[1]; 
			Y2 = beam->ce[k_point2]->mp[2]-beam->ce[k_point1]->mp[2];
			A1 = sqrt(sqr(X2)+sqr(Y2));

			k_circ = beam->ce[0]->graph[5];
			rr = fabs(beam->ce[k_circ]->mp[7]);

			if (beam->ce[beam->ce[0]->graph[2]]->topo_id(k_point2)) { iSWAP(k_point1, k_point2); f0 = -M_PI; }
			map_iso(beam->ce[0]->mp, NULL, f0 += arg0(comp(X2, Y2))+M_PI*.75);
			map_iso(beam->ce[1]->mp, NULL, f0);
			beam->ce[0]->mp[1] = beam->ce[1]->mp[1] = beam->ce[k_point1]->mp[1]+B1*cos(f0+M_PI*.75)-A1*sin(f0+M_PI*.75);
			beam->ce[0]->mp[2] = beam->ce[1]->mp[2] = beam->ce[k_point1]->mp[2]+B1*cos(f0+M_PI*.75)+A1*sin(f0+M_PI*.75);

			if (B[k].type & ERR_CODE) dSWAP(A1, B1);
			if (add_new_maps(beam->ce[0]->mp, (m = size_of_map(beam->ce[0]->mp))+1, 7)) {
				beam->ce[0]->mp[  m] = (CMap)UGOLOK_CELL;
				beam->ce[0]->mp[++m] = (CMap)A1;
				beam->ce[0]->mp[++m] = (CMap)B1;
				beam->ce[0]->mp[++m] = (CMap)B1;
				beam->ce[0]->mp[++m] = (CMap)A1;
				beam->ce[0]->mp[++m] = (CMap)rr;
				beam->ce[0]->mp[++m] = (CMap)0.;
				beam->ce[0]->mp[++m] = (CMap)90.;
			}
			if (add_new_maps(beam->ce[1]->mp, (m = size_of_map(beam->ce[1]->mp))+1, 7)) {
				beam->ce[1]->mp[  m] = (CMap)UGOLOK_CELL;
				beam->ce[1]->mp[++m] = (CMap)B1;
				beam->ce[1]->mp[++m] = (CMap)A1;
				beam->ce[1]->mp[++m] = (CMap)A1;
				beam->ce[1]->mp[++m] = (CMap)B1;
				beam->ce[1]->mp[++m] = (CMap)rr;
				beam->ce[1]->mp[++m] = (CMap)0.;
				beam->ce[1]->mp[++m] = (CMap)90.;
			}
		}

//////////////////////////////////////////
//...преобразуем начальный блок в образце;
		B[k].bar = beam;
		set_block(B[k], SPHERE_GENUS, 1.);

//////////////////////////////////
//...линкование блочной структуры;
		int * link = B[k].link; B[k].link = NULL;
		set_link (B[k], max(NUM_PHASE, link[0]));

		if (m_cell) {
			for (m = min(NUM_PHASE-2, N_arcs), i = 0; i < m; i++)
				if (link[i+1] <= SRF_STATE) {
					B[k].link[i+3] = link[i+1]-shift;	
					for (l = NUM_PHASE; l < link[0]; l++) //...коррекция многофазных сред;
						if (link[l+1] == SRF_STATE-link[i+1]) {
							B[k].link[i+3] += shift; break;
						}
				}
				else B[k].link[i+3] = link[i+1];
		}
		else {
			for (m = min(NUM_PHASE-2, N_arcs), i = 0; i < m; i++)
				if (link[i+1] <= SRF_STATE) {
					B[k].link[m-i+2] = link[i+1]-shift;
					for (l = NUM_PHASE; l < link[0]; l++) //...коррекция многофазных сред;
						if (link[l+1] == SRF_STATE-link[i+1]) {
							B[k].link[m-i+2] += shift; break;
						}
				}
				else B[k].link[m-i+2] = link[i+1];
		}
		for (i = NUM_PHASE; i < link[0]; i++)	B[k].link[i+1] = link[i+1];
		if  (		NUM_PHASE > link[0])				B[k].link[0] = 2+link[0];
															B[k].link[NUM_PHASE] = link[NUM_PHASE];
		if (NN[0] > 0 && NN[1] > 1 || NN[0] > 1) 
		B[k].link[2] = k+N0; else B[k].link[2] = ERR_STATE;
		B[k].link[1] = ERR_STATE;
		delete_struct(link);
	}
	if (bar) {
		CCeBasic * beam = (CCeBasic *)CreateCells(BASIC_CELLS);
		beam->get_beam(m_cell ? bar->bar_cpy(m_cell) : bar->bar_cpy(), leng);
		beam->segms_bar(); 

		delete bar; bar = beam;
	}
	SkeletonRadius(OK_STATE); return(Max_phase);
}

/////////////////////////////////////////////////////
//...block structure on the base of the plane sample;
void CDraft::ConvertToCylStruct(int NN, double axis_dist, double fi0, int m_close)
{
	double pp_cond[] = {axis_dist, 0., 0.}, delta_fi = (fi0 = m_close ? 360. : fi0)/(NN = m_close ? max(NN, 2) : NN);
	int k, i, j, N0 = N;

//////////////////////////////////////////
//...construction of the blocks structure;
	for (k = 0; k < N0; k++) { B[k].bar->cells_iso(pp_cond); B[k].bar->bar_invers();}
	for (j = 1; j < NN; j++)
	for (k = 0; k < N0; k++) {

/////////////////////////////////
//...образование геометрии блока;
		CCeBasic * beam = (CCeBasic *)CreateCells(BASIC_CELLS);

		beam->get_cyl_beam(B[k].bar->bar_cpy(), delta_fi);
		beam->segms_bar();
		beam->cells_iso(NULL, 0., delta_fi*(j+.5)*M_PI/180.);

////////////////////////////////////
//...достраиваем торцевые сегменты;

/////////////////////////////////
//...образование блока в образце;
		add_block(B[k].type);
		B[k+j*N0].bar = beam;
		set_block(B[k+j*N0], SPHERE_GENUS, 1.);

//////////////////////////////////
//...линкование блочной структуры;
		set_link (B[k+j*N0], max(NUM_PHASE, B[k].link[0]));
		if (j > 1)
		B[k+(j-1)*N0].link[2] = k+(j+0)*N0; 
		B[k+(j+0)*N0].link[1] = k+(j-1)*N0; 
		for (i = 0; i < B[k].link[0]; i++) { 
			if (B[k].link[i+1] >= 0)			B[k+j*N0].link[i+3] = B[k].link[i+1]+j*N0; else
			if (B[k].link[i+1] <= SRF_STATE) B[k+j*N0].link[i+3] = B[k].link[i+1]-2;	
			else										B[k+j*N0].link[i+3] = B[k].link[i+1];
		}
	}

///////////////////////////////////////////////////////////////
//...преобразование блоков в нулевом слое и окаймляющего блока;
	for (k = 0; k < N0; k++) {
		CCeBasic * beam = (CCeBasic *)CreateCells(BASIC_CELLS);

		beam->get_cyl_beam(B[k].bar->bar_cpy(), delta_fi);
		beam->segms_bar(); 
		beam->cells_iso(NULL, 0., delta_fi*.5*M_PI/180.); delete B[k].bar;

		B[k].bar = beam;
		set_block(B[k], SPHERE_GENUS, 1.);
		for (i = B[k].link[0]; i > 0 ; i--) 
		B[k].link[i+2] = B[k].link[i];
		if (NN > 1) {
			B[k].link[2] = k+N0; 
			if (m_close) {
				B[k].link[1] = k+(NN-1)*N0; 
				B[k+(NN-1)*N0].link[2] = k; 
			}
		}
	}
	CCeBasic * beam = (CCeBasic *)CreateCells(BASIC_CELLS);
	beam->get_cyl_beam(bar->bar_cpy(), fi0);
	beam->segms_bar(); 
	beam->cells_iso(NULL, 0., fi0*.5*M_PI/180.);

	delete bar; bar = beam;
	SkeletonRadius();
}

/////////////////////////////////////////
//...block structure for cylinder sample;
int CDraft::GetCylinderStruct(double rad, double R, double length, int NX, int NY, int NZ)
{
	double delta_fi = M_PI*2./(NX = max(NX, 2)), delta_rr = (R-rad)/(NY = max(NY, 1));
	CCeBasic * ce; 

//////////////////////
//...обнуляем образец;
   init_blocks(NULL);

//////////////////////////////////////
//...добавлям описание плоского торца;
	ce = new CCeBasic;
	ce->get_ring(rad, R);
	bar = new CCells;	
	bar->bar_add(ce);

///////////////////////////////////////////////
//...образуем блочную структуру плоского торца;
	for (int j = 0; j < NY;  j++)
	for (int i = 0; i < NX;  i++) {
		add_block (NULL_BLOCK);
		ce = new CCeBasic;
		ce->get_ring_segment(i*delta_fi, (i+1)*delta_fi, rad+j*delta_rr, rad+(j+1)*delta_rr);
		B[i+j*NX].bar = new CCells;
		B[i+j*NX].bar->bar_add(ce);
		set_block(B[i+j*NX], SPHERE_GENUS, 1.);
		set_link (B[i+j*NX], 4);
		B[i+j*NX].link[1] = (i-1+NX)%NX+j*NX;
		B[i+j*NX].link[2] = (j+1)%NY ? i+(j+1)*NX : SRF_STATE-1;
		B[i+j*NX].link[3] = (i+1)%NX+j*NX;
		B[i+j*NX].link[4] = j ? i+(j-1)*NX : SRF_STATE;
	}
	SkeletonRadius(OK_STATE);
///////////////////////////////////
//...преобразуем геометрию образца;
	int N0 = N, NN[2] = {1, NZ}, 
		 N_phase = ConvertToBeamStruct(NN, &length, OK_STATE);
	return(N0);
}

/////////////////////////////////////////
//...block structure for cylinder sample;
int CDraft::GetCylinderChannelStruct(double rad, double R, double length, double r0, int NX, int NY, int N2, int NZ)
{
	if (rad >= R) return(0); 

	double delta_rr = (R-rad )/(NY = max(NY, 1)), delta_fi = M_PI*2./(NX = max(NX, 2)), 
			 delta_r0 = (rad-r0)/(N2 = max(N2, 1)), RR = rad;
	CCeBasic * ce;

//////////////////////
//...обнуляем образец;
	init_blocks(NULL);

////////////////////////////////////////////
//...добавлям описание структуры излучателя;
	ce = new CCeBasic;
	ce->get_ring(rad, R);
	bar = new CCells;
	bar->bar_add(ce);

//	ce = new CCeBasic;
//	ce->get_disk(rad);
//	bar = new CCells;
//	bar->bar_add(ce);

//////////////////////////////////////////////////////
//...образуем блочную структуру стенки плоского торца;
	int i, j;
	for (j = 0; j < NY;  j++)
	for (i = 0; i < NX;  i++) {
		add_block (NULL_BLOCK);
		ce = new CCeBasic;
		ce->get_ring_segment(i*delta_fi, (i+1)*delta_fi, rad+j*delta_rr, rad+(j+1)*delta_rr);
		B[i+j*NX].bar = new CCells;
		B[i+j*NX].bar->bar_add(ce);
		set_block(B[i+j*NX], SPHERE_GENUS, 1.);
		set_link (B[i+j*NX], 4);
		B[i+j*NX].link[1] = (i-1+NX)%NX+j*NX;
		B[i+j*NX].link[2] = (j+1)%NY ? i+(j+1)*NX : BOUND3_STATE;
		B[i+j*NX].link[3] = (i+1)%NX+j*NX;
		B[i+j*NX].link[4] = j ? i+(j-1)*NX : SRF_STATE;
	}

//////////////////////////////////////////////////////////
//...образуем блочную структуру излучателя плоского торца;
	int shift = 0, m0 = NX*NY;
	if (r0 < rad) {
		shift = m0; m0 += NX*N2; RR = r0;
		for (j = 0; j < N2;  j++)
		for (i = 0; i < NX;  i++) {
			add_block (NULL_BLOCK);
			ce = new CCeBasic;
			ce->get_ring_segment(i*delta_fi, (i+1)*delta_fi, r0+j*delta_r0, r0+(j+1)*delta_r0);
			B[shift+i+j*NX].bar = new CCells;
			B[shift+i+j*NX].bar->bar_add(ce);
			set_block(B[shift+i+j*NX], SPHERE_GENUS, 1.);
			set_link (B[shift+i+j*NX], 4);
			B[shift+i+j*NX].link[1] = shift+(i-1+NX)%NX+j*NX;
			if ((j+1)%N2) B[shift+i+j*NX].link[2] = shift+i+(j+1)*NX;
			else {
				B[shift+i+j*NX].link[2] = i;
				B[i].link[4] = shift+i+j*NX;
			}
			B[shift+i+j*NX].link[3] = shift+(i+1)%NX+j*NX;
			B[shift+i+j*NX].link[4] = j ? shift+i+(j-1)*NX : SRF_STATE;
		}
	}

////////////////////////////////////////////////
//...образуем центральный блок на плоском торце;
	add_block (NULL_BLOCK);
	ce = new CCeBasic;
	for (i = 0; i < NX;  i++) {
		CCeBasic * arc = new CCeBasic;
		arc->get_arc(RR, i*delta_fi, (i+1)*delta_fi);
		ce->bar_add(arc);
	}
	ce->bar_invers();
	ce->bar_generate();
	CMap * mp = get_map(2, NULL_GENUS, RING_SEGMENT);
	if (mp) {
		mp[++(i = size_of_map(2, NULL_GENUS))] = (CMap)M_PI*2.;
		mp[++i] = (CMap)0.;
		mp[++i] = (CMap)RR;
	}
	ce->bar_span(mp);

	B[m0].bar = new CCells;
	B[m0].bar->bar_add(ce);
	set_block(B[m0], SPHERE_GENUS, 1.);
	set_link (B[m0], NX);

	for (i = 0; i < NX; i++) {
		B[m0].link[i+1] = shift+i;
		B[shift+i].link[4] = m0;
	}
	SkeletonRadius(OK_STATE);

///////////////////////////////////
//...преобразуем геометрию образца;
	int N0 = N, NN[2] = {1, NZ}, 
		 N_phase = ConvertToBeamStruct(NN, &length, OK_STATE);

//////////////////////////////////////////////////////////
//...раставляем граничные условия на нижнем плоском торце;
	for (j = 0; j < NY;  j++)
	for (i = 0; i < NX;  i++) 
		B[i+j*NX].link[1] = BOUND2_STATE;

	if (r0 < rad)
	for (j = 0; j < N2;  j++)
	for (i = 0; i < NX;  i++) 
		B[shift+i+j*NX].link[1] = BOUND1_STATE;

	B[m0].link[1] = BOUND1_STATE;

	return(N0);
}

/////////////////////////////////////////////////////////////////////////
//...block structure for cylinder sample on the base nastran description;
int CDraft::GetCylinderChannelStruct(double rad, double R, double length, char * name, int NZ, int id_curve)
{
	if (rad >= R) return(0); 

////////////////////////////////////////
//...зачитываем блочную структуру торца;
	stru.nodes_in(name);
   bar_condit_in(name);
	LinkUniStruct();
	BlockActivate(NULL_STATE);
//	B[0].bar->cells_out("bar");

////////////////////////////////////////////
//...корректируем криволинейные дуги блоков;
	if (id_curve) {
		int k, j, m, l;
		for (k = 0; k < N; k++)
		for (j = B[k].link[0]; j > 0; j--)
		if ((m = B[k].link[j]) >= 0 && B[k].link[NUM_PHASE] != B[m].link[NUM_PHASE] && B[k].bar && B[k].bar->ce[0] && B[k].bar->ce[0]->graph) {
			for (l = 0; l < NUM_PHASE; l++)
				if (B[k].link[l+1] == -m+SRF_STATE) break;
				if (l < B[k].bar->ce[0]->graph[1]) {
					B[k].bar->ce[B[k].bar->ce[0]->graph[l+2]]->cells_to(SPHERE_GENUS);
					B[k].bar->ce[B[k].bar->ce[0]->graph[l+2]]->circ_correct(rad, B[k].link[NUM_PHASE] == -1 ? -1 : 1, NULL_STATE);
				}
		}
		for (k = 0; k < N; k++) if (B[k].bar && B[k].bar->ce[0] && B[k].bar->ce[0]->graph)
		for (j = min(B[k].bar->ce[0]->graph[1], NUM_PHASE); j > 0; j--)
		if (B[k].link[j] == ERR_STATE) {
			B[k].bar->ce[B[k].bar->ce[0]->graph[j+1]]->cells_to(SPHERE_GENUS);
			B[k].bar->ce[B[k].bar->ce[0]->graph[j+1]]->circ_correct(R, 1, NULL_STATE);
		}
	}

////////////////////////////////////////////////////////////////////////////
//...преобразуем геометрию образца и добавлям описание структуры излучателя;
	int N0 = N, NN[2] = {1, NZ}, 
		 N_phase = ConvertToBeamStruct(NN, &length, OK_STATE);
	add_cyl_surface(rad, 0., 0., 0.);

///////////////////////////////////////////////////////////////////////////////////////////
//...раставляем граничные условия на нижнем торце и корректируем связи (убираем включение);
	for (int k = 0; k < N; k++) {
		if (B[k].link[NUM_PHASE] == -1 && B[k].link[1] <= ERR_STATE) B[k].link[1] = BOUND2_STATE;
		if (B[k].link[NUM_PHASE] == -2 && B[k].link[1] <= ERR_STATE) B[k].link[1] = BOUND1_STATE;

		for (int m = min(NUM_PHASE, B[k].link[0]), i = 0; i < m; i++) if (B[k].link[i+1] <= SRF_STATE)
		for (int l = NUM_PHASE; l < B[k].link[0]; l++) //...коррекция многофазных сред;
		if (B[k].link[i+1] == SRF_STATE-B[k].link[l+1]) {
			 B[k].link[i+1] = B[k].link[l+1]; break;
		}
		B[k].link[0] = B[k].bar->arcs_number();
	}
	return(N0);
}

////////////////////////////////////
//...block structure for box sample;
int CDraft::GetBoxStruct(double AA, double BB, double length, int NX, int NY, int NZ)
{
	double delta_A = AA/(NX = max(NX, 1)), delta_B = BB/(NY = max(NY, 1)), P0[3];
	CCeBasic * ce;

//////////////////////
//...обнуляем образец;
	init_blocks(NULL);

//////////////////////////////////////
//...добавлям описание плоского торца;
	ce = new CCeBasic;
	ce->get_sheet(AA, BB);
	bar = new CCells;	
	bar->bar_add(ce);

///////////////////////////////////////////////
//...образуем блочную структуру плоского торца;
	for (int j = 0; j < NY;  j++)
	for (int i = 0; i < NX;  i++) {
		add_block (NULL_BLOCK);
		ce = new CCeBasic;
		ce->get_sheet(delta_A, delta_B);
		P0[0] = (i+.5)*delta_A-AA*.5;
		P0[1] = (j+.5)*delta_B-BB*.5;
		P0[2] = 0.;
		ce->cells_iso(P0);
		B[i+j*NX].bar = new CCells;
		B[i+j*NX].bar->bar_add(ce);
		set_block(B[i+j*NX], SPHERE_GENUS, 1.);
		set_link (B[i+j*NX], 4);
		B[i+j*NX].link[1] = j ? i+(j-1)*NX : BOUND4_STATE;
		B[i+j*NX].link[2] = (i+1)%NX ? (i+1)+j*NX : BOUND5_STATE;
		B[i+j*NX].link[3] = (j+1)%NY ? i+(j+1)*NX : BOUND6_STATE;
		B[i+j*NX].link[4] = i ? (i-1)+j*NX : BOUND7_STATE;
	}
	SkeletonRadius(OK_STATE);
///////////////////////////////////
//...преобразуем геометрию образца;
	int N0 = N, NN[2] = {1, NZ}, 
		 N_phase = ConvertToBeamStruct(NN, &length, OK_STATE);
	return(N0);
}

////////////////////////////////////////////////////////////////////////
//...block structure for Ono Box sample on the base nastran description;
int CDraft::GetOnoBoxStruct(double AA, double BB, double aa, double bb, double rad, int * NN, double * LL, char * name, int id_curve, int id_phase)
{
////////////////////////////////////////
//...зачитываем блочную структуру торца;
	stru.nodes_in(name);
   bar_condit_in(name);
	LinkUniStruct();
	BlockActivate(NULL_STATE);
//	B[0].bar->cells_out("bar");

////////////////////////////////////////////
//...корректируем криволинейные дуги блоков;
	if (id_curve) {
		double eps = 1e-6/*EE*/;
		for (int j, m, l, k = 0; k < N; k++) 
		for (j = B[k].link[0]; j > 0; j--) 
		if ((m = B[k].link[j]) >= 0 && B[k].link[NUM_PHASE] != B[m].link[NUM_PHASE] && B[k].bar && B[k].bar->ce[0] && B[k].bar->ce[0]->graph) {
			for (l = 0; l < NUM_PHASE; l++) 
				if (B[k].link[l+1] == -m+SRF_STATE) break;
				if (l < B[k].bar->ce[0]->graph[1]) {
					int m1 = get_num(B[k].bar->ce[B[k].bar->ce[0]->graph[l+2]]->graph, 0),
						 m2 = get_num(B[k].bar->ce[B[k].bar->ce[0]->graph[l+2]]->graph, 1);
					double X1 = B[k].bar->ce[m1]->mp[1], Y1 = B[k].bar->ce[m1]->mp[2],
							 X2 = B[k].bar->ce[m2]->mp[1], Y2 = B[k].bar->ce[m2]->mp[2], 
							 X0 = (AA-aa)*.5+rad, Y0 = (BB-bb)*.5+rad, 
							 XX = (AA+aa)*.5-rad, YY = (BB+bb)*.5-rad;
					if (fabs(sqr(X1-X0)+sqr(Y1-Y0)-sqr(rad)) < eps && fabs(sqr(X2-X0)+sqr(Y2-Y0)-sqr(rad)) < eps ||
						 fabs(sqr(X1-X0)+sqr(Y1-YY)-sqr(rad)) < eps && fabs(sqr(X2-X0)+sqr(Y2-YY)-sqr(rad)) < eps ||
						 fabs(sqr(X1-XX)+sqr(Y1-Y0)-sqr(rad)) < eps && fabs(sqr(X2-XX)+sqr(Y2-Y0)-sqr(rad)) < eps ||
						 fabs(sqr(X1-XX)+sqr(Y1-YY)-sqr(rad)) < eps && fabs(sqr(X2-XX)+sqr(Y2-YY)-sqr(rad)) < eps) {
						B[k].bar->ce[B[k].bar->ce[0]->graph[l+2]]->cells_to(SPHERE_GENUS);
						B[k].bar->ce[B[k].bar->ce[0]->graph[l+2]]->circ_correct(rad, B[k].link[NUM_PHASE] == -1 ? -1 : 1, NULL_STATE);
					}
				}
		}
	}

////////////////////////////////////////////////////////////////////////////
//...преобразуем геометрию образца и добавлям описание структуры излучателя;
	int N0 = N, 
		 N_phase = ConvertToBeamStruct(NN, LL, OK_STATE, id_phase);
	if (id_curve) add_cyl_surface(rad, 0., 0., 0.);

///////////////////////////////////////////////////////////////////////////////////////////
//...раставляем граничные условия на нижнем торце и корректируем связи (убираем включение);
	int  k;
	for (k = 0; k < N; k++) {
		if (B[k].link[1] <= ERR_STATE && B[k].link[NUM_PHASE] == -2) B[k].link[1] = BOUND1_STATE;
		if (B[k].link[2] <= ERR_STATE)	B[k].link[2] = BOUND3_STATE;

		for (int m = min(NUM_PHASE, B[k].link[0]), i = 0; i < m; i++) if (B[k].link[i+1] <= SRF_STATE)
		for (int l = NUM_PHASE; l < B[k].link[0]; l++) //...коррекция многофазных сред;
		if (B[k].link[i+1] == SRF_STATE-B[k].link[l+1]) {
			 B[k].link[i+1] = B[k].link[l+1]; break;
		}
		B[k].link[0] = NUM_PHASE;
		B[k].link[NUM_PHASE] = B[k].link[NUM_PHASE]/(1-N_phase)-1;
	}

////////////////////////////////////////////////////////////////////////////////
//...перемещение линков в конец списка (криволинейную границу не устанавливаем);
	for (k = 0; k < N; k++)
	for (int i, j = B[k].link[0]; j > 0; j--) 
	if ((i = B[k].link[j]) >= 0 && B[k].link[NUM_PHASE] != B[i].link[NUM_PHASE]) {
		B[k].link[j] =-i+SRF_STATE;
		link_add(B[k], i);
		link_add(B[k], ERR_STATE, NULL_STATE);
	}

////////////////////////////////////
//...коррекция нескольких включений;
	if (id_phase == OK_STATE)
	for (k = 0; k < N; k++)
	  if (B[k].link[NUM_PHASE] <= -3)
			B[k].link[NUM_PHASE]  = -1;

	return(N0);
}

/////////////////////////////////////////////////////////////
//...block structure for box sample with spherical inclusion;
int CDraft::GetSphBoxStruct(double AA, double BB, double CC, double rad, double ll, int id_bound)
{
	if (rad > AA*.5 || rad > BB*.5 || rad > CC*.5) return(0);
	double P1[3] = {-.5*AA, -.5*BB, -.5*CC},
			 P2[3] = {0.5*AA, -.5*BB, -.5*CC},
			 P3[3] = {0.5*AA, 0.5*BB, -.5*CC},
			 P4[3] = {-.5*AA, 0.5*BB, -.5*CC},
			 P5[3] = {-.5*AA, -.5*BB, 0.5*CC},
			 P6[3] = {0.5*AA, -.5*BB, 0.5*CC},
			 P7[3] = {0.5*AA, 0.5*BB, 0.5*CC},
			 P8[3] = {-.5*AA, 0.5*BB, 0.5*CC};
	CCeBasic * ce;

//////////////////////
//...обнуляем образец;
	int k = -1;
	init_blocks(NULL);

//////////////////////////////////////////////////////////////
//...образуем описание первого блока -- сферическое включение;
	if (rad > 0 && id_bound != OK_STATE) {
		add_block (NULL_BLOCK);
		B[++k].bar = new CCells;
		B[k].type  = POLY_BLOCK;

		ce = new CCeBasic;
		ce->get_sphere(rad);
		B[k].bar->bar_add(ce);
		B[k].bar->bar_ord();
		
		set_block(B[k], SPHERE_GENUS, rad);
		set_link (B[k], NUM_PHASE);
		B[k].link[1] = ll > 0. ? k+2 : k+1;
		B[k].link[NUM_PHASE] = -2;
	}

////////////////////////////////////////////////////////////////////////////////////////////
//...образуем описание второго блока -- параллелепипедная матрица со сферическим включением;
	add_block (NULL_BLOCK);
	B[++k].bar = new CCells;
	B[k].type  = POLY_BLOCK;

	if (rad > 0) {
		ce = new CCeBasic;
		ce->get_sphere(-rad-ll);
		B[k].bar->bar_add(ce);
		B[k].type = ZOOM_BLOCK;
	}

	ce = new CCeBasic;
	ce->get_quad_facet(P1, P5, P8, P4);
	B[k].bar->bar_add(ce);

	ce = new CCeBasic;
	ce->get_quad_facet(P2, P3, P7, P6);
	B[k].bar->bar_add(ce);

	ce = new CCeBasic;
	ce->get_quad_facet(P1, P2, P6, P5);
	B[k].bar->bar_add(ce);

	ce = new CCeBasic;
	ce->get_quad_facet(P4, P8, P7, P3);
	B[k].bar->bar_add(ce);

	ce = new CCeBasic;
	ce->get_quad_facet(P1, P4, P3, P2);
	B[k].bar->bar_add(ce);

	ce = new CCeBasic;
	ce->get_quad_facet(P5, P6, P7, P8);
	B[k].bar->bar_add(ce);
	B[k].bar->bar_ord();

	set_block(B[k], SPHEROID_GENUS, sqrt(AA*AA+BB*BB+CC*CC)*.5*get_param(NUM_MPLS+1));
	B[k].mp[8] = rad+ll;
	set_link (B[k], NUM_PHASE);
	B[k].link[1] = rad > 0. ? (ll > 0. ? k+1 : k-1) : ERR_STATE;
	B[k].link[NUM_PHASE] = -1;

////////////////////////////////////////////////
//...добавляем при необходимости межфазный слой;
	if (rad > 0 && ll > 0.) {
		add_block (NULL_BLOCK);
		B[++k].bar = new CCells;
		B[k].type = DOUBLE_ZOOM_BLOCK;

		ce = new CCeBasic;
		ce->get_sphere(rad+ll);
		B[k].bar->bar_add(ce);

		ce = new CCeBasic;
		ce->get_sphere(-rad);
		B[k].bar->bar_add(ce);
		B[k].bar->bar_ord();
		
		set_block(B[k], SPHEROID_GENUS, rad+ll);
		B[k].mp[8] = rad;
		set_link (B[k], NUM_PHASE);
		B[k].link[1] = k-1;
		B[k].link[2] = k-2;
		B[k].link[NUM_PHASE] = id_bound == OK_STATE ? -2 : -3;
	}
	return(N);
}

///////////////////////////////////////////////////////////////////
//...block structure for spherical sample with spherical inclusion;
int CDraft::GetSph2BodyStruct(double ff_vol, double rad, int id_bound)
{
	if (ff_vol < EE || ff_vol > 1.-EE) return(0);
	double R1 = pow(1./(1.-ff_vol), 1./3.)*rad;
	CCeBasic * ce;

//////////////////////
//...обнуляем образец;
	int k = -1;
	init_blocks(NULL);

//////////////////////////////////////////////////////////////
//...образуем описание первого блока -- сферическое включение;
	if (rad > 0 && id_bound != OK_STATE) {
		add_block (NULL_BLOCK);
		B[++k].bar = new CCells;
		B[k].type  = POLY_BLOCK;

		ce = new CCeBasic;
		ce->get_sphere(rad);
		B[k].bar->bar_add(ce);
		B[k].bar->bar_ord();
		
		set_block(B[k], SPHERE_GENUS, rad);
		set_link (B[k], NUM_PHASE);
		B[k].link[1] = k+1;
		B[k].link[NUM_PHASE] = -2;
	}

//////////////////////////////////////////////////////////////////////////////////////
//...образуем описание второго блока -- сферическая матрица со сферическим включением;
	add_block (NULL_BLOCK);
	B[++k].bar = new CCells;
	B[k].type  = POLY_BLOCK;

	if (rad > 0) {
		ce = new CCeBasic;
		ce->get_sphere(-rad);
		B[k].bar->bar_add(ce);
		B[k].type = ZOOM_BLOCK;
	}
	ce = new CCeBasic;
	ce->get_sphere(R1);
	B[k].bar->bar_add(ce);
	B[k].bar->bar_ord();

	set_block(B[k], SPHEROID_GENUS, R1);
	B[k].mp[8] = rad;
	set_link (B[k], NUM_PHASE);
	B[k].link[1] = rad > 0. ? k-1 : ERR_STATE;
	B[k].link[NUM_PHASE] = -1;

	return(N);
}

///////////////////////////////////////////////////////////////////
//...block structure for cylinder sample with spherical inclusion;
int CDraft::GetCylSphStruct(double ff_vol, double rad, double L, int id_bound)
{
	if (ff_vol < EE || ff_vol > 1.-EE) return(0);
	double R1 = sqrt(rad/(3.*L*(1.-ff_vol)))*2.*rad;
	CCeBasic * ce;

//////////////////////
//...обнуляем образец;
	int k = -1;
	init_blocks(NULL);

//////////////////////////////////////////////////////////////
//...образуем описание первого блока -- сферическое включение;
	if (rad > 0 && id_bound != OK_STATE) {
		add_block (NULL_BLOCK);
		B[++k].bar = new CCells;
		B[k].type  = POLY_BLOCK;

		ce = new CCeBasic;
		ce->get_sphere(rad);
		B[k].bar->bar_add(ce);
		B[k].bar->bar_ord();
		
		set_block(B[k], SPHERE_GENUS, rad);
		set_link (B[k], NUM_PHASE);
		B[k].link[1] = k+1;
		B[k].link[NUM_PHASE] = -2;
	}

/////////////////////////////////////////////////////////////////////////////////////////
//...образуем описание второго блока -- цилиндрическая матрица со сферическим включением;
	add_block (NULL_BLOCK);
	B[++k].bar = new CCells;
	B[k].type  = POLY_BLOCK;

	if (rad > 0) {
		ce = new CCeBasic;
		ce->get_sphere(-rad);
		B[k].bar->bar_add(ce);
		B[k].type = ZOOM_BLOCK;
	}

///////////////////////////////////////////////
//...добавляем два торца и боковую поверхность;
	ce = new CCeBasic;
	ce->get_ring_segment(-M_PI, M_PI, R1, 0.);
	ce->mp[3] = -L*.5;
	B[k].bar->bar_add(ce);

	ce = new CCeBasic;
	ce->get_ring_segment(-M_PI, M_PI, 0., R1);
	ce->mp[3] = L*.5;
	B[k].bar->bar_add(ce);

	ce = new CCeBasic;
	//ce->get_cyl_segment(R1, -M_PI, M_PI, -L*.5, L*.5);
	ce->get_cylinder(R1);
	B[k].bar->bar_add(ce);
	B[k].bar->bar_ord();

	set_block(B[k], SPHEROID_GENUS, sqrt(R1*R1+L*L*.25)*get_param(NUM_MPLS+1));
	B[k].mp[8] = rad;
	set_link (B[k], NUM_PHASE);
	B[k].link[1] = rad > 0. ? k-1 : ERR_STATE;
	B[k].link[NUM_PHASE] = -1;

	return(N);
}


///////////////////////////////////////////////////////////////////////////////////
//...block structure for cylinder sample with spherical inclusion (equal distance);
int CDraft::GetCylSphStruct(double ff_vol, double rad, int id_bound)
{
	if (ff_vol < EE || ff_vol > 1.-EE) return(0);
	double R1 = pow(4./(3.*(1.-ff_vol)), 1./3.)*rad;
	CCeBasic * ce;

//////////////////////
//...обнуляем образец;
	int k = -1;
	init_blocks(NULL);

//////////////////////////////////////////////////////////////
//...образуем описание первого блока -- сферическое включение;
	if (rad > 0 && id_bound != OK_STATE) {
		add_block (NULL_BLOCK);
		B[++k].bar = new CCells;
		B[k].type  = POLY_BLOCK;

		ce = new CCeBasic;
		ce->get_sphere(rad);
		B[k].bar->bar_add(ce);
		B[k].bar->bar_ord();
		
		set_block(B[k], SPHERE_GENUS, rad);
		set_link (B[k], NUM_PHASE);
		B[k].link[1] = k+1;
		B[k].link[NUM_PHASE] = -2;
	}

/////////////////////////////////////////////////////////////////////////////////////////
//...образуем описание второго блока -- цилиндрическая матрица со сферическим включением;
	add_block (NULL_BLOCK);
	B[++k].bar = new CCells;
	B[k].type  = POLY_BLOCK;

	if (rad > 0) {
		ce = new CCeBasic;
		ce->get_sphere(-rad);
		B[k].bar->bar_add(ce);
		B[k].type = ZOOM_BLOCK;
	}

///////////////////////////////////////////////
//...добавляем два торца и боковую поверхность;
	ce = new CCeBasic;
	ce->get_ring_segment(-M_PI, M_PI, R1, 0.);
	ce->mp[3] = -R1*.5;
	B[k].bar->bar_add(ce);

	ce = new CCeBasic;
	ce->get_ring_segment(-M_PI, M_PI, 0., R1);
	ce->mp[3] = R1*.5;
	B[k].bar->bar_add(ce);

	ce = new CCeBasic;
	//ce->get_cyl_segment(R1, -M_PI, M_PI, -R1, R1);
	ce->get_cylinder(R1);
	B[k].bar->bar_add(ce);
	B[k].bar->bar_ord();

	set_block(B[k], SPHEROID_GENUS, M_SQRT2*R1*get_param(NUM_MPLS+1));
	B[k].mp[8] = rad;
	set_link (B[k], NUM_PHASE);
	B[k].link[1] = rad > 0. ? k-1 : ERR_STATE;
	B[k].link[NUM_PHASE] = -1;

	return(N);
}

///////////////////////////////////////////////////////
//...block structure for penetrate spherical inclusion;
int CDraft::GetPenetrateSphere(double rad, double L)
{
	if (rad < L || rad > M_SQRT2*L) return(0);
	double P1[3] = {-L, 0., 0.},
			 P2[3] = {0., -L, 0.},
			 P3[3] = {0., 0., -L}, RR = sqrt(rad*rad-L*L);
	CCeBasic * ce;

//////////////////////
//...обнуляем образец;
	int k = -1;
	init_blocks(NULL);

/////////////////////////////////////////////////////////////////////////////////////
//...образуем описание одного блока -- сферическое включение с проникновением в кубе;
	add_block (NULL_BLOCK);
	B[++k].bar = new CCells;
	B[k].type  = DOUBLE_ZOOM_BLOCK;

////////////////////////////////////////////////////////////////
//...образуем геометрию сферического включения с проникновением;
	ce = new CCeBasic;
	ce->get_sph_intrusion(-rad, L);
	B[k].bar->bar_add(ce);

/////////////////////////////////////////////////
//...образуем геометрию боковых граней с дырками;
	ce = new CCeBasic;
	ce->get_sheet_intrusion(2.*L, 2.*L, RR); ce->cells_iso(P1, 0., -M_PI_2); P1[0] = -P1[0];
	B[k].bar->bar_add(ce);

	ce = new CCeBasic;
	ce->get_sheet_intrusion(2.*L, 2.*L, RR); ce->cells_iso(P1, 0., M_PI_2);
	B[k].bar->bar_add(ce);

	ce = new CCeBasic;
	ce->get_sheet_intrusion(2.*L, 2.*L, RR); ce->cells_iso(P2, M_PI_2, -M_PI_2); P2[1] = -P2[1];
	B[k].bar->bar_add(ce);

	ce = new CCeBasic;
	ce->get_sheet_intrusion(2.*L, 2.*L, RR); ce->cells_iso(P2, M_PI_2, M_PI_2);
	B[k].bar->bar_add(ce);

	ce = new CCeBasic;
	ce->get_sheet_intrusion(2.*L, 2.*L, RR); ce->cells_iso(P3, 0., M_PI); P3[2] = -P3[2];
	B[k].bar->bar_add(ce);

	ce = new CCeBasic;
	ce->get_sheet_intrusion(2.*L, 2.*L, RR); ce->cells_iso(P3, 0., 0.);
	B[k].bar->bar_add(ce);
	B[k].bar->bar_ord();

	set_block(B[k], SPHEROID_GENUS, sqrt(3.)*L*get_param(1));
	B[k].mp[8] = rad;
	set_link (B[k], NUM_PHASE);
	B[k].link[NUM_PHASE] = -1;

	return(N);
}

///////////////////////////////////////
//...block structure for Bars25 sample;
int CDraft::GetBars25Struct(CDraft * bars25, double length, int NZ)
{
//////////////////////
//...обнуляем образец;
   init_blocks(NULL);

//////////////////////////////////////
//...добавлям описание плоского торца;
	bar = bars25->bar; bars25->bar = NULL;

///////////////////////////////////////////////
//...образуем блочную структуру плоского торца;
	for (int k = 0; k < bars25->N;  k++) {
		add_block (bars25->B[k].type);
		B[k].bar = bars25->B[k].bar; bars25->B[k].bar = NULL;
		B[k].bar->cells_common(bars25->B[k].mp);
		set_block(B[k], SPHERE_GENUS, 1.);
		set_link (B[k], bars25->B[k].link[0]);
		memcpy(B[k].link, bars25->B[k].link, (bars25->B[k].link[0]+1)*sizeof(Topo)); 
	}
	SkeletonRadius(OK_STATE);

///////////////////////////////////
//...преобразуем геометрию образца;
	int N0 = N, NN[2] = {1, NZ}, 
		 N_phase = ConvertToBeamStruct(NN, &length, OK_STATE);
	return(N0);
}

////////////////////////////////////////////////////////////////
//...добавление сферической поверхности для коррекции квадратур;
void CDraft::add_sph_surface(double R, double X0, double Y0, double Z0)
{
	if (! bar) bar = new(CCells);
	CCells	* ce  = new(CCells);
	int l = size_of_map(2, SPHERE_GENUS);
	ce->cells_new(1, 2, l+1);
	ce->mp[0] = (CMap)ID_MAP(2, SPHERE_GENUS);
	ce->mp[1] = (CMap)X0;
	ce->mp[2] = (CMap)Y0;
	ce->mp[3] = (CMap)Z0;
	ce->mp[7] = (CMap)R;
	ce->mp[l] = (CMap)NULL_CELL;
	bar->bar_add(ce);
	bar->bar_ord();
}

///////////////////////////////////////////////////////////////////
//...добавление цилиндрической поверхности для коррекции квадратур;
void CDraft::add_cyl_surface(double R, double X0, double Y0, double Z0, int axis, double fi, double theta)
{
	if (! bar) bar = new(CCells);
	CCells	* ce  = new(CCells);
	int l = size_of_map(2, CYL_GENUS);
	ce->cells_new(1, 2, l+1);
	ce->mp[0] = (CMap)ID_MAP(2, CYL_GENUS);
	ce->mp[4] = (CMap)(axis == AXIS_Y ? M_PI_2 : 0.);
	ce->mp[5] = (CMap)(axis == AXIS_Y || axis == AXIS_X ? M_PI_2 : 0.);
	ce->cells_iso(NULL, fi, theta);
	ce->mp[1] = (CMap)X0;
	ce->mp[2] = (CMap)Y0;
	ce->mp[3] = (CMap)Z0;
	ce->mp[7] = (CMap)R;
	ce->mp[l] = (CMap)NULL_CELL;
	bar->bar_add(ce);
	bar->bar_ord();
}

/////////////////////////////////////////
//...hit block structure for beam sample;
void CDraft::hit_beam_struct(CGrid * nd, int N0, int * NN, double * LL, int axis)
{
	if (! NN || ! LL || NN[0] <= 0) return;

	double dd, leng;

	if (axis == AXIS_X) {
		for (int i = 0; i < nd->N;  i++)
		for (int j = 0; j < nd->N1; j++) {
			int k, m, n, l;
			for (leng = 0., l = 0,			  n = 0; n < NN[0]; leng = LL[n], n++)
			for (dd = (LL[n]-leng)/NN[n+1], m = 0; m < NN[n+1]; l++, m++)
			if (leng+m*dd <= nd->Y[j] && nd->Y[j] < leng+(m+1.)*dd+EE_ker) {
				n = NN[0]; break;
			}
			for (k = 0; k < N0; k++)
				if (B[k].bar->ce[0]->in_bar_cell(nd->Z[0], nd->X[i])) break;

			if (k == N0) nd->hit[i+j*nd->N] = ERR_STATE;
			else			 nd->hit[i+j*nd->N] = k+l*N0;
		}
	}
	else
	if (axis == AXIS_Y) {
		for (int i = 0; i < nd->N;  i++)
		for (int j = 0; j < nd->N1; j++) {
			int k, m, n, l;
			for (leng = 0., l = 0,			  n = 0; n < NN[0]; leng = LL[n], n++)
			for (dd = (LL[n]-leng)/NN[n+1], m = 0; m < NN[n+1]; l++, m++)
			if (leng+m*dd <= nd->X[i] && nd->X[i] < leng+(m+1.)*dd+EE_ker) {
				n = NN[0]; break;
			}
			for (k = 0; k < N0; k++)
				if (B[k].bar->ce[0]->in_bar_cell(nd->Y[j], nd->Z[0])) break;

			if (k == N0) nd->hit[i+j*nd->N] = ERR_STATE;
			else			 nd->hit[i+j*nd->N] = k+l*N0;
		}
	}
	else
	if (axis == AXIS_Z) {
		int k, m, n, l;
		for (leng = 0., l = 0,			  n = 0; n < NN[0]; leng = LL[n], n++)
		for (dd = (LL[n]-leng)/NN[n+1], m = 0; m < NN[n+1]; l++, m++)
		if (leng+m*dd <= nd->Z[0] && nd->Z[0] < leng+(m+1.)*dd+EE_ker) {
			n = NN[0]; break;
		}
		for (int i = 0; i < nd->N;  i++)
		for (int j = 0; j < nd->N1; j++) {
			for (k = 0; k < N0; k++)
				if (B[k].bar->ce[0]->in_bar_cell(nd->X[i], nd->Y[j])) break;

			if (k == N0) nd->hit[i+j*nd->N] = ERR_STATE;
			else			 nd->hit[i+j*nd->N] = k+l*N0;
		}
	}
}

////////////////////////////////////////////
//...корректируем криволинейные дуги блоков;
void CDraft::CurveCorrect(double X0, double Y0, double rad)
{
	double eps = 1e-5;
	for (int j, m, l, k = 0; k < N; k++) 
	for (j = B[k].link[0]; j > 0; j--) 
	if ((m = B[k].link[j]) >= 0 && B[k].link[NUM_PHASE] != B[m].link[NUM_PHASE] && B[k].bar && B[k].bar->ce[0] && B[k].bar->ce[0]->graph) {
		for (l = 0; l < NUM_PHASE; l++) 
			if (B[k].link[l+1] == -m+SRF_STATE) break;
			if (l < B[k].bar->ce[0]->graph[1]) {
				int m1 = get_num(B[k].bar->ce[B[k].bar->ce[0]->graph[l+2]]->graph, 0),
					 m2 = get_num(B[k].bar->ce[B[k].bar->ce[0]->graph[l+2]]->graph, 1);
				double X1 = B[k].bar->ce[m1]->mp[1], Y1 = B[k].bar->ce[m1]->mp[2],
						 X2 = B[k].bar->ce[m2]->mp[1], Y2 = B[k].bar->ce[m2]->mp[2];
				if (fabs(sqr(X1-X0)+sqr(Y1-Y0)-sqr(rad)) < eps && fabs(sqr(X2-X0)+sqr(Y2-Y0)-sqr(rad)) < eps) {
					B[k].bar->ce[B[k].bar->ce[0]->graph[l+2]]->cells_to(SPHERE_GENUS);
					B[k].bar->ce[B[k].bar->ce[0]->graph[l+2]]->circ_correct(rad, B[k].link[NUM_PHASE] == -1 ? -1 : 1, NULL_STATE);
				}
			}
	}
}

//////////////////////////////////////////////
//...корректируем номера криволинейных блоков;
void CDraft::StructEllCorrect(int & hit, double X, double Y)
{
	if (hit >= 0 && B[hit].link[NUM_PHASE] == -1 && bar && bar->ce[0] && ((CGrid_el *)&stru)->in_ellipse(bar->ce[0]->mp, X, Y, -1.)) {
		if (B[hit].link[0] > NUM_PHASE && B[hit].link[NUM_PHASE+1] >= 0) hit = B[hit].link[NUM_PHASE+1];
		else { //..проверяем еще один уровень связей;
			int l, k;
			for (l = 0; l < B[hit].link[0]; l++)
			if ((k = B[hit].link[l+1]) >= 0 && B[k].link[NUM_PHASE] == -1 && 
						B[k].link[0] > NUM_PHASE && B[k].link[NUM_PHASE+1] >= 0) break;

			if (l < B[hit].link[0]) hit = B[k].link[NUM_PHASE+1];	
			else							hit = -1;
		}
	}
}
#endif
