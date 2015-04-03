/*===========================================*/
/*                   LIB_TB                  */
/*===========================================*/
#ifndef ___LIB_TB___
#define ___LIB_TB___

#include "cprofile.h"

/*=========================================================================*/
/*                  TABLES of SAMPLES from LIBRARY BARS25                  */
/*=========================================================================*/
Table * circle				  (int id_static_char = NULL_STATE);
Table * ellipse			  (int id_static_char = NULL_STATE);
Table * ring				  (int id_static_char = NULL_STATE);
Table * open_ring			  (int id_static_char = NULL_STATE);
Table * ring_sector		  (int id_static_char = NULL_STATE);
Table * sector				  (int id_static_char = NULL_STATE);
Table * circ_segment		  (int id_static_char = NULL_STATE);
Table * shaft_profile	  (int id_static_char = NULL_STATE);
Table * hexa				  (int id_static_char = NULL_STATE);
Table * hexa_tube			  (int id_static_char = NULL_STATE);
Table * equi_triangle	  (int id_static_char = NULL_STATE);
Table * rounded_triangle  (int id_static_char = NULL_STATE);
Table * rounded_wedge	  (int id_static_char = NULL_STATE);
Table * square				  (int id_static_char = NULL_STATE);
Table * rectangle			  (int id_static_char = NULL_STATE);
Table * rounded_quadrangle(int id_static_char = NULL_STATE);
Table * wedge				  (int id_static_char = NULL_STATE);
Table * right_corner		  (int id_static_char = NULL_STATE);
Table * corner				  (int id_static_char = NULL_STATE);
Table * bubble_corner	  (int id_static_char = NULL_STATE);
Table * frame				  (int id_static_char = NULL_STATE);
Table * triangle_frame	  (int id_static_char = NULL_STATE);
Table * quadrangle_frame  (int id_static_char = NULL_STATE);
Table * arbitrary_frame	  (int id_static_char = NULL_STATE);
Table * un_frame			  (int id_static_char = NULL_STATE);
Table * trap_tube			  (int id_static_char = NULL_STATE);
Table * right_zed			  (int id_static_char = NULL_STATE);
Table * zed					  (int id_static_char = NULL_STATE);
Table * bordered_zed		  (int id_static_char = NULL_STATE);
Table * T_section			  (int id_static_char = NULL_STATE);
Table * T_bubbles			  (int id_static_char = NULL_STATE);
Table * T_top_bulba		  (int id_static_char = NULL_STATE);
Table * T_bordered		  (int id_static_char = NULL_STATE);
Table * I_section			  (int id_static_char = NULL_STATE);
Table * I_common			  (int id_static_char = NULL_STATE);
Table * schweller			  (int id_static_char = NULL_STATE);
Table * schweller_ext	  (int id_static_char = NULL_STATE);
Table * schweller_int	  (int id_static_char = NULL_STATE);
Table * T_shweller_bubbles(int id_static_char = NULL_STATE);
Table * T_shweller		  (int id_static_char = NULL_STATE);
Table * right_schweller	  (int id_static_char = NULL_STATE);
Table * trap_schweller	  (int id_static_char = NULL_STATE);
Table * special_schweller (int id_static_char = NULL_STATE);
Table * zed_schweller	  (int id_static_char = NULL_STATE);
Table * h_section			  (int id_static_char = NULL_STATE);
Table * triangle_star	  (int id_static_char = NULL_STATE);
Table * bench_profile	  (int id_static_char = NULL_STATE);
Table * cross				  (int id_static_char = NULL_STATE);
Table * spangout			  (int id_static_char = NULL_STATE);
Table * gear_profile		  (int id_static_char = NULL_STATE);
Table * trouph_profile	  (int id_static_char = NULL_STATE);
Table * cuffs_profile	  (int id_static_char = NULL_STATE);
Table * hook				  (int id_static_char = NULL_STATE);

///////////////////////////////////////
//...таблицы обpазцов для 3D уpавнений;
void crank_shaft(int k, double & alpha0, double & L0, double & A10, double & A20,
                                                      double & B10, double & B20, double & R0);
void intr_conus (int k, double & r10, double & r20, double & L10, double & L20);
void L_type_room(int k, double & A10, double & A20, double & B10, double & B20,
                        double & C0,  double & R0);
void hole_box   (int k, double & A0,  double & B0,  double & C0,  double & R0,  double & L0);
void poly_box   (int k, double ** P0);

/*===============================================================================================*/
/*                 ТАБЛИЦЫ ОБРАЗЦОВ ИЗ СТАТИЧЕСКОЙ БИБЛИОТЕКИ LIB_CH и LIB_MU                    */
/*===============================================================================================*/
///////////////////////
//...таблица для кpуга;
Table * circle(int id_static_char)
{
	static double r [] = { 0., .10};
	static char * S0[] = { "Section 1" };
	Shablon  records[] = { {DLENGTH_TYPE_RECORD, "r"},
	};
	Table * table = get_shablon_table(1, 1, records, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = r;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = r[k];
	}
	return table;
}

/////////////////////////
//...таблица для эллипса;
Table * ellipse(int id_static_char)
{
	static double a [] = { 0., .20, .50}, 
					  b [] = { 0., .10, .10};
	static char * S0[] = { "Section 1", "Section 2" };
	Shablon  records[] = {	{DLENGTH_TYPE_RECORD, "a"},
									{DLENGTH_TYPE_RECORD, "b"},
	};
	Table * table = get_shablon_table(2, 2, records, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = a;
		table->table[1].parm_values = b;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = a[k];
		((double *)(table->table[1].parm_values))[k] = b[k];
	}
	return table;
}

////////////////////////
//...таблица для кольца;
Table * ring(int id_static_char)
{
	static double r1[] = { 0., .08},
					  r0[] = { 0., .10};
	static char * S0[] = { "Section 1" };
	Shablon  records[] = {	{DLENGTH_TYPE_RECORD, "r1"},
									{DLENGTH_TYPE_RECORD, "r0"},
	};
	Table * table = get_shablon_table(2, 1, records, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = r1;
		table->table[1].parm_values = r0;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = r1[k];
		((double *)(table->table[1].parm_values))[k] = r0[k];
	}
	return table;
}

//////////////////////////////////
//...таблица для открытого кольца;
Table * open_ring(int id_static_char)
{
	static double r1[] = { 0., .01},
					  r0[] = { 0., .02};
	static char * S0[] = { "Section 1" };
	Shablon  records[] = {	{DLENGTH_TYPE_RECORD, "r1"},
									{DLENGTH_TYPE_RECORD, "r0"},
	};
	Table * table = get_shablon_table(2, 1, records, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = r1;
		table->table[1].parm_values = r0;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = r1[k];
		((double *)(table->table[1].parm_values))[k] = r0[k];
	}
	return table;
}

////////////////////////////////////
//...таблица для кольцевого сектора;
Table * ring_sector(int id_static_char)
{
	static double alpha[] = { 0., 60. },
					  r1[]    = { 0., .01 },
					  r0[]    = { 0., .02 };
	static char * S0[] = { "Section 1" };
	Shablon  records[] = {	{ DANGLE_TYPE_RECORD, "alpha"},
									{DLENGTH_TYPE_RECORD, "r1"},
									{DLENGTH_TYPE_RECORD, "r0"},
	};
	Table * table = get_shablon_table(3, 1, records, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = alpha;
		table->table[1].parm_values = r1;
		table->table[2].parm_values = r0;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = alpha[k];
		((double *)(table->table[1].parm_values))[k] = r1[k];
		((double *)(table->table[2].parm_values))[k] = r0[k];
	}
	return table;
}

/////////////////////////
//...таблица для сектора;
Table * sector(int id_static_char)
{
	static double alpha[] = { 0., 60. },
					  r []    = { 0., .01 };
	static char * S0[] = { "Section 1" };
	Shablon  records[] = {	{ DANGLE_TYPE_RECORD, "alpha"},
									{DLENGTH_TYPE_RECORD, "r"},
	};
	Table * table = get_shablon_table(2, 1, records, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = alpha;
		table->table[1].parm_values = r;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = alpha[k];
		((double *)(table->table[1].parm_values))[k] = r[k];
	}
	return table;
}

////////////////////////////////////
//...таблица для кругового сегмента;
Table * circ_segment(int id_static_char)
{
	static double alpha[] = { 0., 60. },
					  r []    = { 0., .01 };
	static char * S0[] = { "Section 1" };
	Shablon  records[] = {	{ DANGLE_TYPE_RECORD, "alpha"},
									{DLENGTH_TYPE_RECORD, "r"},
	};
	Table * table = get_shablon_table(2, 1, records, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = alpha;
		table->table[1].parm_values = r;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = alpha[k];
		((double *)(table->table[1].parm_values))[k] = r[k];
	}
	return table;
}

/////////////////////////////////
//...таблица для вала с выточкой;
Table * shaft_profile(int id_static_char)
{
	static double L [] = { 0., .02},
					  H [] = { 0., .01},
					  R [] = { 0., .02},
					  r [] = { 0.,.005};
	static char * S0[] = { "Section 1" };
	Shablon  records[] = {	{DLENGTH_TYPE_RECORD, "L"},
									{DLENGTH_TYPE_RECORD, "H"},
									{DLENGTH_TYPE_RECORD, "R"},
									{DLENGTH_TYPE_RECORD, "r"},
	};
	Table * table = get_shablon_table(4, 1, records, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = L;
		table->table[1].parm_values = H;
		table->table[2].parm_values = R;
		table->table[3].parm_values = r;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = L[k];
		((double *)(table->table[1].parm_values))[k] = H[k];
		((double *)(table->table[2].parm_values))[k] = R[k];
		((double *)(table->table[3].parm_values))[k] = r[k];
	}
	return table;
}

///////////////////////////////////////
//...таблица для шестигранного стержня;
Table * hexa(int id_static_char)
{
	static double R [] = { 0., .10};
	static char * S0[] = { "Section 1" };
	Shablon  records[] = { {DLENGTH_TYPE_RECORD, "R"},
	};
	Table * table = get_shablon_table(1, 1, records, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = R;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = R[k];
	}
	return table;
}

/////////////////////////////////////
//...таблица для шестигранной трубки;
Table * hexa_tube(int id_static_char)
{
	static double r [] = { 0.,  .01},
					  S [] = { 0., .002};
	static char * S0[] = { "Section 1" };
	Shablon  records[] = {	{DLENGTH_TYPE_RECORD, "R"},
									{DLENGTH_TYPE_RECORD, "S"},
	};
	Table * table = get_shablon_table(2, 1, records, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = r;
		table->table[1].parm_values = S;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = r[k];
		((double *)(table->table[1].parm_values))[k] = S[k];
	}
	return table;
}

//////////////////////////////////////////////
//...таблица для равностороннего треугольника;
Table * equi_triangle(int id_static_char)
{
	static double A [] = { 0., .01};
	static char * S0[] = { "Section 1" };
	Shablon  records[] = { {DLENGTH_TYPE_RECORD, "A"},
	};
	Table * table = get_shablon_table(1, 1, records, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = A;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = A[k];
	}
	return table;
}

///////////////////////////////////////////
//...таблица для скругленного тpеугольника;
Table * rounded_triangle(int id_static_char)
{
	static double A1[] = { 0.,  .07},
					  A2[] = { 0.,  .03},
					  A3[] = { 0.,  .08},
					  r1[] = { 0., .007},
					  r2[] = { 0., .007},
					  r3[] = { 0., .007};
	static char * S0[] = { "Section 1" };
	Shablon  records[] = {	{DLENGTH_TYPE_RECORD, "A1"},
									{DLENGTH_TYPE_RECORD, "A2"},
									{DLENGTH_TYPE_RECORD, "A3"},
									{DLENGTH_TYPE_RECORD, "r1"},
									{DLENGTH_TYPE_RECORD, "r2"},
									{DLENGTH_TYPE_RECORD, "r3"},
	};
	Table * table = get_shablon_table(6, 1, records, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = A1;
		table->table[1].parm_values = A2;
		table->table[2].parm_values = A3;
		table->table[3].parm_values = r1;
		table->table[4].parm_values = r2;
		table->table[5].parm_values = r3;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = A1[k];
		((double *)(table->table[1].parm_values))[k] = A2[k];
		((double *)(table->table[2].parm_values))[k] = A3[k];
		((double *)(table->table[3].parm_values))[k] = r1[k];
		((double *)(table->table[4].parm_values))[k] = r2[k];
		((double *)(table->table[5].parm_values))[k] = r3[k];
	}
	return table;
}

////////////////////////////////////
//...таблица для скругленного клина;
Table * rounded_wedge(int id_static_char)
{
	static double H [] = { 0., .05, .07},
					  B [] = { 0., .13, .10},
					  r [] = { 0., .01, .02};
	static char * S0[] = { "Section 1", "Section 2" };
	Shablon  records[] = {	{DLENGTH_TYPE_RECORD, "H"},
									{DLENGTH_TYPE_RECORD, "B"},
									{DLENGTH_TYPE_RECORD, "r"},
	};
	Table * table = get_shablon_table(3, 2, records, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = H;
		table->table[1].parm_values = B;
		table->table[2].parm_values = r;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = H[k];
		((double *)(table->table[1].parm_values))[k] = B[k];
		((double *)(table->table[2].parm_values))[k] = r[k];
	}
	return table;
}

//////////////////////////
//...таблица для квадрата;
Table * square(int id_static_char)
{
	static double A [] = {0., .01};
	static char * S0[] = { "Section 1" };
	Shablon  records[] = { {DLENGTH_TYPE_RECORD, "A"},
	};
	Table * table = get_shablon_table(1, 1, records, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = A;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = A[k];
	}
	return table;
}

////////////////////////////////
//...таблица для прямоугольника;
Table * rectangle(int id_static_char)
{
	static double A [] = { 0., .20, .50}, 
					  B [] = { 0., .10, .10};
	static char * S0[] = { "Section 1" };
	Shablon  records[] = {	{DLENGTH_TYPE_RECORD, "A"},
									{DLENGTH_TYPE_RECORD, "B"},
	};
	Table * table = get_shablon_table(2, 1, records, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = A;
		table->table[1].parm_values = B;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = A[k];
		((double *)(table->table[1].parm_values))[k] = B[k];
	}
	return table;
}

///////////////////////////////////////////////
//...таблица для скругленного четырехугольника;
Table * rounded_quadrangle(int id_static_char)
{
	static double alpha[] = { 0., 120.},
					  A1[]    = { 0.,  .07},
					  A2[]    = { 0.,  .03},
					  A3[]    = { 0.,  .10},
					  A4[]    = { 0.,  .07},
					  r1[]    = { 0., .007},
					  r2[]    = { 0., .015},
					  r3[]    = { 0., .010},
					  r4[]    = { 0., .007};
	static char * S0[] = { "Section 1" };
	Shablon  records[] = {	{ DANGLE_TYPE_RECORD, "alpha"},
									{DLENGTH_TYPE_RECORD, "A1"},
									{DLENGTH_TYPE_RECORD, "A2"},
									{DLENGTH_TYPE_RECORD, "A3"},
									{DLENGTH_TYPE_RECORD, "A4"},
									{DLENGTH_TYPE_RECORD, "r1"},
									{DLENGTH_TYPE_RECORD, "r2"},
									{DLENGTH_TYPE_RECORD, "r3"},
									{DLENGTH_TYPE_RECORD, "r4"},
	};
	Table * table = get_shablon_table(9, 1, records, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = alpha;
		table->table[1].parm_values = A1;
		table->table[2].parm_values = A2;
		table->table[3].parm_values = A3;
		table->table[4].parm_values = A4;
		table->table[5].parm_values = r1;
		table->table[6].parm_values = r2;
		table->table[7].parm_values = r3;
		table->table[8].parm_values = r4;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = alpha[k];
		((double *)(table->table[1].parm_values))[k] = A1[k];
		((double *)(table->table[2].parm_values))[k] = A2[k];
		((double *)(table->table[3].parm_values))[k] = A3[k];
		((double *)(table->table[4].parm_values))[k] = A4[k];
		((double *)(table->table[5].parm_values))[k] = r1[k];
		((double *)(table->table[6].parm_values))[k] = r2[k];
		((double *)(table->table[7].parm_values))[k] = r3[k];
		((double *)(table->table[8].parm_values))[k] = r4[k];
	}
	return table;
}

///////////////////////
//...таблица для клина;
Table * wedge(int id_static_char)
{
	static double H [] = { 0., .07},
					  B [] = { 0., .12},
					  h [] = { 0., .04},
					  r [] = { 0.,.005};
	static char * S0[] = { "Section 1" };
	Shablon  records[] = {	{DLENGTH_TYPE_RECORD, "H"},
									{DLENGTH_TYPE_RECORD, "B"},
									{DLENGTH_TYPE_RECORD, "h"},
									{DLENGTH_TYPE_RECORD, "r"},
	};
	Table * table = get_shablon_table(4, 1, records, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = H;
		table->table[1].parm_values = B;
		table->table[2].parm_values = h;
		table->table[3].parm_values = r;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = H[k];
		((double *)(table->table[1].parm_values))[k] = B[k];
		((double *)(table->table[2].parm_values))[k] = h[k];
		((double *)(table->table[3].parm_values))[k] = r[k];
	}
	return table;
}

///////////////////////////////////
//...таблица для прямого угольника;
Table * right_corner(int id_static_char)
{
	static double A1[] = { 0.,  .12,  .10},
					  A2[] = { 0.,  .12,  .20},
					  S1[] = { 0.,  .01,  .02},
					  S2[] = { 0.,  .01,  .02},
					  r [] = { 0., .015, .035},
					  r1[] = { 0.,.0075, .012},
					  r2[] = { 0.,.0075, .012},
					  r3[] = { 0., .002,   0.};
	static char * S0[] = { "Section 1", "Section 2" };
	Shablon  records[] = {	{DLENGTH_TYPE_RECORD, "A1"},
									{DLENGTH_TYPE_RECORD, "A2"},
									{DLENGTH_TYPE_RECORD, "S1"},
									{DLENGTH_TYPE_RECORD, "S2"},
									{DLENGTH_TYPE_RECORD, "r"},
									{DLENGTH_TYPE_RECORD, "r1"},
									{DLENGTH_TYPE_RECORD, "r2"},
									{DLENGTH_TYPE_RECORD, "r3"},
	};
	Table * table = get_shablon_table(8, 2, records, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = A1;
		table->table[1].parm_values = A2;
		table->table[2].parm_values = S1;
		table->table[3].parm_values = S2;
		table->table[4].parm_values = r;
		table->table[5].parm_values = r1;
		table->table[6].parm_values = r2;
		table->table[7].parm_values = r3;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = A1[k];
		((double *)(table->table[1].parm_values))[k] = A2[k];
		((double *)(table->table[2].parm_values))[k] = S1[k];
		((double *)(table->table[3].parm_values))[k] = S2[k];
		((double *)(table->table[4].parm_values))[k] = r [k];
		((double *)(table->table[5].parm_values))[k] = r1[k];
		((double *)(table->table[6].parm_values))[k] = r2[k];
		((double *)(table->table[7].parm_values))[k] = r3[k];
	}
	return table;
}

/////////////////////////////////////////
//...таблица для пpоизвольного угольника;
Table * corner(int id_static_char)
{
	static double alpha[] = { 0.,120.},
					  A1[]	 = { 0., .05},
					  A2[]	 = { 0., .05},
					  S1[]	 = { 0.,.012},
					  S2[]	 = { 0.,.012},
					  r []	 = { 0.,.025},
					  r1[]	 = { 0.,.005},
					  r2[]	 = { 0.,.005},
					  r3[]	 = { 0.,.005};
	static char * S0[] = { "Section 1" };
	Shablon  records[] = {	{ DANGLE_TYPE_RECORD, "alpha"},
									{DLENGTH_TYPE_RECORD, "A1"},
									{DLENGTH_TYPE_RECORD, "A2"},
									{DLENGTH_TYPE_RECORD, "S1"},
									{DLENGTH_TYPE_RECORD, "S2"},
									{DLENGTH_TYPE_RECORD, "r"},
									{DLENGTH_TYPE_RECORD, "r1"},
									{DLENGTH_TYPE_RECORD, "r2"},
									{DLENGTH_TYPE_RECORD, "r3"},
	};
	Table * table = get_shablon_table(9, 1, records, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = alpha;
		table->table[1].parm_values = A1;
		table->table[2].parm_values = A2;
		table->table[3].parm_values = S1;
		table->table[4].parm_values = S2;
		table->table[5].parm_values = r;
		table->table[6].parm_values = r1;
		table->table[7].parm_values = r2;
		table->table[8].parm_values = r3;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = alpha[k];
		((double *)(table->table[1].parm_values))[k] = A1[k];
		((double *)(table->table[2].parm_values))[k] = A2[k];
		((double *)(table->table[3].parm_values))[k] = S1[k];
		((double *)(table->table[4].parm_values))[k] = S2[k];
		((double *)(table->table[5].parm_values))[k] = r[k];
		((double *)(table->table[6].parm_values))[k] = r1[k];
		((double *)(table->table[7].parm_values))[k] = r2[k];
		((double *)(table->table[8].parm_values))[k] = r3[k];
	}
	return table;
}

/////////////////////////////////
//...таблица для бульбоугольника;
Table * bubble_corner(int id_static_char)
{
	static double A1[] = { 0., .12},
					  A2[] = { 0., .13},
					  S1[] = { 0., .01},
					  S2[] = { 0., .01},
					  R [] = { 0.,.015},
					  r [] = { 0.,.015},
					  r1[] = { 0.,.005},
					  r2[] = { 0.,.015},
					  r3[] = { 0.,.007};
	static char * S0[] = { "Section 1" };
	Shablon  records[] = {	{DLENGTH_TYPE_RECORD, "A1"},
									{DLENGTH_TYPE_RECORD, "A2"},
									{DLENGTH_TYPE_RECORD, "S1"},
									{DLENGTH_TYPE_RECORD, "S2"},
									{DLENGTH_TYPE_RECORD, "R"},
									{DLENGTH_TYPE_RECORD, "r"},
									{DLENGTH_TYPE_RECORD, "r1"},
									{DLENGTH_TYPE_RECORD, "r2"},
									{DLENGTH_TYPE_RECORD, "r3"},
	};
	Table * table = get_shablon_table(9, 1, records, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = A1;
		table->table[1].parm_values = A2;
		table->table[2].parm_values = S1;
		table->table[3].parm_values = S2;
		table->table[4].parm_values = R;
		table->table[5].parm_values = r;
		table->table[6].parm_values = r1;
		table->table[7].parm_values = r2;
		table->table[8].parm_values = r3;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = A1[k];
		((double *)(table->table[1].parm_values))[k] = A2[k];
		((double *)(table->table[2].parm_values))[k] = S1[k];
		((double *)(table->table[3].parm_values))[k] = S2[k];
		((double *)(table->table[4].parm_values))[k] = R [k];
		((double *)(table->table[5].parm_values))[k] = r [k];
		((double *)(table->table[6].parm_values))[k] = r1[k];
		((double *)(table->table[7].parm_values))[k] = r2[k];
		((double *)(table->table[8].parm_values))[k] = r3[k];
	}
	return table;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//...таблица для пpямоугольной pамы
Table * frame(int id_static_char)
{
	static double A [] = { 0., .10, .10},
					  B [] = { 0., .07, .20},
					  t1[] = { 0., .03, .02},
					  t2[] = { 0., .05, .02},
					  t3[] = { 0., .01, .02},
					  t4[] = { 0., .02, .02},
					  r1[] = { 0.,.007,  0.},
					  r2[] = { 0.,.007,  0.},
					  r3[] = { 0.,.007,  0.},
					  r4[] = { 0.,.007,  0.},
					  R1[] = { 0.,.010,  0.},
					  R2[] = { 0.,.010,  0.},
					  R3[] = { 0.,.010,  0.},
					  R4[] = { 0.,.010,  0.};
	static char * S0[] = { "Section 1" };
	Shablon  records[] = {	{DLENGTH_TYPE_RECORD, "A"},
									{DLENGTH_TYPE_RECORD, "B"},
									{DLENGTH_TYPE_RECORD, "t1"},
									{DLENGTH_TYPE_RECORD, "t2"},
									{DLENGTH_TYPE_RECORD, "t3"},
									{DLENGTH_TYPE_RECORD, "t4"},
									{DLENGTH_TYPE_RECORD, "r1"},
									{DLENGTH_TYPE_RECORD, "r2"},
									{DLENGTH_TYPE_RECORD, "r3"},
									{DLENGTH_TYPE_RECORD, "r4"},
									{DLENGTH_TYPE_RECORD, "R1"},
									{DLENGTH_TYPE_RECORD, "R2"},
									{DLENGTH_TYPE_RECORD, "R3"},
									{DLENGTH_TYPE_RECORD, "R4"},
	};
	Table * table = get_shablon_table(14, 1, records, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = A;
		table->table[1].parm_values = B;
		table->table[2].parm_values = t1;
		table->table[3].parm_values = t2;
		table->table[4].parm_values = t3;
		table->table[5].parm_values = t4;
		table->table[6].parm_values = r1;
		table->table[7].parm_values = r2;
		table->table[8].parm_values = r3;
		table->table[9].parm_values = r4;
		table->table[10].parm_values = R1;
		table->table[11].parm_values = R2;
		table->table[12].parm_values = R3;
		table->table[13].parm_values = R4;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = A[k];
		((double *)(table->table[1].parm_values))[k] = B[k];
		((double *)(table->table[2].parm_values))[k] = t1[k];
		((double *)(table->table[3].parm_values))[k] = t2[k];
		((double *)(table->table[4].parm_values))[k] = t3[k];
		((double *)(table->table[5].parm_values))[k] = t4[k];
		((double *)(table->table[6].parm_values))[k] = r1[k];
		((double *)(table->table[7].parm_values))[k] = r2[k];
		((double *)(table->table[8].parm_values))[k] = r3[k];
		((double *)(table->table[9].parm_values))[k] = r4[k];
		((double *)(table->table[10].parm_values))[k] = R1[k];
		((double *)(table->table[11].parm_values))[k] = R2[k];
		((double *)(table->table[12].parm_values))[k] = R3[k];
		((double *)(table->table[13].parm_values))[k] = R4[k];
	}
	return table;
}

///////////////////////////////////////////
//...таблица для открытой треугольной pамы;
Table * triangle_frame(int id_static_char)
{
	static double H [] = { 0., .10},
					  B [] = { 0., .07},
					  d [] = { 0., .01},
					  h1[] = { 0., .03},
					  h2[] = { 0.,.015},
					  r [] = { 0., .01},
					  r1[] = { 0., .01},
					  r2[] = { 0., .01};
	static char * S0[] = { "Section 1" };
	Shablon  records[] = {	{DLENGTH_TYPE_RECORD, "H"},
									{DLENGTH_TYPE_RECORD, "B"},
									{DLENGTH_TYPE_RECORD, "d"},
									{DLENGTH_TYPE_RECORD, "h1"},
									{DLENGTH_TYPE_RECORD, "h2"},
									{DLENGTH_TYPE_RECORD, "r"},
									{DLENGTH_TYPE_RECORD, "r1"},
									{DLENGTH_TYPE_RECORD, "r2"},
	};
	Table * table = get_shablon_table(8, 1, records, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = H;
		table->table[1].parm_values = B;
		table->table[2].parm_values = d;
		table->table[3].parm_values = h1;
		table->table[4].parm_values = h2;
		table->table[5].parm_values = r;
		table->table[6].parm_values = r1;
		table->table[7].parm_values = r2;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = H[k];
		((double *)(table->table[1].parm_values))[k] = B[k];
		((double *)(table->table[2].parm_values))[k] = d[k];
		((double *)(table->table[3].parm_values))[k] = h1[k];
		((double *)(table->table[4].parm_values))[k] = h2[k];
		((double *)(table->table[5].parm_values))[k] = r[k];
		((double *)(table->table[6].parm_values))[k] = r1[k];
		((double *)(table->table[7].parm_values))[k] = r2[k];
	}
	return table;
}

/////////////////////////////////////////////
//...таблица для открытой пpямоугольной pамы;
Table * quadrangle_frame(int id_static_char)
{
	static double H [] = { 0., .07},
					  B [] = { 0., .10},
					  l1[] = { 0., .03},
					  l2[] = { 0., .06},
					  b1[] = { 0., .01},
					  b2[] = { 0., .01},
					  h1[] = { 0.,.008},
					  h2[] = { 0.,.008},
					  r [] = { 0.,.007},
					  r1[] = { 0.,.006},
					  r2[] = { 0.,.003},
					  r3[] = { 0.,.003};
	static char * S0[] = { "Section 1" };
	Shablon  records[] = {	{DLENGTH_TYPE_RECORD, "H"},
									{DLENGTH_TYPE_RECORD, "B"},
									{DLENGTH_TYPE_RECORD, "l1"},
									{DLENGTH_TYPE_RECORD, "l2"},
									{DLENGTH_TYPE_RECORD, "b1"},
									{DLENGTH_TYPE_RECORD, "b2"},
									{DLENGTH_TYPE_RECORD, "h1"},
									{DLENGTH_TYPE_RECORD, "h2"},
									{DLENGTH_TYPE_RECORD, "r"},
									{DLENGTH_TYPE_RECORD, "r1"},
									{DLENGTH_TYPE_RECORD, "r2"},
									{DLENGTH_TYPE_RECORD, "r3"},
	};
	Table * table = get_shablon_table(12, 1, records, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = H;
		table->table[1].parm_values = B;
		table->table[2].parm_values = l1;
		table->table[3].parm_values = l2;
		table->table[4].parm_values = b1;
		table->table[5].parm_values = b2;
		table->table[6].parm_values = h1;
		table->table[7].parm_values = h2;
		table->table[8].parm_values = r;
		table->table[9].parm_values = r1;
		table->table[10].parm_values = r2;
		table->table[11].parm_values = r3;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = H[k];
		((double *)(table->table[1].parm_values))[k] = B[k];
		((double *)(table->table[2].parm_values))[k] = l1[k];
		((double *)(table->table[3].parm_values))[k] = l2[k];
		((double *)(table->table[4].parm_values))[k] = b1[k];
		((double *)(table->table[5].parm_values))[k] = b2[k];
		((double *)(table->table[6].parm_values))[k] = h1[k];
		((double *)(table->table[7].parm_values))[k] = h2[k];
		((double *)(table->table[8].parm_values))[k] = r[k];
		((double *)(table->table[9].parm_values))[k] = r1[k];
		((double *)(table->table[10].parm_values))[k] = r2[k];
		((double *)(table->table[11].parm_values))[k] = r3[k];
	}
	return table;
}

/////////////////////////////////////////////
//...таблица для профиля "произвольная рама";
Table * arbitrary_frame(int id_static_char)
{
	static double al[] = { 0.,120.},
					  A1[] = { 0., .07},
					  A2[] = { 0., .10},
					  A3[] = { 0., .03},
					  A4[] = { 0., .07},
					  S1[] = { 0.,.015},
					  S2[] = { 0., .02},
					  S3[] = { 0., .03},
					  S4[] = { 0., .01},
					  r1[] = { 0.,.007},
					  r2[] = { 0.,.005},
					  r3[] = { 0.,.003},
					  r4[] = { 0.,.007},
					  R1[] = { 0., .01},
					  R2[] = { 0., .02},
					  R3[] = { 0., .01},
					  R4[] = { 0., .02};
	static char * S0[] = { "Section 1" };
	Shablon  records[] = {	{ DANGLE_TYPE_RECORD, "alpha"},
									{DLENGTH_TYPE_RECORD, "A1"},
									{DLENGTH_TYPE_RECORD, "A2"},
									{DLENGTH_TYPE_RECORD, "A3"},
									{DLENGTH_TYPE_RECORD, "A4"},
									{DLENGTH_TYPE_RECORD, "S1"},
									{DLENGTH_TYPE_RECORD, "S2"},
									{DLENGTH_TYPE_RECORD, "S3"},
									{DLENGTH_TYPE_RECORD, "S4"},
									{DLENGTH_TYPE_RECORD, "r1"},
									{DLENGTH_TYPE_RECORD, "r2"},
									{DLENGTH_TYPE_RECORD, "r3"},
									{DLENGTH_TYPE_RECORD, "r4"},
									{DLENGTH_TYPE_RECORD, "R1"},
									{DLENGTH_TYPE_RECORD, "R2"},
									{DLENGTH_TYPE_RECORD, "R3"},
									{DLENGTH_TYPE_RECORD, "R4"},
	};
	Table * table = get_shablon_table(17, 1, records, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = al;
		table->table[1].parm_values = A1;
		table->table[2].parm_values = A2;
		table->table[3].parm_values = A3;
		table->table[4].parm_values = A4;
		table->table[5].parm_values = S1;
		table->table[6].parm_values = S2;
		table->table[7].parm_values = S3;
		table->table[8].parm_values = S4;
		table->table[9].parm_values = r1;
		table->table[10].parm_values = r2;
		table->table[11].parm_values = r3;
		table->table[12].parm_values = r4;
		table->table[13].parm_values = R1;
		table->table[14].parm_values = R2;
		table->table[15].parm_values = R3;
		table->table[16].parm_values = R4;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = al[k];
		((double *)(table->table[1].parm_values))[k] = A1[k];
		((double *)(table->table[2].parm_values))[k] = A2[k];
		((double *)(table->table[3].parm_values))[k] = A3[k];
		((double *)(table->table[4].parm_values))[k] = A4[k];
		((double *)(table->table[5].parm_values))[k] = S1[k];
		((double *)(table->table[6].parm_values))[k] = S2[k];
		((double *)(table->table[7].parm_values))[k] = S3[k];
		((double *)(table->table[8].parm_values))[k] = S4[k];
		((double *)(table->table[9].parm_values))[k] = r1[k];
		((double *)(table->table[10].parm_values))[k] = r2[k];
		((double *)(table->table[11].parm_values))[k] = r3[k];
		((double *)(table->table[12].parm_values))[k] = r4[k];
		((double *)(table->table[13].parm_values))[k] = R1[k];
		((double *)(table->table[14].parm_values))[k] = R2[k];
		((double *)(table->table[15].parm_values))[k] = R3[k];
		((double *)(table->table[16].parm_values))[k] = R4[k];
	}
	return table;
}

/////////////////////////////////////////////////////////
//...таблица для профиля "незамкнутая произвольная рама";
Table * un_frame(int id_static_char)
{
	static double alpha1[] = { 0.,100.},
					  alpha2[] = { 0.,120.},
					  alpha3[] = { 0., 90.},
					  A1    [] = { 0., .05},
					  A2    [] = { 0., .10},
					  A3    [] = { 0., .10},
					  A4    [] = { 0., .05},
					  S1    [] = { 0., .02},
					  S2    [] = { 0., .02},
					  S3    [] = { 0., .03},
					  S4    [] = { 0., .01},
					  r1    [] = { 0.,.005},
					  r2    [] = { 0., .02},
					  r3    [] = { 0.,.015},
					  R1    [] = { 0., .02},
					  R2    [] = { 0., .02},
					  R3    [] = { 0., .02};
	static char * S0[] = { "Section 1" };
	Shablon  records[] = {	{ DANGLE_TYPE_RECORD, "alpha1"},
									{ DANGLE_TYPE_RECORD, "alpha2"},
									{ DANGLE_TYPE_RECORD, "alpha3"},
									{DLENGTH_TYPE_RECORD, "A1"},
									{DLENGTH_TYPE_RECORD, "A2"},
									{DLENGTH_TYPE_RECORD, "A3"},
									{DLENGTH_TYPE_RECORD, "A4"},
									{DLENGTH_TYPE_RECORD, "S1"},
									{DLENGTH_TYPE_RECORD, "S2"},
									{DLENGTH_TYPE_RECORD, "S3"},
									{DLENGTH_TYPE_RECORD, "S4"},
									{DLENGTH_TYPE_RECORD, "r1"},
									{DLENGTH_TYPE_RECORD, "r2"},
									{DLENGTH_TYPE_RECORD, "r3"},
									{DLENGTH_TYPE_RECORD, "R1"},
									{DLENGTH_TYPE_RECORD, "R2"},
									{DLENGTH_TYPE_RECORD, "R3"},
	};
	Table * table = get_shablon_table(17, 1, records, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = alpha1;
		table->table[1].parm_values = alpha2;
		table->table[2].parm_values = alpha3;
		table->table[3].parm_values = A1;
		table->table[4].parm_values = A2;
		table->table[5].parm_values = A3;
		table->table[6].parm_values = A4;
		table->table[7].parm_values = S1;
		table->table[8].parm_values = S2;
		table->table[9].parm_values = S3;
		table->table[10].parm_values = S4;
		table->table[11].parm_values = r1;
		table->table[12].parm_values = r2;
		table->table[13].parm_values = r3;
		table->table[14].parm_values = R1;
		table->table[15].parm_values = R2;
		table->table[16].parm_values = R3;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = alpha1[k];
		((double *)(table->table[1].parm_values))[k] = alpha2[k];
		((double *)(table->table[2].parm_values))[k] = alpha3[k];
		((double *)(table->table[3].parm_values))[k] = A1[k];
		((double *)(table->table[4].parm_values))[k] = A2[k];
		((double *)(table->table[5].parm_values))[k] = A3[k];
		((double *)(table->table[6].parm_values))[k] = A4[k];
		((double *)(table->table[7].parm_values))[k] = S1[k];
		((double *)(table->table[8].parm_values))[k] = S2[k];
		((double *)(table->table[9].parm_values))[k] = S3[k];
		((double *)(table->table[10].parm_values))[k] = S4[k];
		((double *)(table->table[11].parm_values))[k] = r1[k];
		((double *)(table->table[12].parm_values))[k] = r2[k];
		((double *)(table->table[13].parm_values))[k] = r3[k];
		((double *)(table->table[14].parm_values))[k] = R1[k];
		((double *)(table->table[15].parm_values))[k] = R2[k];
		((double *)(table->table[16].parm_values))[k] = R3[k];
	}
	return table;
}

////////////////////////////////////////
//...таблица для трапецеобразной трубки;
Table * trap_tube(int id_static_char)
{
	static double H [] = { 0.,.007},
					  A [] = { 0., .01},
					  B [] = { 0.,.005},
					  S [] = { 0.,.001};
	static char * S0[] = { "Section 1" };
	Shablon  records[] = {	{DLENGTH_TYPE_RECORD, "H"},
									{DLENGTH_TYPE_RECORD, "A"},
									{DLENGTH_TYPE_RECORD, "B"},
									{DLENGTH_TYPE_RECORD, "S"},
	};
	Table * table = get_shablon_table(4, 1, records, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = H;
		table->table[1].parm_values = A;
		table->table[2].parm_values = B;
		table->table[3].parm_values = S;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = H[k];
		((double *)(table->table[1].parm_values))[k] = A[k];
		((double *)(table->table[2].parm_values))[k] = B[k];
		((double *)(table->table[3].parm_values))[k] = S[k];
	}
	return table;
}

//////////////////////////////////////////
//...таблица для пpофиля "пpавильный зет";
Table * right_zed(int id_static_char)
{
	static double H [] = { 0., .07,  .20},
					  B1[] = { 0., .07, .095},
					  B2[] = { 0., .10, .055},
					  S1[] = { 0., .02,  .06},
					  S2[] = { 0., .03,  .04},
					  S3[] = { 0., .01,  .04},
					  R1[] = { 0., .02,   0.},
					  R2[] = { 0., .01,   0.},
					  r1[] = { 0.,.015,   0.},
					  r2[] = { 0.,.003,   0.},
					  r3[] = { 0.,.007,   0.},
					  r4[] = { 0.,.003,   0.};
	static char * S0[] = { "Section 1", "Section 2" };
	Shablon  records[] = {	{DLENGTH_TYPE_RECORD, "H"},
									{DLENGTH_TYPE_RECORD, "B1"},
									{DLENGTH_TYPE_RECORD, "B2"},
									{DLENGTH_TYPE_RECORD, "S1"},
									{DLENGTH_TYPE_RECORD, "S2"},
									{DLENGTH_TYPE_RECORD, "S3"},
									{DLENGTH_TYPE_RECORD, "R1"},
									{DLENGTH_TYPE_RECORD, "R2"},
									{DLENGTH_TYPE_RECORD, "r1"},
									{DLENGTH_TYPE_RECORD, "r2"},
									{DLENGTH_TYPE_RECORD, "r3"},
									{DLENGTH_TYPE_RECORD, "r4"},
	};
	Table * table = get_shablon_table(12, 2, records, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = H;
		table->table[1].parm_values = B1;
		table->table[2].parm_values = B2;
		table->table[3].parm_values = S1;
		table->table[4].parm_values = S2;
		table->table[5].parm_values = S3;
		table->table[6].parm_values = R1;
		table->table[7].parm_values = R2;
		table->table[8].parm_values = r1;
		table->table[9].parm_values = r2;
		table->table[10].parm_values = r3;
		table->table[11].parm_values = r4;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = H[k];
		((double *)(table->table[1].parm_values))[k] = B1[k];
		((double *)(table->table[2].parm_values))[k] = B2[k];
		((double *)(table->table[3].parm_values))[k] = S1[k];
		((double *)(table->table[4].parm_values))[k] = S2[k];
		((double *)(table->table[5].parm_values))[k] = S3[k];
		((double *)(table->table[6].parm_values))[k] = R1[k];
		((double *)(table->table[7].parm_values))[k] = R2[k];
		((double *)(table->table[8].parm_values))[k] = r1[k];
		((double *)(table->table[9].parm_values))[k] = r2[k];
		((double *)(table->table[10].parm_values))[k] = r3[k];
		((double *)(table->table[11].parm_values))[k] = r4[k];
	}
	return table;
}

////////////////////////////////////////////
//...таблица для пpофиля "пpоизвольный зет";
Table * zed(int id_static_char)
{
	static double alpha1[] = { 0., 60.},
					  alpha2[] = { 0., 90.},
					  H     [] = { 0., .12},
					  B1	  [] = { 0., .08},
					  B2	  [] = { 0., .10},
					  S1	  [] = { 0., .02},
					  S2	  [] = { 0., .03},
					  S3	  [] = { 0., .01},
					  R1	  [] = { 0.,.003},
					  R2	  [] = { 0., .01},
					  r1	  [] = { 0.,.015},
					  r2	  [] = { 0.,.003},
					  r3	  [] = { 0.,.007},
					  r4	  [] = { 0.,.003};
	static char * S0[] = { "Section 1" };
	Shablon  records[] = {	{ DANGLE_TYPE_RECORD, "alpha1"},
									{ DANGLE_TYPE_RECORD, "alpha2"},
									{DLENGTH_TYPE_RECORD, "H"},
									{DLENGTH_TYPE_RECORD, "B1"},
									{DLENGTH_TYPE_RECORD, "B2"},
									{DLENGTH_TYPE_RECORD, "S1"},
									{DLENGTH_TYPE_RECORD, "S2"},
									{DLENGTH_TYPE_RECORD, "S3"},
									{DLENGTH_TYPE_RECORD, "R1"},
									{DLENGTH_TYPE_RECORD, "R2"},
									{DLENGTH_TYPE_RECORD, "r1"},
									{DLENGTH_TYPE_RECORD, "r2"},
									{DLENGTH_TYPE_RECORD, "r3"},
									{DLENGTH_TYPE_RECORD, "r4"},
	};
	Table * table = get_shablon_table(14, 1, records, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = alpha1;
		table->table[1].parm_values = alpha2;
		table->table[2].parm_values = H;
		table->table[3].parm_values = B1;
		table->table[4].parm_values = B2;
		table->table[5].parm_values = S1;
		table->table[6].parm_values = S2;
		table->table[7].parm_values = S3;
		table->table[8].parm_values = R1;
		table->table[9].parm_values = R2;
		table->table[10].parm_values = r1;
		table->table[11].parm_values = r2;
		table->table[12].parm_values = r3;
		table->table[13].parm_values = r4;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = alpha1[k];
		((double *)(table->table[1].parm_values))[k] = alpha2[k];
		((double *)(table->table[2].parm_values))[k] = H[k];
		((double *)(table->table[3].parm_values))[k] = B1[k];
		((double *)(table->table[4].parm_values))[k] = B2[k];
		((double *)(table->table[5].parm_values))[k] = S1[k];
		((double *)(table->table[6].parm_values))[k] = S2[k];
		((double *)(table->table[7].parm_values))[k] = S3[k];
		((double *)(table->table[8].parm_values))[k] = R1[k];
		((double *)(table->table[9].parm_values))[k] = R2[k];
		((double *)(table->table[10].parm_values))[k] = r1[k];
		((double *)(table->table[11].parm_values))[k] = r2[k];
		((double *)(table->table[12].parm_values))[k] = r3[k];
		((double *)(table->table[13].parm_values))[k] = r4[k];
	}
	return table;
}

////////////////////////////////////////
//...таблица для профиля "зет фасонный";
Table * bordered_zed(int id_static_char)
{
	static double H [] = { 0., .35},
					  B [] = { 0., .30},
					  S [] = { 0., .02},
					  h [] = { 0., .08},
					  R1[] = { 0.,.015},
					  R2[] = { 0.,.002},
					  r [] = { 0.,.025};
	static char * S0[] = { "Section 1" };
	Shablon  records[] = {	{DLENGTH_TYPE_RECORD, "H"},
									{DLENGTH_TYPE_RECORD, "B"},
									{DLENGTH_TYPE_RECORD, "S"},
									{DLENGTH_TYPE_RECORD, "h"},
									{DLENGTH_TYPE_RECORD, "R1"},
									{DLENGTH_TYPE_RECORD, "R2"},
									{DLENGTH_TYPE_RECORD, "r"},
	};
	Table * table = get_shablon_table(7, 1, records, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = H;
		table->table[1].parm_values = B;
		table->table[2].parm_values = S;
		table->table[3].parm_values = h;
		table->table[4].parm_values = R1;
		table->table[5].parm_values = R2;
		table->table[6].parm_values = r;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = H[k];
		((double *)(table->table[1].parm_values))[k] = B[k];
		((double *)(table->table[2].parm_values))[k] = S[k];
		((double *)(table->table[3].parm_values))[k] = h[k];
		((double *)(table->table[4].parm_values))[k] = R1[k];
		((double *)(table->table[5].parm_values))[k] = R2[k];
		((double *)(table->table[6].parm_values))[k] = r[k];
	}
	return table;
}

/////////////////////////////////////
//...таблица для T-образного профиля;
Table * T_section(int id_static_char)
{
	static double H [] = { 0., .07,  .2,   .2},
					  B [] = { 0., .15, .15,  .15},
					  S [] = { 0., .02, .04,  .02},
					  S1[] = { 0., .03, .04,  .02},
					  r [] = { 0., .01,  .0, .035},
					  r1[] = { 0.,.007,  .0,   0.},
					  r2[] = { 0.,.005,  .0,   0.};
	static char * S0[] = { "Section 1", "Section 2", "Section 3" };
	Shablon  records[] = {	{DLENGTH_TYPE_RECORD, "H"},
									{DLENGTH_TYPE_RECORD, "B"},
									{DLENGTH_TYPE_RECORD, "S"},
									{DLENGTH_TYPE_RECORD, "S1"},
									{DLENGTH_TYPE_RECORD, "r"},
									{DLENGTH_TYPE_RECORD, "r1"},
									{DLENGTH_TYPE_RECORD, "r2"},
	};
	Table * table = get_shablon_table(7, 3, records, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = H;
		table->table[1].parm_values = B;
		table->table[2].parm_values = S;
		table->table[3].parm_values = S1;
		table->table[4].parm_values = r;
		table->table[5].parm_values = r1;
		table->table[6].parm_values = r2;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = H[k];
		((double *)(table->table[1].parm_values))[k] = B[k];
		((double *)(table->table[2].parm_values))[k] = S[k];
		((double *)(table->table[3].parm_values))[k] = S1[k];
		((double *)(table->table[4].parm_values))[k] = r[k];
		((double *)(table->table[5].parm_values))[k] = r1[k];
		((double *)(table->table[6].parm_values))[k] = r2[k];
	}
	return table;
}

////////////////////////////////////////////////
//...таблица для T-образного профиля с бульбами;
Table * T_bubbles(int id_static_char)
{
	static double H [] = { 0., .12},
					  B [] = { 0., .30},
					  S [] = { 0.,.015},
					  S1[] = { 0., .03},
					  R [] = { 0.,.015},
					  r [] = { 0., .01},
					  r1[] = { 0., .01},
					  r2[] = { 0.,.007};
	static char * S0[] = { "Section 1" };
	Shablon  records[] = {	{DLENGTH_TYPE_RECORD, "H"},
									{DLENGTH_TYPE_RECORD, "B"},
									{DLENGTH_TYPE_RECORD, "S"},
									{DLENGTH_TYPE_RECORD, "S1"},
									{DLENGTH_TYPE_RECORD, "R"},
									{DLENGTH_TYPE_RECORD, "r"},
									{DLENGTH_TYPE_RECORD, "r1"},
									{DLENGTH_TYPE_RECORD, "r2"},
	};
	Table * table = get_shablon_table(8, 1, records, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = H;
		table->table[1].parm_values = B;
		table->table[2].parm_values = S;
		table->table[3].parm_values = S1;
		table->table[4].parm_values = R;
		table->table[5].parm_values = r;
		table->table[6].parm_values = r1;
		table->table[7].parm_values = r2;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = H[k];
		((double *)(table->table[1].parm_values))[k] = B[k];
		((double *)(table->table[2].parm_values))[k] = S[k];
		((double *)(table->table[3].parm_values))[k] = S1[k];
		((double *)(table->table[4].parm_values))[k] = R[k];
		((double *)(table->table[5].parm_values))[k] = r[k];
		((double *)(table->table[6].parm_values))[k] = r1[k];
		((double *)(table->table[7].parm_values))[k] = r2[k];
	}
	return table;
}

///////////////////////////////////////////////////////
//...таблица для T-образного профиля с веpхней бульбой;
Table * T_top_bulba(int id_static_char)
{
	static double H [] = { 0., .14},
					  B [] = { 0., .25},
					  S [] = { 0., .02},
					  S1[] = { 0.,.025},
					  R [] = { 0.,.028},
					  r [] = { 0.,.012},
					  r1[] = { 0.,.007},
					  r2[] = { 0.,.010};
	static char * S0[] = { "Section 1" };
	Shablon  records[] = {	{DLENGTH_TYPE_RECORD, "H"},
									{DLENGTH_TYPE_RECORD, "B"},
									{DLENGTH_TYPE_RECORD, "S"},
									{DLENGTH_TYPE_RECORD, "S1"},
									{DLENGTH_TYPE_RECORD, "R"},
									{DLENGTH_TYPE_RECORD, "r"},
									{DLENGTH_TYPE_RECORD, "r1"},
									{DLENGTH_TYPE_RECORD, "r2"},
	};
	Table * table = get_shablon_table(8, 1, records, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = H;
		table->table[1].parm_values = B;
		table->table[2].parm_values = S;
		table->table[3].parm_values = S1;
		table->table[4].parm_values = R;
		table->table[5].parm_values = r;
		table->table[6].parm_values = r1;
		table->table[7].parm_values = r2;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = H[k];
		((double *)(table->table[1].parm_values))[k] = B[k];
		((double *)(table->table[2].parm_values))[k] = S[k];
		((double *)(table->table[3].parm_values))[k] = S1[k];
		((double *)(table->table[4].parm_values))[k] = R[k];
		((double *)(table->table[5].parm_values))[k] = r[k];
		((double *)(table->table[6].parm_values))[k] = r1[k];
		((double *)(table->table[7].parm_values))[k] = r2[k];
	}
	return table;
}

/////////////////////////////////////////////////////////////
//...таблица для T-образного профиля с концевыми "бордерами";
Table * T_bordered(int id_static_char)
{
	static double H [] = { 0., .17},
					  B [] = { 0., .20},
					  h [] = { 0., .03},
					  S [] = { 0.,.015},
					  S1[] = { 0.,.015},
					  S2[] = { 0., .01},
					  r [] = { 0., .01},
					  r1[] = { 0., .02};
	static char * S0[] = { "Section 1" };
	Shablon  records[] = {	{DLENGTH_TYPE_RECORD, "H"},
									{DLENGTH_TYPE_RECORD, "B"},
									{DLENGTH_TYPE_RECORD, "h"},
									{DLENGTH_TYPE_RECORD, "S"},
									{DLENGTH_TYPE_RECORD, "S1"},
									{DLENGTH_TYPE_RECORD, "S2"},
									{DLENGTH_TYPE_RECORD, "r"},
									{DLENGTH_TYPE_RECORD, "r1"},
	};
	Table * table = get_shablon_table(8, 1, records, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = H;
		table->table[1].parm_values = B;
		table->table[2].parm_values = h;
		table->table[3].parm_values = S;
		table->table[4].parm_values = S1;
		table->table[5].parm_values = S2;
		table->table[6].parm_values = r;
		table->table[7].parm_values = r1;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = H[k];
		((double *)(table->table[1].parm_values))[k] = B[k];
		((double *)(table->table[2].parm_values))[k] = h[k];
		((double *)(table->table[3].parm_values))[k] = S[k];
		((double *)(table->table[4].parm_values))[k] = S1[k];
		((double *)(table->table[5].parm_values))[k] = S2[k];
		((double *)(table->table[6].parm_values))[k] = r[k];
		((double *)(table->table[7].parm_values))[k] = r1[k];
	}
	return table;
}

////////////////////////////////////////////////////
//...таблица для I-обpазного пpофиля (для "pельса");
Table * I_section(int id_static_char)
{
	static double H [] = { 0., .10,   .2},
					  B1[] = { 0., .16,  .16},
					  B2[] = { 0., .12,  .16},
					  S [] = { 0., .03,  .02},
					  S1[] = { 0., .02,  .02},
					  S2[] = { 0., .02,  .02},
					  R1[] = { 0., .01,   0.},
					  R2[] = { 0., .01,   0.},
					  r1[] = { 0., .015,.035},
					  r2[] = { 0., .01, .035};
	static char * S0[] = { "Section 1", "Section 2" };
	Shablon  records[] = {	{DLENGTH_TYPE_RECORD, "H"},
									{DLENGTH_TYPE_RECORD, "B1"},
									{DLENGTH_TYPE_RECORD, "B2"},
									{DLENGTH_TYPE_RECORD, "S"},
									{DLENGTH_TYPE_RECORD, "S1"},
									{DLENGTH_TYPE_RECORD, "S2"},
									{DLENGTH_TYPE_RECORD, "R1"},
									{DLENGTH_TYPE_RECORD, "R2"},
									{DLENGTH_TYPE_RECORD, "r1"},
									{DLENGTH_TYPE_RECORD, "r2"},
	};
	Table * table = get_shablon_table(10, 2, records, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = H;
		table->table[1].parm_values = B1;
		table->table[2].parm_values = B2;
		table->table[3].parm_values = S;
		table->table[4].parm_values = S1;
		table->table[5].parm_values = S2;
		table->table[6].parm_values = R1;
		table->table[7].parm_values = R2;
		table->table[8].parm_values = r1;
		table->table[9].parm_values = r2;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = H[k];
		((double *)(table->table[1].parm_values))[k] = B1[k];
		((double *)(table->table[2].parm_values))[k] = B2[k];
		((double *)(table->table[3].parm_values))[k] = S[k];
		((double *)(table->table[4].parm_values))[k] = S1[k];
		((double *)(table->table[5].parm_values))[k] = S2[k];
		((double *)(table->table[6].parm_values))[k] = R1[k];
		((double *)(table->table[7].parm_values))[k] = R2[k];
		((double *)(table->table[8].parm_values))[k] = r1[k];
		((double *)(table->table[9].parm_values))[k] = r2[k];
	}
	return table;
}

///////////////////////////////////////////////////////////////////
//...таблица для общего I-обpазного пpофиля (для фирмы Build_Soft);
Table * I_common(int id_static_char)
{
	static double B[] = { 0., .06},
					  a[] = { 0., .06},
					  b[] = { 0., .06},
					  c[] = { 0., .05},
					  d[] = { 0., .04},
					  e[] = { 0., .06},
					  f[] = { 0., .07},
					  g[] = { 0., .06},
					  h[] = { 0., .08},
					  H[] = { 0., .20};
	static char * S0[] = { "Section 1" };
	Shablon  records[] = {	{DLENGTH_TYPE_RECORD, "B"},
									{DLENGTH_TYPE_RECORD, "a"},
									{DLENGTH_TYPE_RECORD, "b"},
									{DLENGTH_TYPE_RECORD, "c"},
									{DLENGTH_TYPE_RECORD, "d"},
									{DLENGTH_TYPE_RECORD, "e"},
									{DLENGTH_TYPE_RECORD, "f"},
									{DLENGTH_TYPE_RECORD, "g"},
									{DLENGTH_TYPE_RECORD, "h"},
									{DLENGTH_TYPE_RECORD, "H"},
	};
	Table * table = get_shablon_table(10, 1, records, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = B;
		table->table[1].parm_values = a;
		table->table[2].parm_values = b;
		table->table[3].parm_values = c;
		table->table[4].parm_values = d;
		table->table[5].parm_values = e;
		table->table[6].parm_values = f;
		table->table[7].parm_values = g;
		table->table[8].parm_values = h;
		table->table[9].parm_values = H;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = B[k];
		((double *)(table->table[1].parm_values))[k] = a[k];
		((double *)(table->table[2].parm_values))[k] = b[k];
		((double *)(table->table[3].parm_values))[k] = c[k];
		((double *)(table->table[4].parm_values))[k] = d[k];
		((double *)(table->table[5].parm_values))[k] = e[k];
		((double *)(table->table[6].parm_values))[k] = f[k];
		((double *)(table->table[7].parm_values))[k] = g[k];
		((double *)(table->table[8].parm_values))[k] = h[k];
		((double *)(table->table[9].parm_values))[k] = H[k];
	}
	return table;
}

//////////////////////////
//...таблица для швеллеpа;
Table * schweller(int id_static_char)
{
	static double H [] = { 0., .07,  .20},
					  B1[] = { 0., .10,  .10},
					  B2[] = { 0., .13,  .10},
					  S [] = { 0.,.025,  .02},
					  S1[] = { 0., .01,  .02},
					  S2[] = { 0., .01,  .02},
					  r1[] = { 0.,.015, .035},
					  r2[] = { 0.,.015, .035},
					  R1[] = { 0.,.007,   0.},
					  R2[] = { 0.,.007,   0.};
	static char * S0[] = { "Section 1", "Section 2" };
	Shablon  records[] = {	{DLENGTH_TYPE_RECORD, "H"},
									{DLENGTH_TYPE_RECORD, "B1"},
									{DLENGTH_TYPE_RECORD, "B2"},
									{DLENGTH_TYPE_RECORD, "S"},
									{DLENGTH_TYPE_RECORD, "S1"},
									{DLENGTH_TYPE_RECORD, "S2"},
									{DLENGTH_TYPE_RECORD, "r1"},
									{DLENGTH_TYPE_RECORD, "r2"},
									{DLENGTH_TYPE_RECORD, "R1"},
									{DLENGTH_TYPE_RECORD, "R2"},
	};
	Table * table = get_shablon_table(10, 2, records, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = H;
		table->table[1].parm_values = B1;
		table->table[2].parm_values = B2;
		table->table[3].parm_values = S;
		table->table[4].parm_values = S1;
		table->table[5].parm_values = S2;
		table->table[6].parm_values = r1;
		table->table[7].parm_values = r2;
		table->table[8].parm_values = R1;
		table->table[9].parm_values = R2;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = H[k];
		((double *)(table->table[1].parm_values))[k] = B1[k];
		((double *)(table->table[2].parm_values))[k] = B2[k];
		((double *)(table->table[3].parm_values))[k] = S[k];
		((double *)(table->table[4].parm_values))[k] = S1[k];
		((double *)(table->table[5].parm_values))[k] = S2[k];
		((double *)(table->table[6].parm_values))[k] = r1[k];
		((double *)(table->table[7].parm_values))[k] = r2[k];
		((double *)(table->table[8].parm_values))[k] = R1[k];
		((double *)(table->table[9].parm_values))[k] = R2[k];
	}
	return table;
}

//////////////////////////////////////////
//...таблица для открытого бульбошвеллера;
Table * schweller_ext(int id_static_char)
{
	static double H [] = { 0., .09},
					  B [] = { 0., .10},
					  S [] = { 0.,.015},
					  d [] = { 0., .03},
					  r [] = { 0., .02},
					  R [] = { 0., .02};
	static char * S0[] = { "Section 1" };
	Shablon  records[] = {	{DLENGTH_TYPE_RECORD, "H"},
									{DLENGTH_TYPE_RECORD, "B"},
									{DLENGTH_TYPE_RECORD, "S"},
									{DLENGTH_TYPE_RECORD, "d"},
									{DLENGTH_TYPE_RECORD, "r"},
									{DLENGTH_TYPE_RECORD, "R"},
	};
	Table * table = get_shablon_table(6, 1, records, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = H;
		table->table[1].parm_values = B;
		table->table[2].parm_values = S;
		table->table[3].parm_values = d;
		table->table[4].parm_values = r;
		table->table[5].parm_values = R;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = H[k];
		((double *)(table->table[1].parm_values))[k] = B[k];
		((double *)(table->table[2].parm_values))[k] = S[k];
		((double *)(table->table[3].parm_values))[k] = d[k];
		((double *)(table->table[4].parm_values))[k] = r[k];
		((double *)(table->table[5].parm_values))[k] = R[k];
	}
	return table;
}

//////////////////////////////////////////
//...таблица для закрытого бульбошвеллера;
Table * schweller_int(int id_static_char)
{
	static double H [] = { 0., .09},
					  B [] = { 0., .10},
					  S [] = { 0.,.015},
					  d [] = { 0., .03},
					  r [] = { 0., .02},
					  R [] = { 0., .02};
	static char * S0[] = { "Section 1" };
	Shablon  records[] = {	{DLENGTH_TYPE_RECORD, "H"},
									{DLENGTH_TYPE_RECORD, "B"},
									{DLENGTH_TYPE_RECORD, "S"},
									{DLENGTH_TYPE_RECORD, "d"},
									{DLENGTH_TYPE_RECORD, "r"},
									{DLENGTH_TYPE_RECORD, "R"},
	};
	Table * table = get_shablon_table(6, 1, records, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = H;
		table->table[1].parm_values = B;
		table->table[2].parm_values = S;
		table->table[3].parm_values = d;
		table->table[4].parm_values = r;
		table->table[5].parm_values = R;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = H[k];
		((double *)(table->table[1].parm_values))[k] = B[k];
		((double *)(table->table[2].parm_values))[k] = S[k];
		((double *)(table->table[3].parm_values))[k] = d[k];
		((double *)(table->table[4].parm_values))[k] = r[k];
		((double *)(table->table[5].parm_values))[k] = R[k];
	}
	return table;
}

/////////////////////////////////////////////////
//...таблица для T-обpазного швеллеpа с бульбами;
Table * T_shweller_bubbles(int id_static_char)
{
	static double H [] = { 0., .17},
					  B [] = { 0., .34},
					  L [] = { 0., .17},
					  b [] = { 0., .03},
					  S [] = { 0., .02},
					  S1[] = { 0.,.015},
					  R [] = { 0.,.007},
					  R1[] = { 0., .01},
					  r [] = { 0., .02};
	static char * S0[] = { "Section 1" };
	Shablon  records[] = {	{DLENGTH_TYPE_RECORD, "H"},
									{DLENGTH_TYPE_RECORD, "B"},
									{DLENGTH_TYPE_RECORD, "L"},
									{DLENGTH_TYPE_RECORD, "b"},
									{DLENGTH_TYPE_RECORD, "S"},
									{DLENGTH_TYPE_RECORD, "S1"},
									{DLENGTH_TYPE_RECORD, "R"},
									{DLENGTH_TYPE_RECORD, "R1"},
									{DLENGTH_TYPE_RECORD, "r"},
	};
	Table * table = get_shablon_table(9, 1, records, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = H;
		table->table[1].parm_values = B;
		table->table[2].parm_values = L;
		table->table[3].parm_values = b;
		table->table[4].parm_values = S;
		table->table[5].parm_values = S1;
		table->table[6].parm_values = R;
		table->table[7].parm_values = R1;
		table->table[8].parm_values = r;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = H[k];
		((double *)(table->table[1].parm_values))[k] = B[k];
		((double *)(table->table[2].parm_values))[k] = L[k];
		((double *)(table->table[3].parm_values))[k] = b[k];
		((double *)(table->table[4].parm_values))[k] = S[k];
		((double *)(table->table[5].parm_values))[k] = S1[k];
		((double *)(table->table[6].parm_values))[k] = R[k];
		((double *)(table->table[7].parm_values))[k] = R1[k];
		((double *)(table->table[8].parm_values))[k] = r[k];
	}
	return table;
}

///////////////////////////////////////////////
//...таблица для T-обpазного швеллеpа обычного;
Table * T_shweller(int id_static_char)
{
	static double H [] = { 0., .17},
					  B [] = { 0., .34},
					  a [] = { 0., .19},
					  S [] = { 0., .02},
					  S1[] = { 0.,.015},
					  R [] = { 0., .01},
					  r [] = { 0., .02};
	static char * S0[] = { "Section 1" };
	Shablon  records[] = {	{DLENGTH_TYPE_RECORD, "H"},
									{DLENGTH_TYPE_RECORD, "B"},
									{DLENGTH_TYPE_RECORD, "a"},
									{DLENGTH_TYPE_RECORD, "S"},
									{DLENGTH_TYPE_RECORD, "S1"},
									{DLENGTH_TYPE_RECORD, "R"},
									{DLENGTH_TYPE_RECORD, "r"},
	};
	Table * table = get_shablon_table(7, 1, records, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = H;
		table->table[1].parm_values = B;
		table->table[2].parm_values = a;
		table->table[3].parm_values = S;
		table->table[4].parm_values = S1;
		table->table[5].parm_values = R;
		table->table[6].parm_values = r;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = H[k];
		((double *)(table->table[1].parm_values))[k] = B[k];
		((double *)(table->table[2].parm_values))[k] = a[k];
		((double *)(table->table[3].parm_values))[k] = S[k];
		((double *)(table->table[4].parm_values))[k] = S1[k];
		((double *)(table->table[5].parm_values))[k] = R[k];
		((double *)(table->table[6].parm_values))[k] = r[k];
	}
	return table;
}

/////////////////////////////////////////
//...таблица для швеллера отбортовочного;
Table * right_schweller(int id_static_char)
{
	static double H [] = { 0., .07},
					  B [] = { 0., .10},
					  b [] = { 0., .05},
					  S [] = { 0., .02},
					  S1[] = { 0., .01},
					  S2[] = { 0., .01},
					  r1[] = { 0., .01},
					  r2[] = { 0.,.015},
					  R [] = { 0.,.007},
					  R1[] = { 0.,.007},
					  R2[] = { 0., .00};
	static char * S0[] = { "Section 1" };
	Shablon  records[] = {	{DLENGTH_TYPE_RECORD, "H"},
									{DLENGTH_TYPE_RECORD, "B"},
									{DLENGTH_TYPE_RECORD, "b"},
									{DLENGTH_TYPE_RECORD, "S"},
									{DLENGTH_TYPE_RECORD, "S1"},
									{DLENGTH_TYPE_RECORD, "S2"},
									{DLENGTH_TYPE_RECORD, "r1"},
									{DLENGTH_TYPE_RECORD, "r2"},
									{DLENGTH_TYPE_RECORD, "R"},
									{DLENGTH_TYPE_RECORD, "R1"},
									{DLENGTH_TYPE_RECORD, "R2"},
	};
	Table * table = get_shablon_table(11, 1, records, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = H;
		table->table[1].parm_values = B;
		table->table[2].parm_values = b;
		table->table[3].parm_values = S;
		table->table[4].parm_values = S1;
		table->table[5].parm_values = S2;
		table->table[6].parm_values = r1;
		table->table[7].parm_values = r2;
		table->table[8].parm_values = R;
		table->table[9].parm_values = R1;
		table->table[10].parm_values = R2;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = H[k];
		((double *)(table->table[1].parm_values))[k] = B[k];
		((double *)(table->table[2].parm_values))[k] = b[k];
		((double *)(table->table[3].parm_values))[k] = S[k];
		((double *)(table->table[4].parm_values))[k] = S1[k];
		((double *)(table->table[5].parm_values))[k] = S2[k];
		((double *)(table->table[6].parm_values))[k] = r1[k];
		((double *)(table->table[7].parm_values))[k] = r2[k];
		((double *)(table->table[8].parm_values))[k] = R[k];
		((double *)(table->table[9].parm_values))[k] = R1[k];
		((double *)(table->table[10].parm_values))[k] = R2[k];
	}
	return table;
}

////////////////////////////////////////////////////////
//...таблица для швеллера трапецевидного отбортовочного;
Table * trap_schweller(int id_static_char)
{
	static double H [] = { 0., .17},
					  B [] = { 0., .66},
					  A [] = { 0., .20},
					  C [] = { 0., .18},
					  S [] = { 0., .03},
					  S1[] = { 0.,.025},
					  S2[] = { 0., .02},
					  r1[] = { 0., .02},
					  r2[] = { 0., .05},
					  R [] = { 0.,.007},
					  R1[] = { 0.,.007},
					  R2[] = { 0., .08};
	static char * S0[] = { "Section 1" };
	Shablon  records[] = {	{DLENGTH_TYPE_RECORD, "H"},
									{DLENGTH_TYPE_RECORD, "B"},
									{DLENGTH_TYPE_RECORD, "A"},
									{DLENGTH_TYPE_RECORD, "C"},
									{DLENGTH_TYPE_RECORD, "S"},
									{DLENGTH_TYPE_RECORD, "S1"},
									{DLENGTH_TYPE_RECORD, "S2"},
									{DLENGTH_TYPE_RECORD, "r1"},
									{DLENGTH_TYPE_RECORD, "r2"},
									{DLENGTH_TYPE_RECORD, "R"},
									{DLENGTH_TYPE_RECORD, "R1"},
									{DLENGTH_TYPE_RECORD, "R2"},
	};
	Table * table = get_shablon_table(12, 1, records, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = H;
		table->table[1].parm_values = B;
		table->table[2].parm_values = A;
		table->table[3].parm_values = C;
		table->table[4].parm_values = S;
		table->table[5].parm_values = S1;
		table->table[6].parm_values = S;
		table->table[7].parm_values = r1;
		table->table[8].parm_values = r2;
		table->table[9].parm_values = R;
		table->table[10].parm_values = R1;
		table->table[11].parm_values = R2;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = H[k];
		((double *)(table->table[1].parm_values))[k] = B[k];
		((double *)(table->table[2].parm_values))[k] = A[k];
		((double *)(table->table[3].parm_values))[k] = C[k];
		((double *)(table->table[4].parm_values))[k] = S[k];
		((double *)(table->table[5].parm_values))[k] = S1[k];
		((double *)(table->table[6].parm_values))[k] = S2[k];
		((double *)(table->table[7].parm_values))[k] = r1[k];
		((double *)(table->table[8].parm_values))[k] = r2[k];
		((double *)(table->table[9].parm_values))[k] = R[k];
		((double *)(table->table[10].parm_values))[k] = R1[k];
		((double *)(table->table[11].parm_values))[k] = R2[k];
	}
	return table;
}

///////////////////////////////////////
//...таблица для специального швеллеpа;
Table * special_schweller(int id_static_char)
{
	static double H [] = { 0.,  .25},
					  B [] = { 0.,  .20},
					  h [] = { 0.,  .10},
					  S [] = { 0., .035},
					  S1[] = { 0.,.0175},
					  r [] = { 0.,  .05},
					  r1[] = { 0., .015},
					  R1[] = { 0., .025},
					  R2[] = { 0., .025};
	static char * S0[] = { "Section 1" };
	Shablon  records[] = {	{DLENGTH_TYPE_RECORD, "H"},
									{DLENGTH_TYPE_RECORD, "B"},
									{DLENGTH_TYPE_RECORD, "h"},
									{DLENGTH_TYPE_RECORD, "S"},
									{DLENGTH_TYPE_RECORD, "S1"},
									{DLENGTH_TYPE_RECORD, "r"},
									{DLENGTH_TYPE_RECORD, "r1"},
									{DLENGTH_TYPE_RECORD, "R1"},
									{DLENGTH_TYPE_RECORD, "R2"},
	};
	Table * table = get_shablon_table(9, 1, records, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = H;
		table->table[1].parm_values = B;
		table->table[2].parm_values = h;
		table->table[3].parm_values = S;
		table->table[4].parm_values = S1;
		table->table[5].parm_values = r;
		table->table[6].parm_values = r1;
		table->table[7].parm_values = R1;
		table->table[8].parm_values = R2;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = H[k];
		((double *)(table->table[1].parm_values))[k] = B[k];
		((double *)(table->table[2].parm_values))[k] = h[k];
		((double *)(table->table[3].parm_values))[k] = S[k];
		((double *)(table->table[4].parm_values))[k] = S1[k];
		((double *)(table->table[5].parm_values))[k] = r[k];
		((double *)(table->table[6].parm_values))[k] = r1[k];
		((double *)(table->table[7].parm_values))[k] = R1[k];
		((double *)(table->table[8].parm_values))[k] = R2[k];
	}
	return table;
}

///////////////////////////////////////
//...таблица для профиля "зет-швеллер";
Table * zed_schweller(int id_static_char)
{
	static double alpha1[] = { 0, 100.},
					  alpha2[] = { 0, 100.},
					  alpha3[] = { 0, 120.},
					  A1    [] = { 0,  .07},
					  A2    [] = { 0,  .12},
					  A3    [] = { 0,  .10},
					  A4    [] = { 0,  .05},
					  S1    [] = { 0,  .02},
					  S2    [] = { 0,  .01},
					  S3    [] = { 0,  .05},
					  S4    [] = { 0,  .03},
					  r1    [] = { 0, .015},
					  r2    [] = { 0, .003},
					  r3    [] = { 0, .015},
					  R1    [] = { 0,  .02},
					  R2    [] = { 0,  .02},
					  R3    [] = { 0,  .02};
	static char * S0[] = { "Section 1" };
	Shablon  records[] = {	{ DANGLE_TYPE_RECORD, "alpha1"},
									{ DANGLE_TYPE_RECORD, "alpha2"},
									{ DANGLE_TYPE_RECORD, "alpha3"},
									{DLENGTH_TYPE_RECORD, "A1"},
									{DLENGTH_TYPE_RECORD, "A2"},
									{DLENGTH_TYPE_RECORD, "A3"},
									{DLENGTH_TYPE_RECORD, "A4"},
									{DLENGTH_TYPE_RECORD, "S1"},
									{DLENGTH_TYPE_RECORD, "S2"},
									{DLENGTH_TYPE_RECORD, "S3"},
									{DLENGTH_TYPE_RECORD, "S4"},
									{DLENGTH_TYPE_RECORD, "r1"},
									{DLENGTH_TYPE_RECORD, "r2"},
									{DLENGTH_TYPE_RECORD, "r3"},
									{DLENGTH_TYPE_RECORD, "R1"},
									{DLENGTH_TYPE_RECORD, "R2"},
									{DLENGTH_TYPE_RECORD, "R3"},
	};
	Table * table = get_shablon_table(17, 1, records, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = alpha1;
		table->table[1].parm_values = alpha2;
		table->table[2].parm_values = alpha3;
		table->table[3].parm_values = A1;
		table->table[4].parm_values = A2;
		table->table[5].parm_values = A3;
		table->table[6].parm_values = A4;
		table->table[7].parm_values = S1;
		table->table[8].parm_values = S2;
		table->table[9].parm_values = S3;
		table->table[10].parm_values = S4;
		table->table[11].parm_values = r1;
		table->table[12].parm_values = r2;
		table->table[13].parm_values = r3;
		table->table[14].parm_values = R1;
		table->table[15].parm_values = R2;
		table->table[16].parm_values = R3;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = alpha1[k];
		((double *)(table->table[1].parm_values))[k] = alpha2[k];
		((double *)(table->table[2].parm_values))[k] = alpha3[k];
		((double *)(table->table[3].parm_values))[k] = A1[k];
		((double *)(table->table[4].parm_values))[k] = A2[k];
		((double *)(table->table[5].parm_values))[k] = A3[k];
		((double *)(table->table[6].parm_values))[k] = A4[k];
		((double *)(table->table[7].parm_values))[k] = S1[k];
		((double *)(table->table[8].parm_values))[k] = S2[k];
		((double *)(table->table[9].parm_values))[k] = S3[k];
		((double *)(table->table[10].parm_values))[k] = S4[k];
		((double *)(table->table[11].parm_values))[k] = r1[k];
		((double *)(table->table[12].parm_values))[k] = r2[k];
		((double *)(table->table[13].parm_values))[k] = r3[k];
		((double *)(table->table[14].parm_values))[k] = R1[k];
		((double *)(table->table[15].parm_values))[k] = R2[k];
		((double *)(table->table[16].parm_values))[k] = R3[k];
	}
	return table;
}

/////////////////////////////////////
//...таблица для h-обpазного пpофиля;
Table * h_section(int id_static_char)
{
	static double H [] = { 0., .27},
					  B [] = { 0., .19},
					  h [] = { 0., .10},
					  S [] = { 0., .02},
					  S1[] = { 0., .03},
					  S2[] = { 0., .02},
					  R [] = { 0.,.015},
					  r1[] = { 0., .02},
					  r2[] = { 0., .01},
					  r3[] = { 0., .01};
	static char * S0[] = { "Section 1" };
	Shablon  records[] = {	{DLENGTH_TYPE_RECORD, "H"},
									{DLENGTH_TYPE_RECORD, "B"},
									{DLENGTH_TYPE_RECORD, "h"},
									{DLENGTH_TYPE_RECORD, "S"},
									{DLENGTH_TYPE_RECORD, "S1"},
									{DLENGTH_TYPE_RECORD, "S2"},
									{DLENGTH_TYPE_RECORD, "R"},
									{DLENGTH_TYPE_RECORD, "r1"},
									{DLENGTH_TYPE_RECORD, "r2"},
									{DLENGTH_TYPE_RECORD, "r3"},
	};
	Table * table = get_shablon_table(10, 1, records, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = H;
		table->table[1].parm_values = B;
		table->table[2].parm_values = h;
		table->table[3].parm_values = S;
		table->table[4].parm_values = S1;
		table->table[5].parm_values = S2;
		table->table[6].parm_values = R;
		table->table[7].parm_values = r1;
		table->table[8].parm_values = r2;
		table->table[9].parm_values = r3;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = H[k];
		((double *)(table->table[1].parm_values))[k] = B[k];
		((double *)(table->table[2].parm_values))[k] = h[k];
		((double *)(table->table[3].parm_values))[k] = S[k];
		((double *)(table->table[4].parm_values))[k] = S1[k];
		((double *)(table->table[5].parm_values))[k] = S2[k];
		((double *)(table->table[6].parm_values))[k] = R[k];
		((double *)(table->table[7].parm_values))[k] = r1[k];
		((double *)(table->table[8].parm_values))[k] = r2[k];
		((double *)(table->table[9].parm_values))[k] = r3[k];
	}
	return table;
}

////////////////////////////////////
//...таблица для тркеугольной звезды;
Table * triangle_star(int id_static_char)
{
	static double alpha[] = { 0., 60.},
					  H    [] = { 0., .20},
					  L    [] = { 0., .10},
					  S1   [] = { 0., .03},
					  S2   [] = { 0., .02},
					  r1   [] = { 0.,.015},
					  r2   [] = { 0., .05};
	static char * S0[] = { "Section 1" };
	Shablon  records[] = {	{ DANGLE_TYPE_RECORD, "alpha"},
									{DLENGTH_TYPE_RECORD, "H"},
									{DLENGTH_TYPE_RECORD, "L"},
									{DLENGTH_TYPE_RECORD, "S1"},
									{DLENGTH_TYPE_RECORD, "S2"},
									{DLENGTH_TYPE_RECORD, "r1"},
									{DLENGTH_TYPE_RECORD, "r2"},
	};
	Table * table = get_shablon_table(7, 1, records, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = alpha;
		table->table[1].parm_values = H;
		table->table[2].parm_values = L;
		table->table[3].parm_values = S1;
		table->table[4].parm_values = S2;
		table->table[5].parm_values = r1;
		table->table[6].parm_values = r2;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = alpha[k];
		((double *)(table->table[1].parm_values))[k] = H[k];
		((double *)(table->table[2].parm_values))[k] = L[k];
		((double *)(table->table[3].parm_values))[k] = S1[k];
		((double *)(table->table[4].parm_values))[k] = S2[k];
		((double *)(table->table[5].parm_values))[k] = r1[k];
		((double *)(table->table[6].parm_values))[k] = r2[k];
	}
	return table;
}

////////////////////////////////
//...таблица для двойной звезды;
Table * bench_profile(int id_static_char)
{
	static double alpha[] = { 0., 60.},
					  L    [] = { 0., .25},
					  l    [] = { 0., .07},
					  S    [] = { 0., .03},
					  S1   [] = { 0., .02},
					  r    [] = { 0., .05},
					  r1   [] = { 0.,.015};
	static char * S0[] = { "Section 1" };
	Shablon  records[] = {	{ DANGLE_TYPE_RECORD, "alpha"},
									{DLENGTH_TYPE_RECORD, "L"},
									{DLENGTH_TYPE_RECORD, "l"},
									{DLENGTH_TYPE_RECORD, "S"},
									{DLENGTH_TYPE_RECORD, "S1"},
									{DLENGTH_TYPE_RECORD, "r"},
									{DLENGTH_TYPE_RECORD, "r1"},
	};
	Table * table = get_shablon_table(7, 1, records, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = alpha;
		table->table[1].parm_values = L;
		table->table[2].parm_values = l;
		table->table[3].parm_values = S;
		table->table[4].parm_values = S1;
		table->table[5].parm_values = r;
		table->table[6].parm_values = r1;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = alpha[k];
		((double *)(table->table[1].parm_values))[k] = L[k];
		((double *)(table->table[2].parm_values))[k] = l[k];
		((double *)(table->table[3].parm_values))[k] = S[k];
		((double *)(table->table[4].parm_values))[k] = S1[k];
		((double *)(table->table[5].parm_values))[k] = r[k];
		((double *)(table->table[6].parm_values))[k] = r1[k];
	}
	return table;
}

////////////////////////
//...таблица для кpеста;
Table * cross(int id_static_char)
{
	static double H [] = { 0., .20},
					  L [] = { 0., .15},
					  S [] = { 0., .04},
					  S1[] = { 0., .03},
					  S2[] = { 0.,.015},
					  r1[] = { 0., .03},
					  r2[] = { 0.,.015};
	static char * S0[] = { "Section 1" };
	Shablon  records[] = {	{DLENGTH_TYPE_RECORD, "H"},
									{DLENGTH_TYPE_RECORD, "L"},
									{DLENGTH_TYPE_RECORD, "S"},
									{DLENGTH_TYPE_RECORD, "S1"},
									{DLENGTH_TYPE_RECORD, "S2"},
									{DLENGTH_TYPE_RECORD, "r1"},
									{DLENGTH_TYPE_RECORD, "r2"},
	};
	Table * table = get_shablon_table(7, 1, records, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = H;
		table->table[1].parm_values = L;
		table->table[2].parm_values = S;
		table->table[3].parm_values = S1;
		table->table[4].parm_values = S2;
		table->table[5].parm_values = r1;
		table->table[6].parm_values = r2;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = H[k];
		((double *)(table->table[1].parm_values))[k] = L[k];
		((double *)(table->table[2].parm_values))[k] = S[k];
		((double *)(table->table[3].parm_values))[k] = S1[k];
		((double *)(table->table[4].parm_values))[k] = S2[k];
		((double *)(table->table[5].parm_values))[k] = r1[k];
		((double *)(table->table[6].parm_values))[k] = r2[k];
	}
	return table;
}

//////////////////////////////////
//...таблица для общего шпангоута;
Table * spangout(int id_static_char)
{
	static double alpha[] = { 0., 30., 60.},
					  L    [] = { 0., .60, .80},
					  A    [] = { 0., .18, .20},
					  B    [] = { 0., .50, .50},
					  C    [] = { 0., .30, .30},
					  l1   [] = { 0., .23, .20},
					  l2   [] = { 0., .18, .30},
					  h1   [] = { 0., .08, .08},
					  h2   [] = { 0., .09, .10},
					  h3   [] = { 0.,-.13,  0.},
					  S1   [] = { 0., .05, .05},
					  S2   [] = { 0., .05, .05},
					  S3   [] = { 0., .04, .05},
					  S4   [] = { 0., .05, .05},
					  r1   [] = { 0., .10, .08},
					  r2   [] = { 0., .04, .03},
					  r3   [] = { 0., .03, .03};
		static char * S0[] = { "Initial geometry for complicated-1", "Section2" };
		Shablon  records[] = {	{ DANGLE_TYPE_RECORD, "	alpha"},
										{DLENGTH_TYPE_RECORD, "L"},
										{DLENGTH_TYPE_RECORD, "A"},
										{DLENGTH_TYPE_RECORD, "B"},
										{DLENGTH_TYPE_RECORD, "C"},
										{DLENGTH_TYPE_RECORD, "l1"},
										{DLENGTH_TYPE_RECORD, "l2"},
										{DLENGTH_TYPE_RECORD, "h1"},
										{DLENGTH_TYPE_RECORD, "h2"},
										{DLENGTH_TYPE_RECORD, "h3"},
										{DLENGTH_TYPE_RECORD, "S1"},
										{DLENGTH_TYPE_RECORD, "S2"},
										{DLENGTH_TYPE_RECORD, "S3"},
										{DLENGTH_TYPE_RECORD, "S4"},
										{DLENGTH_TYPE_RECORD, "r1"},
										{DLENGTH_TYPE_RECORD, "r2"},
										{DLENGTH_TYPE_RECORD, "r3"},
		};
	Table * table = get_shablon_table(17, 2, records, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = alpha;
		table->table[1].parm_values = L;
		table->table[2].parm_values = A;
		table->table[3].parm_values = B;
		table->table[4].parm_values = C;
		table->table[5].parm_values = l1;
		table->table[6].parm_values = l2;
		table->table[7].parm_values = h1;
		table->table[8].parm_values = h2;
		table->table[9].parm_values = h3;
		table->table[10].parm_values = S1;
		table->table[11].parm_values = S2;
		table->table[12].parm_values = S3;
		table->table[13].parm_values = S4;
		table->table[14].parm_values = r1;
		table->table[15].parm_values = r2;
		table->table[16].parm_values = r3;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = alpha[k];
		((double *)(table->table[1].parm_values))[k] = L[k];
		((double *)(table->table[2].parm_values))[k] = A[k];
		((double *)(table->table[3].parm_values))[k] = B[k];
		((double *)(table->table[4].parm_values))[k] = C[k];
		((double *)(table->table[5].parm_values))[k] = l1[k];
		((double *)(table->table[6].parm_values))[k] = l2[k];
		((double *)(table->table[7].parm_values))[k] = h1[k];
		((double *)(table->table[8].parm_values))[k] = h2[k];
		((double *)(table->table[9].parm_values))[k] = h3[k];
		((double *)(table->table[10].parm_values))[k] = S1[k];
		((double *)(table->table[11].parm_values))[k] = S2[k];
		((double *)(table->table[12].parm_values))[k] = S3[k];
		((double *)(table->table[13].parm_values))[k] = S4[k];
		((double *)(table->table[14].parm_values))[k] = r1[k];
		((double *)(table->table[15].parm_values))[k] = r2[k];
		((double *)(table->table[16].parm_values))[k] = r3[k];
	}
	return table;
}

////////////////////////////////
//...таблица для зубчатого вала;
Table * gear_profile(int id_static_char)
{
	static int    N [] = { 0,   4,   2};
	static double L [] = { 0.,.02,.015},
					  R1[] = { 0.,.03, .01},
					  R2[] = { 0.,.01, .03},
					  r [] = { 0.,.01,.002};
	static char * S0[] = { "Section 1", "Section 2" };
	Shablon  records[] = {	{	  INT_TYPE_RECORD, "N"},
									{DLENGTH_TYPE_RECORD, "L"},
									{DLENGTH_TYPE_RECORD, "R1"},
									{DLENGTH_TYPE_RECORD, "R2"},
									{DLENGTH_TYPE_RECORD, "r"},
	};
	Table * table = get_shablon_table(5, 2, records, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = N;
		table->table[1].parm_values = L;
		table->table[2].parm_values = R1;
		table->table[3].parm_values = R2;
		table->table[4].parm_values = r;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((int    *)(table->table[0].parm_values))[k] = N[k];
		((double *)(table->table[1].parm_values))[k] = L[k];
		((double *)(table->table[2].parm_values))[k] = R1[k];
		((double *)(table->table[3].parm_values))[k] = R2[k];
		((double *)(table->table[4].parm_values))[k] = r[k];
	}
	return table;
}

/////////////////////////////////////////////////
//...таблица для профиля корытообразного профиля;
Table * trouph_profile(int id_static_char)
{
	static double alpha1[] = { 0.,120.},
					  alpha2[] = { 0.,120.},
					  alpha3[] = { 0., 90.},
					  alpha4[] = { 0., 60.},
					  A1    [] = { 0., .05},
					  A2    [] = { 0., .15},
					  A3    [] = { 0., .10},
					  A4    [] = { 0., .15},
					  A5    [] = { 0., .07},
					  S1    [] = { 0., .02},
					  S2    [] = { 0., .02},
					  S3    [] = { 0., .03},
					  S4    [] = { 0., .01},
					  S5    [] = { 0., .02},
					  r1    [] = { 0.,.015},
					  r2    [] = { 0.,.003},
					  r3    [] = { 0.,.003},
					  r4    [] = { 0.,.003},
					  R1    [] = { 0., .02},
					  R2    [] = { 0., .02},
					  R3    [] = { 0., .02},
					  R4    [] = { 0., .02};
	static char * S0[] = { "Section 1" };
	Shablon  records[] = {	{ DANGLE_TYPE_RECORD, "alpha1"},
									{ DANGLE_TYPE_RECORD, "alpha2"},
									{ DANGLE_TYPE_RECORD, "alpha3"},
									{ DANGLE_TYPE_RECORD, "alpha4"},
									{DLENGTH_TYPE_RECORD, "A1"},
									{DLENGTH_TYPE_RECORD, "A2"},
									{DLENGTH_TYPE_RECORD, "A3"},
									{DLENGTH_TYPE_RECORD, "A4"},
									{DLENGTH_TYPE_RECORD, "A5"},
									{DLENGTH_TYPE_RECORD, "S1"},
									{DLENGTH_TYPE_RECORD, "S2"},
									{DLENGTH_TYPE_RECORD, "S3"},
									{DLENGTH_TYPE_RECORD, "S4"},
									{DLENGTH_TYPE_RECORD, "S5"},
									{DLENGTH_TYPE_RECORD, "r1"},
									{DLENGTH_TYPE_RECORD, "r2"},
									{DLENGTH_TYPE_RECORD, "r3"},
									{DLENGTH_TYPE_RECORD, "r4"},
									{DLENGTH_TYPE_RECORD, "R1"},
									{DLENGTH_TYPE_RECORD, "R2"},
									{DLENGTH_TYPE_RECORD, "R3"},
									{DLENGTH_TYPE_RECORD, "R4"},
	};
	Table * table = get_shablon_table(22, 1, records, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = alpha1;
		table->table[1].parm_values = alpha2;
		table->table[2].parm_values = alpha3;
		table->table[3].parm_values = alpha4;
		table->table[4].parm_values = A1;
		table->table[5].parm_values = A2;
		table->table[6].parm_values = A3;
		table->table[7].parm_values = A4;
		table->table[8].parm_values = A5;
		table->table[9].parm_values = S1;
		table->table[10].parm_values = S2;
		table->table[11].parm_values = S3;
		table->table[12].parm_values = S4;
		table->table[13].parm_values = S5;
		table->table[14].parm_values = r1;
		table->table[15].parm_values = r2;
		table->table[16].parm_values = r3;
		table->table[17].parm_values = r4;
		table->table[18].parm_values = R1;
		table->table[19].parm_values = R2;
		table->table[20].parm_values = R3;
		table->table[21].parm_values = R4;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = alpha1[k];
		((double *)(table->table[1].parm_values))[k] = alpha2[k];
		((double *)(table->table[2].parm_values))[k] = alpha3[k];
		((double *)(table->table[3].parm_values))[k] = alpha4[k];
		((double *)(table->table[4].parm_values))[k] = A1[k];
		((double *)(table->table[5].parm_values))[k] = A2[k];
		((double *)(table->table[6].parm_values))[k] = A3[k];
		((double *)(table->table[7].parm_values))[k] = A4[k];
		((double *)(table->table[8].parm_values))[k] = A5[k];
		((double *)(table->table[9].parm_values))[k] = S1[k];
		((double *)(table->table[10].parm_values))[k] = S2[k];
		((double *)(table->table[11].parm_values))[k] = S3[k];
		((double *)(table->table[12].parm_values))[k] = S4[k];
		((double *)(table->table[13].parm_values))[k] = S5[k];
		((double *)(table->table[14].parm_values))[k] = r1[k];
		((double *)(table->table[15].parm_values))[k] = r2[k];
		((double *)(table->table[16].parm_values))[k] = r3[k];
		((double *)(table->table[17].parm_values))[k] = r4[k];
		((double *)(table->table[18].parm_values))[k] = R1[k];
		((double *)(table->table[19].parm_values))[k] = R2[k];
		((double *)(table->table[20].parm_values))[k] = R3[k];
		((double *)(table->table[21].parm_values))[k] = R4[k];
	}
	return table;
}

/////////////////////////////////////////////////////
//...таблица для корытообразного профиля с манжетами;
Table * cuffs_profile(int id_static_char)
{
	static double alpha1[] = { 0.,100.},
					  alpha2[] = { 0.,120.},
					  alpha3[] = { 0.,120.},
					  alpha4[] = { 0., 90.},
					  alpha5[] = { 0., 60.},
					  alpha6[] = { 0.,120.},
					  A1    [] = { 0., .05},
					  A2    [] = { 0., .10},
					  A3    [] = { 0., .12},
					  A4    [] = { 0., .07},
					  A5    [] = { 0., .15},
					  A6    [] = { 0., .10},
					  A7    [] = { 0., .04},
					  S1    [] = { 0.,.015},
					  S2    [] = { 0., .02},
					  S3    [] = { 0.,.015},
					  S4    [] = { 0., .02},
					  S5    [] = { 0.,.005},
					  S6    [] = { 0., .01},
					  S7    [] = { 0.,.015},
					  r1    [] = { 0.,.005},
					  r2    [] = { 0.,.015},
					  r3    [] = { 0., .01},
					  r4    [] = { 0.,.007},
					  r5    [] = { 0.,.003},
					  r6    [] = { 0.,.007},
					  R1    [] = { 0.,.007},
					  R2    [] = { 0.,.007},
					  R3    [] = { 0.,.007},
					  R4    [] = { 0.,.007},
					  R5    [] = { 0.,.007},
					  R6    [] = { 0.,.007};
	static char * S0[] = { "Section 1" };
	Shablon  records[] = {	{ DANGLE_TYPE_RECORD, "alpha1"},
									{ DANGLE_TYPE_RECORD, "alpha2"},
									{ DANGLE_TYPE_RECORD, "alpha3"},
									{ DANGLE_TYPE_RECORD, "alpha4"},
									{ DANGLE_TYPE_RECORD, "alpha5"},
									{ DANGLE_TYPE_RECORD, "alpha6"},
									{DLENGTH_TYPE_RECORD, "A1"},
									{DLENGTH_TYPE_RECORD, "A2"},
									{DLENGTH_TYPE_RECORD, "A3"},
									{DLENGTH_TYPE_RECORD, "A4"},
									{DLENGTH_TYPE_RECORD, "A5"},
									{DLENGTH_TYPE_RECORD, "A6"},
									{DLENGTH_TYPE_RECORD, "A7"},
									{DLENGTH_TYPE_RECORD, "S1"},
									{DLENGTH_TYPE_RECORD, "S2"},
									{DLENGTH_TYPE_RECORD, "S3"},
									{DLENGTH_TYPE_RECORD, "S4"},
									{DLENGTH_TYPE_RECORD, "S5"},
									{DLENGTH_TYPE_RECORD, "S6"},
									{DLENGTH_TYPE_RECORD, "S7"},
									{DLENGTH_TYPE_RECORD, "r1"},
									{DLENGTH_TYPE_RECORD, "r2"},
									{DLENGTH_TYPE_RECORD, "r3"},
									{DLENGTH_TYPE_RECORD, "r4"},
									{DLENGTH_TYPE_RECORD, "r5"},
									{DLENGTH_TYPE_RECORD, "r6"},
									{DLENGTH_TYPE_RECORD, "R1"},
									{DLENGTH_TYPE_RECORD, "R2"},
									{DLENGTH_TYPE_RECORD, "R3"},
									{DLENGTH_TYPE_RECORD, "R4"},
									{DLENGTH_TYPE_RECORD, "R5"},
									{DLENGTH_TYPE_RECORD, "R6"},
	};
	Table * table = get_shablon_table(32, 1, records, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = alpha1;
		table->table[1].parm_values = alpha2;
		table->table[2].parm_values = alpha3;
		table->table[3].parm_values = alpha4;
		table->table[4].parm_values = alpha5;
		table->table[5].parm_values = alpha6;
		table->table[6].parm_values = A1;
		table->table[7].parm_values = A2;
		table->table[8].parm_values = A3;
		table->table[9].parm_values = A4;
		table->table[10].parm_values = A5;
		table->table[11].parm_values = A6;
		table->table[12].parm_values = A7;
		table->table[13].parm_values = S1;
		table->table[14].parm_values = S2;
		table->table[15].parm_values = S3;
		table->table[16].parm_values = S4;
		table->table[17].parm_values = S5;
		table->table[18].parm_values = S6;
		table->table[19].parm_values = S7;
		table->table[20].parm_values = r1;
		table->table[21].parm_values = r2;
		table->table[22].parm_values = r3;
		table->table[23].parm_values = r4;
		table->table[24].parm_values = r5;
		table->table[25].parm_values = r6;
		table->table[26].parm_values = R1;
		table->table[27].parm_values = R2;
		table->table[28].parm_values = R3;
		table->table[29].parm_values = R4;
		table->table[30].parm_values = R5;
		table->table[31].parm_values = R6;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = alpha1[k];
		((double *)(table->table[1].parm_values))[k] = alpha2[k];
		((double *)(table->table[2].parm_values))[k] = alpha3[k];
		((double *)(table->table[3].parm_values))[k] = alpha4[k];
		((double *)(table->table[4].parm_values))[k] = alpha5[k];
		((double *)(table->table[5].parm_values))[k] = alpha6[k];
		((double *)(table->table[6].parm_values))[k] = A1[k];
		((double *)(table->table[7].parm_values))[k] = A2[k];
		((double *)(table->table[8].parm_values))[k] = A3[k];
		((double *)(table->table[9].parm_values))[k] = A4[k];
		((double *)(table->table[10].parm_values))[k] = A5[k];
		((double *)(table->table[11].parm_values))[k] = A6[k];
		((double *)(table->table[12].parm_values))[k] = A7[k];
		((double *)(table->table[13].parm_values))[k] = S1[k];
		((double *)(table->table[14].parm_values))[k] = S2[k];
		((double *)(table->table[15].parm_values))[k] = S3[k];
		((double *)(table->table[16].parm_values))[k] = S4[k];
		((double *)(table->table[17].parm_values))[k] = S5[k];
		((double *)(table->table[18].parm_values))[k] = S6[k];
		((double *)(table->table[19].parm_values))[k] = S7[k];
		((double *)(table->table[20].parm_values))[k] = r1[k];
		((double *)(table->table[21].parm_values))[k] = r2[k];
		((double *)(table->table[22].parm_values))[k] = r3[k];
		((double *)(table->table[23].parm_values))[k] = r4[k];
		((double *)(table->table[24].parm_values))[k] = r5[k];
		((double *)(table->table[25].parm_values))[k] = r6[k];
		((double *)(table->table[26].parm_values))[k] = R1[k];
		((double *)(table->table[27].parm_values))[k] = R2[k];
		((double *)(table->table[28].parm_values))[k] = R3[k];
		((double *)(table->table[29].parm_values))[k] = R4[k];
		((double *)(table->table[30].parm_values))[k] = R5[k];
		((double *)(table->table[31].parm_values))[k] = R6[k];
	}
	return table;
}

////////////////////////////////
//...таблица для профиля "крюк";
Table * hook(int id_static_char)
{
	static double alpha1[] = { 0.,150.},
					  alpha2[] = { 0., 60.},
					  alpha3[] = { 0., 90.},
					  A1    [] = { 0., .06},
					  A2    [] = { 0., .10},
					  A3    [] = { 0., .07},
					  A4    [] = { 0., .03},
					  S1    [] = { 0., .02},
					  S2    [] = { 0.,.007},
					  S3    [] = { 0., .01},
					  S4    [] = { 0.,.003},
					  r1    [] = { 0., .01},
					  r2    [] = { 0.,.003},
					  r3    [] = { 0., .01},
					  R1    [] = { 0., .02},
					  R2    [] = { 0., .01},
					  R3    [] = { 0., .01};
	static char * S0[] = { "Section 1" };
	Shablon  records[] = {	{ DANGLE_TYPE_RECORD, "alpha1"},
									{ DANGLE_TYPE_RECORD, "alpha2"},
									{ DANGLE_TYPE_RECORD, "alpha3"},
									{DLENGTH_TYPE_RECORD, "A1"},
									{DLENGTH_TYPE_RECORD, "A2"},
									{DLENGTH_TYPE_RECORD, "A3"},
									{DLENGTH_TYPE_RECORD, "A4"},
									{DLENGTH_TYPE_RECORD, "S1"},
									{DLENGTH_TYPE_RECORD, "S2"},
									{DLENGTH_TYPE_RECORD, "S3"},
									{DLENGTH_TYPE_RECORD, "S4"},
									{DLENGTH_TYPE_RECORD, "r1"},
									{DLENGTH_TYPE_RECORD, "r2"},
									{DLENGTH_TYPE_RECORD, "r3"},
									{DLENGTH_TYPE_RECORD, "R1"},
									{DLENGTH_TYPE_RECORD, "R2"},
									{DLENGTH_TYPE_RECORD, "R3"},
	};
	Table * table = get_shablon_table(17, 1, records, S0, id_static_char);
	if (id_static_char == SPECIAL_STATE) {
		table->table[0].parm_values = alpha1;
		table->table[1].parm_values = alpha2;
		table->table[2].parm_values = alpha3;
		table->table[3].parm_values = A1;
		table->table[4].parm_values = A2;
		table->table[5].parm_values = A3;
		table->table[6].parm_values = A4;
		table->table[7].parm_values = S1;
		table->table[8].parm_values = S2;
		table->table[9].parm_values = S3;
		table->table[10].parm_values = S4;
		table->table[11].parm_values = r1;
		table->table[12].parm_values = r2;
		table->table[13].parm_values = r3;
		table->table[14].parm_values = R1;
		table->table[15].parm_values = R2;
		table->table[16].parm_values = R3;
	}
	else
	for (int k = 1; k <= table->N_group; k++) {
		((double *)(table->table[0].parm_values))[k] = alpha1[k];
		((double *)(table->table[1].parm_values))[k] = alpha2[k];
		((double *)(table->table[2].parm_values))[k] = alpha3[k];
		((double *)(table->table[3].parm_values))[k] = A1[k];
		((double *)(table->table[4].parm_values))[k] = A2[k];
		((double *)(table->table[5].parm_values))[k] = A3[k];
		((double *)(table->table[6].parm_values))[k] = A4[k];
		((double *)(table->table[7].parm_values))[k] = S1[k];
		((double *)(table->table[8].parm_values))[k] = S2[k];
		((double *)(table->table[9].parm_values))[k] = S3[k];
		((double *)(table->table[10].parm_values))[k] = S4[k];
		((double *)(table->table[11].parm_values))[k] = r1[k];
		((double *)(table->table[12].parm_values))[k] = r2[k];
		((double *)(table->table[13].parm_values))[k] = r3[k];
		((double *)(table->table[14].parm_values))[k] = R1[k];
		((double *)(table->table[15].parm_values))[k] = R2[k];
		((double *)(table->table[16].parm_values))[k] = R3[k];
	}
	return table;
}

/*===============================================================================================*/
/*                       ТАБЛИЦЫ ОБРАЗЦОВ ИЗ 3D БИБЛИОТЕКИ LIB_SM                                */
/*===============================================================================================*/
/////////////////////////////////////////////////
//...таблица для crank_shaft (геометpия области);
void crank_shaft(int k, double & alpha0, double & L0, double & A10, double & A20,
                                                      double & B10, double & B20, double & R0)
{
  double alpha[3] = {50.,  5.,  5.},
             L[3] = {25., 25., 25.},
            A1[3] = {15., 10., 12.},
            A2[3] = {10., 10., 10.},
            B1[3] = { 5.,  5.,  5.},
            B2[3] = { 7.,  5.,  7.},
             R[3] = { 2.,  1.,  .5};
  if (k < 1 || k > 3) k = 1;

  alpha0 = alpha[k-1]; L0 = L[k-1]; A10 = A1[k-1]; A20 = A2[k-1];
  B10 = B1[k-1]; B20 = B2[k-1]; R0 = R[k-1]; return;
}

///////////////////////////////////////////////////////
//...таблица для усеченного конуса (геометpия области);
void intr_conus(int k, double & r10, double & r20, double & L10, double & L20)
{
  double r1[3] = {5.,  10.,  5.},
         r2[3] = {10.,  5.,  5.},
         L1[3] = {0.,   0.,  0.},
         L2[3] = {15., 15., 25.};
  if (k < 1 || k > 3) k = 1;

  r10 = r1[k-1]; r20 = r2[k-1]; L10 = L1[k-1]; L20 = L2[k-1]; return;
}

////////////////////////////////////////////////////////
//...таблица для L-образной комнаты (геометpия области);
void L_type_room(int k, double & A10, double & A20, double & B10, double & B20,
                        double & C0,  double & R0)
{
  double A1[3] = {1.5,  1.,  1.0001},
         A2[3] = {1.5,  1.5, 1.5},
         B1[3] = {1.,   1.5, 1.0},
         B2[3] = {1.,   1.,  1.0},
         C [3] = {1.,   1.,  1.0},
         R [3] = {0.5,  0.,  0.0};
  if (k < 1 || k > 3) k = 1;

  A10 = A1[k-1]; A20 = A2[k-1]; B10 = B1[k-1]; B20 = B2[k-1]; C0 = C[k-1]; R0 = R[k-1]; return;
}

//////////////////////////////////////////////
//...таблица для hole_box (геометpия области);
void hole_box(int k, double & A0, double & B0, double & C0, double & R0, double & L0)
{
  double A[2] = {2.,  2.},
         B[2] = {2.,  4.},
         C[2] = {2.,  5.},
         R[2] = {.5, .75},
         L[2] = {1.,  3.};
  if (k < 1 || k > 2) k = 1;

  A0 = A[k-1]; B0 = B[k-1]; C0 = C[k-1]; R0 = R[k-1]; L0 = L[k-1]; return;
}

//////////////////////////////////////////////
//...таблица для poly_box (геометpия области);
void poly_box(int k, double ** P0)
{
  double P[2][24] = {{ 0., 0., 0.,  2., 0., 0.,  2., 2., 0.,  0., 2., 0.,
                       0., 0., 2.,  1., 0., 2.,  1., 1., 2.,  0., 1., 2.},
                     { 0., 0., 0.,  2., 0., 0.,  2., 2., 0.,  0., 2., 0.,
                       0., 0., 2.,  2., 0., 2.,  2., 2., 2.,  0., 2., 2.}};
  if (k < 1 || k > 2) k = 1;

  for (int m = 0; m < 24; m++) P0[m][k] = P[k-1][m]; return;
}
#endif

