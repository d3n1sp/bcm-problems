#include "stdafx.h"

#include "cprofile.h"
#include "cpro_skin.h"

#include "ccebasic.h"
#include "ccedemo.h"

#include "clame2d.h"

#include "ccohes2d.h"
#include "cgrid_el.h"

#include "unit_mes.h"

#ifdef ___WINDOWS_LOG_MESSAGE___
#include "..\vc2d\externs.h"
#define  Message(Msg)    theMainFrame->Message(Msg)
#else
#define  Message(Msg)    printf(Msg);  printf("\n")
#endif

#ifdef ___PRO_SKIN2D_cpp___
///////////////////////////////////////////////////////////////////////////////////////////////////
//...число 2D статических образцов дл€ скин-задач и их имена
int    GetSK2DSampleCount(void) { return(NUM_SK2D_SAMPLES); }
char * GetSK2DSampleName (int N_sm)
{
  static char * S[] = { "Cell with ellipse  ", //_SK1 
                        "Extern cell  "        //_SK2
  };
  return S[N_sm];
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//...функци€ образовани€ контекста дл€ скин-задачки
void * CreateSK2DContext(int N_sm)
{
  if (N_sm < 0 || N_sm >= NUM_SK2D_SAMPLES) N_sm = 0;

  Context * cont = (Context *)new_struct(sizeof(Context));
  if (! cont) return(NULL);

  cont->N           = N_sm+SHIFT_SK2D_SAMPLES;
  cont->static_char = 1;
  cont->sample_name = GetSK2DSampleName(N_sm);
  cont->units       = UNIT_SI;

  switch (cont->N-SHIFT_SK2D_SAMPLES) {
		case _SK1:  { 
						static char  * S[] = { "Double plane problem", "Classic model", "Cohesion model" };
						Shablon  records[] = {	{DLENGTH_TYPE_RECORD, "Semi-length"},//(0) -- L
														{DLENGTH_TYPE_RECORD, "Height"},     //(1) -- H
														{DLENGTH_TYPE_RECORD, "Major axis"}, //(2) -- a
														{DLENGTH_TYPE_RECORD, "Minor axis"}, //(3) -- b
														{ DANGLE_TYPE_RECORD, "Angle"},      //(4) -- theta
														{ DOUBLE_TYPE_RECORD, "eps"},        //(5) -- eps
														{ DOUBLE_TYPE_RECORD, "nju1"},       //(6) -- nju1
														{ DOUBLE_TYPE_RECORD, "nju2"},       //(7) -- nju2
														{ DOUBLE_TYPE_RECORD, "Ef/Em"},      //(8) -- E_ratio
														{ DOUBLE_TYPE_RECORD, "C0"},         //(9) -- C0
														{	  INT_TYPE_RECORD, "N1"},         //(10) -- N1
														{	  INT_TYPE_RECORD, "N2"},         //(11) -- N2
														{	  INT_TYPE_RECORD, "N"},          //(12) -- N
														{	  INT_TYPE_RECORD, "index"},      //(13) -- index
						};
						if ((cont->table = get_shablon_table(14, 3, records, S)) == NULL) {
							delete[] cont; return(NULL);
						}
						int k;
						for (k = 1; k <= cont->table->N_group; k++)
							((double *)(cont->table->table[0].parm_values))[k] = 1.;
							((double *)(cont->table->table[1].parm_values))[k] = 1.2;
							((double *)(cont->table->table[2].parm_values))[k] = 1.4;
							((double *)(cont->table->table[3].parm_values))[k] = 0.5;
							((double *)(cont->table->table[4].parm_values))[k] = 0.;
							((double *)(cont->table->table[5].parm_values))[k] = 0.9;
							((double *)(cont->table->table[6].parm_values))[k] = 0.3;
							((double *)(cont->table->table[7].parm_values))[k] = 0.3;
							((double *)(cont->table->table[8].parm_values))[k] = 2.;
							((double *)(cont->table->table[9].parm_values))[k] = 10.;
							((int    *)(cont->table->table[10].parm_values))[k] = 20;
							((int    *)(cont->table->table[11].parm_values))[k] = 20;
							((int    *)(cont->table->table[12].parm_values))[k] = 5;
							((int    *)(cont->table->table[13].parm_values))[k] = k;
						cont->left   = 0.; 
						cont->right  = 0.; 
						cont->bottom = 0.; 
						cont->top    = 0.;
		}				break;
		case _NUM_SK2D_SAMPLES-1:  { 
						static char  * S[] = { "Double plane problem", "Classic model", "Cohesion model" };
						Shablon  records[] = {	{DOUBLE_TYPE_RECORD, "nju1"}, //(0) -- nju1
														{DOUBLE_TYPE_RECORD, "nju2"}, //(1) -- nju2
														{DOUBLE_TYPE_RECORD, "Ef/Em"},//(2) -- E_ratio
														{DOUBLE_TYPE_RECORD, "C0"},   //(3) -- C0
														{	INT_TYPE_RECORD,  "N"},    //(4) -- N
														{	INT_TYPE_RECORD,  "index"},//(5) -- index
						};
						if ((cont->table = get_shablon_table(6, 3, records, S)) == NULL) {
							delete[] cont; return(NULL);
						}
						int k;
						for (k = 1; k <= cont->table->N_group; k++)
							((double *)(cont->table->table[0].parm_values))[k] = 0.3;
							((double *)(cont->table->table[1].parm_values))[k] = 0.3;
							((double *)(cont->table->table[2].parm_values))[k] = 2.;
							((double *)(cont->table->table[3].parm_values))[k] = 10.;
							((int    *)(cont->table->table[4].parm_values))[k] = 2;
							((int    *)(cont->table->table[5].parm_values))[k] = k;
						cont->left   = 0.; 
						cont->right  = 0.; 
						cont->bottom = 0.; 
						cont->top    = 0.;
		}				break;
  }
  set_default    (cont);
  set_table_units(cont->table, UNIT_SI); return (void *)cont;
}

/*===============================================================================================*/
/*                                 »Ќ»÷»јЋ»«ј÷»я ќЅ–ј«÷ј                                         */
/*===============================================================================================*/
///////////////////////////////////////////////////////////////////////////////////////////////////
//...пpедваpительна€ инициализаци€ 2D статического обpазца дл€ скин-задач
int skin2D_init(void * context, CGrid * block_nd)
{
  Context * cont = (Context *)context;
  if (! is_Skin2D(cont)  || ! cont->table) return(0);
  if (cont->table->index > 0) cont->units = cont->table->table_units[cont->table->index-1];

  double  L, H, eps, a, b, x0 = 0., y0 = 0., theta, nju1, nju2, E_ratio, C0;
  int     N1, N2, N, i, index;
  Table * tab = cont->table;

  cont->left = cont->right = cont->bottom = cont->top = cont->back = cont->front = 0.;

  if (cont->N != _NUM_SK2D_SAMPLES-1+SHIFT_SK2D_SAMPLES) 
      DeleteSample(cont->sm);

  switch (cont->N-SHIFT_SK2D_SAMPLES) {
     case _SK1: { L       = ((double *)tab->table[0].parm_values)[tab->index];
                  H       = ((double *)tab->table[1].parm_values)[tab->index];
                  eps     = ((double *)tab->table[5].parm_values)[tab->index];
                  a       = ((double *)tab->table[2].parm_values)[tab->index]*eps;
                  b       = ((double *)tab->table[3].parm_values)[tab->index]*eps;
                  theta   = ((double *)tab->table[4].parm_values)[tab->index];
                  nju1    = ((double *)tab->table[6].parm_values)[tab->index];
                  nju2    = ((double *)tab->table[7].parm_values)[tab->index];
                  E_ratio = ((double *)tab->table[8].parm_values)[tab->index];
                  C0      = ((double *)tab->table[9].parm_values)[tab->index];
                  N1      = ((int    *)tab->table[10].parm_values)[tab->index];
                  N2      = ((int    *)tab->table[11].parm_values)[tab->index];
                  N       = ((int    *)tab->table[12].parm_values)[tab->index];
                  index   = ((int    *)tab->table[13].parm_values)[tab->index];
                  if (index == 2 && (cont->sm = CreateSample(LAME2D_SAMPLE)) != 0 ||
                      index == 3 && (cont->sm = CreateSample(COHES2D_SAMPLE)) != 0) {

                      CCeDemo * ce = new CCeDemo;

                      ce->get_test_inclusion2D(L*2., H, a, b, x0, y0, theta*M_PI/180.);
                      cont->sm->init_blocks(ce);

///////////////////////////
//...parameters of problem;
                      cont->sm->set_param(0, (double)N); //...number of multipoles;
                      cont->sm->set_param(1, .92);
                      cont->sm->set_param(2, 100.);
                      switch (index) {
                      case 2 : ((CLame2D *)cont->sm)->set_fasa_hmg (nju1, nju2, 1., E_ratio); break;
                      case 3 : ((CCohes2D *)cont->sm)->set_fasa_hmg(nju1, nju2, 1., E_ratio, C0, C0*E_ratio); break;
                      }
                  }
/////////////////////////////
//...переустанавливаем рамку;
                  cont->left   = 0.;
                  cont->right  = L*2.;
                  cont->bottom = 0.;
                  cont->top    = H; 
     }            break;
     case _NUM_SK2D_SAMPLES-1: { 
                  nju1    = ((double *)tab->table[0].parm_values)[tab->index];
                  nju2    = ((double *)tab->table[1].parm_values)[tab->index];
                  E_ratio = ((double *)tab->table[2].parm_values)[tab->index];
                  C0      = ((double *)tab->table[3].parm_values)[tab->index];
                  N       = ((int    *)tab->table[4].parm_values)[tab->index];
                  index   = ((int    *)tab->table[5].parm_values)[tab->index];
                  if ( cont->sm || 
                       index == 2 && (cont->sm = CreateSample(LAME2D_SAMPLE)) != 0 ||
                       index == 3 && (cont->sm = CreateSample(COHES2D_SAMPLE)) != 0) {

/////////////////////////////////////////
//...перегрузка геометрии в новую задачу;
                       if (cont->sm->N) {
                          CSample * sm = 0;
                          if (index == 2) sm = CreateSample(LAME2D_SAMPLE); else
                          if (index == 3) sm = CreateSample(COHES2D_SAMPLE);

                          if (sm) {
                              sm->init_blocks(NULL);
                              for (int k = 0; k < cont->sm->N; k++) {
                                   sm->add_block(POLY_BLOCK);
                                   sm->B[k].bar   = cont->sm->B[k].bar;   cont->sm->B[k].bar = 0;
                                   sm->B[k].link  = cont->sm->B[k].link;  cont->sm->B[k].link = 0;
                                   sm->B[k].mp    = cont->sm->B[k].mp;    cont->sm->B[k].mp = 0;
                                   sm->B[k].shape = cont->sm->B[k].shape; cont->sm->B[k].shape = 0;
                              }
                              SetSample(cont, sm);
                          }
                       }

///////////////////////////
//...parameters of problem;
                       cont->sm->set_param(0, (double)N); //...number of multipoles;
                       cont->sm->set_param(1, .92);
                       cont->sm->set_param(2, 100.);
                       switch (index) {
                       case 2 : ((CLame2D *)cont->sm)->set_fasa_hmg (nju1, nju2, 1., E_ratio); break;
                       case 3 : ((CCohes2D *)cont->sm)->set_fasa_hmg(nju1, nju2, 1., E_ratio, C0, C0*E_ratio); break;
                       }
                       if (! cont->sm->N) {
                           double HX = 2, HY = 1.2, X = 1., Y = 0.6, R = sqrt(HX*HX+HY*HY);
                           if (cont->sm->add_block(POLY_BLOCK)                                     &&
                               cont->sm->set_block(cont->sm->B[cont->sm->N-1], SPHERE_GENUS, R, 0, X, Y, 0.) &&
                               cont->sm->set_link (cont->sm->B[cont->sm->N-1], 1)) {
                               cont->sm->B[cont->sm->N-1].bar = new CCells;
                               double P[12] = {0., 0., 0., 2., 0., 0., 2., 1.2, 0., 0., 1.2, 0.}, 
                                      pp[4] = {0., 0., 0., 0.}, A, B, C, D, E, R, Q;

                               CCeBasic * ce = new CCeBasic;
                               ce->get_facet_directly(P, 4, 0, 2, 0);

                               i = cont->sm->B[cont->sm->N-1].bar->bar_add(ce);
                               A = abs_point(P,   P+3);
                               B = abs_point(P+3, P+6);
                               E = abs_point(P,   P+6); R = (A+B+E)*.5;
                               C = abs_point(P+6, P+9);
                               D = abs_point(P+9, P);   Q = (C+D+E)*.5;

                               ((CCeBasic *)cont->sm->B[cont->sm->N-1].bar)->SetFacetParam(i-1, sqrt(R*(R-A)*(R-B)*(R-E))+sqrt(Q*(Q-C)*(Q-D)*(Q-E)), pp);
                           }
                       }
                  }
/////////////////////////////
//...переустанавливаем рамку;
                  double par[6] = {MAX_HIT, MIN_HIT, MAX_HIT, MIN_HIT, MAX_HIT, MIN_HIT};
                  for (i = 0; i < cont->sm->N; i++) 
                      cont->sm->SkeletonBounding(cont->sm->B[i], par);

                  if (cont->sm->N) {
                      cont->left   = par[0];
                      cont->right  = par[1];
                      cont->bottom = par[2];
                      cont->top    = par[3]; 
                  }
     }            break;
  }
  return(cont->sm != NULL);
}

//////////////////////////////////////////////////////////////////////
//...pешатель функций фоpмы дл€ 2D образца (задачи пограничного сло€);
int skin2D_solv(void * context, int & sec, int & hund, int id_solv)
{
    Context  * cont = (Context *)context; sec = hund = 0;
    if ((is_Skin2D(cont) || is_SkGp3D(cont) || is_SkIb3D(cont)) && cont->sm && cont->sm->N > 1) {
        CSample * sm = cont->sm;

///////////////////////////
//...solving of the probem;
        clock_t start = clock();
        int  counting = (cont->N == _NUM_SK2D_SAMPLES-1+SHIFT_SK2D_SAMPLES ? 
                         ANALYTICAL_COUNTING : 
                             BASIC_COUNTING);

        if (sm->counting_kernel(counting) != OK_COUNTING) {
            Message("Error in sample counting...");
            return(0);
        }
        Message("O'K");

////////////////////////////////////////////
//...фиксируем врем€ и выходим из программы;
        timing_process(start, & hund, & sec);
    }
    return(1);
}

//////////////////////////////////////////////////////////
//...результаты дл€ 2D образца (задачи пограничного сло€);
char * GetFuncSKName(int num)
{ static char *L[]={ "Full displacement Rx",     //_FK1
                     "Classic displacement Ux",  //_FK2
                     "Full displacement Ry",     //_FK3
                     "Classic displacement Uy",  //_FK4
                     "Full normal stress",       //_FK5
                     "Classic normal stress",    //_FK6
                     "Full shear stress",        //_FK7
                     "Classic shear stress",     //_FK8
                     "Full transverse stress",   //_FK9
                     "Classic transverse stress" //_FK10
                   };
  return L[num];
}

void solver_SK(void * context, CGrid * nd, double * F, double * par, int id_F)
{
  Context * cont = (Context *)context;
  if (cont && cont->sm && nd && F) {
      int id_block, l;

////////////////////////////////////
//...обpабатываем все входные точки;
      if (nd->X && nd->Y && cont->sm->B)
      for (l = 0; l < nd->N; l++) {
           if (nd->hit) id_block = nd->hit[l];
           else         id_block = -1;
//////////////////////////////////
//...search for appropriate block;
           if (id_block < 0 && Sph3D_struc_in(cont->sm, id_block, nd->X[l], nd->Y[l], 0.) &&
               cont->sm->B[id_block].link && par  &&
               cont->sm->B[id_block].link[0] == 6 && cont->sm->bar && 
             ((CCeDemo *) cont->sm->bar)->in_test_inclusion2D(nd->X[l], nd->Y[l], par) != -cont->sm->B[id_block].link[5])
               id_block = cont->sm->B[id_block].link[6];
//////////////////////
//...data calculation;
           if (id_block >= 0 && id_block < cont->sm->N) {
               double FF[3] = {0., 0., 0.};
               switch (id_F) {
                      case _FK1: cont->sm->GetFuncAllValues(nd->X[l], nd->Y[l], 0., FF, id_block, 0);
                            F[l] = FF[0];
                            break;
                      case _FK2: cont->sm->GetFuncAllValues(nd->X[l], nd->Y[l], 0., FF, id_block, 1);
                            F[l] = FF[0];
                            break;
                      case _FK3: cont->sm->GetFuncAllValues(nd->X[l], nd->Y[l], 0., FF, id_block, 0);
                            F[l] = FF[1];
                            break;
                      case _FK4: cont->sm->GetFuncAllValues(nd->X[l], nd->Y[l], 0., FF, id_block, 1);
                            F[l] = FF[1];
                            break;
                      case _FK5: cont->sm->GetFuncAllValues(nd->X[l], nd->Y[l], 0., FF, id_block, 2);
                            F[l] = FF[0];
                            break;
                      case _FK6: cont->sm->GetFuncAllValues(nd->X[l], nd->Y[l], 0., FF, id_block, 4);
                            F[l] = FF[0];
                            break;
                      case _FK7: cont->sm->GetFuncAllValues(nd->X[l], nd->Y[l], 0., FF, id_block, 2);
                            F[l] = FF[1];
                            break;
                      case _FK8: cont->sm->GetFuncAllValues(nd->X[l], nd->Y[l], 0., FF, id_block, 4);
                            F[l] = FF[1];
                            break;
                      case _FK9: cont->sm->GetFuncAllValues(nd->X[l], nd->Y[l], 0., FF, id_block, 3);
                            F[l] = FF[1];
                            break;
                      case _FK10: cont->sm->GetFuncAllValues(nd->X[l], nd->Y[l], 0., FF, id_block, 5);
                            F[l] = FF[1];
                            break;
                      case _NUM_SK_FUNCTIONS: 
                            F[l] = id_block;
                      default:
                            F[l] = 0.;
               }
           }
           if (nd->hit) nd->hit[l] = id_block;
      }
  }
}
#endif
#undef  Message
