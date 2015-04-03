#include "stdafx.h"
#include "ccebasic.h"

///////////////////////////////////////////////////////////
//         CONSTRUCTION of the ALL CELLS TYPES           //
///////////////////////////////////////////////////////////
//////////////////////////////////////////////////////
//...global function for construction all cells types;
CCells * CreateCells(Num_Cells id_CELLS)
{
	switch (id_CELLS) {
	  case  NULL_CELLS: return new CCells;
	  case BASIC_CELLS: return new CCeBasic;
	}
   return NULL;
}
