#include "stdafx.h"

#include "cgrid_el.h"
#include "cgrid_qg.h"

///////////////////////////////////////////////////////////
//         CONSTRUCTION of the ALL GRID NODES            //
///////////////////////////////////////////////////////////
/////////////////////////////////////////////////////
//...global function for construction all grid nodes;
CGrid * CreateNodes(Num_Nodes id_NODES)
{
	switch (id_NODES) {
     case    NULL_NODES: return new CGrid;
     case GRID_EL_NODES: return new CGrid_el;
     case GRID_QG_NODES: return new CGrid_QG;
	}
   return NULL;
}
