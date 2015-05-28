#include "stdafx.h"

#include "drafts.h"

//////////////////////////////
//...construction real drafts;
CDraft<double> * CreateDraftR(Num_Draft id_DRAFT, int id_dop)
{
	return CreateDraft<double>(id_DRAFT, id_dop);
}

/////////////////////////////////
//...construction complex drafts;
CDraft<complex> * CreateDraftC(Num_Draft id_DRAFT, int id_dop)
{
	return CreateDraft<complex>(id_DRAFT, id_dop);
}

///////////////////////////////////////
//...construction double-double drafts;
CDraft<dd_real> * CreateDraftD(Num_Draft id_DRAFT, int id_dop)
{
	return CreateDraft<dd_real>(id_DRAFT, id_dop);
}

/////////////////////////////////////
//...construction quad-double drafts;
CDraft<qd_real> * CreateDraftQ(Num_Draft id_DRAFT, int id_dop)
{
	return CreateDraft<qd_real>(id_DRAFT, id_dop);
}
