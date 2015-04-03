#include "stdafx.h"

#include "cbeam.h"
#include "cexpp.h"
#include "cpoly.h"
#include "cmapi.h"
#include "cheat.h"
#include "cwave.h"
#include "cacou.h"

////////////////////////////////////////////////////////////////////////////
//            CONSTRUCTION of the AANALYTICAL SHAPE FUNCTIONS             //
////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//...global functions for construction all analytical shape functions;
template <typename T>
CShape<T> * CreateShape(Num_Shape id_SHAPE, int id_dop)
{
	switch (id_SHAPE) {
      case			 NULL_SHAPE: return new CShape<T>;
      case	  MP2D_POLY_SHAPE: return new CMapi2DPoly<T>;
      case	  MP3D_POLY_SHAPE: return new CMapi3DPoly<T>;
      case	  MP3D_ZOOM_SHAPE: return new CMapi3DZoom<T>;
      case	  MP3D_ELLI_SHAPE: return new CMapi3DEll<T>;
      case	  BEAM_POLY_SHAPE: return new CBeamShape<T>;
      case	MP2D_CORNER_SHAPE: return new CMapi2DCorner<T>;
      case	MP2D_CLAYER_SHAPE: return new CMapi2DClayer<T>;
      case	MP2D_CIRCLE_SHAPE: return new CMapi2DCircle<T>;
      case	MP3D_CLAYER_SHAPE: return new CMapi3DClayer<T>;
      case	MP3D_SPHERE_SHAPE: return new CMapi3DSphere<T>;
      case	  HT1D_POLY_SHAPE: return new CHeat1DPoly<T>;
      case	  HT2D_POLY_SHAPE: return new CHeat2DPoly<T>;
      case	  HT3D_POLY_SHAPE: return new CHeat3DPoly<T>;
      case	HT2D_STREEP_SHAPE: return new CHeat2DStreep<T>;
      case	HT2D_CIRCLE_SHAPE: return new CHeat2DCircle<T>;
      case	HT3D_SPHERE_SHAPE: return new CHeat3DSphere<T>;
      case    SK2D_BEAMZSHAPE: return new CSkin2DBeamZ<T>;
      case    SK3D_BEAM_SHAPE: return new CSkin3DBeam<T>;
      case    SK2D_ELLI_SHAPE: return new CSkin2DEll<T>;
      case	  SK2D_POLY_SHAPE: return new CSkin2DPoly<T>;
      case	  SK3D_POLY_SHAPE: return new CSkin3DPoly<T>;
      case	  SK3D_ZOOM_SHAPE: return new CSkin3DZoom<T>;
      case	  SK3D_EXPP_SHAPE: return new CSkin3DExpp<T>(id_dop);
	}
   return NULL;
}
CShape <double> * CreateShapeD(Num_Shape id_SHAPE, int id_dop) {return CreateShape <double>(id_SHAPE, id_dop);}
CShape<dd_real> * CreateShapeR(Num_Shape id_SHAPE, int id_dop) {return CreateShape<dd_real>(id_SHAPE, id_dop);}
CShape<qd_real> * CreateShapeQ(Num_Shape id_SHAPE, int id_dop) {return CreateShape<qd_real>(id_SHAPE, id_dop);}
CShape<complex> * CreateShapeC(Num_Shape id_SHAPE, int id_dop) {return CreateShape<complex>(id_SHAPE, id_dop);}

////////////////////////////////////////////////////////////////
//...construction additional list of the double shape functions;
CShape<double> * CreateShapeS(Num_Shape id_SHAPE, int id_dop)
{
	switch (id_SHAPE) {
      case		 MP3D_SPHEROID_SHAPE: return new CMapi3DSpheroid;
      case MP3D_SPHEROID_FULL_SHAPE: return new CMapi3DSpheroidFull;
      case			  WV2D_POLY_SHAPE: return new CWave2DPoly;
      case			  WV3D_POLY_SHAPE: return new CWave3DPoly;
      case			  AU3D_ELLI_SHAPE: return new CAcou3DEll;
      case			  AU3D_POLY_SHAPE: return new CAcou3DPoly;
      case			  AU3D_ZOOM_SHAPE: return new CAcou3DZoom;
      case			  AU3D_WAVE_SHAPE: return new CAcou3DWave;
      case			  AU3D_BEAM_SHAPE: return new CAcou3DBeam;
      case			  AU3D_BEAMZSHAPE: return new CAcou3DBeamZ;
	}
   return NULL;
}
