/*=========================================*/
/*                  SHAPES                 */
/*=========================================*/
#ifndef ___SHAPES___
#define ___SHAPES___

#include <typeinfo>

#include "cbeam.h"
#include "cexpp.h"
#include "cpoly.h"
#include "cmapi.h"
#include "cheat.h"
#include "cwave.h"
#include "cacou.h"

////////////////////////////////////////////////////////////////////////////////
//                    ALL EXISTING TYPES OF SHAPE FUNCTIONS                   //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//...global functions for construction of all existing types of shape functions;
template <typename T>
CShape<T> * CreateShape(Num_Shape id_SHAPE = NULL_SHAPE, int id_dop = 0)
{
	switch (id_SHAPE) {
      case			  MP2D_POLY_SHAPE: return new CMapi2DPoly<T>;
      case			  MP3D_POLY_SHAPE: return new CMapi3DPoly<T>;
      case			  MP3D_ZOOM_SHAPE: return new CMapi3DZoom<T>;
      case			  MP3D_ELLI_SHAPE: return new CMapi3DEll<T>;
      case			  BEAM_POLY_SHAPE: return new CBeamShape<T>;
      case			MP2D_CORNER_SHAPE: return new CMapi2DCorner<T>;
      case			MP2D_CLAYER_SHAPE: return new CMapi2DClayer<T>;
      case			MP2D_CIRCLE_SHAPE: return new CMapi2DCircle<T>;
      case			MP3D_CLAYER_SHAPE: return new CMapi3DClayer<T>;
      case			MP3D_SPHERE_SHAPE: return new CMapi3DSphere<T>;
      case		 MP3D_SPHEROID_SHAPE: if (typeid(T) == typeid(double)) return (CShape<T> *)(new CMapi3DSpheroid);
      case MP3D_SPHEROID_FULL_SHAPE: if (typeid(T) == typeid(double)) return (CShape<T> *)(new CMapi3DSpheroidFull);
      case			  HT1D_POLY_SHAPE: return new CHeat1DPoly<T>;
      case			  HT2D_POLY_SHAPE: return new CHeat2DPoly<T>;
      case			  HT3D_POLY_SHAPE: return new CHeat3DPoly<T>;
      case			HT2D_STREEP_SHAPE: return new CHeat2DStreep<T>;
      case			HT2D_CIRCLE_SHAPE: return new CHeat2DCircle<T>;
      case			HT3D_SPHERE_SHAPE: return new CHeat3DSphere<T>;
      case			  SK2D_BEAMZSHAPE: return new CSkin2DBeamZ<T>;
      case			  SK3D_BEAM_SHAPE: return new CSkin3DBeam<T>;
      case			  SK2D_ELLI_SHAPE: return new CSkin2DEll<T>;
      case			  SK2D_POLY_SHAPE: return new CSkin2DPoly<T>;
      case			  SK3D_POLY_SHAPE: return new CSkin3DPoly<T>;
      case			  SK3D_ZOOM_SHAPE: return new CSkin3DZoom<T>;
      case			  SK3D_EXPP_SHAPE: return new CSkin3DExpp<T>(id_dop);
      case			  WV2D_POLY_SHAPE: if (typeid(T) == typeid(double)) return (CShape<T> *)(new CWave2DPoly);
      case			  WV3D_POLY_SHAPE: if (typeid(T) == typeid(double)) return (CShape<T> *)(new CWave3DPoly);
      case			  AU3D_ELLI_SHAPE: if (typeid(T) == typeid(double)) return (CShape<T> *)(new CAcou3DEll);
      case			  AU3D_POLY_SHAPE: if (typeid(T) == typeid(double)) return (CShape<T> *)(new CAcou3DPoly);
      case			  AU3D_ZOOM_SHAPE: if (typeid(T) == typeid(double)) return (CShape<T> *)(new CAcou3DZoom);
      case			  AU3D_WAVE_SHAPE: if (typeid(T) == typeid(double)) return (CShape<T> *)(new CAcou3DWave);
      case			  AU3D_BEAM_SHAPE: if (typeid(T) == typeid(double)) return (CShape<T> *)(new CAcou3DBeam);
      case			  AU3D_BEAMZSHAPE: if (typeid(T) == typeid(double)) return (CShape<T> *)(new CAcou3DBeamZ);
	}
   return new CShape<T>;
}
#endif

