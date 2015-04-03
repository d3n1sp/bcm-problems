//------------------------------------------------------------------------------------------------
// File: SolverPar.h
// Creator: Kharchenko S.A.
//------------------------------------------------------------------------------------------------
#include<cstdlib>
#include<cstring>

#include "globals.h"
#include "slvparam.h"

#ifdef __FV_2004__
namespace flowvision {
#endif

#ifndef __SolverPar
#define __SolverPar

/////////////////////////////////////////////////////////////////////////////
// Classes declared in this file

#ifdef __FV_2004__
class EqnSolver::CSlvParam;
#else
class CSlvParam;
#endif

class CSolverPar
{
// Data
#ifdef __FV_2004__
	EqnSolver::CSlvParam slvparam; // slvparam structure specifies parameters of the solver
#else
	CSlvParam slvparam; // slvparam structure specifies parameters of the solver
#endif
public:
// Get/set functions
#ifdef __FV_2004__
	EqnSolver::CSlvParam * GetParam () {return &slvparam;}; // Get slvparam parameter
#else
	CSlvParam * GetParam () {return &slvparam;}; // Get slvparam parameter
#endif
	int GetNcpuext () {return slvparam.ncpuext;}; // Get ncpuext parameter
	int GetNcpuextsch () {return slvparam.ncpuextsch;}; // Get ncpuextsch parameter
	void SetDecomp (int _nchilds, int _weightcellvolume, int _weightcellboundary) { // Set decomposition parameters
		slvparam.nchilds = _nchilds;
		slvparam.weightcellvolume   = _weightcellvolume;
		slvparam.weightcellboundary = _weightcellboundary;
	};
	void SetNcycle (int _ncycle) { // Set ncycle parameter
		slvparam.ncycle = _ncycle;
	};
	void SetNcpupart (int _ncpupart) { // Set ncpupart parameter
		slvparam.ncpupart = _ncpupart;
	};
	void SetNcpuext (int _ncpuext) { // Set ncpuext parameter
		slvparam.ncpuext = _ncpuext;
	};
	void SetNcpuextsch (int _ncpuextsch) { // Set ncpuext parameter
		slvparam.ncpuextsch = _ncpuextsch;
	};
	void SetSchurtype (int _schurtype) { // Set schurtype parameter
		slvparam.schurtype = _schurtype;
	};
	void SetNlevschur (int _nlevschur) { // Set nlevschur parameter
		slvparam.nlevschur = _nlevschur;
	};
#ifdef __FV_2004__
	void SetParam (EqnSolver::CSlvParam &_slvparam) { // Set ncpuext parameter
#else
	void SetParam (CSlvParam &_slvparam) { // slvparam structure specifies parameters of the solver
#endif
		slvparam = _slvparam;
	};
// Friend classes
#ifdef __FV_2004__
	friend class EqnSolver::CDecomp;
#else
	friend class CDecomp;
#endif
};

#endif

#ifdef __FV_2004__
} // namespace flowvision
#endif


