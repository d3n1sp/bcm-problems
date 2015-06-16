//------------------------------------------------------------------------------------------------
// File: slvparam.h
// Creator: Kharchenko S.A.
//------------------------------------------------------------------------------------------------
#include<cstdlib>
#include<cstring>

#ifndef __SlvParam
#define __SlvParam

/////////////////////////////////////////////////////////////////////////////
// Classes declared in this file

#ifdef __FlowVisionOut
	class CRegion;    // FlowVision class
	class sMethodSet; // FlowVision class
#endif

class CSlvParam
{
public:
// Data
	int msglev;       // msglev is the output level
	int nchilds;      // nchilds is the number of childs in the processors tree
	int ncpupart;     // ncpupart is the number cpus used for partitioning computations
	int ncpuext;      // ncpuext is the number additional blocks in the extended binary tree
	int ncpuextsch;   // ncpuextsch is the number additional blocks in the extended binary tree for Schur complement blocks
	int decomptype;   // decomptype is the type of decomposition
	int schurtype;    // schurtype is the type of ND decomposition (with or without Schur blocks distribution)
	int nsupmax;      // nsupmax is the number of supernodes per block in the processors tree
	int weightcellvolume;   // weightcellvolume   is the total allowable weight of the volume   hypercell
	int weightcellboundary; // weightcellboundary is the total allowable weight of the boundary hypercell
	int fcttyp;       // fcttyp is the type of the factorization (0 - ilu0, 1 - second order, 2 - second order with special thresh for bordering, 3 or 4 - second order with limited number of elements in a row, 5 - combined element level and thresholds are used for filtering)
	int ooctyp;       // ooctyp is the out-of-core type
	int ordtype;      // ordtype is the type of the ordering
	int ordtypeblk;   // ordtypeblk is the type of the block ordering
	int collap;       // collap is the collapse parameter
	int ncycle;       // ncycle parameter specifies the number of cycles used to determine extended sparsity
	int nprecini;     // nprecini is the number of iterations with preconditioner recomputation at each iteration
	int nprec;        // nprec is the number of iterations per preconditioner recomputation
	int sclpttype;    // sclpttype is the type of point explicit scaling
	double alpha;     // alpha is the Tikhonov regularization parameter
	double dshift;    // dshift parameter specifies the diagonal shift after scaling
	double dshiftim;  // dshiftim parameter specifies the imaginary part of the diagonal shift after scaling
	double memory;    // memory parameter specifies the relative amount of memory available for ich2 fct arrays
	double tau1;      // tau1 is the first  order factorization threshold
	double tau2;      // tau2 is the second order factorization threshold
	double tau1save;  // tau1save is the saved first  order factorization threshold for future use
	double tau2save;  // tau2save is the saved second order factorization threshold for future use
	int ndiagfct;     // ndiagfct  is the maximal number of first  order elements per row in diagonal blocks
	int ndiagfct2;    // ndiagfct2 is the maximal number of second order elements per row in diagonal blocks
	int nlevfct;      // nlevfct  is the maximal level of the element to be considered as the first  order element in a row
	int nlevfct2;     // nlevfct2 is the maximal level of the element to be considered as the second order element in a row
	int nlevschur;    // nlevschur is the maximal level of the element to be considered in updates as the Schur element
	double tau1bord;  // tau1bord is the first  order factorization threshold in the bordering
	double tau2bord;  // tau2bord is the second order factorization threshold in the bordering
	double theta;     // theta is the diagonal correction parameter
	double sclmin;    // sclmin specifies the minimal allowable singular value when scaling
	double pivmin;    // pivmin specifies the minimal allowable pivot singular value when performing incomplete factorization
	int ichkfct;      // ichkfct specifies the number of supernode row per histogram output in fct
	double omega;     // omega is the parameter of the SOR algorithm
	double permismin; // permismin specifies the minimal allowable value for the solution vector component
	int jobitr;       // jobitr is the job to be performed
	int oocitr;       // oocitr is the out-of-core type of the iterative scheme
	int nrestarts;    // nrestarts is the maximal number of restarts
	int ncorrmax;     // ncorrmax is the maximal number of correctors
	int nrhsrestart;  // nrhsrestart is the maximal number of right hand sides per restart
	int ittype;       // ittype is the type iterative scheme used
	int sttype;       // sttype is the type of the stopping criterion: 0 - rel., 1 - abs., 2 - norm vector
	int niter;        // niter is the maximal number of iterations
	int niter2;       // niter2 is the maximal number of outer iterations
	int niterlocal;   // niterlocal is the local number of iterations
	int nblkmax;      // nblkmax is the maximal number of blocks used to store data
	int ichk;         // ichk specifies the number of iterations per residual checks in the iterative scheme
	int ncolbuf;      // ncolbuf is the number of vectors used as exchange buffer in the iterative scheme
	int nsteps;       // nsteps is the number of nonlinear shifts used to control nonlinear convergence
	double thrfltmtr; // thrfltmtr is the matrix filtering parameter (<0 - no filtering)
	double thrfltrhs; // thrfltrhs is the right hand side filtering parameter (<0 - no filtering)
	double eps;       // eps is the stopping criterion
	double eps2;      // eps2 is the outer stopping criterion
	double valuemin;  // valuemin parameter specifies the minimal allowable value for each component of the nonlinear solver
	double xmin;      // xmin is the threshold value of the iterative scheme
	double xmax;      // xmax is the threshold value of the iterative scheme
	char path[512];   // path is the path to the out-of-core data
	double density;   // density is the preconditioner density
	double density2;  // density2 is the preconditioner density
	int nitperf;      // nitperf is the number of performed iterations
	int nitperf2;     // nitperf2 is the number of performed secondary iterations
	double resfin;    // resfin is the final residual
	int nfreer;       // nfreer is the number of R part freeings per incomplete factorization
	double timefct;   // timefct is the factorization time
	double timeitr;   // timeitr is the iterative time
	double timefctout;// timefctout is the factorization output time
	int dirslvtype;   // dirslvtype is the type of the direct problem solver
	int iuseteupl;    // iuseteupl is the parameter that specifies is it necessary to use teuplitz structure or not
	int ngroups;      // ngroups is the number of groups on which all processors are divided
	int nlev;         // nlev is the number of levels in the 3D partitioning of the domain
	int collapini;    // collapini  is the initial 3D condensation of the elements
	int collapxini;   // collapxini  is the initial X condensation of the elements
	int collapyini;   // collapyini  is the initial Y condensation of the elements
	int collapzini;   // collapzini  is the initial Z condensation of the elements
	int collapmult;   // collapmult is the 3D condensation of the elements per each level
	int npt0;         // npt0 is one directional number of intergration points in domain 0 (closest region)
	int npt1;         // npt1 is one directional number of intergration points in domain 1
	int npt2;         // npt2 is one directional number of intergration points in domain 2
	double rmin;      // rmin is the radius that describes the closest region
	double rmax;      // rmax is the radius that describes the other regions
	int ndiag;        // ndiag is the number of diagonals in the sparsified matrix
	int nneib;        // nneib is the number of neibouring nodes in the sparsified matrix
	int ndiagsave;    // ndiagsave is the saved number of diagonals in the sparsified matrix for future use
	int nneibsave;    // nneibsave is the saved number of neibouring nodes in the sparsified matrix for future use
	int cmptype;      // cmptype is the compression type
	double epscmp;    // epscmp is the compression accuracy
#ifdef __FlowVisionOut
	CRegion* pregion;   // FlowVision parameter pregion
	sMethodSet*  meth;  // FlowVision parameter meth
	int numpar;         // FlowVision parameter numPar
	sMethodSet** methj; // FlowVision parameter methj
#endif
// Functions
	CSlvParam () { // Constructor
		msglev      = 0;
		nchilds     = 3;
		ncpupart    = 1;
		ncpuext     = 1;
		ncpuextsch  = 1;
		decomptype  = 0;
		schurtype   = 0;
		nsupmax     = 100;
		weightcellvolume   = 1000;
		weightcellboundary = 1000;
		fcttyp      = 1;
		ooctyp      = 0;
		ordtype     = 1;
		ordtypeblk  = 0;
//		collap      = 5;
		collap      = 0;
		ncycle      = 0;
		nprecini    = 1;
		nprec       = 1;
		sclpttype   = 0;
		alpha       = 0.0e0;
		dshift      = 0.0e0;
		dshiftim    = 0.0e0;
		memory      = 1.0e0;
		tau1        = 1.0e-2;
		tau2        = 1.0e-4;
		tau1save    = 1.0e-2;
		tau2save    = 1.0e-4;
		ndiagfct    = 20;
		ndiagfct2   = 100;
		nlevfct     = 200000;
		nlevfct2    = 1000000;
		nlevschur   = -1;
		tau1bord    = 1.0e-2;
		tau2bord    = 1.0e-4;
		theta       = 1.0e-1;
		sclmin      = 1.0e-8;
		pivmin      = 1.0e-8;
		ichkfct     = 50000;
		omega       = 0.9e0;
		permismin   = 1.0e-50;
		jobitr      = 0;
		oocitr      = 0;
		nrestarts   = -1;
		ncorrmax    = 0;
		nrhsrestart = 1;
		ittype      = 2;
		sttype      = 2;
		niter       = 1000;
		niter2      = 10;
		niterlocal  = 20;
		nblkmax     = 50;
		ichk        = 5;
		ncolbuf     = 10;
		nsteps      = 5;
		thrfltmtr   = 1.0e-2;
		thrfltrhs   = -1.0e0;
		eps         = 1.0e-9;
		eps2        = 1.0e-9;
		valuemin    = -1.0e100;
		xmin        = 0.2e0;
		xmax        = 2.0e0;
		strcpy (path,"");
		density     = 0.0e0;
		density2    = 0.0e0;
		nitperf     = 0;
		nitperf2    = 0;
		resfin      = 1.0e10;
		nfreer      = 0;
		timefct     = 0.0e0;
		timeitr     = 0.0e0;
		timefctout  = 10.0e0;
		dirslvtype  = 0;
		iuseteupl   = 0;
		ngroups     = 1;
		nlev        = 3;
		collapini   = 3;
		collapxini  = 3;
		collapyini  = 3;
		collapzini  = 3;
		collapmult  = 2;
		npt0        = 3;
		npt1        = 2;
		npt2        = 1;
		rmin        = 2.0e0;
		rmax        = 10.0e0;
		ndiag       = 100;
		nneib       = 1;
		ndiagsave   = 100;
		nneibsave   = 1;
		cmptype     = 0;
		epscmp      = 1.0e-8;
	};
	void SetShift (double _dshift, double _dshiftim) { // Set shift parameters
		dshift    = _dshift;
		dshiftim  = _dshiftim;
	};
	void SetFct (int _fcttyp, double _memory, double _tau1, double _tau2, double _theta) { // Set factorization parameters
		fcttyp    = _fcttyp;
		memory    = _memory;
		tau1      = _tau1;
		tau2      = _tau2;
		theta     = _theta;
	};
	void SetFctBord (double _tau1bord, double _tau2bord) { // Set factorization parameters for bordering
		fcttyp   = 2;
		tau1bord = _tau1bord;
		tau2bord = _tau2bord;
	};
	void SetFctNDiag (int _ndiagfct, int _ndiagfct2) { // Set factorization parameters for diagonal blocks
		fcttyp   = 3;
		ndiagfct  = _ndiagfct;
		ndiagfct2 = _ndiagfct2;
	};
	void SetFctNLev (int _nlevfct, int _nlevfct2) { // Set factorization parameters
		fcttyp   = 5;
		nlevfct  = _nlevfct;
		nlevfct2 = _nlevfct2;
	};
	void SetOrd (int _ordtype) { // Set the ordering type
		ordtype   = _ordtype;
	};
	void SetCollap (int _collap) { // Set the ordering type
		collap    = _collap;
	};
	void SetPiv (double _sclmin, double _pivmin) { // Set pivoting factorization parameters
		sclmin    = _sclmin;
		pivmin    = _pivmin;
	};
	void SetSor (double _omega, double _permismin) { // Set SOR parameters
		omega     = _omega;
		permismin = _permismin;
	};
	void SetFltMtr (double _thrfltmtr) { // Set matrix filtering parameter
		thrfltmtr = _thrfltmtr;
	};
	void SetFltRhs (double _thrfltrhs) { // Set right hand side filtering parameter
		thrfltrhs = _thrfltrhs;
	};
	void SetIter (int _ittype, int _niter, double _eps) { // Set iterative scheme parameters
		ittype    = _ittype;
		niter     = _niter;
		eps       = _eps;
	};
	void SetIchkFct (int _ichkfct) { // Set ichkfct
		ichkfct   = _ichkfct;
	};
	void SetIchk (int _ichk) { // Set ichk
		ichk      = _ichk;
	};
	void SetNcolbuf (int _ncolbuf) { // Set ncolbuf
		ncolbuf   = _ncolbuf;
	};
	void SetOoc (int _nfiles, char *_path) { // Set out-of-core control data
		ooctyp    = _nfiles;
		strcpy (path,_path);
	};
	void SetNcycle (int _ncycle) { // Set ncycle
		ncycle    = _ncycle;
	};
	void SetNewton (int _niter2, double _eps2, int _nprec, int _nsteps) { // Set Newton parameters
		niter2    = _niter2;
		eps2      = _eps2;
		nprec     = _nprec;
		nsteps    = _nsteps;
	};
	void SetCorr (double _xmin, double _xmax) { // Set corrector parameters
		xmin      = _xmin;
		xmax      = _xmax;
	};
	void SetJobItr (int _jobitr) { // Set iterative scheme job
		jobitr    = _jobitr;
	};
	void SetRestarts (int _nrestarts, int _ncorrmax, int _nrhsrestart) { // Set restarts parameters
		nrestarts   = _nrestarts;
		ncorrmax    = _ncorrmax;
		nrhsrestart = _nrhsrestart;
	};
	void SetPtScale (int _sclpttype) { // Set sclpttype parameters
		sclpttype   = _sclpttype;
	};
	void SetOocItr (int _oocitr) { // Set out-of-core type of the iterative scheme
		oocitr    = _oocitr;
	};
	friend class CSVector;
	friend class CSVectorC;
	friend class CSMatrix;
	friend class CSMatrixR;
	friend class CSMatrixRB;
	friend class CSMatrixRS;
	friend class CSMatrixC;
	friend class CSMatrixCS;
	friend class CGSMatrixCS;
	friend class CFctC;
};

#endif
