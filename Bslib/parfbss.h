
//------------------------------------------------------------------------------------------------
// File: parfbss.h
// Creator: Kharchenko S.A. 
// Goal: Stores system of equations data and solves system of linear equations in parallel mode
//------------------------------------------------------------------------------------------------

#ifndef __PARFBSS_SOLVER__
#define __PARFBSS_SOLVER__

int ParFBSS_Create (void *_Comm, int _NI, int _NZAI, int *_ObjId);
int ParFBSS_ReCreate (int _ObjId, int _NI, int _NZAI);
int ParFBSS_LinearEquationsData (int _ObjId, 
											int **_LIST, int **_IA, int **_JA, double **_A, 
											double **_RightHandSide, double **_InitialGuess);
int ParFBSS_NormalizationData (int _ObjId, double **_Norm);
int ParFBSS_NonLinearEquationsData (int _ObjId, double **_Poly2);
int ParFBSS_SetSolverParameters (int _ObjId, 
											int _NpDistr,
											double _Tau1, double _Tau2, double _Theta, 
											int _NIterMax, double _Eps);
int ParFBSS_SetMsglevAndFilename (int _ObjId, int _Msglev, char *_FileName);
int ParFBSS_SetSchurtype (int _ObjId, int _Schurtype, int _Ncycle);
int ParFBSS_SetNcpuext (int _ObjId, int _NcpuExt);
int ParFBSS_SolveSymm (int _ObjId, int _Job, int *_NIterPerf);
int ParFBSS_Solve (int _ObjId, int _Job, int *_NIterPerf);
int ParFBSS_Solution (int _ObjId, double **_Solution);
int ParFBSS_Distribution (int _ObjId, int *_NIOpt, int **_LISTOpt);
int ParFBSS_Destroy (int _ObjId);

#endif
