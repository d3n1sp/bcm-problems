#include "stdafx.h"
#include "unit_mes.h"

///////////////////////////////////////////////////////////////////////////////////////////////////
//...ñïèñîê ñèñòåì åäèíèö èçìåðåíèÿ âåëè÷èí
  int  GetUnitsCount  (void) { return(NUM_UNITS);}
  int  GetUnitsDefault(void) { return(UNIT_SI);  }

//////////////////////////////////////////////////////////////////////////////////////////
//...ñïèñêè íàçâàíèé óñòàíîâëåííîé åäèíèöû è íîðìèðóþùèé ìíîæèòåëü äëÿ âûâîäèìîé âåëè÷èíû;
const char * GetUnitsNameR(int N_ui)
{ const char *L[]={ "Ñèñòåìà SI",
                     "Ìîäèôèöèðîâàííàÿ SI",
                     "Cèñòåìà CGS", 
                     "Ñïåöèàëüíàÿ òåõíè÷åñêàÿ",
                     "Àíãëèéñêàÿ ñèñòåìà IPS",
                     "Àíãëèéñêàÿ ñèñòåìà FPS"};
  return L[N_ui];
}

const char * GetAngleNameR(int N_ui)
{ const char *L[]={ "ðàä",
                     "ðàä",
                     "ðàä",
                     "ãðàä",
                     "ãðàä", 
                     "ãðàä"};
  return L[N_ui];
}

const char * GetLengthNameR(int N_ui)
{ const char *L[]={ "ì",
                     "ìì",
                     "ñì", 
                     "ìì",
                     "äþéì",
                     "ôóò"};
  return L[N_ui];
}

const char * GetTimeNameR(int N_ui)
{ const char *L[]={ "ñåê",
                     "ñåê",
                     "ñåê", 
                     "ñåê",
                     "ñåê",
                     "ñåê"};
  return L[N_ui];
}

const char * GetVelocityNameR(int N_ui)
{ const char *L[]={ "ì/ñåê",
                     "ìì/ñåê",
                     "ñì/ñåê", 
                     "ìì/ñåê",
                     "äþéì/ñåê",
                     "ôóò/ñåê"};
  return L[N_ui];
}

const char * GetMassNameR(int N_ui)
{ const char *L[]={ "êã",
                     "êã",
                     "êã",
                     "êã",
                     "ôñ*ñåê^2/ôóò",
                     "ôñ*ñåê^2/ôóò"};
  return L[N_ui];
}

const char * GetTemperatureNameR(int N_ui)
{ const char *L[]={ "K = C-273.15",
                     "K = C-273.15",
                     "C",
                     "C",
                     "F = (C-32)*5/9",
                     "F = (C-32)*5/9"};
  return L[N_ui];
}
  
const char * GetLineWeightNameR(int N_ui)
{ const char *L[]={ "Í/ì",
                     "ìÍ/ìì",
                     "êãc/ñì",
                     "êãc/êì",
                     "ôñ/ôóò (ôóíò-ñèëà/ôóò)",
                     "ôñ/ôóò (ôóíò-ñèëà/ôóò)"};
  return L[N_ui];
}

const char * GetMassDensityNameR(int N_ui)
{ const char *L[]={ "êã/ì^3",
                     "êã/ìì^3",
                     "ãðàì/ñì^3",
                     "êã/ìì^3",
                     "ôñ*ñåê^2/äþéì^4",
                     "ôñ*ñåê^2/ôóò^4"};
  return L[N_ui];
}

const char * GetStressNameR(int N_ui)
{ const char *L[]={ "Ïà = Í/ì^2",
                     "êÏà = ìÍ/ìì^2",
                     "êãñ/ñì^2",
                     "ÌÏà",
                     "ôñ/äþéì^2 (ôóíò-ñèëà/äþéì^2)",
                     "ôñ/ôóò^2 (ôóíò-ñèëà/ôóò^2)"};
  return L[N_ui];
}

const char * GetStrainNameR(int N_ui)
{ const char *L[]={ "",
                     "",
                     "",
                     "ìèêðî-ñòðåéí",
                     "",
                     ""};
  return L[N_ui];
}

const char * GetPointForceNameR(int N_ui)
{ const char *L[]={ "Í = êã*ì/ñåê^2",
                     "ìÍ = êã*ìì/ñåê^2",
                     "êãñ = êã*981ñì/ñåê^2",
                     "Í = êã*ì/ñåê^2",
                     "ôñ = ôóíò-ñèëà",
                     "ôñ = ôóíò-ñèëà"};
  return L[N_ui];
}

const char * GetLineForceNameR(int N_ui)
{ const char *L[]={ "Í/ì",
                     "ìÍ/ìì",
                     "êãñ/ñì",
                     "Í/ìì",
                     "ôñ/äþéì (ôóíò-ñèëà/äþéì)",
                     "ôñ/ôóò (ôóíò-ñèëà/ôóò)"};
  return L[N_ui];
}

const char * GetPointForceMomentNameR(int N_ui)
{ const char *L[]={ "Í*ì",
                     "ìÍ*ìì",
                     "êãñ*ñì",
                     "Í*ìì",
                     "ôñ*äþéì",
                     "ôñ*ôóò"};
  return L[N_ui];
}

const char * GetBeamStiffnessNameR(int N_ui)
{ const char *L[]={ "Í*ì^2",
                     "ìÍ*ìì^2",
                     "êãñ*ñì^2",
                     "Í*ìì^2",
                     "ôñ*äþéì^2",
                     "ôñ*ôóò^2"};
  return L[N_ui];
}
  
const char * GetUnitsNameA(int N_ui)
{ const char *L[]={ "SI",
                     "Modified SI",
                     "Russian technical", 
                     "Wire special",
                     "English IPS",
                     "English FPS"};
  return L[N_ui];
}

const char * GetAngleNameA(int N_ui)
{ const char *L[]={ "rad",
                     "rad",
                     "rad",
                     "deg",
                     "deg",
                     "deg"};
  return L[N_ui];
}
  
const char * GetLengthNameA(int N_ui)
{ const char *L[]={ "m",
                     "mm",
                     "cm",
                     "mm",
                     "in",
                     "ft"};
  return L[N_ui];
}

const char * GetTimeNameA(int N_ui)
{ const char *L[]={ "sec",
                     "sec",
                     "sec", 
                     "sec",
                     "sec",
                     "sec"};
  return L[N_ui];
}

const char * GetVelocityNameA(int N_ui)
{ const char *L[]={ "m/sec",
                     "mm/sec",
                     "cm/sec", 
                     "mm/sec",
                     "in/sec",
                     "ft/sec"};
  return L[N_ui];
}
  
const char * GetMassNameA(int N_ui)
{ const char *L[]={ "kg",
                     "kg",
                     "kg",
                     "kg",
                     "lbf*sec^2/ft",
                     "lbf*sec^2/ft"};
  return L[N_ui];
}

const char * GetTemperatureNameA(int N_ui)
{ const char *L[]={ "K = C-273.15",
                     "K = C-273.15",
                     "C",
                     "C",
                     "F = (C-32)*5/9",
                     "F = (C-32)*5/9"};
  return L[N_ui];
}
  
const char * GetLineWeightNameA(int N_ui)
{ const char *L[]={ "N/m",
                     "mN/mm",
                     "kgf/cm",
                     "kgf/km",
                     "lbf*sec^2/ft^2",
                     "lbf*sec^2/ft^2"};
  return L[N_ui];
}

const char * GetMassDensityNameA(int N_ui)
{ const char *L[]={ "kg/m^3",
                     "kg/mm^3",
                     "g/cm^3",
                     "kg/mm^3",
                     "lbf*sec^2/in^4",
                     "lbf*sec^2/ft^4"};
  return L[N_ui];
}

const char * GetStressNameA(int N_ui)
{ const char *L[]={ "Pa = N/m^2",
                     "kPa = mN/mm^2",
                     "kgf/cm^2",
                     "MPa",
                     "lbf/in^2",
                     "lbf/ft^2"};
  return L[N_ui];
}

const char * GetStrainNameA(int N_ui)
{ const char *L[]={ "",
                     "",
                     "",
                     "micro-strain",
                     "",
                     ""};
  return L[N_ui];
}

const char * GetPointForceNameA(int N_ui)
{ const char *L[]={ "N = kg*m/sec^2",
                     "mN = kg*mm/sec^2",
                     "kgf = kg*981cm/sec^2",
                     "N",
                     "lbf",
                     "lbf"};
  return L[N_ui];
}

const char * GetLineForceNameA(int N_ui)
{ const char *L[]={ "N/m",
                     "mN/mm",
                     "kgf/cm",
                     "N/mm",
                     "lbf/in",
                     "lbf/ft"};
  return L[N_ui];
}

const char * GetPointForceMomentNameA(int N_ui)
{ const char *L[]={ "N*m",
                     "mN*mm",
                     "kgf*cm",
                     "N*mm",
                     "lbf*in",
                     "lbf*ft"};
  return L[N_ui];
}

const char * GetBeamStiffnessNameA(int N_ui)
{ const char *L[]={ "N*m^2",
                     "mN*mm^2",
                     "kgf*cm^2",
                     "N*mm^2",
                     "lbf*in^2",
                     "lbf*ft^2"};
  return L[N_ui];
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//...metric coefficients for translating one system of units to enother
double GetMetricAngle(int units_to, int units_from)
{
  if (units_to < 0 || units_to >= NUM_UNITS
                   || units_from < 0 || units_from >= NUM_UNITS) return(1.);
  double metric[NUM_UNITS][NUM_UNITS] = {
         {1.       ,1.       ,1.       ,M_PI/180.,M_PI/180.,M_PI/180.},
         {1.       ,1.       ,1.       ,M_PI/180.,M_PI/180.,M_PI/180.},
         {1.       ,1.       ,1.       ,M_PI/180.,M_PI/180.,M_PI/180.},
         {180./M_PI,180./M_PI,180./M_PI,1.       ,1.       ,1.       },
         {180./M_PI,180./M_PI,180./M_PI,1.       ,1.       ,1.       },
         {180./M_PI,180./M_PI,180./M_PI,1.       ,1.       ,1.       },
  };
  return metric[units_to][units_from];
}

double GetMetricLength(int units_to, int units_from)
{
  if (units_to < 0 || units_to >= NUM_UNITS
                   || units_from < 0 || units_from >= NUM_UNITS) return(1.);
  double metric[NUM_UNITS][NUM_UNITS] = {
         {1.       ,1e-3    ,1e-2    ,1e-3    ,0.0254,0.3048},
         {1e3      ,1.      ,10.     ,1.      ,25.4  ,304.8 },
         {1e2      ,1e-1    ,1.      ,1e-1    ,2.54  ,30.48 },
         {1e3      ,1.      ,10.     ,1.      ,25.4  ,304.8 },
         {1./0.0254,1./25.4 ,1./2.54 ,1./25.4 ,1.    ,12.   },
         {1./0.3048,1./304.8,1./30.48,1./304.8,1./12.,1.    },
  };
  return metric[units_to][units_from];
}

double GetMetricPointForce(int units_to, int units_from)
{
  if (units_to < 0 || units_to >= NUM_UNITS
                   || units_from < 0 || units_from >= NUM_UNITS) return(1.);
  double metric[NUM_UNITS][NUM_UNITS] = {
         {1.                  ,1e-3               ,9.81         ,1.                  ,9.81*0.45359237,9.81*0.45359237},
         {1e3                 ,1.                 ,9810.        ,1e3                 ,981.*4.5359237 ,981.*4.5359237 },
         {1./9.81             ,1./9810.           ,1.           ,1./9.81             ,0.45359237     ,0.45359237     },
         {1.                  ,1e-3               ,9.81         ,1.                  ,9.81*0.45359237,9.81*0.45359237},
         {10./(9.81*4.5359237),1./(981.*4.5359237),1./0.45359237,10./(9.81*4.5359237),1.             ,1.             },
         {10./(9.81*4.5359237),1./(981.*4.5359237),1./0.45359237,10./(9.81*4.5359237),1.             ,1.             },
  };
  return metric[units_to][units_from];
}

double GetMetricDensity(int units_to, int units_from)
{
  if (units_to < 0 || units_to >= NUM_UNITS
                   || units_from < 0 || units_from >= NUM_UNITS) return(1.);
  double metric[NUM_UNITS][NUM_UNITS] = {
         {1.                ,1e9                ,1e3               ,1e9                ,10695515.18   ,10695515.18/20736.   },
         {1e-9              ,1.                 ,1e-6              ,1.                 ,1.069551518e-2,1.069551518e-2/20736.},
         {1e-3              ,1e6                ,1.                ,1e6                ,10695.51518   ,10695.51518/20736.   },
         {1e-9              ,1.                 ,1e-6              ,1.                 ,1.069551518e-2,1.069551518e-2/20736.},
         {1./10695515.18    ,1e2/1.069551518    ,1./10695.51518    ,1e2/1.069551518    ,1.            ,1./20736.            },
         {20736./10695515.18,20736e2/1.069551518,20736./10695.51518,20736e2/1.069551518,20736.        ,1.                   },
  };
  return metric[units_to][units_from];
}

double GetMetricStrain(int units_to, int units_from)
{
  if (units_to < 0 || units_to >= NUM_UNITS
                   || units_from < 0 || units_from >= NUM_UNITS) return(1.);
  double metric[NUM_UNITS][NUM_UNITS] = {
         {1. ,1. ,1. ,1e-6,1. ,1. },
         {1. ,1. ,1. ,1e-6,1. ,1. },
         {1. ,1. ,1. ,1e-6,1. ,1. },
         {1e6,1e6,1e6,1.  ,1e6,1e6},
         {1. ,1. ,1. ,1e-6,1. ,1. },
         {1. ,1. ,1. ,1e-6,1. ,1. },
  };
  return metric[units_to][units_from];
}

double GetMetricTemperature(int units_to, int units_from)
{
  if (units_to < 0 || units_to >= NUM_UNITS
                   || units_from < 0 || units_from >= NUM_UNITS) return(1.);
  double metric[NUM_UNITS][NUM_UNITS] = {
         {1.   ,1.   ,1.   ,1.   ,1.8,1.8},
         {1.   ,1.   ,1.   ,1.   ,1.8,1.8},
         {1.   ,1.   ,1.   ,1.   ,1.8,1.8},
         {1.   ,1.   ,1.   ,1.   ,1.8,1.8},
         {5./9.,5./9.,5./9.,5./9.,1. ,1. },
         {5./9.,5./9.,5./9.,5./9.,1. ,1. },
  };
  return metric[units_to][units_from];
}

double GetMetricShiftTemperature(int units_to, int units_from)
{
  if (units_to < 0 || units_to >= NUM_UNITS
                   || units_from < 0 || units_from >= NUM_UNITS) return(1.);
  double metric[NUM_UNITS][NUM_UNITS] = {
         {0.    ,0.    ,-273.15,-273.15,-241.15 ,-241.15 },
         {0.    ,0.    ,-273.15,-273.15,-241.15 ,-241.15 },
         {273.15,273.15,0.     ,0.     ,-160./9.,-160./9.},
         {273.15,273.15,0.     ,0.     ,32.     ,32.     },
         {305.15,305.15,32.    ,32.    ,0.      ,0.      },
         {305.15,305.15,32.    ,32.    ,0.      ,0.      },
  };
  return metric[units_to][units_from];
}
double GetMetricLineForce(int units_to, int units_from)
{
  return GetMetricPointForce(units_to, units_from)*GetMetricLength(units_from, units_to);
}
double GetMetricPointForceMoment(int units_to, int units_from)
{
  return GetMetricPointForce(units_to, units_from)*GetMetricLength(units_to, units_from);
}
double GetMetricLineForceMoment(int units_to, int units_from)
{
  return GetMetricPointForce(units_to, units_from);
}
double GetMetricStress(int units_to, int units_from)
{
  return GetMetricLineForce(units_to, units_from)*GetMetricLength(units_from, units_to);
}
double GetMetricStiffness(int units_to, int units_from)
{
  return GetMetricPointForceMoment(units_to, units_from)*GetMetricLength(units_to, units_from);
}
double GetMetricLineWeight(int units_to, int units_from)
{
  double metric = GetMetricLineForce(units_to, units_from), f = 1.;
  if (units_to == UNIT_WIRE) f = GetMetricPointForce(UNIT_CGS, UNIT_SI)*GetMetricLength(UNIT_MSI, UNIT_SI)*1e-3;
  else                       f = GetMetricPointForce(UNIT_SI, UNIT_CGS)*GetMetricLength(UNIT_SI, UNIT_MSI)*1e3;
  return metric*f;
}


