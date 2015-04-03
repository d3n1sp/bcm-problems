/*============================================*/
/*                  UNIT_MES                  */
/*============================================*/
#ifndef ___UNIT_MES___
#define ___UNIT_MES___

#include "utils.h"

////////////////////////////////////////////////////
//...абсолютна€ нумераци€ систем единиц в программе;
enum Units_Type {
     UNIT_SI, UNIT_MSI, UNIT_CGS, UNIT_WIRE, UNIT_IPS, UNIT_FPS,
     NUM_UNITS
};

///////////////////////////////////////////////////////////////////////////////////////////////////
//                     ќрганизаци€ диалога по выбору системы единиц.                             //
///////////////////////////////////////////////////////////////////////////////////////////////////
  int  GetUnitsCount  (void);
  int  GetUnitsDefault(void);
  const char * GetUnitsNameR           (int N_ui);
  const char * GetAngleNameR           (int N_ui);
  const char * GetLengthNameR          (int N_ui);
  const char * GetTimeNameR            (int N_ui);
  const char * GetVelocityNameR        (int N_ui);
  const char * GetMassNameR            (int N_ui);
  const char * GetTemperatureNameR     (int N_ui);
  const char * GetLineWeightNameR      (int N_ui);
  const char * GetMassDensityNameR     (int N_ui);
  const char * GetStressNameR          (int N_ui);
  const char * GetStrainNameR          (int N_ui);
  const char * GetPointForceNameR      (int N_ui);
  const char * GetLineForceNameR       (int N_ui);
  const char * GetLineForceMomentNameR (int N_ui);
  const char * GetBeamStiffnessNameR   (int N_ui);
  
  const char * GetUnitsNameA           (int N_ui);
  const char * GetAngleNameA           (int N_ui);
  const char * GetLengthNameA          (int N_ui);
  const char * GetTimeNameA            (int N_ui);
  const char * GetVelocityNameA        (int N_ui);
  const char * GetMassNameA            (int N_ui);
  const char * GetTemperatureNameA     (int N_ui);
  const char * GetLineWeightNameA      (int N_ui);
  const char * GetMassDensityNameA     (int N_ui);
  const char * GetStressNameA          (int N_ui);
  const char * GetStrainNameA          (int N_ui);
  const char * GetPointForceNameA      (int N_ui);
  const char * GetLineForceNameA       (int N_ui);
  const char * GetPointForceMomentNameA(int N_ui);
  const char * GetBeamStiffnessNameA   (int N_ui);

  double GetMetricAngle           (int units_to, int units_from = UNIT_SI);
  double GetMetricLength          (int units_to, int units_from = UNIT_SI);
  double GetMetricPointForce      (int units_to, int units_from = UNIT_SI);
  double GetMetricDensity         (int units_to, int units_from = UNIT_SI);
  double GetMetricStrain          (int units_to, int units_from = UNIT_SI);
  double GetMetricTemperature     (int units_to, int units_from = UNIT_SI);
  double GetMetricShiftTemperature(int units_to, int units_from = UNIT_SI);

  double GetMetricLineForce       (int units_to, int units_from = UNIT_SI);
  double GetMetricPointForceMoment(int units_to, int units_from = UNIT_SI);
  double GetMetricLineForceMoment (int units_to, int units_from = UNIT_SI);
  double GetMetricStress          (int units_to, int units_from = UNIT_SI);
  double GetMetricStiffness       (int units_to, int units_from = UNIT_SI);
  double GetMetricLineWeight      (int units_to, int units_from = UNIT_SI);

#endif
