/*===========================================*/
/*                  PROFILE                  */
/*===========================================*/
#ifndef ___PROFILE___
#define ___PROFILE___

extern const char * COLONSTR;
extern const char * SEMICOLONSTR;
extern const char * STRSEPARATOR;
extern const char * PARSEPARATOR;
extern const char * TABLESTR;
extern const char * VERSIONSTR;

////////////////////////////////////////////////////////
//...включение урезанного режима работы контекста задач;
#define ___ABRIDGE_PROFILE_MODE___
#ifndef ___ABRIDGE_PROFILE_MODE___
#include "csample.h"
#else
#include "cbase.h"
#include "kernel.h"
#endif
#define SHIFT_LM3D_SAMPLES  0
enum Samples_Lame3D {
    _LM1, _LM2, _LM3,
    _NUM_LM3D_SAMPLES
};
#define NUM_LM3D_SAMPLES _NUM_LM3D_SAMPLES
#define SHIFT_MX3D_SAMPLES (SHIFT_LM3D_SAMPLES+NUM_LM3D_SAMPLES)
enum Samples_Maxw3D {
    _MX1, _MX2, _MX3,
    _NUM_MX3D_SAMPLES
};
#define NUM_MX3D_SAMPLES _NUM_MX3D_SAMPLES
#define SHIFT_AU3D_SAMPLES (SHIFT_MX3D_SAMPLES+NUM_MX3D_SAMPLES)
enum Samples_Acou3D {
    _AU1, _AU2, _AU3, 
    _NUM_AU3D_SAMPLES
};
#define NUM_AU3D_SAMPLES _NUM_AU3D_SAMPLES
#define SHIFT_MP3D_SAMPLES (SHIFT_AU3D_SAMPLES+NUM_AU3D_SAMPLES)
enum Samples_Mapi3D {
    _MP1, _MP2, _MP3,
    _NUM_MP3D_SAMPLES
};
#define NUM_MP3D_SAMPLES _NUM_MP3D_SAMPLES
#define SHIFT_AU2D_SAMPLES (SHIFT_MP3D_SAMPLES+NUM_MP3D_SAMPLES)
enum Samples_Acou2D {
    _AU4, _AU5, _AU6, _AU7, _AU8, _AU9,
    _NUM_AU2D_SAMPLES
};
#define NUM_AU2D_SAMPLES _NUM_AU2D_SAMPLES
#define SHIFT_SK2D_SAMPLES (SHIFT_AU2D_SAMPLES+NUM_AU2D_SAMPLES)
enum Samples_Skin2D {
    _SK1, _SK2, 
    _NUM_SK2D_SAMPLES
};
#define NUM_SK2D_SAMPLES _NUM_SK2D_SAMPLES
#define SHIFT_MP2D_SAMPLES (SHIFT_SK2D_SAMPLES+NUM_SK2D_SAMPLES)
enum Samples_Mapi2D {
    _MP10, _MP11, _MP12, _MP13, _MP14, _MP15, _MP16, _MP17, _MP18, _MP19, _MP20,//35!!!
    _MP21, _MP22, _MP23, _MP24, _MP25, _MP26, _MP27, _MP28, _MP29, _MP30, _MP31, _MP32, _MP33,
    _NUM_MP2D_SAMPLES
};
#define NUM_MP2D_SAMPLES _NUM_MP2D_SAMPLES
#define SHIFT_BARS_GOST_SAMPLES (SHIFT_MP2D_SAMPLES+NUM_MP2D_SAMPLES)
enum Samples_BarsGost2D {
    _LU_GOST, _LN_GOST, _IU_GOST, _SU_GOST, _LGUGOST, _LGNGOST, _SGUGOST, _RC_GOST, 
	 _NUM_BARS_GOST_SAMPLES
};
#define NUM_BARS_GOST_SAMPLES _NUM_BARS_GOST_SAMPLES-4
#define SHIFT_BARS_SAMPLES (SHIFT_BARS_GOST_SAMPLES+NUM_BARS_GOST_SAMPLES)
enum Samples_Bars2D {
    _BARS1,  _BARS2,  _BARS3,  _BARS4,  _BARS5,  _BARS6,  _BARS7,  _BARS8,  _BARS9,  _BARS10, 
	 _BARS11, _BARS12, _BARS13, _BARS14, _BARS15, _BARS16, _BARS17, _BARS18, _BARS19, _BARS20, 
	 _BARS21, _BARS22, _BARS23, _BARS24, _BARS25, _BARS26, _BARS27, _BARS28, _BARS29, _BARS30, 
	 _BARS31, _BARS32, _BARS33, _BARS34, _BARS35, _BARS36, _BARS37, _BARS38, _BARS39, _BARS40, 
	 _BARS41, _BARS42, _BARS43, _BARS44, _BARS45, _BARS46, _BARS47, _BARS48, _BARS49, _BARS50, 
	 _BARS51, _BARS52, _BARS53, _CHAIN0,
	 _NUM_BARS_SAMPLES
};
#define NUM_BARS_SAMPLES _NUM_BARS_SAMPLES
#define SHIFT_WIRE_SAMPLES (SHIFT_BARS_SAMPLES+NUM_BARS_SAMPLES)
enum Samples_Wire2D {
    _WIRE1, _WIRE2, _WIRE3, _WIRE4, _WIRE5, _WIRE6, _WIRE7, _WIRE8, _WIRE9, _WIRE10,
    _WIRE11,_WIRE12,_WIRE13,_WIRE14,_WIRE15,_WIRE16,_WIRE17,
    _NUM_WIRE_SAMPLES
};
#define NUM_WIRE_SAMPLES _NUM_WIRE_SAMPLES
#define SHIFT_CABLE_SAMPLES (SHIFT_WIRE_SAMPLES+NUM_WIRE_SAMPLES)
enum Samples_Cable {
    _CABLE1, _CABLE2, _CABLE3,
    _NUM_CABLE_SAMPLES
};
#define NUM_CABLE_SAMPLES _NUM_CABLE_SAMPLES
#define SHIFT_REDUCEDWIRE_SAMPLES (SHIFT_CABLE_SAMPLES+NUM_CABLE_SAMPLES)
enum Samples_ReducedWire2D {
    _REDUCEDWIRE1, _REDUCEDWIRE2, _REDUCEDWIRE3, _REDUCEDWIRE4, _REDUCEDWIRE5, _REDUCEDWIRE6,
    _NUM_REDUCEDWIRE_SAMPLES
};
#define NUM_REDUCEDWIRE_SAMPLES _NUM_REDUCEDWIRE_SAMPLES-3
#define SHIFT_TOWER_SAMPLES (SHIFT_REDUCEDWIRE_SAMPLES+NUM_REDUCEDWIRE_SAMPLES)
enum Samples_Tower {
    _TOWER1, _TOWER2, _TOWER3, _TOWER4, _TOWER5, _TOWER6,
    _NUM_TOWER_SAMPLES
};
#define NUM_TOWER_SAMPLES _NUM_TOWER_SAMPLES-1
#define SHIFT_FITTING_SAMPLES (SHIFT_TOWER_SAMPLES+NUM_TOWER_SAMPLES)
enum Samples_Fitting {
    _FITTING1, _FITTING2, _FITTING3, _FITTING4,
    _NUM_FITTING_SAMPLES
};
#define NUM_FITTING_SAMPLES _NUM_FITTING_SAMPLES
#define NUM_SAMPLE_AUGPF3D (SHIFT_FITTING_SAMPLES+NUM_FITTING_SAMPLES)

#define NUM_SAMPLE_ERR      -1
#define NUM_SAMPLE_LMGPF3D (NUM_SAMPLE_AUGPF3D+2)
#define NUM_SAMPLE_MPGPF3D (NUM_SAMPLE_LMGPF3D+2)
#define NUM_SAMPLE_MXGPF3D (NUM_SAMPLE_MPGPF3D+2)
#define NUM_SAMPLE_SKGPF3D (NUM_SAMPLE_MXGPF3D+2)
#define NUM_SAMPLE_CHAIN   (NUM_SAMPLE_SKGPF3D+2)

//////////////////////////////////
//...typical examples in acoustic;
#define ACOU_L_DOMAIN      SHIFT_AU3D_SAMPLES+_AU1
#define ACOU_BOX_DOMAIN    SHIFT_AU3D_SAMPLES+_AU2
#define ACOU_USTUP_DOMAIN  SHIFT_AU3D_SAMPLES+_AU3

#define MAXW_L_DOMAIN      SHIFT_MX3D_SAMPLES+_MX1
#define MAXW_BOX_DOMAIN    SHIFT_MX3D_SAMPLES+_MX2
#define MAXW_USTUP_DOMAIN  SHIFT_MX3D_SAMPLES+_MX3

///////////////////////////////////////////////////////////////////////
//...структура одной записи для организации работы с таблицей образцов;
struct Record {
	int  type;
	char * parm_name;		// Parameter name
	void * parm_values;	// Parameter values
};

enum Num_type_Record {
					 ERR_TYPE_RECORD = -1,
				   CHAR_TYPE_RECORD,
		 DLENGTH_CHAR_TYPE_RECORD,
		 DSTRESS_CHAR_TYPE_RECORD,
		 DTEMPER_CHAR_TYPE_RECORD,
DTEMPER_EXTENT_CHAR_TYPE_RECORD,
		 DFREQUE_CHAR_TYPE_RECORD,
		 COMMENT_CHAR_TYPE_RECORD,
					 INT_TYPE_RECORD,
				  FLOAT_TYPE_RECORD,
				 DOUBLE_TYPE_RECORD,
			   DLENGTH_TYPE_RECORD,
				 DANGLE_TYPE_RECORD,
  			   DTEMPER_TYPE_RECORD,
	  DTEMPER_EXTENT_TYPE_RECORD,
			 DVELOCITY_TYPE_RECORD,
			   DSQUARE_TYPE_RECORD,
			  DINERTIA_TYPE_RECORD,
		  DRESISTANCE_TYPE_RECORD,
				  DMASS_TYPE_RECORD,
	   DMASS_DENSITY_TYPE_RECORD,
			   DWEIGHT_TYPE_RECORD,
				 DFORCE_TYPE_RECORD,
		  DFORCE_LINE_TYPE_RECORD,
			   DSTRESS_TYPE_RECORD,
			  DPERCENT_TYPE_RECORD, 
	 DDIMLESS_FREQUE_TYPE_RECORD, 
				 NUM_OF_TYPE_RECORD, 
	};

struct Shablon {
    Num_type_Record parm_type;
    char * parm_name;
};

/////////////////////////////////////////////////////////////////////
//...структура таблицы для организации работы с библиотекой образцов;
struct Table {
    int  N;              // Parameter count
    int  N_group;        // Group count
    int  index;          // Current group index (0-custom group)
    int  *  table_units; // Group units array
    char ** table_names; // Group names array
    Record  table[1];
};

////////////////////////////////////////////////////////////
//...структура контекста для организации работы с образцами;
struct Context {
    int     N, units, concat, method, Do, id_struc, static_char;
    double  E, nju, Ro, Hz, velocity, Om, Eps, h, h_min, h_med;
    double  left, right, bottom, top, back, front;
    Table   * table;
    Param   * param;
    char * sample_name;
    char * GOST_name;
    char    * sample_comment;
    char    * sample_description;
#ifndef ___ABRIDGE_PROFILE_MODE___
    CSample * sm;
#else
    void    * sm;
#endif
    int     * err_code;
    char    * MSG;
};

/////////////////////////////////////////////////////////////
//...вспомогательные функции для работы с таблицами образцов;
void delete_table   (Table *& tab, int id_static_char = NULL_STATE);
void set_table_units(Table *  tab, int id_units);
void set_table_param(Table *  tab, int k, Num_type_Record type, char * parm_name, void * parm_values);
void table_cpy      (Table *  tab, Table * tab_source, int id_static_char = NULL_STATE);
Table * get_table        (int N, int N_group, char * table_names[]);
Table * get_shablon_table(int N, int N_group, Shablon * records, char * table_names[], int id_static_char = ADDITIONAL_STATE);
Table * get_shablon_table(Table * table, int id_static_char = NULL_STATE);

////////////////////////////////////////////////
//...добавление и удаление параметра из таблицы;
void add_table_param(Table *& tab, int N, Shablon * records, int N_ini = -1);
void del_table_param(Table *& tab, int N, int N_ini = -1, int id_static_char = NULL_STATE);
void add_table_group(Table *  tab, int N_group, char ** table_names, int N_ini = -1, int id_static_char = NULL_STATE);
void del_table_group(Table *  tab, int N_group, int N_ini = -1, int id_static_char = NULL_STATE);

////////////////////////////////////////////////////
//...перемена порядка группы параметров на обратный;
void inverse_table_group(Table * tab);

///////////////////////////////////////////////////////////////////
//...проверка попадания точки в пространственный размерный обpазец;
int Sample3DIn(void * context, int & k, double * P);

//////////////////////////////////////
//...абсолютная идентификация образца;
int is_Lame3D(void * context);
int is_Maxw3D(void * context);
int is_Acou3D(void * context);
int is_Acou2D(void * context);
int is_Skin2D(void * context);
int is_Mapi3D(void * context);
int is_Mapi2D(void * context);
int is_BarsGOST2D(void * context);
int is_Bars2D(void * context);
int is_Wire2D(void * context);
int is_Cable2D(void * context);
int is_Fitting(void * context);
int is_Tower3D(void * context);
int is_LmGp3D(void * context);
int is_LmIb3D(void * context);
int is_MpGp3D(void * context);
int is_MpIb3D(void * context);
int is_MxGp3D(void * context);
int is_MxIb3D(void * context);
int is_AuGp3D(void * context);
int is_AuIb3D(void * context);
int is_SkGp3D(void * context);
int is_SkIb3D(void * context);
int is_Chain (void * context);

//////////////////////////////////////////////////////////////////////////////////
//...уничтожение всего контекста образца и операции над образцом внутри контекста;
void DeleteContext(void *& context);
void DeleteSample (void *  context);
void *  GetSample (void *  context);
void    SetSample (void *  context, void * sample);
void	  OutSample (void *  context, FILE * device);
#ifndef ___ABRIDGE_PROFILE_MODE___
CSample * GetContextSample(void * context);
#endif

///////////////////////////////////////////////////////////////////
//...вспомогательные функции, устанавливающие параметры в контекст;
void set_default (Context * cont);
void SetContextParam(void * context, Param * param,  int id_free = 1);
void SetSampleNumber(void * context, int Num);
int  SetMaterial    (void * context, double E, double nju, double Ro);
int  SetVelocity    (void * context, double velocity);
void SetTableIndex  (void * context, int index);
void SetTableParam  (Table * table,  int index, int N, double f, int id_static_char = NULL_STATE);
void SetTableParam  (Table * table,  int N, double f, int id_static_char = NULL_STATE);
void SetTableParam  (void * context, int index, int N, double f);
void SetTableParam  (void * context, int N, double f);
double GetTableParam(Table * table,  int index, int N);
double GetTableParam(Table * table,  int N);
double GetTableParam(void * context, int index, int N);
double GetTableParam(void * context, int N);
int    GetTableIndex(void * context);
int    GetNParam    (void * context);
int    GetNGroup    (void * context);
void SetUnits       (void * context, int units);
void SetConcat      (void * context, int concat);
void SetIdStruc     (void * context, int id_struc);
void SetMethod      (void * context, int method, double * param = NULL);
void SetGridRect    (void * context, double h = 0., double h_min = 0., double h_med = 0., int id_struc = 1);

void str_reinst(char *& str, char * str_concat, int & buf_count, int id_reinst = 0, int buf_length = 200);
void str_reinst(char *& str, const char * str_concat, int & buf_count, int id_reinst = 0, int buf_length = 200);

///////////////////////////////////////////////////
//...запись выделенной строки параметров в таблицу;
char * SetTableParamsAsString(Table * table, int index, char * pchar, char *& p_end, const char * ColonStr = COLONSTR, const char * SemiColonStr = SEMICOLONSTR);
const char * SetTableParamsAsString(Table * table, int index, const char * strParams, char *& strWork, int & buf_count, const char * ColonStr = COLONSTR, const char * SemiColonStr = SEMICOLONSTR);

////////////////////////////////////////////////////////////////////////
//...преобразование строки параметров таблицы в текстовое представление;
void GetTableParamsAsString(Table * table, int index, char *& str, int & buf_count, const char * ColonStr = COLONSTR, const char * SemiColonStr = SEMICOLONSTR);

//////////////////////////////////////////////////////////////////////
//...загрузка и обнуление параметров таблицы в позицию "по-умолчанию";
void LoadingParameters(Table * table,  int index);
void LoadingParameters(void * context, int index);
void ZeroCustomParameters(Table * table);
void ZeroCustomParameters(void * context);

///////////////////////////////////////////////////////////////
//...копирование параметров "по-умолчанию" в группу параметров;
void CopyTableParameters(Table * table,  int index, int index_ini = 0, int N_group = 1);
void CopyTableParameters(void * context, int index, int index_ini = 0, int N_group = 1);

///////////////////////////////////////////////////////////////////////////////////////////////////
//...добавление параметров "по-умолчанию" в группу параметров и преобразование статической таблицы;
void AddCustomParameters(Table * table,  char * table_name, int index = -1, int id_static_char = NULL_STATE);
void AddCustomParameters(void * context, char * table_name, int index = -1);
void AddTableParameters (void * context, int N_group, char ** table_names = NULL, int index = -1);
void ConvertStaticTable (void * context, int id_static_char = NULL_STATE);

/////////////////////////////////////////////
//...сохранение таблицы из контекста в файле;
void version_out(FILE * device, int wrongVersion = NULL_STATE);
void table_out	 (FILE * device, void * context, int head_writing = NULL_STATE, int id_csv_state = 0);
int  TableOut(char * ch_TABLE,  void * context, int head_writing = NULL_STATE);
int  TableOut(const char * ch_TABLE, void * context, int head_writing = NULL_STATE);

/////////////////////////////////////////////////////////////////////////////
//...шаблон CSV таблицы и шаблон таблиц для описания образцов из модуля BARs;
Table * GetCSVTable (int N_group, int N);
Table * GetBARsTable(int N_group, int N);

////////////////////////////////////////////////////
//...вспомогательные функции для зачитывания таблиц;
unsigned long VersionReading(char * pchar);

/////////////////////////////////////////////////
//...зачитывание csv представления любой таблицы;
void TableCSVReading(void * context, char *& pchar, char * p_end);

/////////////////////////////////////////////////////////////////
//...вспомогательные функции, извлекающие параметры из контекста;
void GetMaterial    (void * context, double & E, double & nju, double & Ro);
int  GetIdStaticChar(void * context);
int  GetUnits       (void * context);
int  GetConcat      (void * context);

///////////////////////////////////////////////////////////////////////////////////////////////////
//...извлечение границ рамки образца
double  Get2DLeft   (void * context);
double  Get2DRight  (void * context);
double  Get2DBottom (void * context);
double  Get2DTop    (void * context);
double  Get3DNear   (void * context);
double  Get3DFar    (void * context);
void SetFrameContext(void * context, double X1, double X2, 
                                     double Y1, double Y2, double Z1, double Z2);

////////////////////////////////////////////////////////////////////////////////////////
//...вспомогательные функции, описывающие одномерную связную компоненту контуpа образца;
int  Sample2DContourLinks (void * context);
int  Sample2DContourPoints(void * context, int id_link, int Max_knee);
void Sample2DContour      (void * context, int id_link, int Max_knee, double * X, double * Y);

///////////////////////////////////////////////////////////////////////////////////////////////
//...извлечение имени обpазца, его номера, таблицы, описания и степени решенности из контекста;
Param * GetContextParam     (void * context);
char  * GetSampleName (void * context);
char  * GetGOSTName   (void * context);
char  * GetProfileName(void * context, int index);
int     GetSampleNumber     (void * context);
int     GetSampleGpfNumber  (int id_problem = 0);
int     GetSampleIbrNumber  (int id_problem = 0);
int     GetSampleChainNumber();
void  * GetTable            (void * context);
char  * GetSampleDescription(void * context);
char  * GetSampleComment    (void * context);
int     GetSampleBlockNum   (void * context);

///////////////////////////////////////////////////////////////
//...извлечение и установка геометрии обpазца по его контексту;
#ifndef ___ABRIDGE_PROFILE_MODE___
CCells * get_bar (void * context, int id_block = -1);
void * GetContextBar (void * context, int id_block = -1);
void   SetContextBar (void * sample, void * bar, int id_block = -1);
CMap * GetContextMap (void * context, int id_block);
#endif

///////////////////////////////////////////////////////////////////////////////////////
//...функции образования контекста для образца по его абсолютному номеру и по описанию;
void * CreateContext  (int N_sm);
void * CreateContext	 (void * context);
void * CreateGPContext(char * description, int id_problem = 0);
void * CreateIBContext(char * description, int id_problem = 0);
void * CreateDSContext(char * description);
void	SetSampleDescription	(void * context, char * description);
void	SetSampleDescription	(void * context, const char * description);
void  SetSampleComment(void * context, char * comment);
void  SetSampleComment(void * context, const char * comment);
void	SetTable		 (void * context, Table * table, int id_static_char = NULL_STATE);

////////////////////////////////////////
//...list of several problems operating;
void ** CreateNewList(int iSamplesCount);
void    DestroyList(void **& list);
void ** CreateSKList();
void ** CreateAUList();
void ** Create3DList();
void ** CreateBARSGOSTList();
void ** CreateBARSList();
void ** CreateWIREList();
void ** CreateReducedWIREList();
void ** CreateCABLEList();
void ** CreateTOWERList();
void ** CreateFITTINGList();

///////////////////////////////////////////////////////////////
//...предварительная инициализация для статического 3D образца;
#ifndef ___ABRIDGE_PROFILE_MODE___
int sample3D_init(void * context, CGrid * block_nd = NULL, int type = NULL_BLOCK);
#endif

#endif
