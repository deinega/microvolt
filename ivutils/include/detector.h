#ifndef _DETECTOR_H
#define _DETECTOR_H

/** @file detector.h
  \en \brief 
  Classes for detectors that record field values during simulation. 
  Field values are E and H for electromagnetic simulation; n, p and psi for transport simulation etc.
  Detectors do not correspond to any objects in space and do not influence the field values. 
  Detectors can be at arbitrary points inside calculated volume and shouldn't necessary be aligned with mesh. 
  Field values at that points are interpolated by mesh.
  \ru \brief 
  Классы, реализующие детекторы, снимающие поля на сетке в течение численного эксперимента. 
  Детекторам не отвечают какие-то реальные тела в пространстве, и на значения полей они не влияют. 
  Детекторы могут размещаться в произвольном месте вычислительного объема и не привязаны к сетке. 
  Значения полей на них получаются путем интерполяции по соседним сеточным узлам.
*/

#include "loggerio.h"
#include "refobj.h"
#include "component.h"
#include "data_flags.h"
#include "cvector_3.h"

using namespace std;

/**\en
  Class for recording data.
  There are following phases to work with this class:
  1. Initialization using SetSize and SetRecordingFileNode
  2. record
  3. StartRecord
  4. NextRecord wiht array of current field values as an argument
  5. EndRecord
*/
/**\ru 
  Базовый класс, реализующий запись данных.
  Работа с ним предполагает следующие стадии:
  1. Конструирование, вызов SetSize, SetRecordingFileNode
  2. Последовательный вызов функции record (используется в случае Фурье на лету)
  3. Вызов StartRecord.
  4. Последовательный вызов функций NextRecord, аргументом которой является массив текущих значений полей
  5. Вызов EndRecord.
*/
class FieldsRecorder: public apComponent, public restorer{
protected:
  using apComponent::comm;
  using apComponent::config;

  ///\en file name
  ///\ru имя файла
  string fname;
  ///\en if first bit is turned on, then recorded fields values are assumed to be zero
  ///\ru если включен первый бит, то происходит холостое обновление нулями
  int wf;
  size_t rec_size; ///\en number of fields per point
  size_t bytes_size; ///en sizeof(vec_type)*rec_size
  ///\en number of points where fields are recorded
  ///\ru количество записываемых точек
  size_t vnum; 
public:
  FieldsRecorder(const string &fn):fname(fn),wf(0),rec_size(0),bytes_size(0),vnum(0){}
  void SetWorkFlag(int wf_){wf=wf_;}
  void SetFileName(const string &fn){fname=fn;}
  ///<\en size of one record in bytes
  ///<\ru размер единичной записи в байтах
  virtual void SetBytesSize(size_t ds){
    bytes_size=ds;
    rec_size=ds/sizeof(vec_type);
  }
  ///<\en number of points where field is recorded
  ///<\ru количество записываемых точек
  void SetSize(size_t sz){
    vnum=sz;
  }
  ///\en set node which records file (in parallel regime file can be recorded only from one node)
  ///\ru запись домена, который записывает итоговый файл (при параллельной работе файл записывается только на одном домене)
  virtual void SetRecordingFileNode(int node)=0;
  ///\en set node which process data for current point (used for Fourier-on-fly)
  ///\ru запись домена, который обрабатывает данные для текущей точки (используется в случае Фурье на лету)
  virtual void record(int node){}
  virtual void SetVectorSet(void *vs){}
  ///\en set using space coordinates, output fields and names for them
  virtual void SetTypes(int argtype, int outtype, vector<string> *columns=NULL){}
  ///\en set frequency range for recorders that perform Fourieron-fly transformation
  virtual void SetFrequencyRange(vec_type fmin, vec_type fmax, vec_type df){}
  ///\en returns false if this recorder is not working
  virtual bool active(){return true;}

  ///\en if data is recorded to file, flush data from buffer to disc
  ///\ru если происходит запись в файл, сбросить данные из буфера на жесткий диск
  virtual void Flush(){} 
  ///\en start record
  ///\en \param iteration from which recording is started (not zero if record is continuing after restarting the numerical experiment).
  ///\ru начать работу
  ///\ru \param с какой итерации начинается запись (не равно нулю, если запись продолжается после перезапуска эксперимента).
  virtual int StartRecord(int it)=0;
  virtual int NextRecord(void *filebuf)=0;
  virtual int EndRecord()=0;
  ///\en set extra string which will be written in the text file, used in text recorders only
  virtual void SetExtra(const string &extra_){}
  virtual size_t data_size() const{return 0;}
};

/**\en
  This class records fields to the binary file.
  After calculation this binary file can be used to plot fields history.
*/
/**\ru
  Последовательно записывает данные в файл.
  По окончании расчета на основании этого файла можно восстановить историю сигнала.
*/
class FileRecorder: public FieldsRecorder{
  using FieldsRecorder::comm;
  using FieldsRecorder::config;
  using FieldsRecorder::wf;
  using FieldsRecorder::fname;

  ///\en only node, which has if_record=true, records fields to the file.
  /// other nodes only collect field and send it to this node
  ///\ru записывает данные в файл тот домен, у которого установлен флаг if_record
  bool if_record;

  ///<\en 1 is initialization is successful, otherwise -1
  ///<\ru равен 1 при успешной инициализации, в противном случае равен -1
  int inifile; 
  SeqRecord SR;

public:
  FileRecorder(const string &fname_):
    if_record(false),inifile(-1),FieldsRecorder(fname_){}
  void SetBytesSize(size_t bs){
    FieldsRecorder::SetBytesSize(bs);
    SR.SetDataSize(bs);
  }
  void SetRecordingFileNode(int node){
    if(comm->get_myrank()==node)
      if_record=true;
  }
  void Flush();
  int StartRecord(int it);
  int NextRecord(void *filebuf);
  int EndRecord();
};

/**\en
  This class records fields to the text file.
  After calculation this text file can be used to plot fields history.
*/
template<class vset_t>
class TextRecorder: public FieldsRecorder{
protected:
  using FieldsRecorder::comm;
  using FieldsRecorder::config;
  using FieldsRecorder::wf;
  using FieldsRecorder::vnum;
  using FieldsRecorder::rec_size;
  
  ///\en only node, which has if_record=true, records fields to the file.
  /// other nodes only collect field and send it to this node
  bool if_record;

  ///<\en 1 is initialization is successful, otherwise -1
  int inifile;
  FILE *f; ///\en text file descriptor
  string format; ///\en data fprintf format

  vset_t *vset; ///\en points where fields are recorded

  int argtype; ///\en recorded space coordinates
  vector<string> columns; ///\en columns names which will be printed as a header of the file

public:
  TextRecorder(const string &fname_=""):
    if_record(false),inifile(-1),f(NULL),format("%g"),vset(NULL),argtype(0),FieldsRecorder(fname_){}

  void SetFormat(const string format_){
    format=format_;
  }

  void SetVectorSet(void *vs){
    vset=(vset_t *)vs;
  }

  void SetTypes(int argtype_, int outtype, vector<string> *columns_=NULL){
    argtype=argtype_;
    if(columns_)
      columns=*columns_;
  }

  void SetRecordingFileNode(int node){
    if(comm->get_myrank()==node)
      if_record=true;
  }
  void Flush();
  int StartRecord(int it);
  int NextRecord(void *filebuf);
  int EndRecord();
};

/**\en
  This class records fields fluxes through the surface defined by v_set to the text file.
*/
template<class vset_t>
class TextFluxRecorder: public TextRecorder<vset_t>{
  using TextRecorder<vset_t>::comm;
  using TextRecorder<vset_t>::config;
  using TextRecorder<vset_t>::wf;

  using TextRecorder<vset_t>::fname;
  using TextRecorder<vset_t>::vset;
  using TextRecorder<vset_t>::if_record;
  using TextRecorder<vset_t>::inifile;
  using TextRecorder<vset_t>::f;

  using TextRecorder<vset_t>::FieldsRecorder::vnum;
  using TextRecorder<vset_t>::FieldsRecorder::rec_size;

  using TextRecorder<vset_t>::argtype;
  using TextRecorder<vset_t>::columns;

  string extra; // extra column which will be placed at the text file
  
public:
  TextFluxRecorder(const string &fname=""):TextRecorder<vset_t>(fname){}

  void SetExtra(const string &extra_){
    if(extra_.length())
      extra=extra_+"\t";
  }

  int StartRecord(int it);
  int NextRecord(void *filebuf);
};

/**\en
  Perform Fourier-on-fly transformation.
  Forier-on-fly is performed at nodes having nonzero bufer size bufsize.
  Node which records binary file at the end of the simulation is defined by SetRecordingFileNode.
*/
/**\ru
  Осуществляет преобразование Фурье на лету.
  Делают Фурье на лету домены, у которых ненулевой размер буфера bufsize.
  Домен, который записывает итоговый файл, задается функцией SetRecordingFileNode.
*/
class emFourierRecorder: public FieldsRecorder{
  using FieldsRecorder::comm;
  using FieldsRecorder::config;
  using FieldsRecorder::wf;

  ///\en node which records binary file at the end of the simulation
  ///\ru домен, записывающий итоговый файл
  int wr_rank; 

  ///\en time step
  ///\ru шаг по времени
  vec_type dt;
  ///\en calculated frequencies
  ///\ru рассчитваемые частоты 
  vec_type fmin, fmax, df;
  ///\en frequencies number
  ///\ru количество частот
  int fnum; 
  ///\en representation of wave in time domain
  ///\ru представление плоской волны
  int exp_dir; 
  /// \en update frequency in steps (1 by default)
  int ufreq;
  ///\en current time iteration
  ///\ru текущий шаг по времени
  int ti; 
  ///\en number of points which are processed at the current node
  ///\ru количество точек, для которых делается Фурье на лету на текущем домене
  size_t bufsize; 
  ///\en nodes which make Fourier-on-fly for given points sequence.
  ///\en used in parallel realization of FileRecord
  ///\ru пакер доменов, которые делают Фурье на лету для заданной последовательности точек
  ///\ru используется в параллельной реализации FileRecord
//  int_pack nodes;
  vector<int> nodes;
  ///\en buffer for intermediate results
  ///\ru буфер для результатов
  vec_type *F;
  ///\en tables with some values of exponenta to optimize Fourier transform
  cvec_type *exp_table_dt, *exp_table_t;
  ///\en how often this table should be used (or zero, if this table is not used)
  int exp_upd;

public:

  emFourierRecorder(vec_type dt_, vec_type fmin_, vec_type fmax_, vec_type df_, int exp_dir_, int ufreq_, const string &fname_):
    wr_rank(0), FieldsRecorder(fname_), dt(dt_), fmin(fmin_), fmax(fmax_), df(df_), ufreq(ufreq_),
    exp_dir(exp_dir_==TF_BACKWARD?1:-1), ti(-1), fnum((df_ ? int(acdiv(fmax_-fmin_,df_)): 0)+1), 
    bufsize(0), F(), exp_table_dt(), exp_table_t(), exp_upd(0){}

  typedef size_t write_read_f(void *p, size_t size, size_t n, FILE *f);

  void SetRecordingFileNode(int node){
    wr_rank=node;
  }
  void record(int node);
  void SetFrequencyRange(vec_type fmin_, vec_type fmax_, vec_type df_){
    fmin=fmin_;
    fmax=fmax_;
    df=df_;
    if(fmax>=fmin)
      fnum=(df_ ? int(acdiv(fmax_-fmin_,df_)): 0)+1;
    else
      fnum=0;
  }
  virtual bool active(){
    return fnum>0;
  }
  ///<\en calls FileRecord
  ///<\ru вызывает FileRecord
  void Flush(); 
  int StartRecord(int it);
  int NextRecord(void *filebuf);
  int EndRecord();
  ///\en dump results to the file
  ///\ru скинуть результаты в файл
  void FileRecord(const string &fname); 
  int write_read(FILE *f, write_read_f *fun);
  size_t data_size() const{
    return fnum*bufsize*2*bytes_size;
  }
};

/**\en
  Base class for detectors.

  Work with this class assumes three phases:
  1. Intialization using Init and StartRecord.
  2. StartNextRecord and CompleteNextRecord.
  3. EndRecord.
*/
/**\ru 
  Базовый класс, реализующий детекторы.

  Работа с классом предполагает следующие стадии:
  1. Конструирование, инициализация с помощью Init и вызов StartRecord.
  2. Последовательный вызов функций StartNextRecord и CompleteNextRecord.
  3. Вызов EndRecord.
*/
class BaseDetectorSet: public apComponent, public restorer{
public:
  ///\en set bit flag
  ///\ru установка битового флага
  virtual void SetWorkFlag(int wf)=0; 
  virtual void SetFileName(const string &fn){}
  virtual void SetFrequencyRange(vec_type fmin, vec_type fmax, vec_type df){}
  virtual int Init()=0;
  virtual void clear(){}
  ///\en flush data from memory to the disc
  ///\ru сбросить данные на диск
  virtual void Flush()=0; 
  ///\en start record
  ///\en \param iteration from which recording is started (not zero if record is continuing after restarting the numerical experiment).
  ///\ru начать работу
  ///\ru \param it с какой итерации начинается запись (не равно нулю, если запись продолжается после перезапуска эксперимента).
  virtual int StartRecord(int it=0)=0;
  ///\en start record current field values
  ///\ru начать запись текущих значений
  virtual int StartNextRecord()=0; 
  ///\en end record current field values
  ///\ru закончить запись текущих значений
  virtual int CompleteNextRecord()=0; 
  ///\ru завершить работу
  ///\en finish work
  virtual int EndRecord()=0; 
  ///\en set extra string which will be written in the text file, used in text recorders only
  virtual void SetExtra(const string &extra){};

  virtual base_vset_it *begin(){return NULL;}
  virtual size_t vec_size(){return 0;}
};

/**\en
  Detectors recording field values in points sequence vset_t.
  Field values are given by container container_t.
*/
/**\ru 
  Детекторы, считывающие предоставляемые контейнером container_t значения 
  в совокупности точек, задаваемой vset_t
*/
template<class vset_tt, class xfer_t>
class DetectorSet: public BaseDetectorSet{
public:
  typedef vset_tt vset_t;
  typedef typename vset_t::iterator vec_it;

protected:

  ///\en points where fields are recorded
  ///\ru точки, из которых записываются данные
  mngptr<vset_t> vsetp; 
  mngptr<VecTransform> transform; ///<\en vector transformation
  ///\en class which is used to get field values
  ///\ru класс, используемый для получения значений полей
  xfer_t xfer;
  int rec_out; // record values in points outside container as a zero
  ///\en this object directly records fields
  ///\ru непосредственно отвечает за запись данных
  mngptr<FieldsRecorder> R; 
  ///\en buffer for fields values from container
  ///\ru буфер для данных
  char *filebuf; 

  int bufnum; ///<\en id of buffer filebuf in xfer
  int bufsize; ///<\en detectors number

  ///\en 1 if initialization is successfull, otherwise -1
  ///\ru равен 1 при успешной инициализации, в противном случае равен -1
  int inirec; 
  ///\en if first bit is turned on, then recorded fields are assumed to be zero
  ///\ru битовый флаг. если поднят первый бит, то происходит холостое обновление нулями
  int wf; 
  /**\en
  Flag to control how files are recorded in parallel regime.
  First bit is off if zeroth node records the file, 
  otherwise the file is recorded by node that have more points.
  For the case of Fourier-on-fly (R==emFourierRecorder), one can use second bit. 
  Second bit is on if then points field values are Fouirier transformed at the node that have these points, 
  otherwise they are sent to the recording node, which is defined by first bit.
  */
  /**\ru 
  Битовый флаг, регулирующий запись файлов процессорами в случае параллельного исполения.
  Первый бит равен единице, если пишет файл нулевой процессор,
  в противном случае пишет тот процессор, на который больше всего попадает точек.
  Если выполняется Фурье на лету (R==emFourierRecorder), то можно использовать второй бит. 
  Если он включен, то данные детекторов подвергаются Фурье на лету на своем процессоре, 
  в противном случае, они пересылаются на пишущий процессор, который зависит от значения первого бита.
  */
  int local_write;
  ///\en auxiliary флаг, that is turned on at StartNextRecord and turned off at CompleteNextRecord
  ///\en used to control if these function are called in the right order
  ///\ru вспомогательный флаг, включаемый при StartNextRecord и выключаемый при CompleteNextRecord
  ///\ru нужен для контроля вызова этих функций в правильной последовательности
  bool if_sent;

  ///\en timers numbers
  ///\ru номера таймеров
  struct{
    int det_pre, det_rec;
  } tid;

  ///\en gets the domain having maximal number of detectors
  int find_major_rank();

public:

  DetectorSet(mngarg<FieldsRecorder> R_, typename xfer_t::container_t *container_, int lwr,
  mngarg<vset_t> vset, mngarg<VecTransform> tr=make_mngarg(new VecTransform(),1), int rec_out_=0):
    inirec(-1),xfer(container_),rec_out(rec_out_),filebuf(NULL),bufsize(0),wf(0),R(R_),
    local_write(lwr),vsetp(vset),transform(tr){

    R->SetVectorSet(vsetp.ptr());
    R->SetBytesSize(xfer.get_full_data_size());
  }

  virtual ~DetectorSet() {
    if(filebuf)delete []filebuf;
  }

  void SetFrequencyRange(vec_type fmin, vec_type fmax, vec_type df){
    R->SetFrequencyRange(fmin, fmax, df);
  }

  const vset_t *GetVectorSet() const{
    return vsetp.ptr();
  }

  const VecTransform *GetVecTransform() const{
    return transform.ptr();
  }

  void AddDefaultTimers(const vector<int> &p_id);

  virtual void SetWorkFlag(int wf);

  virtual void SetFileName(const string &fn){
    R->SetFileName(fn);
  }

  void set_rec_out(int rec_out_){
    rec_out=rec_out_;
  }

  void Dump(bool lim=false){
//    if(dump&dmpDETPOINTS)
      apComponent::Dump(lim);
  }

  ///\en set using space coordinates, output fields and names for them
  void SetTypes(int argtype, int outtype, vector<string> *columns=NULL){
    xfer.set_full_data(outtype);
    R->SetTypes(argtype, outtype, columns);
    R->SetBytesSize(xfer.get_full_data_size());
  }

  ///\en prepare parallel transfers, allocate buffer memory
  ///\ru подготавливает параллельные пересылки, резервирует буфер
  virtual int Init();

  virtual void clear(){
     inirec=-1;
     xfer.clear(); 
  }

  ///\en calls R->Flush
  ///\ru вызывает R->Flush
  void Flush();

  virtual int StartRecord(int it=0);

  ///\en start parallel transfers in xfer
  ///\ru стартует параллельные пересылки xfer
  virtual int StartNextRecord();

  ///\en records data from xfer to buffer and calls R->NextRecord
  ///\ru записывает данные xfer в буфер и вызывает R->NextRecord
  virtual int CompleteNextRecord();

  ///\en calls R->EndRecord
  /// вызывает R->EndRecord
  virtual int EndRecord();

  void SetExtra(const string &extra){
    R->SetExtra(extra);
  }

  int write_read(FILE *f, write_read_f *fun){
    return R->write_read(f,fun);
  }

  virtual base_vset_it *begin(){return new virt_vset_it<typename vset_t::iterator>(vsetp->begin());}
  virtual size_t vec_size(){return vsetp->size();}
};


///\en This macro provides all instantiations with given transfer type
#define INSTANTIATE_DETECTORS(transfer_t) \
  template class DetectorSet<UniformGrid<Vector_3>, transfer_t>; \
  template class DetectorSet<BoxSurfaceSet, transfer_t>;\
  template class DetectorSet<SphereSurfaceSet, transfer_t>;\
  template class DetectorSet<vector<Vector_3>, transfer_t>;\
  template class DetectorSet<SpaceVectorSet, transfer_t>;


#endif
