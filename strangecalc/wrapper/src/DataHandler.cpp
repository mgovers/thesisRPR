#include "DataHandler.h"
#include <TDataset.h>
#include <iostream>
#include <stdexcept>
#include <sstream>

using std::string;


/* DataHandler class *********************************************************************/

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 * Constructor
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
DataHandler::DataHandler()
  : m_VarEntry(NULL), m_ObsEntry(NULL)
{}


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 * Destructor
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
DataHandler::~DataHandler()
{
  clear();
}


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 * Public static function that constructs a static
 * DataHandler object (featuring resource-managing facilities
 * through its destructor) which reads and imports the data
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
void DataHandler::importData(Data dataPoints[], int* dataCount, const std::string& dataFile)
{
  static DataHandler dataHandler;

  try {
    dataHandler.readData(dataPoints, dataCount, dataFile);
  }
  catch(const std::exception& e)
  {
    std::cerr << e.what() << "\n";
  }
}


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 * Public static function that constructs a DataHandler
 * object (featuring resource-managing facilities through its
 * destructor) which verifies the integrity of the data file
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
bool DataHandler::isValidDataFile(const std::string& dataFile)
{
  try {
    DataHandler dataHandler;
    dataHandler.validate(dataFile);
  }
  catch(const std::exception& e)
  {
    std::cerr << e.what() << "\n";
    return false;
  }

  return true;
}


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 * Peek function that ignores single and tab spaces
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
char DataHandler::peek()
{
  while(m_Input.peek() == ' ' || m_Input.peek() == '\t')
    m_Input.ignore();

  return m_Input.peek();
}


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 * Verifies whether 'entry' is a valid number and stores the
 * result in 'number'
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
bool DataHandler::isValidNumber(const string& entry, double& number)
{
  std::istringstream iss(entry);
  iss >> number;

  return iss.eof() && !iss.fail();
}


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 * Print error and exit
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
void DataHandler::abort(const string& message)
{
  throw std::logic_error("\nError: " + message + "!\n" + "(data file: " + m_DataFile + ")");
}


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 * Ignore comment and empty lines
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
void DataHandler::ignoreComments()
{
  while(peek() == '#' || peek() == '\n')
    getline(m_Input, m_Entry);
}



/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 * Add variable entry
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
void DataHandler::add(VariableEntry* varEntry)
{
  m_DataEntries.push_back(varEntry);
  ++m_NumberOfVariables;
}


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 * Add observable entry
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
void DataHandler::add(ObservableEntry* obsEntry)
{
  m_DataEntries.push_back(obsEntry);
  m_ObservableEntries.push_back(obsEntry);
}


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 * Add dump entry
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
void DataHandler::add(DataEntry* dumpEntry)
{
  m_DataEntries.push_back(dumpEntry);
}


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 * Close acquired resources and release dynamically
 * allocated memory
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
void DataHandler::clear()
{
  if(m_Input.is_open())
    m_Input.close();

  if(m_VarEntry)
    delete m_VarEntry;

  for(it = m_DataEntries.begin(); it != m_DataEntries.end(); ++it)
    if(*it)
      delete *it;
}


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 * Open data file and intialize variables, observables,
 * and photo or electro struct
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
void DataHandler::init()
{
  clear();
  m_NumberOfVariables = 0;
  m_DataEntries.clear();
  m_ObservableEntries.clear();

  const string type = m_DataFile.substr(0, m_DataFile.find('.'));

  if(!(!type.compare("photo") || !type.compare("electro")))
    abort("Unknown reaction type '" + type + "'");

  m_IsPhoto = !type.compare("photo");

  const size_t pos = m_DataFile.find("iso.") + 4;
  const string isoString = m_DataFile.substr(pos, m_DataFile.find('.', pos) - pos);
  const int iso = atoi(isoString.c_str());

  if(iso < 1 || iso > 10) // excluding radiative capture
    abort("Unknown isospin channel '" + isoString + "'");

  m_Input.open(TDataset::GetDataFolder() + type + "/iso." + isoString + "/" + m_DataFile);

  if(!m_Input.is_open())
    abort("Could not open data file");


  // Initialize photo or electro struct
  // ----------------------------------
  if(m_IsPhoto)
    {
      m_Photo = Photo();
      m_Photo.iso = iso;
      m_Photo.emax = -1;
    }
  else
    {
      m_Electro = Electro();
      m_Electro.iso = iso;
      m_Electro.cos_ang = 1;
      m_Electro.beam_ener_input = 1;  // beam energy as default electron beam variable
    }

  interpretHeader();
}


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 * Interpret data file header
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
void DataHandler::interpretHeader()
{
  ignoreComments(); m_Input >> m_Entry;

  if(!m_Entry.compare("$vars:"))
    while(peek() != '\n')
      {
        m_Input >> m_Entry;

        if(!m_Entry.compare("/"))
          {
            add(new DataEntry());

            while(!(peek() == '&' || peek() == '\n'))
              m_Input >> m_Entry;
          }
        else if(!interpretKinematicVariable())
          interpretObservable();

        if(peek() == '&')
          m_Input >> m_Entry;
      }
  else
    abort("The '$vars:' specifier should be the first entry in the data file; received '" + m_Entry + "'");

  if(m_NumberOfVariables == 0)
    abort("No kinematic variables are specified");
  else if(m_ObservableEntries.size() == 0)
    abort("No observables are specified");

  m_Input.ignore();  // Dump current line
  ignoreComments();
}


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 * Interpret constant variables in data file
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
void DataHandler::checkForConstSpecifier()
{
  if(peek() == '$')
    {
      m_Input >> m_Entry;

      if(!m_Entry.compare("$const:"))
        {
          while(peek() != '\n')
            {
              m_Input >> m_Entry;

              if(!interpretKinematicVariable(true))
                abort("Unknown kinematic variable '" + m_Entry + "'");

              if(peek() == '&')
                m_Input >> m_Entry;
            }

           m_Input.ignore();  // Dump current line
           ignoreComments();
        }
      else
        abort("Only the '$const:' specifier is allowed in between data rows; received '" + m_Entry + "'");
    }
}


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 * Interpret kinematic variables in data file header
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
bool DataHandler::interpretKinematicVariable(bool isConst)
{
  if(m_IsPhoto)
    if(!m_Entry.compare("wlab") || !m_Entry.compare("wlab_min"))
      m_VarEntry = new VariableEntry(m_Entry, &m_Photo.emin);
    else if(!m_Entry.compare("w"))
      {
        m_VarEntry = new VariableEntry(m_Entry, &m_Photo.emin);
        m_Photo.is_w = 1;
      }
    else if(!m_Entry.compare("wlab_max"))
      m_VarEntry = new VariableEntry(m_Entry, &m_Photo.emax);
    else if(!m_Entry.compare("costhkcm") || !m_Entry.compare("costhkcm_min"))
      m_VarEntry = new VariableEntry(m_Entry, &m_Photo.cos);
    else if(!m_Entry.compare("costhkcm_max"))
      {
        m_VarEntry = new VariableEntry(m_Entry, &m_Photo.cos_max);
        m_Photo.is_cos_bin = 1;
      }
    else if(!m_Entry.compare("thetacm"))
      m_VarEntry = new VariableEntry(m_Entry, &m_Photo.cos, cos);
    else if(!m_Entry.compare("-t") || !m_Entry.compare("t"))
      {
        m_VarEntry = new VariableEntry(m_Entry, &m_Photo.t);
        m_Photo.t_ang = 1;
        m_VarEntry->m_Unit = (m_Entry.front() == '-' ? -1 : 1);
      }
    else
      return false;
  else
    if(!m_Entry.compare("w"))
      {
        m_VarEntry = new VariableEntry(m_Entry, &m_Electro.s);
        m_Electro.is_w = 1;
      }
    else if(!m_Entry.compare("s"))
      m_VarEntry = new VariableEntry(m_Entry, &m_Electro.s);
    else if(!m_Entry.compare("costhkcm"))
      m_VarEntry = new VariableEntry(m_Entry, &m_Electro.cos);
    else if(!m_Entry.compare("thetacm"))
      m_VarEntry = new VariableEntry(m_Entry, &m_Electro.cos, cos);
    else if(!m_Entry.compare("-t") || !m_Entry.compare("t"))
      {
        m_VarEntry = new VariableEntry(m_Entry, &m_Electro.t);
        m_Electro.cos_ang = 0;
        m_VarEntry->m_Unit = (m_Entry.front() == '-' ? -1 : 1);
      }
    else if(!m_Entry.compare("qsquared"))
      m_VarEntry = new VariableEntry(m_Entry, &m_Electro.qsquared);
    else if(!m_Entry.compare("E_beam"))
      m_VarEntry = new VariableEntry(m_Entry, &m_Electro.e_beam_ener);
    else if(!m_Entry.compare("epsilon"))
      {
        m_VarEntry = new VariableEntry(m_Entry, &m_Electro.eps);
        m_Electro.beam_ener_input = 0;
      }
    else if(!m_Entry.compare("phi"))
      m_VarEntry = new VariableEntry(m_Entry, &m_Electro.phi);
    else
      return false;


  if(isConst)
    m_Input >> *m_VarEntry->m_Value;


  // Units and options
  // -----------------
  while(!(peek() == '&' || peek() == '\n'))
    {
      m_Input >> m_Entry;

      if(!m_Entry.compare("(GeV)"))
        m_VarEntry->m_Unit *= 1e3;
      else if(!m_Entry.compare("(GeV^2)"))
        m_VarEntry->m_Unit *= 1e6;
      else if(!m_Entry.compare("(deg)"))
        m_VarEntry->m_Unit = PI/180.;
      else
        abort("Option '" + m_Entry + "' unknown for variable '" + m_VarEntry->m_Name + "'");
    }

  if(isConst)
    {
      *m_VarEntry->m_Value *= m_VarEntry->m_Unit;

      if(m_VarEntry->m_Function)
        *m_VarEntry->m_Value = (*m_VarEntry->m_Function)(*m_VarEntry->m_Value);

      delete m_VarEntry;
    }
  else
    add(m_VarEntry);

  m_VarEntry = NULL;

  return true;
}


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 * Interpret observables in data file header
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
void DataHandler::interpretObservable()
{
  add(m_ObsEntry = new ObservableEntry(m_Entry));

  // Units and options
  // -----------------
  while(!(peek() == '&' || peek() == '\n'))
    {
      m_Input >> m_Entry;

      if(!m_Entry.compare("-sysErr") || !m_Entry.compare("-sysErr:abs"))
        m_ObsEntry->m_HasSysError = true;
      else if(!m_Entry.compare("-sysErr:rel"))
        {
          m_ObsEntry->m_HasSysError = true;
          m_ObsEntry->m_HasRelativeSysError = true;
        }
      else if(m_Entry.find("-sysErr:rel:") != string::npos)
        {
          m_ObsEntry->m_HasSysError = true;
          m_ObsEntry->m_HasRelativeSysError = true;
          m_ObsEntry->m_HasPredefinedSysError = true;

          if(!isValidNumber(m_Entry.substr(12), m_ObsEntry->m_SysError))
            abort("Option '" + m_Entry + "' does not contain a valid number");
        }
      else if(m_Entry.find("-sysErr:abs:") != string::npos)
        {
          m_ObsEntry->m_HasSysError = true;
          m_ObsEntry->m_HasPredefinedSysError = true;

          if(!isValidNumber(m_Entry.substr(12), m_ObsEntry->m_SysError))
            abort("Option '" + m_Entry + "' does not contain a valid number");
        }
      else if(!m_Entry.compare("-norm:1"))
        m_ObsEntry->m_HasAltNormConvention = true;
      else if(!m_Entry.compare("-dt"))
        m_ObsEntry->m_IsDcsDt = true;
      else if(!m_Entry.compare("-dcos"))
        m_ObsEntry->m_Unit = 1./2./PI;
      else if(!m_Entry.compare("(nb)"))
        m_ObsEntry->m_Unit *= 1e-3;
      else if(!m_Entry.compare("(mb)"))
        m_ObsEntry->m_Unit *= 1e3;
      else
        abort("Option '" + m_Entry + "' unknown for observable '" + m_ObsEntry->m_Name + "'");
    }
}


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 * Copy Photo struct
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
void DataHandler::copyPhotoStruct(Photo& target)
{
  target.iso = m_Photo.iso;
  target.emin = m_Photo.emin;
  target.emax = m_Photo.emax > 0 ? m_Photo.emax : m_Photo.emin;
  target.ds_dt = m_Photo.ds_dt;
  target.t_ang = m_Photo.t_ang;
  target.is_w = m_Photo.is_w;
  target.ampli = m_Photo.ampli;
  target.error = m_Photo.error;

  if(m_Photo.t_ang)
    target.t = m_Photo.t;
  else
    {
      target.cos = m_Photo.cos;
      target.is_cos_bin = m_Photo.is_cos_bin;

      if(m_Photo.is_cos_bin)
        target.cos_max = m_Photo.cos_max;
    }

  strcpy(target.observable, m_ObsEntry->m_Name);
}


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 * Copy Electro struct
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
void DataHandler::copyElectroStruct(Electro& target)
{
  target.iso = m_Electro.iso;
  target.s = m_Electro.s;
  target.qsquared = m_Electro.qsquared;
  target.phi = m_Electro.phi;
  target.ds_dt = m_Electro.ds_dt;
  target.cos_ang = m_Electro.cos_ang;
  target.is_w = m_Electro.is_w;
  target.beam_ener_input = m_Electro.beam_ener_input;
  target.cs_convention = m_Electro.cs_convention;
  target.ampli = m_Electro.ampli;
  target.error = m_Electro.error;

  if(m_Electro.cos_ang)
    target.cos = m_Electro.cos;
  else
    target.t = m_Electro.t;

  if(m_Electro.beam_ener_input)
    target.e_beam_ener = m_Electro.e_beam_ener;
  else
    target.eps = m_Electro.eps;

  strcpy(target.observable, m_ObsEntry->m_Name);
}


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 * Read data from input stream and store it in dataPoints[]
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
void DataHandler::readData(Data dataPoints[], int* dataCount, const string& dataFile)
{
  m_DataFile = dataFile;
  init();

  double* const obsValue = m_IsPhoto ? &m_Photo.ampli : &m_Electro.ampli;
  double* const obsError = m_IsPhoto ? &m_Photo.error : &m_Electro.error;
  short* const isDcsDt = m_IsPhoto ? &m_Photo.ds_dt : &m_Electro.ds_dt;
  const int numObs = m_ObservableEntries.size();

  while(peek() != EOF)
    {
      checkForConstSpecifier();

      for(it = m_DataEntries.begin(); it != m_DataEntries.end(); ++it)
        (*it)->read(m_Input);

      for(int i = 0; i < numObs; ++i)
        if((m_ObsEntry = m_ObservableEntries[i])->m_Error > 0)
          {
            dataPoints[*dataCount].photo_prod = m_IsPhoto ? 1 : 0;
            dataPoints[*dataCount].electro_prod = m_IsPhoto ? 0 : 1;
            dataPoints[*dataCount].kaoncapture = 0;

            if(m_ObsEntry->m_IsDcsDt)
              *isDcsDt = 1;

            if(m_ObsEntry->m_HasAltNormConvention)
              m_Electro.cs_convention = 1;

            *obsValue = m_ObsEntry->m_Value;
            *obsError = m_ObsEntry->m_Error;

            if(m_IsPhoto)
              copyPhotoStruct(dataPoints[*dataCount].photo);
            else
              copyElectroStruct(dataPoints[*dataCount].elec);

            ++(*dataCount);
          }

      ignoreComments();
    }
}


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 * Validate data and variables in data file
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
void DataHandler::validate(const string& dataFile)
{
  m_DataFile = dataFile;
  init();

  const int numObs = m_ObservableEntries.size();
  int numCols = m_DataEntries.size();

  for(int i = 0; i < numObs; ++i)
    {
      numCols += (m_ObservableEntries[i]->m_HasSysError &&
                  !m_ObservableEntries[i]->m_HasPredefinedSysError ? 2 : 1);
    }

  double number;

  while(peek() != EOF)
    {
      checkForConstSpecifier();

      for(int i = 0; i < numCols; ++i)
        {
          m_Input >> m_Entry;

          if(!isValidNumber(m_Entry, number))
            abort("Entry '" + m_Entry + "' is not a valid number");

          if(peek() == '\n' && i != numCols - 1)
            abort("Specified number of data columns does not match actual number");
        }

      ignoreComments();
    }
}


/* DataEntry utility class ***************************************************************/

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 * Constructor
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
DataEntry::DataEntry()
{}


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 * Destructor
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
DataEntry::~DataEntry()
{}


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 * Read entry and dump it
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
void DataEntry::read(std::ifstream& m_Input)
{
  static string dump;
  m_Input >> dump;
}


/* VariableEntry utility class ***********************************************************/

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 * Constructor
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
VariableEntry::VariableEntry(const std::string& name, double* value,
                             double (*function)(double))
  : DataEntry(), m_Value(value), m_Unit(1), m_Name(name), m_Function(function)
{}


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 * Read variable entry from m_Input stream
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
void VariableEntry::read(std::ifstream& m_Input)
{
  m_Input >> *m_Value;
  *m_Value *= m_Unit;

  if(m_Function)
    *m_Value = (*m_Function)(*m_Value);
}


/* ObservableEntry utility class *********************************************************/

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 * Constructor
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
ObservableEntry::ObservableEntry(const std::string& name)
  : DataEntry(), m_Unit(1), m_Value(0), m_Error(0), m_SysError(0),
    m_HasSysError(false), m_HasRelativeSysError(false),
    m_HasPredefinedSysError(false), m_IsDcsDt(false),
    m_HasAltNormConvention(false), m_ReadData(true)
{
  strcpy(m_Name, name.c_str());
}


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 * Read observable entry from m_Input stream
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
void ObservableEntry::read(std::ifstream& m_Input)
{
  m_Input >> m_Value >> m_Error;

  if(m_HasSysError && !m_HasPredefinedSysError)
    m_Input >> m_SysError;

  if(m_Error == 0.)  // CAUTION !!!
    return;

  if(m_HasSysError)
    m_Error = sqrt(m_Error*m_Error + m_SysError*m_SysError*
                   (m_HasRelativeSysError ? m_Value*m_Value : 1));

  m_Value *= m_Unit;
  m_Error *= m_Unit;
}
