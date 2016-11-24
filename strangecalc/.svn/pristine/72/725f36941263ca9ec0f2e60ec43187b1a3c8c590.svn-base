#ifndef DATAHANDLER_H
#define DATAHANDLER_H

#include <string>
#include <fstream>
#include <vector>
#include <fitting.h>


/* DataEntry utility class ***************************************************************/
class DataEntry
{
private:
  explicit DataEntry();
  virtual ~DataEntry();

  virtual void read(std::ifstream&);

  friend class DataHandler;
  friend class VariableEntry;
  friend class ObservableEntry;
};


/* VariableEntry utility class ***********************************************************/
class VariableEntry : private DataEntry
{
private:
  explicit VariableEntry(const std::string&, double*, double (*)(double) = NULL);

  double* const m_Value;
  double m_Unit;
  std::string m_Name;
  double (*m_Function)(double);

  void read(std::ifstream&);

  friend class DataHandler;
};


/* ObservableEntry utility class *********************************************************/
class ObservableEntry : private DataEntry
{
private:
  explicit ObservableEntry(const std::string&);

  double m_Unit, m_Value, m_Error, m_SysError;
  bool m_HasSysError, m_HasRelativeSysError, m_HasPredefinedSysError, m_IsDcsDt,
    m_HasAltNormConvention, m_ReadData;
  char m_Name[15];

  void read(std::ifstream&);

  friend class DataHandler;
};


/* Uncopyable utility class **************************************************************/
class Uncopyable
{
protected:
  Uncopyable() {}
  ~Uncopyable() {}

private:
  Uncopyable(const Uncopyable&);
  Uncopyable& operator=(const Uncopyable&);
};


/* DataHandler class *********************************************************************/
class DataHandler : private Uncopyable
{
public:
 static void importData(Data[], int*, const std::string&);
 static bool isValidDataFile(const std::string&);

private:
  explicit DataHandler();
  ~DataHandler();

  bool m_IsPhoto;
  int m_NumberOfVariables;
  std::string m_DataFile, m_Entry;
  std::ifstream m_Input;
  std::vector<DataEntry*> m_DataEntries;
  std::vector<ObservableEntry*> m_ObservableEntries;
  std::vector<DataEntry*>::iterator it;
  VariableEntry* m_VarEntry;
  ObservableEntry* m_ObsEntry;
  Photo m_Photo;
  Electro m_Electro;

  char peek();
  bool isValidNumber(const std::string&, double&);

  void abort(const std::string&);
  void ignoreComments();
  void add(VariableEntry*);
  void add(ObservableEntry*);
  void add(DataEntry*);
  void clear();
  void init();
  void interpretHeader();
  void checkForConstSpecifier();
  bool interpretKinematicVariable(bool = false);
  void interpretObservable();
  void copyPhotoStruct(Photo&);
  void copyElectroStruct(Electro&);
  void readData(Data[], int*, const std::string&);
  void validate(const std::string&);
};


#endif // DATAHANDLER_H
