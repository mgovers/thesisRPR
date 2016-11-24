// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME TDatasetDict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "TDataset.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_TDataset(void *p = 0);
   static void *newArray_TDataset(Long_t size, void *p);
   static void delete_TDataset(void *p);
   static void deleteArray_TDataset(void *p);
   static void destruct_TDataset(void *p);
   static void directoryAutoAdd_TDataset(void *obj, TDirectory *dir);
   static void streamer_TDataset(TBuffer &buf, void *obj);
   static Long64_t merge_TDataset(void *obj, TCollection *coll,TFileMergeInfo *info);
   static void reset_TDataset(void *obj, TFileMergeInfo *info);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TDataset*)
   {
      ::TDataset *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TDataset >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TDataset", ::TDataset::Class_Version(), "TDataset.h", 39,
                  typeid(::TDataset), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TDataset::Dictionary, isa_proxy, 17,
                  sizeof(::TDataset) );
      instance.SetNew(&new_TDataset);
      instance.SetNewArray(&newArray_TDataset);
      instance.SetDelete(&delete_TDataset);
      instance.SetDeleteArray(&deleteArray_TDataset);
      instance.SetDestructor(&destruct_TDataset);
      instance.SetDirectoryAutoAdd(&directoryAutoAdd_TDataset);
      instance.SetStreamerFunc(&streamer_TDataset);
      instance.SetMerge(&merge_TDataset);
      instance.SetResetAfterMerge(&reset_TDataset);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TDataset*)
   {
      return GenerateInitInstanceLocal((::TDataset*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TDataset*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr TDataset::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TDataset::Class_Name()
{
   return "TDataset";
}

//______________________________________________________________________________
const char *TDataset::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TDataset*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TDataset::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TDataset*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TDataset::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TDataset*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TDataset::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TDataset*)0x0)->GetClass(); }
   return fgIsA;
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TDataset(void *p) {
      return  p ? new(p) ::TDataset : new ::TDataset;
   }
   static void *newArray_TDataset(Long_t nElements, void *p) {
      return p ? new(p) ::TDataset[nElements] : new ::TDataset[nElements];
   }
   // Wrapper around operator delete
   static void delete_TDataset(void *p) {
      delete ((::TDataset*)p);
   }
   static void deleteArray_TDataset(void *p) {
      delete [] ((::TDataset*)p);
   }
   static void destruct_TDataset(void *p) {
      typedef ::TDataset current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around the directory auto add.
   static void directoryAutoAdd_TDataset(void *p, TDirectory *dir) {
      ((::TDataset*)p)->DirectoryAutoAdd(dir);
   }
   // Wrapper around a custom streamer member function.
   static void streamer_TDataset(TBuffer &buf, void *obj) {
      ((::TDataset*)obj)->::TDataset::Streamer(buf);
   }
   // Wrapper around the merge function.
   static Long64_t merge_TDataset(void *obj,TCollection *coll,TFileMergeInfo *info) {
      return ((::TDataset*)obj)->Merge(coll,info);
   }
   // Wrapper around the Reset function.
   static void reset_TDataset(void *obj,TFileMergeInfo *info) {
      ((::TDataset*)obj)->ResetAfterMerge(info);
   }
} // end of namespace ROOT for class ::TDataset

namespace {
  void TriggerDictionaryInitialization_TDatasetDict_Impl() {
    static const char* headers[] = {
"TDataset.h",
0
    };
    static const char* includePaths[] = {
"/usr/local/include",
"/usr/include",
"/usr/local/include",
"/home/mgovers/Software/strangecalc/wrapper/src/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "TDatasetDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate(R"ATTRDUMP(Dataset for strangeViewer)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TDataset.h")))  TDataset;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "TDatasetDict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "TDataset.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"TDataset", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("TDatasetDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_TDatasetDict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_TDatasetDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_TDatasetDict() {
  TriggerDictionaryInitialization_TDatasetDict_Impl();
}
