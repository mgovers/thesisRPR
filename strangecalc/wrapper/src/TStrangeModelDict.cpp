// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME TStrangeModelDict

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
#include "TStrangeModel.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_TStrangeModel(void *p = 0);
   static void delete_TStrangeModel(void *p);
   static void deleteArray_TStrangeModel(void *p);
   static void destruct_TStrangeModel(void *p);
   static void streamer_TStrangeModel(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TStrangeModel*)
   {
      ::TStrangeModel *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TStrangeModel >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TStrangeModel", ::TStrangeModel::Class_Version(), "TStrangeModel.h", 37,
                  typeid(::TStrangeModel), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TStrangeModel::Dictionary, isa_proxy, 17,
                  sizeof(::TStrangeModel) );
      instance.SetNew(&new_TStrangeModel);
      instance.SetDelete(&delete_TStrangeModel);
      instance.SetDeleteArray(&deleteArray_TStrangeModel);
      instance.SetDestructor(&destruct_TStrangeModel);
      instance.SetStreamerFunc(&streamer_TStrangeModel);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TStrangeModel*)
   {
      return GenerateInitInstanceLocal((::TStrangeModel*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TStrangeModel*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr TStrangeModel::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TStrangeModel::Class_Name()
{
   return "TStrangeModel";
}

//______________________________________________________________________________
const char *TStrangeModel::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TStrangeModel*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TStrangeModel::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TStrangeModel*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TStrangeModel::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TStrangeModel*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TStrangeModel::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TStrangeModel*)0x0)->GetClass(); }
   return fgIsA;
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TStrangeModel(void *p) {
      return  p ? new(p) ::TStrangeModel( (TRootIOCtor *)nullptr ) : new ::TStrangeModel( (TRootIOCtor *)nullptr );
   }
   // Wrapper around operator delete
   static void delete_TStrangeModel(void *p) {
      delete ((::TStrangeModel*)p);
   }
   static void deleteArray_TStrangeModel(void *p) {
      delete [] ((::TStrangeModel*)p);
   }
   static void destruct_TStrangeModel(void *p) {
      typedef ::TStrangeModel current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_TStrangeModel(TBuffer &buf, void *obj) {
      ((::TStrangeModel*)obj)->::TStrangeModel::Streamer(buf);
   }
} // end of namespace ROOT for class ::TStrangeModel

namespace {
  void TriggerDictionaryInitialization_TStrangeModelDict_Impl() {
    static const char* headers[] = {
"TStrangeModel.h",
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
#line 1 "TStrangeModelDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate(R"ATTRDUMP(strangecalc model)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(strangecalc model)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TStrangeModel.h")))  TStrangeModel;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "TStrangeModelDict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "TStrangeModel.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"TStrangeModel", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("TStrangeModelDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_TStrangeModelDict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_TStrangeModelDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_TStrangeModelDict() {
  TriggerDictionaryInitialization_TStrangeModelDict_Impl();
}
