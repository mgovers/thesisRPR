// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME TMultiModelDict

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
#include "TMultiModel.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_TMultiModel(void *p = 0);
   static void delete_TMultiModel(void *p);
   static void deleteArray_TMultiModel(void *p);
   static void destruct_TMultiModel(void *p);
   static void streamer_TMultiModel(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TMultiModel*)
   {
      ::TMultiModel *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TMultiModel >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TMultiModel", ::TMultiModel::Class_Version(), "TMultiModel.h", 36,
                  typeid(::TMultiModel), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TMultiModel::Dictionary, isa_proxy, 17,
                  sizeof(::TMultiModel) );
      instance.SetNew(&new_TMultiModel);
      instance.SetDelete(&delete_TMultiModel);
      instance.SetDeleteArray(&deleteArray_TMultiModel);
      instance.SetDestructor(&destruct_TMultiModel);
      instance.SetStreamerFunc(&streamer_TMultiModel);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TMultiModel*)
   {
      return GenerateInitInstanceLocal((::TMultiModel*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TMultiModel*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr TMultiModel::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TMultiModel::Class_Name()
{
   return "TMultiModel";
}

//______________________________________________________________________________
const char *TMultiModel::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TMultiModel*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TMultiModel::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TMultiModel*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TMultiModel::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TMultiModel*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TMultiModel::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TMultiModel*)0x0)->GetClass(); }
   return fgIsA;
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TMultiModel(void *p) {
      return  p ? new(p) ::TMultiModel( (TRootIOCtor *)nullptr ) : new ::TMultiModel( (TRootIOCtor *)nullptr );
   }
   // Wrapper around operator delete
   static void delete_TMultiModel(void *p) {
      delete ((::TMultiModel*)p);
   }
   static void deleteArray_TMultiModel(void *p) {
      delete [] ((::TMultiModel*)p);
   }
   static void destruct_TMultiModel(void *p) {
      typedef ::TMultiModel current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_TMultiModel(TBuffer &buf, void *obj) {
      ((::TMultiModel*)obj)->::TMultiModel::Streamer(buf);
   }
} // end of namespace ROOT for class ::TMultiModel

namespace ROOT {
   static TClass *maplEintcOTStrangeCalcmUgR_Dictionary();
   static void maplEintcOTStrangeCalcmUgR_TClassManip(TClass*);
   static void *new_maplEintcOTStrangeCalcmUgR(void *p = 0);
   static void *newArray_maplEintcOTStrangeCalcmUgR(Long_t size, void *p);
   static void delete_maplEintcOTStrangeCalcmUgR(void *p);
   static void deleteArray_maplEintcOTStrangeCalcmUgR(void *p);
   static void destruct_maplEintcOTStrangeCalcmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<int,TStrangeCalc*>*)
   {
      map<int,TStrangeCalc*> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<int,TStrangeCalc*>));
      static ::ROOT::TGenericClassInfo 
         instance("map<int,TStrangeCalc*>", -2, "map", 96,
                  typeid(map<int,TStrangeCalc*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEintcOTStrangeCalcmUgR_Dictionary, isa_proxy, 0,
                  sizeof(map<int,TStrangeCalc*>) );
      instance.SetNew(&new_maplEintcOTStrangeCalcmUgR);
      instance.SetNewArray(&newArray_maplEintcOTStrangeCalcmUgR);
      instance.SetDelete(&delete_maplEintcOTStrangeCalcmUgR);
      instance.SetDeleteArray(&deleteArray_maplEintcOTStrangeCalcmUgR);
      instance.SetDestructor(&destruct_maplEintcOTStrangeCalcmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<int,TStrangeCalc*> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const map<int,TStrangeCalc*>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEintcOTStrangeCalcmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<int,TStrangeCalc*>*)0x0)->GetClass();
      maplEintcOTStrangeCalcmUgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEintcOTStrangeCalcmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEintcOTStrangeCalcmUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,TStrangeCalc*> : new map<int,TStrangeCalc*>;
   }
   static void *newArray_maplEintcOTStrangeCalcmUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,TStrangeCalc*>[nElements] : new map<int,TStrangeCalc*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEintcOTStrangeCalcmUgR(void *p) {
      delete ((map<int,TStrangeCalc*>*)p);
   }
   static void deleteArray_maplEintcOTStrangeCalcmUgR(void *p) {
      delete [] ((map<int,TStrangeCalc*>*)p);
   }
   static void destruct_maplEintcOTStrangeCalcmUgR(void *p) {
      typedef map<int,TStrangeCalc*> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<int,TStrangeCalc*>

namespace {
  void TriggerDictionaryInitialization_TMultiModelDict_Impl() {
    static const char* headers[] = {
"TMultiModel.h",
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
#line 1 "TMultiModelDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate(R"ATTRDUMP(A set of TStrangeModel's)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(A set of TStrangeModel's)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TMultiModel.h")))  TMultiModel;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "TMultiModelDict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "TMultiModel.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"TMultiModel", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("TMultiModelDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_TMultiModelDict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_TMultiModelDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_TMultiModelDict() {
  TriggerDictionaryInitialization_TMultiModelDict_Impl();
}
