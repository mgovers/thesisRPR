// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME TCalcInfoDict

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
#include "TCalcInfo.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_TCalcInfo(void *p = 0);
   static void *newArray_TCalcInfo(Long_t size, void *p);
   static void delete_TCalcInfo(void *p);
   static void deleteArray_TCalcInfo(void *p);
   static void destruct_TCalcInfo(void *p);
   static void streamer_TCalcInfo(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TCalcInfo*)
   {
      ::TCalcInfo *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TCalcInfo >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TCalcInfo", ::TCalcInfo::Class_Version(), "TCalcInfo.h", 29,
                  typeid(::TCalcInfo), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TCalcInfo::Dictionary, isa_proxy, 16,
                  sizeof(::TCalcInfo) );
      instance.SetNew(&new_TCalcInfo);
      instance.SetNewArray(&newArray_TCalcInfo);
      instance.SetDelete(&delete_TCalcInfo);
      instance.SetDeleteArray(&deleteArray_TCalcInfo);
      instance.SetDestructor(&destruct_TCalcInfo);
      instance.SetStreamerFunc(&streamer_TCalcInfo);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TCalcInfo*)
   {
      return GenerateInitInstanceLocal((::TCalcInfo*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TCalcInfo*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr TCalcInfo::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TCalcInfo::Class_Name()
{
   return "TCalcInfo";
}

//______________________________________________________________________________
const char *TCalcInfo::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TCalcInfo*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TCalcInfo::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TCalcInfo*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TCalcInfo::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TCalcInfo*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TCalcInfo::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TCalcInfo*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void TCalcInfo::Streamer(TBuffer &R__b)
{
   // Stream an object of class TCalcInfo.

   TObject::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TCalcInfo(void *p) {
      return  p ? new(p) ::TCalcInfo : new ::TCalcInfo;
   }
   static void *newArray_TCalcInfo(Long_t nElements, void *p) {
      return p ? new(p) ::TCalcInfo[nElements] : new ::TCalcInfo[nElements];
   }
   // Wrapper around operator delete
   static void delete_TCalcInfo(void *p) {
      delete ((::TCalcInfo*)p);
   }
   static void deleteArray_TCalcInfo(void *p) {
      delete [] ((::TCalcInfo*)p);
   }
   static void destruct_TCalcInfo(void *p) {
      typedef ::TCalcInfo current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_TCalcInfo(TBuffer &buf, void *obj) {
      ((::TCalcInfo*)obj)->::TCalcInfo::Streamer(buf);
   }
} // end of namespace ROOT for class ::TCalcInfo

namespace {
  void TriggerDictionaryInitialization_TCalcInfoDict_Impl() {
    static const char* headers[] = {
"TCalcInfo.h",
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
#line 1 "TCalcInfoDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate(R"ATTRDUMP(Data struct wrapper)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TCalcInfo.h")))  TCalcInfo;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "TCalcInfoDict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "TCalcInfo.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"TCalcInfo", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("TCalcInfoDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_TCalcInfoDict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_TCalcInfoDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_TCalcInfoDict() {
  TriggerDictionaryInitialization_TCalcInfoDict_Impl();
}
