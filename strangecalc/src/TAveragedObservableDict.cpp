// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME TAveragedObservableDict

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
#include "TAveragedObservable.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_TAveragedObservable(void *p = 0);
   static void delete_TAveragedObservable(void *p);
   static void deleteArray_TAveragedObservable(void *p);
   static void destruct_TAveragedObservable(void *p);
   static void streamer_TAveragedObservable(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TAveragedObservable*)
   {
      ::TAveragedObservable *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TAveragedObservable >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TAveragedObservable", ::TAveragedObservable::Class_Version(), "TAveragedObservable.h", 39,
                  typeid(::TAveragedObservable), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TAveragedObservable::Dictionary, isa_proxy, 16,
                  sizeof(::TAveragedObservable) );
      instance.SetNew(&new_TAveragedObservable);
      instance.SetDelete(&delete_TAveragedObservable);
      instance.SetDeleteArray(&deleteArray_TAveragedObservable);
      instance.SetDestructor(&destruct_TAveragedObservable);
      instance.SetStreamerFunc(&streamer_TAveragedObservable);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TAveragedObservable*)
   {
      return GenerateInitInstanceLocal((::TAveragedObservable*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TAveragedObservable*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void delete_TAveragedObservablecLcLIntegrand(void *p);
   static void deleteArray_TAveragedObservablecLcLIntegrand(void *p);
   static void destruct_TAveragedObservablecLcLIntegrand(void *p);
   static void streamer_TAveragedObservablecLcLIntegrand(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TAveragedObservable::Integrand*)
   {
      ::TAveragedObservable::Integrand *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TAveragedObservable::Integrand >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TAveragedObservable::Integrand", ::TAveragedObservable::Integrand::Class_Version(), "TAveragedObservable.h", 101,
                  typeid(::TAveragedObservable::Integrand), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TAveragedObservable::Integrand::Dictionary, isa_proxy, 16,
                  sizeof(::TAveragedObservable::Integrand) );
      instance.SetDelete(&delete_TAveragedObservablecLcLIntegrand);
      instance.SetDeleteArray(&deleteArray_TAveragedObservablecLcLIntegrand);
      instance.SetDestructor(&destruct_TAveragedObservablecLcLIntegrand);
      instance.SetStreamerFunc(&streamer_TAveragedObservablecLcLIntegrand);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TAveragedObservable::Integrand*)
   {
      return GenerateInitInstanceLocal((::TAveragedObservable::Integrand*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TAveragedObservable::Integrand*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr TAveragedObservable::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TAveragedObservable::Class_Name()
{
   return "TAveragedObservable";
}

//______________________________________________________________________________
const char *TAveragedObservable::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TAveragedObservable*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TAveragedObservable::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TAveragedObservable*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TAveragedObservable::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TAveragedObservable*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TAveragedObservable::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TAveragedObservable*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TAveragedObservable::Integrand::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TAveragedObservable::Integrand::Class_Name()
{
   return "TAveragedObservable::Integrand";
}

//______________________________________________________________________________
const char *TAveragedObservable::Integrand::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TAveragedObservable::Integrand*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TAveragedObservable::Integrand::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TAveragedObservable::Integrand*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TAveragedObservable::Integrand::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TAveragedObservable::Integrand*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TAveragedObservable::Integrand::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TAveragedObservable::Integrand*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void TAveragedObservable::Streamer(TBuffer &R__b)
{
   // Stream an object of class TAveragedObservable.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TObject::Streamer(R__b);
      R__b.ReadStaticArray((Int_t*)fVarState);
      R__b.ReadStaticArray((double*)fLowerBounds);
      R__b.ReadStaticArray((double*)fUpperBounds);
      R__b >> fCenteredKinematics;
      R__b >> fCalcInfo;
      R__b.CheckByteCount(R__s, R__c, TAveragedObservable::IsA());
   } else {
      R__c = R__b.WriteVersion(TAveragedObservable::IsA(), kTRUE);
      TObject::Streamer(R__b);
      R__b.WriteArray((Int_t*)fVarState, 3);
      R__b.WriteArray(fLowerBounds, 3);
      R__b.WriteArray(fUpperBounds, 3);
      R__b << fCenteredKinematics;
      R__b << fCalcInfo;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TAveragedObservable(void *p) {
      return  p ? new(p) ::TAveragedObservable( (TRootIOCtor *)nullptr ) : new ::TAveragedObservable( (TRootIOCtor *)nullptr );
   }
   // Wrapper around operator delete
   static void delete_TAveragedObservable(void *p) {
      delete ((::TAveragedObservable*)p);
   }
   static void deleteArray_TAveragedObservable(void *p) {
      delete [] ((::TAveragedObservable*)p);
   }
   static void destruct_TAveragedObservable(void *p) {
      typedef ::TAveragedObservable current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_TAveragedObservable(TBuffer &buf, void *obj) {
      ((::TAveragedObservable*)obj)->::TAveragedObservable::Streamer(buf);
   }
} // end of namespace ROOT for class ::TAveragedObservable

//______________________________________________________________________________
void TAveragedObservable::Integrand::Streamer(TBuffer &R__b)
{
   // Stream an object of class TAveragedObservable::Integrand.

   ::Error("TAveragedObservable::Integrand::Streamer", "version id <=0 in ClassDef, dummy Streamer() called"); if (R__b.IsReading()) { }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_TAveragedObservablecLcLIntegrand(void *p) {
      delete ((::TAveragedObservable::Integrand*)p);
   }
   static void deleteArray_TAveragedObservablecLcLIntegrand(void *p) {
      delete [] ((::TAveragedObservable::Integrand*)p);
   }
   static void destruct_TAveragedObservablecLcLIntegrand(void *p) {
      typedef ::TAveragedObservable::Integrand current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_TAveragedObservablecLcLIntegrand(TBuffer &buf, void *obj) {
      ((::TAveragedObservable::Integrand*)obj)->::TAveragedObservable::Integrand::Streamer(buf);
   }
} // end of namespace ROOT for class ::TAveragedObservable::Integrand

namespace {
  void TriggerDictionaryInitialization_TAveragedObservableDict_Impl() {
    static const char* headers[] = {
"TAveragedObservable.h",
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
#line 1 "TAveragedObservableDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate(R"ATTRDUMP(Observables averaged over binned kinematic variables)ATTRDUMP"))) __attribute__((annotate(R"ATTRDUMP(Observables averaged over binned kinematic variables)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TAveragedObservable.h")))  TAveragedObservable;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "TAveragedObservableDict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "TAveragedObservable.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"TAveragedObservable", payloadCode, "@",
"TAveragedObservable::Integrand", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("TAveragedObservableDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_TAveragedObservableDict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_TAveragedObservableDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_TAveragedObservableDict() {
  TriggerDictionaryInitialization_TAveragedObservableDict_Impl();
}
