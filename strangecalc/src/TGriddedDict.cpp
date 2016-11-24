// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME TGriddedDict

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
#include "TGridded.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *TGriddedlETStrangeModelgR_Dictionary();
   static void TGriddedlETStrangeModelgR_TClassManip(TClass*);
   static void *new_TGriddedlETStrangeModelgR(void *p = 0);
   static void delete_TGriddedlETStrangeModelgR(void *p);
   static void deleteArray_TGriddedlETStrangeModelgR(void *p);
   static void destruct_TGriddedlETStrangeModelgR(void *p);
   static void streamer_TGriddedlETStrangeModelgR(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TGridded<TStrangeModel>*)
   {
      ::TGridded<TStrangeModel> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TGridded<TStrangeModel> >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TGridded<TStrangeModel>", ::TGridded<TStrangeModel>::Class_Version(), "TGridded.h", 39,
                  typeid(::TGridded<TStrangeModel>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &TGriddedlETStrangeModelgR_Dictionary, isa_proxy, 16,
                  sizeof(::TGridded<TStrangeModel>) );
      instance.SetNew(&new_TGriddedlETStrangeModelgR);
      instance.SetDelete(&delete_TGriddedlETStrangeModelgR);
      instance.SetDeleteArray(&deleteArray_TGriddedlETStrangeModelgR);
      instance.SetDestructor(&destruct_TGriddedlETStrangeModelgR);
      instance.SetStreamerFunc(&streamer_TGriddedlETStrangeModelgR);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TGridded<TStrangeModel>*)
   {
      return GenerateInitInstanceLocal((::TGridded<TStrangeModel>*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TGridded<TStrangeModel>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *TGriddedlETStrangeModelgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::TGridded<TStrangeModel>*)0x0)->GetClass();
      TGriddedlETStrangeModelgR_TClassManip(theClass);
   return theClass;
   }

   static void TGriddedlETStrangeModelgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *TGriddedlETMultiModelgR_Dictionary();
   static void TGriddedlETMultiModelgR_TClassManip(TClass*);
   static void *new_TGriddedlETMultiModelgR(void *p = 0);
   static void delete_TGriddedlETMultiModelgR(void *p);
   static void deleteArray_TGriddedlETMultiModelgR(void *p);
   static void destruct_TGriddedlETMultiModelgR(void *p);
   static void streamer_TGriddedlETMultiModelgR(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TGridded<TMultiModel>*)
   {
      ::TGridded<TMultiModel> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TGridded<TMultiModel> >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TGridded<TMultiModel>", ::TGridded<TMultiModel>::Class_Version(), "TGridded.h", 39,
                  typeid(::TGridded<TMultiModel>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &TGriddedlETMultiModelgR_Dictionary, isa_proxy, 16,
                  sizeof(::TGridded<TMultiModel>) );
      instance.SetNew(&new_TGriddedlETMultiModelgR);
      instance.SetDelete(&delete_TGriddedlETMultiModelgR);
      instance.SetDeleteArray(&deleteArray_TGriddedlETMultiModelgR);
      instance.SetDestructor(&destruct_TGriddedlETMultiModelgR);
      instance.SetStreamerFunc(&streamer_TGriddedlETMultiModelgR);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TGridded<TMultiModel>*)
   {
      return GenerateInitInstanceLocal((::TGridded<TMultiModel>*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::TGridded<TMultiModel>*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *TGriddedlETMultiModelgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::TGridded<TMultiModel>*)0x0)->GetClass();
      TGriddedlETMultiModelgR_TClassManip(theClass);
   return theClass;
   }

   static void TGriddedlETMultiModelgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

//______________________________________________________________________________
template <> atomic_TClass_ptr TGridded<TStrangeModel>::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
template <> const char *TGridded<TStrangeModel>::Class_Name()
{
   return "TGridded<TStrangeModel>";
}

//______________________________________________________________________________
template <> const char *TGridded<TStrangeModel>::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TGridded<TStrangeModel>*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
template <> int TGridded<TStrangeModel>::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TGridded<TStrangeModel>*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
template <> TClass *TGridded<TStrangeModel>::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TGridded<TStrangeModel>*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
template <> TClass *TGridded<TStrangeModel>::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TGridded<TStrangeModel>*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
template <> atomic_TClass_ptr TGridded<TMultiModel>::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
template <> const char *TGridded<TMultiModel>::Class_Name()
{
   return "TGridded<TMultiModel>";
}

//______________________________________________________________________________
template <> const char *TGridded<TMultiModel>::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TGridded<TMultiModel>*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
template <> int TGridded<TMultiModel>::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TGridded<TMultiModel>*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
template <> TClass *TGridded<TMultiModel>::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TGridded<TMultiModel>*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
template <> TClass *TGridded<TMultiModel>::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TGridded<TMultiModel>*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
template <> void TGridded<TStrangeModel>::Streamer(TBuffer &R__b)
{
   // Stream an object of class TGridded<TStrangeModel>.

   TStrangeModel::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TGriddedlETStrangeModelgR(void *p) {
      return  p ? new(p) ::TGridded<TStrangeModel>( (TRootIOCtor *)nullptr ) : new ::TGridded<TStrangeModel>( (TRootIOCtor *)nullptr );
   }
   // Wrapper around operator delete
   static void delete_TGriddedlETStrangeModelgR(void *p) {
      delete ((::TGridded<TStrangeModel>*)p);
   }
   static void deleteArray_TGriddedlETStrangeModelgR(void *p) {
      delete [] ((::TGridded<TStrangeModel>*)p);
   }
   static void destruct_TGriddedlETStrangeModelgR(void *p) {
      typedef ::TGridded<TStrangeModel> current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_TGriddedlETStrangeModelgR(TBuffer &buf, void *obj) {
      ((::TGridded<TStrangeModel>*)obj)->::TGridded<TStrangeModel>::Streamer(buf);
   }
} // end of namespace ROOT for class ::TGridded<TStrangeModel>

//______________________________________________________________________________
template <> void TGridded<TMultiModel>::Streamer(TBuffer &R__b)
{
   // Stream an object of class TGridded<TMultiModel>.

   TMultiModel::Streamer(R__b);
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TGriddedlETMultiModelgR(void *p) {
      return  p ? new(p) ::TGridded<TMultiModel>( (TRootIOCtor *)nullptr ) : new ::TGridded<TMultiModel>( (TRootIOCtor *)nullptr );
   }
   // Wrapper around operator delete
   static void delete_TGriddedlETMultiModelgR(void *p) {
      delete ((::TGridded<TMultiModel>*)p);
   }
   static void deleteArray_TGriddedlETMultiModelgR(void *p) {
      delete [] ((::TGridded<TMultiModel>*)p);
   }
   static void destruct_TGriddedlETMultiModelgR(void *p) {
      typedef ::TGridded<TMultiModel> current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_TGriddedlETMultiModelgR(TBuffer &buf, void *obj) {
      ((::TGridded<TMultiModel>*)obj)->::TGridded<TMultiModel>::Streamer(buf);
   }
} // end of namespace ROOT for class ::TGridded<TMultiModel>

namespace {
  void TriggerDictionaryInitialization_TGriddedDict_Impl() {
    static const char* headers[] = {
"TGridded.h",
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
#line 1 "TGriddedDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$TGridded.h")))  TStrangeModel;
template <typename Model> class __attribute__((annotate("$clingAutoload$TGridded.h")))  TGridded;

class __attribute__((annotate("$clingAutoload$TGridded.h")))  TMultiModel;
typedef TGridded<TStrangeModel> TGriddedStrangeModel __attribute__((annotate("$clingAutoload$TGridded.h"))) ;
typedef TGridded<TMultiModel> TGriddedMultiModel __attribute__((annotate("$clingAutoload$TGridded.h"))) ;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "TGriddedDict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "TGridded.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"TGridded<TMultiModel>", payloadCode, "@",
"TGridded<TStrangeModel>", payloadCode, "@",
"TGriddedMultiModel", payloadCode, "@",
"TGriddedStrangeModel", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("TGriddedDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_TGriddedDict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_TGriddedDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_TGriddedDict() {
  TriggerDictionaryInitialization_TGriddedDict_Impl();
}
