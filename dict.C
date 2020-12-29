// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
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

// The generated code does not explicitly qualifies STL entities
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "RooDoubleCBFast.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_RooDoubleCBFast(void *p = 0);
   static void *newArray_RooDoubleCBFast(Long_t size, void *p);
   static void delete_RooDoubleCBFast(void *p);
   static void deleteArray_RooDoubleCBFast(void *p);
   static void destruct_RooDoubleCBFast(void *p);
   static void streamer_RooDoubleCBFast(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::RooDoubleCBFast*)
   {
      ::RooDoubleCBFast *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::RooDoubleCBFast >(0);
      static ::ROOT::TGenericClassInfo 
         instance("RooDoubleCBFast", ::RooDoubleCBFast::Class_Version(), "RooDoubleCBFast.h", 8,
                  typeid(::RooDoubleCBFast), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::RooDoubleCBFast::Dictionary, isa_proxy, 16,
                  sizeof(::RooDoubleCBFast) );
      instance.SetNew(&new_RooDoubleCBFast);
      instance.SetNewArray(&newArray_RooDoubleCBFast);
      instance.SetDelete(&delete_RooDoubleCBFast);
      instance.SetDeleteArray(&deleteArray_RooDoubleCBFast);
      instance.SetDestructor(&destruct_RooDoubleCBFast);
      instance.SetStreamerFunc(&streamer_RooDoubleCBFast);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::RooDoubleCBFast*)
   {
      return GenerateInitInstanceLocal((::RooDoubleCBFast*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::RooDoubleCBFast*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr RooDoubleCBFast::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *RooDoubleCBFast::Class_Name()
{
   return "RooDoubleCBFast";
}

//______________________________________________________________________________
const char *RooDoubleCBFast::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RooDoubleCBFast*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int RooDoubleCBFast::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::RooDoubleCBFast*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *RooDoubleCBFast::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RooDoubleCBFast*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *RooDoubleCBFast::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::RooDoubleCBFast*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void RooDoubleCBFast::Streamer(TBuffer &R__b)
{
   // Stream an object of class RooDoubleCBFast.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      RooAbsPdf::Streamer(R__b);
      x.Streamer(R__b);
      mean.Streamer(R__b);
      width.Streamer(R__b);
      alpha1.Streamer(R__b);
      n1.Streamer(R__b);
      alpha2.Streamer(R__b);
      n2.Streamer(R__b);
      R__b.CheckByteCount(R__s, R__c, RooDoubleCBFast::IsA());
   } else {
      R__c = R__b.WriteVersion(RooDoubleCBFast::IsA(), kTRUE);
      RooAbsPdf::Streamer(R__b);
      x.Streamer(R__b);
      mean.Streamer(R__b);
      width.Streamer(R__b);
      alpha1.Streamer(R__b);
      n1.Streamer(R__b);
      alpha2.Streamer(R__b);
      n2.Streamer(R__b);
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_RooDoubleCBFast(void *p) {
      return  p ? new(p) ::RooDoubleCBFast : new ::RooDoubleCBFast;
   }
   static void *newArray_RooDoubleCBFast(Long_t nElements, void *p) {
      return p ? new(p) ::RooDoubleCBFast[nElements] : new ::RooDoubleCBFast[nElements];
   }
   // Wrapper around operator delete
   static void delete_RooDoubleCBFast(void *p) {
      delete ((::RooDoubleCBFast*)p);
   }
   static void deleteArray_RooDoubleCBFast(void *p) {
      delete [] ((::RooDoubleCBFast*)p);
   }
   static void destruct_RooDoubleCBFast(void *p) {
      typedef ::RooDoubleCBFast current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_RooDoubleCBFast(TBuffer &buf, void *obj) {
      ((::RooDoubleCBFast*)obj)->::RooDoubleCBFast::Streamer(buf);
   }
} // end of namespace ROOT for class ::RooDoubleCBFast

namespace {
  void TriggerDictionaryInitialization_dict_Impl() {
    static const char* headers[] = {
"RooDoubleCBFast.h",
0
    };
    static const char* includePaths[] = {
"/cvmfs/atlas.cern.ch/repo/sw/software/21.2/AnalysisBaseExternals/21.2.150/InstallArea/x86_64-centos7-gcc8-opt/include/",
"/afs/cern.ch/work/r/rcoelhol/jet_studies/run/analysis/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$RooDoubleCBFast.h")))  RooDoubleCBFast;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "dict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "RooDoubleCBFast.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"RooDoubleCBFast", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_dict() {
  TriggerDictionaryInitialization_dict_Impl();
}
