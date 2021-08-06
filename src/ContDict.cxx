// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME ContDict
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
#include "vector"
#include "DDRec/Vector3D.h"

// Header files passed via #pragma extra_include

namespace dd4hep {
   namespace ROOT {
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance();
      static TClass *dd4hep_Dictionary();

      // Function generating the singleton type initializer
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance()
      {
         static ::ROOT::TGenericClassInfo 
            instance("dd4hep", 0 /*version*/, "DDRec/Vector3D.h", 21,
                     ::ROOT::Internal::DefineBehavior((void*)0,(void*)0),
                     &dd4hep_Dictionary, 0);
         return &instance;
      }
      // Insure that the inline function is _not_ optimized away by the compiler
      ::ROOT::TGenericClassInfo *(*_R__UNIQUE_DICT_(InitFunctionKeeper))() = &GenerateInitInstance;  
      // Static variable to force the class initialization
      static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstance(); R__UseDummy(_R__UNIQUE_DICT_(Init));

      // Dictionary for non-ClassDef classes
      static TClass *dd4hep_Dictionary() {
         return GenerateInitInstance()->GetClass();
      }

   }
}

namespace dd4hep {
   namespace rec {
   namespace ROOT {
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance();
      static TClass *dd4hepcLcLrec_Dictionary();

      // Function generating the singleton type initializer
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance()
      {
         static ::ROOT::TGenericClassInfo 
            instance("dd4hep::rec", 0 /*version*/, "DDRec/Vector3D.h", 21,
                     ::ROOT::Internal::DefineBehavior((void*)0,(void*)0),
                     &dd4hepcLcLrec_Dictionary, 0);
         return &instance;
      }
      // Insure that the inline function is _not_ optimized away by the compiler
      ::ROOT::TGenericClassInfo *(*_R__UNIQUE_DICT_(InitFunctionKeeper))() = &GenerateInitInstance;  
      // Static variable to force the class initialization
      static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstance(); R__UseDummy(_R__UNIQUE_DICT_(Init));

      // Dictionary for non-ClassDef classes
      static TClass *dd4hepcLcLrec_Dictionary() {
         return GenerateInitInstance()->GetClass();
      }

   }
}
}

namespace ROOT {
   static TClass *dd4hepcLcLreccLcLVector3D_Dictionary();
   static void dd4hepcLcLreccLcLVector3D_TClassManip(TClass*);
   static void *new_dd4hepcLcLreccLcLVector3D(void *p = 0);
   static void *newArray_dd4hepcLcLreccLcLVector3D(Long_t size, void *p);
   static void delete_dd4hepcLcLreccLcLVector3D(void *p);
   static void deleteArray_dd4hepcLcLreccLcLVector3D(void *p);
   static void destruct_dd4hepcLcLreccLcLVector3D(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::dd4hep::rec::Vector3D*)
   {
      ::dd4hep::rec::Vector3D *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::dd4hep::rec::Vector3D));
      static ::ROOT::TGenericClassInfo 
         instance("dd4hep::rec::Vector3D", "DDRec/Vector3D.h", 32,
                  typeid(::dd4hep::rec::Vector3D), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &dd4hepcLcLreccLcLVector3D_Dictionary, isa_proxy, 4,
                  sizeof(::dd4hep::rec::Vector3D) );
      instance.SetNew(&new_dd4hepcLcLreccLcLVector3D);
      instance.SetNewArray(&newArray_dd4hepcLcLreccLcLVector3D);
      instance.SetDelete(&delete_dd4hepcLcLreccLcLVector3D);
      instance.SetDeleteArray(&deleteArray_dd4hepcLcLreccLcLVector3D);
      instance.SetDestructor(&destruct_dd4hepcLcLreccLcLVector3D);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::dd4hep::rec::Vector3D*)
   {
      return GenerateInitInstanceLocal((::dd4hep::rec::Vector3D*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::dd4hep::rec::Vector3D*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *dd4hepcLcLreccLcLVector3D_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::dd4hep::rec::Vector3D*)0x0)->GetClass();
      dd4hepcLcLreccLcLVector3D_TClassManip(theClass);
   return theClass;
   }

   static void dd4hepcLcLreccLcLVector3D_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_dd4hepcLcLreccLcLVector3D(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::dd4hep::rec::Vector3D : new ::dd4hep::rec::Vector3D;
   }
   static void *newArray_dd4hepcLcLreccLcLVector3D(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::dd4hep::rec::Vector3D[nElements] : new ::dd4hep::rec::Vector3D[nElements];
   }
   // Wrapper around operator delete
   static void delete_dd4hepcLcLreccLcLVector3D(void *p) {
      delete ((::dd4hep::rec::Vector3D*)p);
   }
   static void deleteArray_dd4hepcLcLreccLcLVector3D(void *p) {
      delete [] ((::dd4hep::rec::Vector3D*)p);
   }
   static void destruct_dd4hepcLcLreccLcLVector3D(void *p) {
      typedef ::dd4hep::rec::Vector3D current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::dd4hep::rec::Vector3D

namespace ROOT {
   static TClass *vectorlEdd4hepcLcLreccLcLVector3DgR_Dictionary();
   static void vectorlEdd4hepcLcLreccLcLVector3DgR_TClassManip(TClass*);
   static void *new_vectorlEdd4hepcLcLreccLcLVector3DgR(void *p = 0);
   static void *newArray_vectorlEdd4hepcLcLreccLcLVector3DgR(Long_t size, void *p);
   static void delete_vectorlEdd4hepcLcLreccLcLVector3DgR(void *p);
   static void deleteArray_vectorlEdd4hepcLcLreccLcLVector3DgR(void *p);
   static void destruct_vectorlEdd4hepcLcLreccLcLVector3DgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<dd4hep::rec::Vector3D>*)
   {
      vector<dd4hep::rec::Vector3D> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<dd4hep::rec::Vector3D>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<dd4hep::rec::Vector3D>", -2, "vector", 339,
                  typeid(vector<dd4hep::rec::Vector3D>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEdd4hepcLcLreccLcLVector3DgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<dd4hep::rec::Vector3D>) );
      instance.SetNew(&new_vectorlEdd4hepcLcLreccLcLVector3DgR);
      instance.SetNewArray(&newArray_vectorlEdd4hepcLcLreccLcLVector3DgR);
      instance.SetDelete(&delete_vectorlEdd4hepcLcLreccLcLVector3DgR);
      instance.SetDeleteArray(&deleteArray_vectorlEdd4hepcLcLreccLcLVector3DgR);
      instance.SetDestructor(&destruct_vectorlEdd4hepcLcLreccLcLVector3DgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<dd4hep::rec::Vector3D> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<dd4hep::rec::Vector3D>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEdd4hepcLcLreccLcLVector3DgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<dd4hep::rec::Vector3D>*)0x0)->GetClass();
      vectorlEdd4hepcLcLreccLcLVector3DgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEdd4hepcLcLreccLcLVector3DgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEdd4hepcLcLreccLcLVector3DgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<dd4hep::rec::Vector3D> : new vector<dd4hep::rec::Vector3D>;
   }
   static void *newArray_vectorlEdd4hepcLcLreccLcLVector3DgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<dd4hep::rec::Vector3D>[nElements] : new vector<dd4hep::rec::Vector3D>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEdd4hepcLcLreccLcLVector3DgR(void *p) {
      delete ((vector<dd4hep::rec::Vector3D>*)p);
   }
   static void deleteArray_vectorlEdd4hepcLcLreccLcLVector3DgR(void *p) {
      delete [] ((vector<dd4hep::rec::Vector3D>*)p);
   }
   static void destruct_vectorlEdd4hepcLcLreccLcLVector3DgR(void *p) {
      typedef vector<dd4hep::rec::Vector3D> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<dd4hep::rec::Vector3D>

namespace {
  void TriggerDictionaryInitialization_ContDict_Impl() {
    static const char* headers[] = {
"vector",
"DDRec/Vector3D.h",
0
    };
    static const char* includePaths[] = {
"../include",
"/cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/root/6.18.04/include",
"/afs/desy.de/user/d/dudarboh/TOFAnalysis/src/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "ContDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
namespace dd4hep{namespace rec{class __attribute__((annotate("$clingAutoload$DDRec/Vector3D.h")))  Vector3D;}}
namespace std{template <typename _Tp> class __attribute__((annotate("$clingAutoload$bits/allocator.h")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "ContDict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "vector"
#include "DDRec/Vector3D.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"dd4hep::rec::Vector3D", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("ContDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_ContDict_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_ContDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_ContDict() {
  TriggerDictionaryInitialization_ContDict_Impl();
}
