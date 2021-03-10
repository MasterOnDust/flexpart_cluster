/* File: trajectory_clusteringmodule.c
 * This file is auto-generated with f2py (version:2).
 * f2py is a Fortran to Python Interface Generator (FPIG), Second Edition,
 * written by Pearu Peterson <pearu@cens.ioc.ee>.
 * Generation date: Wed Mar 10 14:52:22 2021
 * Do not edit this file directly unless you know what you are doing!!!
 */

#ifdef __cplusplus
extern "C" {
#endif

/*********************** See f2py2e/cfuncs.py: includes ***********************/
#include "Python.h"
#include <stdarg.h>
#include "fortranobject.h"
#include <math.h>

/**************** See f2py2e/rules.py: mod_rules['modulebody'] ****************/
static PyObject *trajectory_clustering_error;
static PyObject *trajectory_clustering_module;

/*********************** See f2py2e/cfuncs.py: typedefs ***********************/
/*need_typedefs*/

/****************** See f2py2e/cfuncs.py: typedefs_generated ******************/
/*need_typedefs_generated*/

/********************** See f2py2e/cfuncs.py: cppmacros **********************/
#if defined(PREPEND_FORTRAN)
#if defined(NO_APPEND_FORTRAN)
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) _##F
#else
#define F_FUNC(f,F) _##f
#endif
#else
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) _##F##_
#else
#define F_FUNC(f,F) _##f##_
#endif
#endif
#else
#if defined(NO_APPEND_FORTRAN)
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) F
#else
#define F_FUNC(f,F) f
#endif
#else
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) F##_
#else
#define F_FUNC(f,F) f##_
#endif
#endif
#endif
#if defined(UNDERSCORE_G77)
#define F_FUNC_US(f,F) F_FUNC(f##_,F##_)
#else
#define F_FUNC_US(f,F) F_FUNC(f,F)
#endif

#define rank(var) var ## _Rank
#define shape(var,dim) var ## _Dims[dim]
#define old_rank(var) (PyArray_NDIM((PyArrayObject *)(capi_ ## var ## _tmp)))
#define old_shape(var,dim) PyArray_DIM(((PyArrayObject *)(capi_ ## var ## _tmp)),dim)
#define fshape(var,dim) shape(var,rank(var)-dim-1)
#define len(var) shape(var,0)
#define flen(var) fshape(var,0)
#define old_size(var) PyArray_SIZE((PyArrayObject *)(capi_ ## var ## _tmp))
/* #define index(i) capi_i ## i */
#define slen(var) capi_ ## var ## _len
#define size(var, ...) f2py_size((PyArrayObject *)(capi_ ## var ## _tmp), ## __VA_ARGS__, -1)

#define CHECKSCALAR(check,tcheck,name,show,var)\
    if (!(check)) {\
        char errstring[256];\
        sprintf(errstring, "%s: "show, "("tcheck") failed for "name, var);\
        PyErr_SetString(trajectory_clustering_error,errstring);\
        /*goto capi_fail;*/\
    } else 
#ifdef DEBUGCFUNCS
#define CFUNCSMESS(mess) fprintf(stderr,"debug-capi:"mess);
#define CFUNCSMESSPY(mess,obj) CFUNCSMESS(mess) \
    PyObject_Print((PyObject *)obj,stderr,Py_PRINT_RAW);\
    fprintf(stderr,"\n");
#else
#define CFUNCSMESS(mess)
#define CFUNCSMESSPY(mess,obj)
#endif

#ifndef max
#define max(a,b) ((a > b) ? (a) : (b))
#endif
#ifndef min
#define min(a,b) ((a < b) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a,b) ((a > b) ? (a) : (b))
#endif
#ifndef MIN
#define MIN(a,b) ((a < b) ? (a) : (b))
#endif

#if defined(PREPEND_FORTRAN)
#if defined(NO_APPEND_FORTRAN)
#if defined(UPPERCASE_FORTRAN)
#define F_WRAPPEDFUNC(f,F) _F2PYWRAP##F
#else
#define F_WRAPPEDFUNC(f,F) _f2pywrap##f
#endif
#else
#if defined(UPPERCASE_FORTRAN)
#define F_WRAPPEDFUNC(f,F) _F2PYWRAP##F##_
#else
#define F_WRAPPEDFUNC(f,F) _f2pywrap##f##_
#endif
#endif
#else
#if defined(NO_APPEND_FORTRAN)
#if defined(UPPERCASE_FORTRAN)
#define F_WRAPPEDFUNC(f,F) F2PYWRAP##F
#else
#define F_WRAPPEDFUNC(f,F) f2pywrap##f
#endif
#else
#if defined(UPPERCASE_FORTRAN)
#define F_WRAPPEDFUNC(f,F) F2PYWRAP##F##_
#else
#define F_WRAPPEDFUNC(f,F) f2pywrap##f##_
#endif
#endif
#endif
#if defined(UNDERSCORE_G77)
#define F_WRAPPEDFUNC_US(f,F) F_WRAPPEDFUNC(f##_,F##_)
#else
#define F_WRAPPEDFUNC_US(f,F) F_WRAPPEDFUNC(f,F)
#endif


/************************ See f2py2e/cfuncs.py: cfuncs ************************/
static int
double_from_pyobj(double* v, PyObject *obj, const char *errmess)
{
    PyObject* tmp = NULL;
    if (PyFloat_Check(obj)) {
        *v = PyFloat_AsDouble(obj);
        return !(*v == -1.0 && PyErr_Occurred());
    }

    tmp = PyNumber_Float(obj);
    if (tmp) {
        *v = PyFloat_AsDouble(tmp);
        Py_DECREF(tmp);
        return !(*v == -1.0 && PyErr_Occurred());
    }
    if (PyComplex_Check(obj))
        tmp = PyObject_GetAttrString(obj,"real");
    else if (PyBytes_Check(obj) || PyUnicode_Check(obj))
        /*pass*/;
    else if (PySequence_Check(obj))
        tmp = PySequence_GetItem(obj,0);
    if (tmp) {
        PyErr_Clear();
        if (double_from_pyobj(v,tmp,errmess)) {Py_DECREF(tmp); return 1;}
        Py_DECREF(tmp);
    }
    {
        PyObject* err = PyErr_Occurred();
        if (err==NULL) err = trajectory_clustering_error;
        PyErr_SetString(err,errmess);
    }
    return 0;
}

static int f2py_size(PyArrayObject* var, ...)
{
  npy_int sz = 0;
  npy_int dim;
  npy_int rank;
  va_list argp;
  va_start(argp, var);
  dim = va_arg(argp, npy_int);
  if (dim==-1)
    {
      sz = PyArray_SIZE(var);
    }
  else
    {
      rank = PyArray_NDIM(var);
      if (dim>=1 && dim<=rank)
        sz = PyArray_DIM(var, dim-1);
      else
        fprintf(stderr, "f2py_size: 2nd argument value=%d fails to satisfy 1<=value<=%d. Result will be 0.\n", dim, rank);
    }
  va_end(argp);
  return sz;
}

static int
int_from_pyobj(int* v, PyObject *obj, const char *errmess)
{
    PyObject* tmp = NULL;

    if (PyLong_Check(obj)) {
        *v = Npy__PyLong_AsInt(obj);
        return !(*v == -1 && PyErr_Occurred());
    }

    tmp = PyNumber_Long(obj);
    if (tmp) {
        *v = Npy__PyLong_AsInt(tmp);
        Py_DECREF(tmp);
        return !(*v == -1 && PyErr_Occurred());
    }

    if (PyComplex_Check(obj))
        tmp = PyObject_GetAttrString(obj,"real");
    else if (PyBytes_Check(obj) || PyUnicode_Check(obj))
        /*pass*/;
    else if (PySequence_Check(obj))
        tmp = PySequence_GetItem(obj, 0);
    if (tmp) {
        PyErr_Clear();
        if (int_from_pyobj(v, tmp, errmess)) {
            Py_DECREF(tmp);
            return 1;
        }
        Py_DECREF(tmp);
    }
    {
        PyObject* err = PyErr_Occurred();
        if (err == NULL) {
            err = trajectory_clustering_error;
        }
        PyErr_SetString(err, errmess);
    }
    return 0;
}

static int
float_from_pyobj(float* v, PyObject *obj, const char *errmess)
{
    double d=0.0;
    if (double_from_pyobj(&d,obj,errmess)) {
        *v = (float)d;
        return 1;
    }
    return 0;
}


/********************* See f2py2e/cfuncs.py: userincludes *********************/
/*need_userincludes*/

/********************* See f2py2e/capi_rules.py: usercode *********************/


/* See f2py2e/rules.py */
extern void F_FUNC(clustering,CLUSTERING)(float*,float*,float*,int*,int*,int*,int*,float*,float*,float*,float*,float*,float*,float*,int*);
extern void F_WRAPPEDFUNC(distance2,DISTANCE2)(float*,float*,float*,float*,float*);
/*eof externroutines*/

/******************** See f2py2e/capi_rules.py: usercode1 ********************/


/******************* See f2py2e/cb_rules.py: buildcallback *******************/
/*need_callbacks*/

/*********************** See f2py2e/rules.py: buildapi ***********************/

/********************************* clustering *********************************/
static char doc_f2py_rout_trajectory_clustering_clustering[] = "\
xclust,yclust,zclust,fclust,rms,rmsclust,zrms,iterations = clustering(xl,yl,zl,[n,niters,ntime,ncluster])\n\nWrapper for ``clustering``.\
\n\nParameters\n----------\n"
"xl : input rank-2 array('f') with bounds (ntime,n)\n"
"yl : input rank-2 array('f') with bounds (ntime,n)\n"
"zl : input rank-2 array('f') with bounds (ntime,n)\n"
"\nOther Parameters\n----------------\n"
"n : input int, optional\n    Default: shape(xl,1)\n"
"niters : input int, optional\n    Default: 100\n"
"ntime : input int, optional\n    Default: shape(xl,0)\n"
"ncluster : input int, optional\n    Default: 5\n"
"\nReturns\n-------\n"
"xclust : rank-2 array('f') with bounds (ntime,ncluster)\n"
"yclust : rank-2 array('f') with bounds (ntime,ncluster)\n"
"zclust : rank-2 array('f') with bounds (ntime,ncluster)\n"
"fclust : rank-1 array('f') with bounds (ncluster)\n"
"rms : float\n"
"rmsclust : rank-1 array('f') with bounds (ncluster)\n"
"zrms : float\n"
"iterations : int";
/* extern void F_FUNC(clustering,CLUSTERING)(float*,float*,float*,int*,int*,int*,int*,float*,float*,float*,float*,float*,float*,float*,int*); */
static PyObject *f2py_rout_trajectory_clustering_clustering(const PyObject *capi_self,
                           PyObject *capi_args,
                           PyObject *capi_keywds,
                           void (*f2py_func)(float*,float*,float*,int*,int*,int*,int*,float*,float*,float*,float*,float*,float*,float*,int*)) {
    PyObject * volatile capi_buildvalue = NULL;
    volatile int f2py_success = 1;
/*decl*/

  float *xl = NULL;
  npy_intp xl_Dims[2] = {-1, -1};
  const int xl_Rank = 2;
  PyArrayObject *capi_xl_tmp = NULL;
  int capi_xl_intent = 0;
  PyObject *xl_capi = Py_None;
  float *yl = NULL;
  npy_intp yl_Dims[2] = {-1, -1};
  const int yl_Rank = 2;
  PyArrayObject *capi_yl_tmp = NULL;
  int capi_yl_intent = 0;
  PyObject *yl_capi = Py_None;
  float *zl = NULL;
  npy_intp zl_Dims[2] = {-1, -1};
  const int zl_Rank = 2;
  PyArrayObject *capi_zl_tmp = NULL;
  int capi_zl_intent = 0;
  PyObject *zl_capi = Py_None;
  int n = 0;
  PyObject *n_capi = Py_None;
  int niters = 0;
  PyObject *niters_capi = Py_None;
  int ntime = 0;
  PyObject *ntime_capi = Py_None;
  int ncluster = 0;
  PyObject *ncluster_capi = Py_None;
  float *xclust = NULL;
  npy_intp xclust_Dims[2] = {-1, -1};
  const int xclust_Rank = 2;
  PyArrayObject *capi_xclust_tmp = NULL;
  int capi_xclust_intent = 0;
  float *yclust = NULL;
  npy_intp yclust_Dims[2] = {-1, -1};
  const int yclust_Rank = 2;
  PyArrayObject *capi_yclust_tmp = NULL;
  int capi_yclust_intent = 0;
  float *zclust = NULL;
  npy_intp zclust_Dims[2] = {-1, -1};
  const int zclust_Rank = 2;
  PyArrayObject *capi_zclust_tmp = NULL;
  int capi_zclust_intent = 0;
  float *fclust = NULL;
  npy_intp fclust_Dims[1] = {-1};
  const int fclust_Rank = 1;
  PyArrayObject *capi_fclust_tmp = NULL;
  int capi_fclust_intent = 0;
  float rms = 0;
  float *rmsclust = NULL;
  npy_intp rmsclust_Dims[1] = {-1};
  const int rmsclust_Rank = 1;
  PyArrayObject *capi_rmsclust_tmp = NULL;
  int capi_rmsclust_intent = 0;
  float zrms = 0;
  int iterations = 0;
    static char *capi_kwlist[] = {"xl","yl","zl","n","niters","ntime","ncluster",NULL};

/*routdebugenter*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_clock();
#endif
    if (!PyArg_ParseTupleAndKeywords(capi_args,capi_keywds,\
        "OOO|OOOO:trajectory_clustering.clustering",\
        capi_kwlist,&xl_capi,&yl_capi,&zl_capi,&n_capi,&niters_capi,&ntime_capi,&ncluster_capi))
        return NULL;
/*frompyobj*/
  /* Processing variable xl */
  ;
  capi_xl_intent |= F2PY_INTENT_IN;
  capi_xl_tmp = array_from_pyobj(NPY_FLOAT,xl_Dims,xl_Rank,capi_xl_intent,xl_capi);
  if (capi_xl_tmp == NULL) {
    PyObject *exc, *val, *tb;
    PyErr_Fetch(&exc, &val, &tb);
    PyErr_SetString(exc ? exc : trajectory_clustering_error,"failed in converting 1st argument `xl' of trajectory_clustering.clustering to C/Fortran array" );
    npy_PyErr_ChainExceptionsCause(exc, val, tb);
  } else {
    xl = (float *)(PyArray_DATA(capi_xl_tmp));

  /* Processing variable ncluster */
  if (ncluster_capi == Py_None) ncluster = 5; else
    f2py_success = int_from_pyobj(&ncluster,ncluster_capi,"trajectory_clustering.clustering() 4th keyword (ncluster) can't be converted to int");
  if (f2py_success) {
  /* Processing variable niters */
  if (niters_capi == Py_None) niters = 100; else
    f2py_success = int_from_pyobj(&niters,niters_capi,"trajectory_clustering.clustering() 2nd keyword (niters) can't be converted to int");
  if (f2py_success) {
  /* Processing variable iterations */
  /* Processing variable rms */
  /* Processing variable zrms */
  /* Processing variable n */
  if (n_capi == Py_None) n = shape(xl,1); else
    f2py_success = int_from_pyobj(&n,n_capi,"trajectory_clustering.clustering() 1st keyword (n) can't be converted to int");
  if (f2py_success) {
  CHECKSCALAR(shape(xl,1)==n,"shape(xl,1)==n","1st keyword n","clustering:n=%d",n) {
  /* Processing variable ntime */
  if (ntime_capi == Py_None) ntime = shape(xl,0); else
    f2py_success = int_from_pyobj(&ntime,ntime_capi,"trajectory_clustering.clustering() 3rd keyword (ntime) can't be converted to int");
  if (f2py_success) {
  CHECKSCALAR(shape(xl,0)==ntime,"shape(xl,0)==ntime","3rd keyword ntime","clustering:ntime=%d",ntime) {
  /* Processing variable xclust */
  xclust_Dims[0]=ntime,xclust_Dims[1]=ncluster;
  capi_xclust_intent |= F2PY_INTENT_OUT|F2PY_INTENT_HIDE;
  capi_xclust_tmp = array_from_pyobj(NPY_FLOAT,xclust_Dims,xclust_Rank,capi_xclust_intent,Py_None);
  if (capi_xclust_tmp == NULL) {
    PyObject *exc, *val, *tb;
    PyErr_Fetch(&exc, &val, &tb);
    PyErr_SetString(exc ? exc : trajectory_clustering_error,"failed in converting hidden `xclust' of trajectory_clustering.clustering to C/Fortran array" );
    npy_PyErr_ChainExceptionsCause(exc, val, tb);
  } else {
    xclust = (float *)(PyArray_DATA(capi_xclust_tmp));

  /* Processing variable yclust */
  yclust_Dims[0]=ntime,yclust_Dims[1]=ncluster;
  capi_yclust_intent |= F2PY_INTENT_OUT|F2PY_INTENT_HIDE;
  capi_yclust_tmp = array_from_pyobj(NPY_FLOAT,yclust_Dims,yclust_Rank,capi_yclust_intent,Py_None);
  if (capi_yclust_tmp == NULL) {
    PyObject *exc, *val, *tb;
    PyErr_Fetch(&exc, &val, &tb);
    PyErr_SetString(exc ? exc : trajectory_clustering_error,"failed in converting hidden `yclust' of trajectory_clustering.clustering to C/Fortran array" );
    npy_PyErr_ChainExceptionsCause(exc, val, tb);
  } else {
    yclust = (float *)(PyArray_DATA(capi_yclust_tmp));

  /* Processing variable zclust */
  zclust_Dims[0]=ntime,zclust_Dims[1]=ncluster;
  capi_zclust_intent |= F2PY_INTENT_OUT|F2PY_INTENT_HIDE;
  capi_zclust_tmp = array_from_pyobj(NPY_FLOAT,zclust_Dims,zclust_Rank,capi_zclust_intent,Py_None);
  if (capi_zclust_tmp == NULL) {
    PyObject *exc, *val, *tb;
    PyErr_Fetch(&exc, &val, &tb);
    PyErr_SetString(exc ? exc : trajectory_clustering_error,"failed in converting hidden `zclust' of trajectory_clustering.clustering to C/Fortran array" );
    npy_PyErr_ChainExceptionsCause(exc, val, tb);
  } else {
    zclust = (float *)(PyArray_DATA(capi_zclust_tmp));

  /* Processing variable fclust */
  fclust_Dims[0]=ncluster;
  capi_fclust_intent |= F2PY_INTENT_OUT|F2PY_INTENT_HIDE;
  capi_fclust_tmp = array_from_pyobj(NPY_FLOAT,fclust_Dims,fclust_Rank,capi_fclust_intent,Py_None);
  if (capi_fclust_tmp == NULL) {
    PyObject *exc, *val, *tb;
    PyErr_Fetch(&exc, &val, &tb);
    PyErr_SetString(exc ? exc : trajectory_clustering_error,"failed in converting hidden `fclust' of trajectory_clustering.clustering to C/Fortran array" );
    npy_PyErr_ChainExceptionsCause(exc, val, tb);
  } else {
    fclust = (float *)(PyArray_DATA(capi_fclust_tmp));

  /* Processing variable rmsclust */
  rmsclust_Dims[0]=ncluster;
  capi_rmsclust_intent |= F2PY_INTENT_OUT|F2PY_INTENT_HIDE;
  capi_rmsclust_tmp = array_from_pyobj(NPY_FLOAT,rmsclust_Dims,rmsclust_Rank,capi_rmsclust_intent,Py_None);
  if (capi_rmsclust_tmp == NULL) {
    PyObject *exc, *val, *tb;
    PyErr_Fetch(&exc, &val, &tb);
    PyErr_SetString(exc ? exc : trajectory_clustering_error,"failed in converting hidden `rmsclust' of trajectory_clustering.clustering to C/Fortran array" );
    npy_PyErr_ChainExceptionsCause(exc, val, tb);
  } else {
    rmsclust = (float *)(PyArray_DATA(capi_rmsclust_tmp));

  /* Processing variable yl */
  yl_Dims[0]=ntime,yl_Dims[1]=n;
  capi_yl_intent |= F2PY_INTENT_IN;
  capi_yl_tmp = array_from_pyobj(NPY_FLOAT,yl_Dims,yl_Rank,capi_yl_intent,yl_capi);
  if (capi_yl_tmp == NULL) {
    PyObject *exc, *val, *tb;
    PyErr_Fetch(&exc, &val, &tb);
    PyErr_SetString(exc ? exc : trajectory_clustering_error,"failed in converting 2nd argument `yl' of trajectory_clustering.clustering to C/Fortran array" );
    npy_PyErr_ChainExceptionsCause(exc, val, tb);
  } else {
    yl = (float *)(PyArray_DATA(capi_yl_tmp));

  /* Processing variable zl */
  zl_Dims[0]=ntime,zl_Dims[1]=n;
  capi_zl_intent |= F2PY_INTENT_IN;
  capi_zl_tmp = array_from_pyobj(NPY_FLOAT,zl_Dims,zl_Rank,capi_zl_intent,zl_capi);
  if (capi_zl_tmp == NULL) {
    PyObject *exc, *val, *tb;
    PyErr_Fetch(&exc, &val, &tb);
    PyErr_SetString(exc ? exc : trajectory_clustering_error,"failed in converting 3rd argument `zl' of trajectory_clustering.clustering to C/Fortran array" );
    npy_PyErr_ChainExceptionsCause(exc, val, tb);
  } else {
    zl = (float *)(PyArray_DATA(capi_zl_tmp));

/*end of frompyobj*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_call_clock();
#endif
/*callfortranroutine*/
        (*f2py_func)(xl,yl,zl,&n,&niters,&ntime,&ncluster,xclust,yclust,zclust,fclust,&rms,rmsclust,&zrms,&iterations);
if (PyErr_Occurred())
  f2py_success = 0;
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_call_clock();
#endif
/*end of callfortranroutine*/
        if (f2py_success) {
/*pyobjfrom*/
/*end of pyobjfrom*/
        CFUNCSMESS("Building return value.\n");
        capi_buildvalue = Py_BuildValue("NNNNfNfi",capi_xclust_tmp,capi_yclust_tmp,capi_zclust_tmp,capi_fclust_tmp,rms,capi_rmsclust_tmp,zrms,iterations);
/*closepyobjfrom*/
/*end of closepyobjfrom*/
        } /*if (f2py_success) after callfortranroutine*/
/*cleanupfrompyobj*/
  if((PyObject *)capi_zl_tmp!=zl_capi) {
    Py_XDECREF(capi_zl_tmp); }
  }  /*if (capi_zl_tmp == NULL) ... else of zl*/
  /* End of cleaning variable zl */
  if((PyObject *)capi_yl_tmp!=yl_capi) {
    Py_XDECREF(capi_yl_tmp); }
  }  /*if (capi_yl_tmp == NULL) ... else of yl*/
  /* End of cleaning variable yl */
  }  /*if (capi_rmsclust_tmp == NULL) ... else of rmsclust*/
  /* End of cleaning variable rmsclust */
  }  /*if (capi_fclust_tmp == NULL) ... else of fclust*/
  /* End of cleaning variable fclust */
  }  /*if (capi_zclust_tmp == NULL) ... else of zclust*/
  /* End of cleaning variable zclust */
  }  /*if (capi_yclust_tmp == NULL) ... else of yclust*/
  /* End of cleaning variable yclust */
  }  /*if (capi_xclust_tmp == NULL) ... else of xclust*/
  /* End of cleaning variable xclust */
  } /*CHECKSCALAR(shape(xl,0)==ntime)*/
  } /*if (f2py_success) of ntime*/
  /* End of cleaning variable ntime */
  } /*CHECKSCALAR(shape(xl,1)==n)*/
  } /*if (f2py_success) of n*/
  /* End of cleaning variable n */
  /* End of cleaning variable zrms */
  /* End of cleaning variable rms */
  /* End of cleaning variable iterations */
  } /*if (f2py_success) of niters*/
  /* End of cleaning variable niters */
  } /*if (f2py_success) of ncluster*/
  /* End of cleaning variable ncluster */
  if((PyObject *)capi_xl_tmp!=xl_capi) {
    Py_XDECREF(capi_xl_tmp); }
  }  /*if (capi_xl_tmp == NULL) ... else of xl*/
  /* End of cleaning variable xl */
/*end of cleanupfrompyobj*/
    if (capi_buildvalue == NULL) {
/*routdebugfailure*/
    } else {
/*routdebugleave*/
    }
    CFUNCSMESS("Freeing memory.\n");
/*freemem*/
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_clock();
#endif
    return capi_buildvalue;
}
/***************************** end of clustering *****************************/

/********************************* distance2 *********************************/
static char doc_f2py_rout_trajectory_clustering_distance2[] = "\
distance2 = distance2(rlat1,rlon1,rlat2,rlon2)\n\nWrapper for ``distance2``.\
\n\nParameters\n----------\n"
"rlat1 : input float\n"
"rlon1 : input float\n"
"rlat2 : input float\n"
"rlon2 : input float\n"
"\nReturns\n-------\n"
"distance2 : float";
/* extern void F_WRAPPEDFUNC(distance2,DISTANCE2)(float*,float*,float*,float*,float*); */
static PyObject *f2py_rout_trajectory_clustering_distance2(const PyObject *capi_self,
                           PyObject *capi_args,
                           PyObject *capi_keywds,
                           void (*f2py_func)(float*,float*,float*,float*,float*)) {
    PyObject * volatile capi_buildvalue = NULL;
    volatile int f2py_success = 1;
/*decl*/

  float distance2 = 0;
  float rlat1 = 0;
  PyObject *rlat1_capi = Py_None;
  float rlon1 = 0;
  PyObject *rlon1_capi = Py_None;
  float rlat2 = 0;
  PyObject *rlat2_capi = Py_None;
  float rlon2 = 0;
  PyObject *rlon2_capi = Py_None;
    static char *capi_kwlist[] = {"rlat1","rlon1","rlat2","rlon2",NULL};

/*routdebugenter*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_clock();
#endif
    if (!PyArg_ParseTupleAndKeywords(capi_args,capi_keywds,\
        "OOOO|:trajectory_clustering.distance2",\
        capi_kwlist,&rlat1_capi,&rlon1_capi,&rlat2_capi,&rlon2_capi))
        return NULL;
/*frompyobj*/
  /* Processing variable rlat1 */
    f2py_success = float_from_pyobj(&rlat1,rlat1_capi,"trajectory_clustering.distance2() 1st argument (rlat1) can't be converted to float");
  if (f2py_success) {
  /* Processing variable rlon1 */
    f2py_success = float_from_pyobj(&rlon1,rlon1_capi,"trajectory_clustering.distance2() 2nd argument (rlon1) can't be converted to float");
  if (f2py_success) {
  /* Processing variable rlat2 */
    f2py_success = float_from_pyobj(&rlat2,rlat2_capi,"trajectory_clustering.distance2() 3rd argument (rlat2) can't be converted to float");
  if (f2py_success) {
  /* Processing variable rlon2 */
    f2py_success = float_from_pyobj(&rlon2,rlon2_capi,"trajectory_clustering.distance2() 4th argument (rlon2) can't be converted to float");
  if (f2py_success) {
  /* Processing variable distance2 */
/*end of frompyobj*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_call_clock();
#endif
/*callfortranroutine*/
  (*f2py_func)(&distance2,&rlat1,&rlon1,&rlat2,&rlon2);
if (PyErr_Occurred())
  f2py_success = 0;
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_call_clock();
#endif
/*end of callfortranroutine*/
        if (f2py_success) {
/*pyobjfrom*/
/*end of pyobjfrom*/
        CFUNCSMESS("Building return value.\n");
        capi_buildvalue = Py_BuildValue("f",distance2);
/*closepyobjfrom*/
/*end of closepyobjfrom*/
        } /*if (f2py_success) after callfortranroutine*/
/*cleanupfrompyobj*/
  /* End of cleaning variable distance2 */
  } /*if (f2py_success) of rlon2*/
  /* End of cleaning variable rlon2 */
  } /*if (f2py_success) of rlat2*/
  /* End of cleaning variable rlat2 */
  } /*if (f2py_success) of rlon1*/
  /* End of cleaning variable rlon1 */
  } /*if (f2py_success) of rlat1*/
  /* End of cleaning variable rlat1 */
/*end of cleanupfrompyobj*/
    if (capi_buildvalue == NULL) {
/*routdebugfailure*/
    } else {
/*routdebugleave*/
    }
    CFUNCSMESS("Freeing memory.\n");
/*freemem*/
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_clock();
#endif
    return capi_buildvalue;
}
/****************************** end of distance2 ******************************/
/*eof body*/

/******************* See f2py2e/f90mod_rules.py: buildhooks *******************/
/*need_f90modhooks*/

/************** See f2py2e/rules.py: module_rules['modulebody'] **************/

/******************* See f2py2e/common_rules.py: buildhooks *******************/

/*need_commonhooks*/

/**************************** See f2py2e/rules.py ****************************/

static FortranDataDef f2py_routine_defs[] = {
  {"clustering",-1,{{-1}},0,(char *)F_FUNC(clustering,CLUSTERING),(f2py_init_func)f2py_rout_trajectory_clustering_clustering,doc_f2py_rout_trajectory_clustering_clustering},
  {"distance2",-1,{{-1}},0,(char *)F_WRAPPEDFUNC(distance2,DISTANCE2),(f2py_init_func)f2py_rout_trajectory_clustering_distance2,doc_f2py_rout_trajectory_clustering_distance2},

/*eof routine_defs*/
  {NULL}
};

static PyMethodDef f2py_module_methods[] = {

  {NULL,NULL}
};

static struct PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  "trajectory_clustering",
  NULL,
  -1,
  f2py_module_methods,
  NULL,
  NULL,
  NULL,
  NULL
};

PyMODINIT_FUNC PyInit_trajectory_clustering(void) {
  int i;
  PyObject *m,*d, *s, *tmp;
  m = trajectory_clustering_module = PyModule_Create(&moduledef);
  Py_SET_TYPE(&PyFortran_Type, &PyType_Type);
  import_array();
  if (PyErr_Occurred())
    {PyErr_SetString(PyExc_ImportError, "can't initialize module trajectory_clustering (failed to import numpy)"); return m;}
  d = PyModule_GetDict(m);
  s = PyUnicode_FromString("$Revision: $");
  PyDict_SetItemString(d, "__version__", s);
  Py_DECREF(s);
  s = PyUnicode_FromString(
    "This module 'trajectory_clustering' is auto-generated with f2py (version:2).\nFunctions:\n"
"  xclust,yclust,zclust,fclust,rms,rmsclust,zrms,iterations = clustering(xl,yl,zl,n=shape(xl,1),niters=100,ntime=shape(xl,0),ncluster=5)\n"
"  distance2 = distance2(rlat1,rlon1,rlat2,rlon2)\n"
".");
  PyDict_SetItemString(d, "__doc__", s);
  Py_DECREF(s);
  s = PyUnicode_FromString("1.20.1");
  PyDict_SetItemString(d, "__f2py_numpy_version__", s);
  Py_DECREF(s);
  trajectory_clustering_error = PyErr_NewException ("trajectory_clustering.error", NULL, NULL);
  /*
   * Store the error object inside the dict, so that it could get deallocated.
   * (in practice, this is a module, so it likely will not and cannot.)
   */
  PyDict_SetItemString(d, "_trajectory_clustering_error", trajectory_clustering_error);
  Py_DECREF(trajectory_clustering_error);
  for(i=0;f2py_routine_defs[i].name!=NULL;i++) {
    tmp = PyFortranObject_NewAsAttr(&f2py_routine_defs[i]);
    PyDict_SetItemString(d, f2py_routine_defs[i].name, tmp);
    Py_DECREF(tmp);
  }


    {
      extern float F_FUNC(distance2,DISTANCE2)(void);
      PyObject* o = PyDict_GetItemString(d,"distance2");
      tmp = F2PyCapsule_FromVoidPtr((void*)F_FUNC(distance2,DISTANCE2),NULL);
      PyObject_SetAttrString(o,"_cpointer", tmp);
      Py_DECREF(tmp);
      s = PyUnicode_FromString("distance2");
      PyObject_SetAttrString(o,"__name__", s);
      Py_DECREF(s);
    }
    
/*eof initf2pywraphooks*/
/*eof initf90modhooks*/

/*eof initcommonhooks*/


#ifdef F2PY_REPORT_ATEXIT
  if (! PyErr_Occurred())
    on_exit(f2py_report_on_exit,(void*)"trajectory_clustering");
#endif
  return m;
}
#ifdef __cplusplus
}
#endif
