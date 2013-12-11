#include <Python.h>
#include <numpy/arrayobject.h>
#include "lightcurve.h"

struct module_state {
    PyObject *error;
};

#if PY_MAJOR_VERSION >= 3
#define GETSTATE(m) ((struct module_state*)PyModule_GetState(m))
#else
#define GETSTATE(m) (&_state)
static struct module_state _state;
#endif

#define PARSE_ARRAY(o) (PyArrayObject*) PyArray_FROM_OTF(o, NPY_DOUBLE, \
        NPY_IN_ARRAY)

PyObject *kois_light_curve (PyObject *self, PyObject *args)
{
    int K;
    double texp, mu1, mu2;
    PyObject *t_obj, *periods_obj, *epochs_obj, *durations_obj,
             *rors_obj, *impacts_obj;
    if (!PyArg_ParseTuple(args, "OidOOOOOdd", &t_obj, &K, &texp, &periods_obj,
                          &epochs_obj, &durations_obj, &rors_obj,
                          &impacts_obj, &mu1, &mu2))
        return NULL;

    PyArrayObject *t_array = PARSE_ARRAY(t_obj),
                  *periods_array = PARSE_ARRAY(periods_obj),
                  *epochs_array = PARSE_ARRAY(epochs_obj),
                  *durations_array = PARSE_ARRAY(durations_obj),
                  *rors_array = PARSE_ARRAY(rors_obj),
                  *impacts_array = PARSE_ARRAY(impacts_obj);

    int n = (int) PyArray_DIM(t_array, 0),
        np = (int) PyArray_DIM(periods_array, 0);

    npy_intp dim[1] = {n};
    PyArrayObject *flux_array =
                        (PyArrayObject*)PyArray_SimpleNew(1, dim, NPY_DOUBLE);

    double *t = PyArray_DATA(t_array),
           *periods = PyArray_DATA(periods_array),
           *epochs = PyArray_DATA(epochs_array),
           *durations = PyArray_DATA(durations_array),
           *rors = PyArray_DATA(rors_array),
           *impacts = PyArray_DATA(impacts_array),
           *f = PyArray_DATA(flux_array);

    lightcurve (n, t, K, texp, np, periods, epochs, durations, rors, impacts,
                mu1, mu2, f);

    Py_DECREF(t_array);
    Py_DECREF(periods_array);
    Py_DECREF(epochs_array);
    Py_DECREF(durations_array);
    Py_DECREF(rors_array);
    Py_DECREF(impacts_array);

    return (PyObject*)flux_array;
}

static PyMethodDef kois_methods[] = {
    {"light_curve", (PyCFunction) kois_light_curve, METH_VARARGS, ""},
    {NULL, NULL, 0, NULL}
};

#if PY_MAJOR_VERSION >= 3

static int kois_traverse(PyObject *m, visitproc visit, void *arg) {
    Py_VISIT(GETSTATE(m)->error);
    return 0;
}

static int kois_clear(PyObject *m) {
    Py_CLEAR(GETSTATE(m)->error);
    return 0;
}

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_kois",
    NULL,
    sizeof(struct module_state),
    kois_methods,
    NULL,
    kois_traverse,
    kois_clear,
    NULL
};

#define INITERROR return NULL

PyObject *PyInit__kois(void)
#else
#define INITERROR return

void init_kois(void)
#endif
{
#if PY_MAJOR_VERSION >= 3
    PyObject *module = PyModule_Create(&moduledef);
#else
    PyObject *module = Py_InitModule("_kois", kois_methods);
#endif

    if (module == NULL)
        INITERROR;
    struct module_state *st = GETSTATE(module);

    st->error = PyErr_NewException("_kois.Error", NULL, NULL);
    if (st->error == NULL) {
        Py_DECREF(module);
        INITERROR;
    }

    import_array();

#if PY_MAJOR_VERSION >= 3
    return module;
#endif
}
