/******************************************
Copyright (C) 2023-2025 Authors of GANAK, see AUTHORS file

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
***********************************************/

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <ganak.hpp>
#include <arjun/arjun.h>
#include <cryptominisat5/solvertypesmini.h>

#include <vector>
#include <set>
#include <string>
#include <memory>
#include <stdexcept>

#ifndef PYGANAK_VERSION
#define PYGANAK_VERSION "1.1.0"
#endif

// -------------------------------------------------------------------------
// Internal state stored per Counter object
// -------------------------------------------------------------------------
struct CounterState {
    int verbose = 0;
    uint64_t seed = 0;
    uint32_t max_var = 0;  // highest 1-indexed variable seen so far
    std::vector<std::vector<CMSat::Lit>> clauses;
    std::vector<uint32_t> sampling_set;   // 0-indexed (CMSat convention)
    bool sampling_set_given = false;
};

typedef struct {
    PyObject_HEAD
    CounterState* state;
} Counter;

// -------------------------------------------------------------------------
// Helpers
// -------------------------------------------------------------------------

// Convert a Python integer literal value (signed, 1-indexed) to CMSat::Lit.
// Returns false and sets a Python exception on error.
static bool lit_from_py_long(PyObject* item, CMSat::Lit& out_lit, uint32_t& max_var) {
    long val = PyLong_AsLong(item);
    if (val == -1L && PyErr_Occurred()) return false;
    if (val == 0) {
        PyErr_SetString(PyExc_ValueError, "Clause contains 0; literals must be nonzero integers");
        return false;
    }
    uint32_t var = (uint32_t)((val < 0) ? -val : val) - 1;  // convert to 0-indexed
    bool sign = (val < 0);                                    // true = negative literal
    out_lit = CMSat::Lit(var, sign);
    if (var + 1 > max_var) max_var = var + 1;
    return true;
}

// Parse an iterable of signed integers into a vector of CMSat::Lit.
static bool parse_clause(PyObject* clause_obj,
                          std::vector<CMSat::Lit>& lits,
                          uint32_t& max_var) {
    PyObject* iter = PyObject_GetIter(clause_obj);
    if (!iter) return false;

    PyObject* item;
    while ((item = PyIter_Next(iter)) != nullptr) {
        CMSat::Lit lit;
        bool ok = lit_from_py_long(item, lit, max_var);
        Py_DECREF(item);
        if (!ok) { Py_DECREF(iter); return false; }
        lits.push_back(lit);
    }
    Py_DECREF(iter);
    return !PyErr_Occurred();
}

// -------------------------------------------------------------------------
// Counter.__new__ / __init__ / __dealloc__
// -------------------------------------------------------------------------

static PyObject* Counter_new(PyTypeObject* type, PyObject* /*args*/, PyObject* /*kwds*/) {
    Counter* self = (Counter*)type->tp_alloc(type, 0);
    if (self) {
        self->state = new (std::nothrow) CounterState();
        if (!self->state) {
            Py_DECREF(self);
            return PyErr_NoMemory();
        }
    }
    return (PyObject*)self;
}

static int Counter_init(Counter* self, PyObject* args, PyObject* kwds) {
    static const char* kwlist[] = {"verbose", "seed", nullptr};
    int verbose = 0;
    unsigned long long seed = 0;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|iK",
                                      const_cast<char**>(kwlist),
                                      &verbose, &seed))
        return -1;
    self->state->verbose = verbose;
    self->state->seed = (uint64_t)seed;
    return 0;
}

static void Counter_dealloc(Counter* self) {
    delete self->state;
    self->state = nullptr;
    Py_TYPE(self)->tp_free((PyObject*)self);
}

// -------------------------------------------------------------------------
// Counter.add_clause(clause)
// -------------------------------------------------------------------------

static PyObject* Counter_add_clause(Counter* self, PyObject* args) {
    PyObject* clause_obj = nullptr;
    if (!PyArg_ParseTuple(args, "O", &clause_obj)) return nullptr;

    std::vector<CMSat::Lit> lits;
    if (!parse_clause(clause_obj, lits, self->state->max_var)) return nullptr;

    self->state->clauses.push_back(std::move(lits));
    Py_RETURN_NONE;
}

// -------------------------------------------------------------------------
// Counter.add_clauses(clauses)   -- iterable of clauses
// -------------------------------------------------------------------------

static PyObject* Counter_add_clauses(Counter* self, PyObject* args) {
    PyObject* clauses_obj = nullptr;
    if (!PyArg_ParseTuple(args, "O", &clauses_obj)) return nullptr;

    PyObject* outer_iter = PyObject_GetIter(clauses_obj);
    if (!outer_iter) return nullptr;

    PyObject* clause_obj;
    while ((clause_obj = PyIter_Next(outer_iter)) != nullptr) {
        std::vector<CMSat::Lit> lits;
        bool ok = parse_clause(clause_obj, lits, self->state->max_var);
        Py_DECREF(clause_obj);
        if (!ok) { Py_DECREF(outer_iter); return nullptr; }
        self->state->clauses.push_back(std::move(lits));
    }
    Py_DECREF(outer_iter);
    if (PyErr_Occurred()) return nullptr;
    Py_RETURN_NONE;
}

// -------------------------------------------------------------------------
// Counter.set_sampling_set(vars)
// vars: iterable of 1-indexed variable numbers
// -------------------------------------------------------------------------

static PyObject* Counter_set_sampling_set(Counter* self, PyObject* args) {
    PyObject* vars_obj = nullptr;
    if (!PyArg_ParseTuple(args, "O", &vars_obj)) return nullptr;

    PyObject* iter = PyObject_GetIter(vars_obj);
    if (!iter) return nullptr;

    std::vector<uint32_t> sampl;
    PyObject* item;
    while ((item = PyIter_Next(iter)) != nullptr) {
        long v = PyLong_AsLong(item);
        Py_DECREF(item);
        if (v == -1L && PyErr_Occurred()) { Py_DECREF(iter); return nullptr; }
        if (v <= 0) {
            PyErr_SetString(PyExc_ValueError, "Sampling set variables must be positive integers (1-indexed)");
            Py_DECREF(iter);
            return nullptr;
        }
        uint32_t var0 = (uint32_t)v - 1;   // 0-indexed
        sampl.push_back(var0);
        if (var0 + 1 > self->state->max_var) self->state->max_var = var0 + 1;
    }
    Py_DECREF(iter);
    if (PyErr_Occurred()) return nullptr;

    self->state->sampling_set = std::move(sampl);
    self->state->sampling_set_given = true;
    Py_RETURN_NONE;
}

// -------------------------------------------------------------------------
// Counter.count()
// Returns a Python int (exact model count).
// Arjun preprocessing is run automatically.
// -------------------------------------------------------------------------

static PyObject* Counter_count(Counter* self, PyObject*) {
    // Snapshot mutable state while we hold the GIL
    const int verbose         = self->state->verbose;
    const uint64_t seed       = self->state->seed;
    const uint32_t max_var    = self->state->max_var;
    const bool ss_given       = self->state->sampling_set_given;

    // Take copies so we can release the GIL safely
    auto clauses_copy    = self->state->clauses;
    auto sampl_set_copy  = self->state->sampling_set;

    // ---- Result / error storage (no Python objects!) ----
    std::unique_ptr<CMSat::Field> result_field;
    std::string error_msg;
    bool success = false;

    Py_BEGIN_ALLOW_THREADS
    try {
        // ---- Field generator: integer (mpz) ----
        std::unique_ptr<CMSat::FieldGen> fg = std::make_unique<ArjunNS::FGenMpz>();

        // ---- Build SimplifiedCNF ----
        ArjunNS::SimplifiedCNF cnf(fg);
        if (max_var > 0) cnf.new_vars(max_var);
        for (const auto& cl : clauses_copy) cnf.add_clause(cl);

        // ---- Sampling set ----
        if (ss_given) {
            cnf.set_sampl_vars(sampl_set_copy);
        } else {
            // All variables are in the sampling set
            std::vector<uint32_t> all_vars;
            all_vars.reserve(max_var);
            for (uint32_t i = 0; i < max_var; i++) all_vars.push_back(i);
            cnf.set_sampl_vars(all_vars);
        }

        // ---- Arjun preprocessing ----
        ArjunNS::Arjun arjun;
        arjun.set_verb(0);
        ArjunNS::Arjun::ElimToFileConf etof_conf;  // default settings
        ArjunNS::SimpConf simp_conf;                // default settings
        arjun.standalone_minimize_indep(cnf, /*all_indep=*/false);
        arjun.standalone_elim_to_file(cnf, etof_conf, simp_conf);

        // ---- Ganak counting ----
        GanakInt::CounterConfiguration conf;
        conf.verb = verbose;
        conf.seed = seed;
        Ganak counter(conf, fg);
        setup_ganak(cnf, counter);

        auto cnt = cnf.get_multiplier_weight()->dup();
        if (!cnf.get_multiplier_weight()->is_zero()) {
            *cnt *= *counter.count();
        }
        result_field = std::move(cnt);
        success = true;
    } catch (const std::exception& e) {
        error_msg = e.what();
    } catch (...) {
        error_msg = "Unknown C++ exception in count()";
    }
    Py_END_ALLOW_THREADS

    if (!success) {
        PyErr_SetString(PyExc_RuntimeError, error_msg.c_str());
        return nullptr;
    }

    // Convert mpz result to Python int via string representation
    auto& mpz_result = static_cast<ArjunNS::FMpz&>(*result_field);
    std::string s = mpz_result.val.get_str(10);
    return PyLong_FromString(s.c_str(), nullptr, 10);
}

// -------------------------------------------------------------------------
// Counter.nof_vars()  -- number of variables declared so far
// -------------------------------------------------------------------------

static PyObject* Counter_nof_vars(Counter* self, PyObject*) {
    return PyLong_FromUnsignedLong((unsigned long)self->state->max_var);
}

// -------------------------------------------------------------------------
// Counter.nof_clauses()
// -------------------------------------------------------------------------

static PyObject* Counter_nof_clauses(Counter* self, PyObject*) {
    return PyLong_FromSize_t(self->state->clauses.size());
}

// -------------------------------------------------------------------------
// Method table
// -------------------------------------------------------------------------

static PyMethodDef counter_methods[] = {
    {"add_clause",       (PyCFunction)Counter_add_clause,       METH_VARARGS,
     "add_clause(clause)\n\n"
     "Add a clause. *clause* is an iterable of nonzero signed integers;\n"
     "positive integers represent positive literals, negative integers\n"
     "represent negative literals.  Variable indices are 1-based.\n\n"
     "Example::\n\n"
     "    c.add_clause([1, -2, 3])   # x1 OR NOT x2 OR x3\n"},

    {"add_clauses",      (PyCFunction)Counter_add_clauses,      METH_VARARGS,
     "add_clauses(clauses)\n\n"
     "Add multiple clauses at once. *clauses* is an iterable of clauses,\n"
     "where each clause is an iterable of nonzero signed integers.\n"},

    {"set_sampling_set", (PyCFunction)Counter_set_sampling_set, METH_VARARGS,
     "set_sampling_set(vars)\n\n"
     "Set the *sampling set* (also called the independent support).\n"
     "*vars* is an iterable of positive integers (1-based variable indices).\n"
     "Only variables in this set are counted; variables not in the set are\n"
     "treated as existentially quantified.  If this method is never called,\n"
     "all variables are included in the sampling set (unweighted projected\n"
     "model counting becomes plain model counting).\n"},

    {"count",            (PyCFunction)Counter_count,            METH_NOARGS,
     "count() -> int\n\n"
     "Run Arjun preprocessing followed by the Ganak model counter and\n"
     "return the exact model count as a Python integer.\n\n"
     "The formula may be called multiple times (e.g. after adding more\n"
     "clauses), but note that a fresh Arjun+Ganak run is performed every\n"
     "time.\n"},

    {"nof_vars",         (PyCFunction)Counter_nof_vars,         METH_NOARGS,
     "nof_vars() -> int\n\nReturn the number of variables seen so far.\n"},

    {"nof_clauses",      (PyCFunction)Counter_nof_clauses,      METH_NOARGS,
     "nof_clauses() -> int\n\nReturn the number of clauses added so far.\n"},

    {nullptr, nullptr, 0, nullptr}
};

// -------------------------------------------------------------------------
// Type object
// -------------------------------------------------------------------------

static PyTypeObject counter_type = {
    PyVarObject_HEAD_INIT(nullptr, 0)
    "pyganak.Counter",                           /* tp_name */
    sizeof(Counter),                             /* tp_basicsize */
    0,                                           /* tp_itemsize */
    (destructor)Counter_dealloc,                 /* tp_dealloc */
    0,                                           /* tp_vectorcall_offset */
    nullptr,                                     /* tp_getattr */
    nullptr,                                     /* tp_setattr */
    nullptr,                                     /* tp_as_async */
    nullptr,                                     /* tp_repr */
    nullptr,                                     /* tp_as_number */
    nullptr,                                     /* tp_as_sequence */
    nullptr,                                     /* tp_as_mapping */
    nullptr,                                     /* tp_hash */
    nullptr,                                     /* tp_call */
    nullptr,                                     /* tp_str */
    nullptr,                                     /* tp_getattro */
    nullptr,                                     /* tp_setattro */
    nullptr,                                     /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,                          /* tp_flags */
    "Counter(verbose=0, seed=0)\n\n"
    "Exact model counter backed by the Ganak solver.\n\n"
    "Usage example::\n\n"
    "    from pyganak import Counter\n"
    "    c = Counter()\n"
    "    c.add_clause([1, 2])    # x1 OR x2\n"
    "    c.add_clause([-1, 2])   # NOT x1 OR x2\n"
    "    print(c.count())        # prints 2\n",
                                                 /* tp_doc */
    nullptr,                                     /* tp_traverse */
    nullptr,                                     /* tp_clear */
    nullptr,                                     /* tp_richcompare */
    0,                                           /* tp_weaklistoffset */
    nullptr,                                     /* tp_iter */
    nullptr,                                     /* tp_iternext */
    counter_methods,                             /* tp_methods */
    nullptr,                                     /* tp_members */
    nullptr,                                     /* tp_getset */
    nullptr,                                     /* tp_base */
    nullptr,                                     /* tp_dict */
    nullptr,                                     /* tp_descr_get */
    nullptr,                                     /* tp_descr_set */
    0,                                           /* tp_dictoffset */
    (initproc)Counter_init,                      /* tp_init */
    nullptr,                                     /* tp_alloc */
    Counter_new,                                 /* tp_new */
};

// -------------------------------------------------------------------------
// Module definition
// -------------------------------------------------------------------------

static PyModuleDef pyganak_module = {
    PyModuleDef_HEAD_INIT,
    "pyganak",
    "Python bindings for Ganak, an exact model counter.\n\n"
    "Ganak supports exact projected model counting with Arjun preprocessing.\n"
    "See Counter for the main API.\n",
    -1,
    nullptr,  /* module-level methods */
};

PyMODINIT_FUNC PyInit_pyganak(void) {
    if (PyType_Ready(&counter_type) < 0) return nullptr;

    PyObject* mod = PyModule_Create(&pyganak_module);
    if (!mod) return nullptr;

    Py_INCREF(&counter_type);
    if (PyModule_AddObject(mod, "Counter", (PyObject*)&counter_type) < 0) {
        Py_DECREF(&counter_type);
        Py_DECREF(mod);
        return nullptr;
    }

    if (PyModule_AddStringConstant(mod, "__version__", PYGANAK_VERSION) < 0) {
        Py_DECREF(mod);
        return nullptr;
    }
    if (PyModule_AddStringConstant(mod, "VERSION", PYGANAK_VERSION) < 0) {
        Py_DECREF(mod);
        return nullptr;
    }

    return mod;
}
