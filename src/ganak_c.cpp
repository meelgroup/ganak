/******************************************
Copyright (C) 2024 Authors of GANAK, see AUTHORS file

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

#include "ganak_c.h"
#include "ganak.hpp"

#include <arjun/arjun_c_priv.hpp>

#include <memory>
#include <string>
#include <exception>

using namespace GanakInt;

namespace {

thread_local std::string g_last_error;

void set_error(const std::string& msg) { g_last_error = msg; }
void clear_error() { g_last_error.clear(); }

}  // anonymous

struct ganak_counter {
    Ganak*           counter;
    arjun_field_kind kind;
    CounterConfiguration conf;
    std::unique_ptr<CMSat::FieldGen> fg;  // kept alive for field outputs
};

extern "C" {

const char* ganak_last_error(void) { return g_last_error.c_str(); }

ganak_counter_t* ganak_counter_new(uint32_t verb, uint64_t seed,
                                   const arjun_fgen_t* fg) {
    try {
        clear_error();
        auto* w = new ganak_counter_t();
        w->kind = fg->kind;
        w->fg   = fg->fg->dup();
        w->conf.verb = (int)verb;
        w->conf.seed = seed;
        w->counter = new Ganak(w->conf, w->fg);
        return w;
    } catch (const std::exception& e) { set_error(e.what()); return nullptr; }
}

void ganak_counter_free(ganak_counter_t* c) {
    if (c == nullptr) return;
    delete c->counter;
    delete c;
}

int ganak_counter_setup_from_simpcnf(ganak_counter_t* c,
                                     const arjun_simpcnf_t* cnf) {
    clear_error();
    try {
        setup_ganak(*cnf->cnf, *c->counter);
        return 0;
    } catch (const std::exception& e) { set_error(e.what()); return -1; }
}

arjun_field_t* ganak_counter_count(ganak_counter_t* c,
                                   uint8_t bits_threads,
                                   int num_threads,
                                   int debug_threads) {
    try {
        clear_error();
        auto cnt = c->counter->count(bits_threads, num_threads, debug_threads != 0);
        auto* w = new arjun_field_t();
        w->f = std::move(cnt);
        w->kind = c->kind;
        return w;
    } catch (const std::exception& e) { set_error(e.what()); return nullptr; }
}

int ganak_counter_is_approximate(const ganak_counter_t* c) {
    return c->counter->get_is_approximate() ? 1 : 0;
}

}  // extern "C"
