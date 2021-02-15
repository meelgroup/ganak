#ifndef GANAK_PACKED_COMPONENT_H
#define GANAK_PACKED_COMPONENT_H


#include "cacheable_component.h"

template<class T>
class BitStuffer {
public:
    BitStuffer(T *data) : data_start_(data), p(data) {
        *p = 0;
    }

    void stuff(const unsigned val, const unsigned num_bits_val) {
        assert(num_bits_val > 0);
        assert((val >> num_bits_val) == 0);
        if (end_of_bits_ == 0)
            *p = 0;
        assert((*p >> end_of_bits_) == 0);
        *p |= val << end_of_bits_;
        end_of_bits_ += num_bits_val;
        if (end_of_bits_ > _bits_per_block) {
            //assert(*p);
            end_of_bits_ -= _bits_per_block;
            *(++p) = val >> (num_bits_val - end_of_bits_);
            assert(!(end_of_bits_ == 0) | (*p == 0));
        } else if (end_of_bits_ == _bits_per_block) {
            end_of_bits_ -= _bits_per_block;
            p++;
        }
    }

    void assert_size(unsigned size) {
        if (end_of_bits_ == 0)
            p--;
        assert(p - data_start_ == size - 1);
    }

private:
    T *data_start_ = nullptr;
    T *p = nullptr;
    // in the current block
    // the bit postion just after the last bit written
    unsigned end_of_bits_ = 0;

    static const unsigned _bits_per_block = (sizeof(T) << 3);

};

class PackedComponent : public BaseComponent {
public:

    virtual ~PackedComponent() = default;

    static unsigned bits_per_variable() {
        return _bits_per_variable;
    }

    static unsigned variable_mask() {
        return _variable_mask;
    }

    static unsigned bits_per_clause() {
        return _bits_per_clause;
    }

    static unsigned bits_per_block() {
        return _bits_per_block;
    }

    static unsigned bits_of_data_size() {
        return _bits_of_data_size;
    }

    static void adjustPackSize(unsigned int maxVarId, unsigned int maxClId);

    static void outbit(unsigned v) {
        for (auto i = 0; i < 32; i++) {
            cout << ((v & 2147483648) ? "1" : "0");
            v &= 2147483648 - 1;
            v <<= 1;
        }
    }

    static unsigned log2(unsigned v) {
        // taken from
        // http://graphics.stanford.edu/~seander/bithacks.html#IntegerLogLookup
        static const char LogTable256[256] =
                {
#define LT(n) n, n, n, n, n, n, n, n, n, n, n, n, n, n, n, n
                        -1, 0, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3,
                        LT(4), LT(5), LT(5), LT(6), LT(6), LT(6), LT(6),
                        LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7), LT(7)
                };

        unsigned r;     // r will be lg(v)
        unsigned int t, tt; // temporaries

        if ((tt = (v >> 16))) {
            r = (t = (tt >> 8)) ? 24 + LogTable256[t] : 16 + LogTable256[tt];
        } else {
            r = (t = (v >> 8)) ? 8 + LogTable256[t] : LogTable256[v];
        }
        return r;
    }

    static unsigned _debug_static_val;

protected:
    static unsigned _data_size_mask;

    static unsigned _bits_per_clause, _bits_per_variable; // bitsperentry
    static unsigned _bits_of_data_size; // number of bits needed to store the data size

    static unsigned _variable_mask, _clause_mask;
    static const unsigned _bits_per_block = (sizeof(unsigned) << 3);
};


#endif //GANAK_PACKED_COMPONENT_H
