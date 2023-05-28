// Debug code

#pragma once

/* #define VERBOSE_DEBUG */
/* #define SLOW_DEBUG */

#define COLRED "\033[31m"
#define COLYEL2 "\033[35m"
#define COLYEL "\033[33m"
#define COLCYN "\033[36m"
#define COLWHT "\033[97m"
#define COLORG "\033[43m"
#define COLBLBACK  "\033[44m"
#define COLORGBG "\033[100m"
#define COLREDBG "\033[41m"

//default
#define COLDEF "\033[0m"

// slow debug
#ifdef SLOW_DEBUG
#define SLOW_DEBUG_DO(x) \
    do { x; } while (0)
#else
#define SLOW_DEBUG_DO(x) do { } while (0)
#endif
/////

#define verb_print(a, b) if (config_.verb >= 1) cout << "c o " << b << endl;

// verbose debug
#ifdef VERBOSE_DEBUG
#define VERBOSE_PRINT(x) \
    do { std::cout << x << std::endl; } while (0)
#define VERBOSE_DEBUG_DO(x) do { x; } while (0)
#else
#define VERBOSE_PRINT(x) do { } while (0)
#define VERBOSE_DEBUG_DO(x) do { } while (0)
#endif

#ifdef VERBOSE_DEBUG
#define print_debug(x) std::cout << COLDEF << x << COLDEF << endl
#define print_debug_noendl(x) std::cout << x
#else
#define print_debug(x) do {} while(0)
#define print_debug_noendl(x) do {} while (0)
#endif
#define print_tmpdebug(x) std::cout << COLDEF << x << COLDEF << endl

#define release_assert(a) \
    do { \
        if (!(a)) {\
            fprintf(stderr, "*** ASSERTION FAILURE in %s() [%s:%d]: %s\n", \
            __FUNCTION__, __FILE__, __LINE__, #a); \
            abort(); \
        } \
    } while (0)
