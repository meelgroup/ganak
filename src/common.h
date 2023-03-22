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

// verbose debug
#ifdef VERBOSE_DEBUG
#define print_debug(x) std::cout << COLDEF << x << endl
#define print_debug_noendl(x) std::cout << x
#else
#define print_debug(x) do {} while(0)
#define print_debug_noendl(x) do {} while (0)
#endif

