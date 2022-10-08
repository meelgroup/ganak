// Debug code

#pragma once

#define VERBOSE_DEBUG

#define COLRED "\033[31m"
#define COLYEL "\033[33m"
#define COLCYN "\033[36m"
#define COLWHT "\033[97m"
#define COLORG "\033[43m"
#define COLGRN "\033[43m"
#define COLORGBG "\033[100m"
#define COLREDBG "\033[41m"

//default
#define COLDEF "\033[0m"


#ifdef VERBOSE_DEBUG
#define print_debug(x) std::cout << COLDEF << x << endl
#else
#define print_debug(x) do {} while(0)
#endif

