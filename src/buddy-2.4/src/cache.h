/*========================================================================
               Copyright (C) 1996-2002 by Jorn Lind-Nielsen
                            All rights reserved

    Permission is hereby granted, without written agreement and without
    license or royalty fees, to use, reproduce, prepare derivative
    works, distribute, and display this software and its documentation
    for any purpose, provided that (1) the above copyright notice and
    the following two paragraphs appear in all copies of the source code
    and (2) redistributions, including without limitation binaries,
    reproduce these notices in the supporting documentation. Substantial
    modifications to this software may be copyrighted by their authors
    and need not follow the licensing terms described here, provided
    that the new terms are clearly indicated in all files where they apply.

    IN NO EVENT SHALL JORN LIND-NIELSEN, OR DISTRIBUTORS OF THIS
    SOFTWARE BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL,
    INCIDENTAL, OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OF THIS
    SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE AUTHORS OR ANY OF THE
    ABOVE PARTIES HAVE BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

    JORN LIND-NIELSEN SPECIFICALLY DISCLAIM ANY WARRANTIES, INCLUDING,
    BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
    FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE PROVIDED HEREUNDER IS
    ON AN "AS IS" BASIS, AND THE AUTHORS AND DISTRIBUTORS HAVE NO
    OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR
    MODIFICATIONS.
========================================================================*/

/*************************************************************************
  $Header: /cvsroot/buddy/buddy/src/cache.h,v 1.1.1.1 2004/06/25 13:22:34 haimcohen Exp $
  FILE:  cache.h
  DESCR: Cache class for caching apply/exist etc. results
  AUTH:  Jorn Lind
  DATE:  (C) june 1997
*************************************************************************/

#ifndef _CACHE_H
#define _CACHE_H

#include <stdint.h>

typedef struct
{
   union
   {
      uint64_t res;
   } r;
   int a,b,c,proj_end;
} BddCacheData;


typedef struct
{
   BddCacheData *table;
   int tablesize;
} BddCache;


extern int  BddCache_init(BddCache *, int);
extern void BddCache_done(BddCache *);
extern int  BddCache_resize(BddCache *, int);
extern void BddCache_reset(BddCache *);

#define BddCache_lookup(cache, hash) (&(cache)->table[hash % (cache)->tablesize])


#endif /* _CACHE_H */


/* EOF */
