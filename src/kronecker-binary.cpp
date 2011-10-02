/* -*- mode: c++; mode: flymake; coding: utf-8-unix -*- */
/*
  Copyright (c) 2011-2011 Tadanori TERUYA (tell) <tadanori.teruya@gmail.com>

  Permission is hereby granted, free of charge, to any person
  obtaining a copy of this software and associated documentation files
  (the "Software"), to deal in the Software without restriction,
  including without limitation the rights to use, copy, modify, merge,
  publish, distribute, sublicense, and/or sell copies of the Software,
  and to permit persons to whom the Software is furnished to do so,
  subject to the following conditions:

  The above copyright notice and this permission notice shall be
  included in all copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
  MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
  BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
  ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
  CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.

  @license: The MIT license <http://opensource.org/licenses/MIT>
*/

#include <cassert>

#include "util.hpp"
#include "kronecker-jacobi.hpp"

namespace integer {

extern int tbl1[];

/*
  Kronecker-binary
*/
int kronecker(int64_t x, int64_t y)
{
  // #1
  if (y == 0) {
    return (x == 1 || x == -1) ? 1 : 0;
  }

  // #2
  if ((x % 2) == 0 && (y % 2) == 0) {
    return 0;
  }
  unsigned int v = 0;
  int k; // return value.
  while ((y % 2) == 0) {
    v += 1;
    y /= 2;
  }
  if ((v % 2) == 0) {
    k = 1;
  } else {
    k = tbl1[x & 0x7];
  }
  if (y < 0) {
    y = -y;
    if (x < 0) {
      k = -k;
    }
  }

  // #3
#if 0
  x %= y;
  // @note: if a < 0, then a % b is negative.
  if (x < 0) {
    x += y;
  }
#else
  // @note: avoid remainder operation.
  if (x < 0) {
    x = -x;
    if ((y & 0x3) == 3) {
      k = -k;
    }
  }
#endif

  assert(x >= 0);
  assert(y >= 0);

  for (;;) {
    // #4
    if (x == 0) {
      if (y > 1) {
        return 0;
      } else {
        return k;
      }
    }

    // #5
    v = 0;
    while ((x % 2) == 0) {
      v += 1;
      x /= 2;
    }
    if ((v % 2) == 1) {
      k = tbl1[y & 0x7] * k;
    }

    // #6
    int64_t r = y - x;
    if (r > 0) {
      if (x & y & 0x2) {
        k = -k;
      }
      y = x;
      x = r;
    } else {
      x = -r;
    }
  }
}

} // namespace integer
