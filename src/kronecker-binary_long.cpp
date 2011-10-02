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
#include <x86intrin.h>

#include "util.hpp"
#include "mpint.hpp"
#include "kronecker-jacobi.hpp"

namespace integer {

extern int tbl1[];

namespace impl {

/*
  Kronecker-binary
*/
int kronecker(const mpint::MPInt& in_x, const mpint::MPInt& in_y)
{
  using namespace mpint;

  MPInt x(in_x), y(in_y);

  // #1
  if (y.isZero()) {
    return x.size() == 1 && x[0] == 1 ? 1 : 0;
  }

  // #2
  if (! ((x[0] | y[0]) & 0x1)) {
    return 0;
  }
  size_t v = y.NTZ();
  y >>= v;
  int k; // return value.
  if (! (v & 0x1)) {
    k = 1;
  } else {
    k = tbl1[x[0] & 0x7];
  }
  if (y.isNeg() && x.isNeg()) {
    k = -k;
  }
  MPInt::absolute(y, y);

  // #3
#if 0
  x %= y;
#else
  // @note: avoid remainder operation, use recipro.
  if (x < 0) {
    if ((y[0] & 0x3) == 0x3) {
      k = -k;
    }
  }
  MPInt::absolute(x, x);
#endif

  assert(x >= 0);
  assert(y >= 0);

  MPInt r;
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
    size_t v = x.NTZ();
    x >>= v;
    if (v & 0x1) {
      k = tbl1[y[0] & 0x7] * k;
    }

    // #6
    r = y - x;
    if (r.isPos()) {
      if (x[0] & y[0] & 0x2) {
        k = -k;
      }
      y.swap(x);
      x.swap(r);
    } else {
      MPInt::absolute(x, r);
    }
  }
}

} // namespace impl

} // namespace integer
