/* -*- coding: utf-8-unix -*- */
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

#ifndef KRONECKER_JACOBI_HPP
#define KRONECKER_JACOBI_HPP

#include <cstdint>

#include "mpint.hpp"

namespace integer {

int kronecker(int64_t, int64_t);

namespace impl {
int kronecker(const mpint::MPInt&, const mpint::MPInt&);
}

} // namespace integer

#endif // KRONECKER_JACOBI_HPP
