/* -*- mode: c++; coding: utf-8-unix -*- */
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

#include <cstdint>
#include <iostream>
#include <string>
#include <sstream>
#include <xbyak/xbyak_util.h>

#include <gmpxx.h>
#define USE_GMP

#include "util.hpp"
#include "mpint.hpp"
#include "kronecker-jacobi.hpp"

using namespace ff_util;

const int N = 10000;

#define BENCHF "%s:\t% 10.2f clk\n"
#define GNUPLOTF " % 15.2f"

#define OUTPUT_GNUPLOT

void test_mpint()
{
  PUTSERR(__func__);

  using namespace mpint;

  {
    MPInt a(0), b(1), c(2);

    TEST_ASSERT(a.isZero());
    TEST_EQ(a, a);
    TEST_EQ(a, 0);
    TEST_EQ(a.size(), 0);
    TEST_ASSERT(! b.isZero());
    TEST_EQ(b, b);
    TEST_EQ(b, 1);
    TEST_ASSERT(! c.isZero());
    TEST_EQ(c, c);
    TEST_EQ(c, 2);
  }

  {
    const uint64_t ary[] = {
      0, 1, 2, 3, 4,
    };
    const uint64_t ary1[] = {
      5, 6, 7,
    };
    MPInt a(1), b(2), c(3), d(ary), e(ary1);
    TEST_EQ(uintptr_t(a.get()) & 0xf, 0);
    TEST_EQ(uintptr_t(b.get()) & 0xf, 0);
    TEST_EQ(uintptr_t(c.get()) & 0xf, 0);
    TEST_EQ(uintptr_t(d.get()) & 0xf, 0);
    TEST_EQ(uintptr_t(e.get()) & 0xf, 0);
  }

  {
    mpint::MPInt a("-0x14");
    TEST_EQ(a.capacity(), 1);
    TEST_EQ(a.size(), 1);
    TEST_EQ(a.sign(), -1);

    mpint::MPInt b;
    TEST_EQ(b.capacity(), 0);
    TEST_EQ(b.size(), 0);
    TEST_EQ(b.sign(), 0);
    b = a;
    TEST_EQ(b.capacity(), a.capacity());
    TEST_EQ(b.size(), a.size());
    TEST_EQ(b.sign(), a.sign());
    TEST_EQ(a, b);

    mpint::MPInt c(20);
    TEST_EQ(a.abs(), c);

    mpint::MPInt d = a.abs();
    TEST_EQ(c, d);

    a.swap(c);
    TEST_EQ(a, d);
  }

  {
    const uint64_t ary[] = {0, 1, 0};
    mpint::MPInt c(ary);
    TEST_EQ(c.capacity(), 3);
    TEST_EQ(c.size(), 2);
    TEST_EQ(c.sign(), 1);
  }

  {
    const uint64_t iary0[] = {0, 1, 0};
    const uint64_t iary1[] = {0, 1, 1};
    const mpint::MPInt ary[] = {
      mpint::MPInt(iary1, true),
      mpint::MPInt(iary0, true),
      mpint::MPInt(-123),
      mpint::MPInt(0),
      mpint::MPInt(1213),
      mpint::MPInt(iary0),
      mpint::MPInt(iary1),
    };

    for (auto& entry : ary) {
      if (entry.isEven()) {
        TEST_ASSERT(! entry.isOdd());
      } else {
        TEST_ASSERT(entry.isOdd());
      }
    }

    for (size_t i = 0; i < (sizeof(ary)/sizeof(ary[0])) - 1; ++i) {
      TEST_ASSERT(ary[i] == ary[i]);
      if (!(ary[i] == ary[i])) {
        PUT(i);
        PUT(ary[i]);
      }
    }

    for (size_t i = 0; i < (sizeof(ary)/sizeof(ary[0])) - 1; ++i) {
      for (size_t j = 0; j < i; ++j) {
        TEST_ASSERT(ary[j] < ary[i]);
        if (!(ary[j] < ary[i])) {
          PUT(j);
          PUT(i);
          PUT(ary[j]);
          PUT(ary[i]);
        }

        TEST_ASSERT(ary[j] <= ary[i]);
        if (!(ary[j] <= ary[i])) {
          PUT(j);
          PUT(i);
          PUT(ary[j]);
          PUT(ary[i]);
        }
      }
    }

    for (size_t i = 0; i < (sizeof(ary)/sizeof(ary[0])) - 1; ++i) {
      for (size_t j = 0; j < i; ++j) {
        TEST_ASSERT(ary[i] > ary[j]);
        if (!(ary[i] > ary[j])) {
          PUT(j);
          PUT(i);
          PUT(ary[j]);
          PUT(ary[i]);
        }

        TEST_ASSERT(ary[i] >= ary[j]);
        if (!(ary[i] >= ary[j])) {
          PUT(j);
          PUT(i);
          PUT(ary[j]);
          PUT(ary[i]);
        }
      }
    }
  }
}

void test_mpint_sub()
{
  PUTSERR(__func__);

  using namespace mpint;

  {
    MPInt a(0), b;

    MPInt::sub(b, a, a);
    TEST_ASSERT(b.isZero());

    MPInt::sub(b, a, -a);
    TEST_ASSERT(b.isZero());

    MPInt::sub(b, -a, -a);
    TEST_ASSERT(b.isZero());

    MPInt::sub(b, -a, a);
    TEST_ASSERT(b.isZero());
  }

  {
    MPInt a(0), b(1), c;

    MPInt::sub(c, a, b);
    TEST_EQ(c, -1);

    MPInt::sub(c, b, a);
    TEST_EQ(c, 1);
  }

  {
    MPInt a(1), b(1), c;

    MPInt::sub(c, a, b);
    TEST_ASSERT(c.isZero());

    MPInt::sub(c, b, a);
    TEST_ASSERT(c.isZero());
  }

  {
    mpint::MPInt e, f, g;

    e.set(0);
    f.set(0);
    mpint::MPInt::sub(g, e, f);
    TEST_ASSERT(g >= 0);
    TEST_ASSERT(g.isZero());
    mpint::MPInt::sub(g, f, e);
    TEST_ASSERT(g <= 0);
    TEST_ASSERT(g.isZero());
  }

  {
    mpint::MPInt e, f, g;

    e.set(10);
    f.set(10);
    mpint::MPInt::sub(g, e, f);
    TEST_ASSERT(g >= 0);
    TEST_ASSERT(g.isZero());
    mpint::MPInt::sub(g, f, e);
    TEST_ASSERT(g <= 0);
    TEST_ASSERT(g.isZero());

    e.set(-10);
    f.set(-10);
    mpint::MPInt::sub(g, e, f);
    TEST_ASSERT(g >= 0);
    TEST_ASSERT(g.isZero());
    mpint::MPInt::sub(g, f, e);
    TEST_ASSERT(g <= 0);
    TEST_ASSERT(g.isZero());
  }

  {
    mpint::MPInt e, f, g;

    e.set(100);
    f.set(5);
    mpint::MPInt::sub(g, e, f);
    TEST_ASSERT(g >= 0);
    mpint::MPInt::sub(g, f, e);
    TEST_ASSERT(g <= 0);

    e.set(-10);
    f.set(-300);
    mpint::MPInt::sub(g, e, f);
    TEST_ASSERT(g >= 0);
    mpint::MPInt::sub(g, f, e);
    TEST_ASSERT(g <= 0);
  }

  {
    const uint64_t arye[] = {0, 0, 1};
    mpint::MPInt e, f, g;

    e.set(arye);
    f.set(10);
    mpint::MPInt::sub(g, e, f);
    TEST_ASSERT(g >= 0);
    TEST_EQ(g.size(), 2);
    mpint::MPInt::sub(g, f, e);
    TEST_ASSERT(g <= 0);
    TEST_EQ(g.size(), 2);

    e.set(-20);
    f.set(arye, true);
    mpint::MPInt::sub(g, e, f);
    TEST_ASSERT(g >= 0);
    mpint::MPInt::sub(g, f, e);
    TEST_ASSERT(g <= 0);
  }

  {
    const uint64_t arye[] = { 5, 0, 100, 2, 3,};
    const uint64_t aryf[] = {10, 0,  10, 0,};
    mpint::MPInt e, f, g;

    e.set(arye);
    f.set(aryf);
    mpint::MPInt::sub(g, e, f);
    TEST_ASSERT(g >= 0);
    mpint::MPInt::sub(g, f, e);
    TEST_ASSERT(g <= 0);

    e.set(aryf, true);
    f.set(arye, true);
    mpint::MPInt::sub(g, e, f);
    TEST_ASSERT(g >= 0);
    mpint::MPInt::sub(g, f, e);
    TEST_ASSERT(g <= 0);
  }

  {
    const uint64_t arye[] = { 5, 0,  0, 0, 0, 0, 0, 1, 1, 0, 1, 0,};
    const uint64_t aryf[] = {10, 0, 10, 0,};
    mpint::MPInt e, f, g;

    e.set(arye);
    f.set(aryf);
    mpint::MPInt::sub(g, e, f);
    TEST_ASSERT(g >= 0);
    mpint::MPInt::sub(g, f, e);
    TEST_ASSERT(g <= 0);

    e.set(aryf, true);
    f.set(arye, true);
    mpint::MPInt::sub(g, e, f);
    TEST_ASSERT(g >= 0);
    mpint::MPInt::sub(g, f, e);
    TEST_ASSERT(g <= 0);
  }
}

void test_mpint_NTZ()
{
  PUTSERR(__func__);

  using namespace mpint;

  {
    MPInt a(1);

    size_t an = a.NTZ();
    TEST_ASSERT(a > 0);
    TEST_EQ(an, 0);
    a = -a;
    TEST_ASSERT(a < 0);
    an = a.NTZ();
    TEST_EQ(an, 0);
  }

  {
    MPInt a(2);

    size_t an = a.NTZ();
    TEST_ASSERT(a > 0);
    TEST_EQ(an, 1);
    a = -a;
    TEST_ASSERT(a < 0);
    an = a.NTZ();
    TEST_EQ(an, 1);
  }

  {
    MPInt a(3);

    size_t an = a.NTZ();
    TEST_ASSERT(a > 0);
    TEST_EQ(an, 0);
    a = -a;
    TEST_ASSERT(a < 0);
    an = a.NTZ();
    TEST_EQ(an, 0);
  }

  {
    MPInt a(4);

    size_t an = a.NTZ();
    TEST_ASSERT(a > 0);
    TEST_EQ(an, 2);
    a = -a;
    TEST_ASSERT(a < 0);
    an = a.NTZ();
    TEST_EQ(an, 2);
  }

  {
    const uint64_t ary[] = {0, 1};
    MPInt a(ary);

    size_t an = a.NTZ();
    TEST_ASSERT(a > 0);
    TEST_EQ(an, 64);
    a = -a;
    TEST_ASSERT(a < 0);
    an = a.NTZ();
    TEST_EQ(an, 64);
  }

  {
    const uint64_t ary[] = {1, 1};
    MPInt a(ary);

    size_t an = a.NTZ();
    TEST_ASSERT(a > 0);
    TEST_EQ(an, 0);
    a = -a;
    TEST_ASSERT(a < 0);
    an = a.NTZ();
    TEST_EQ(an, 0);
  }

  {
    const uint64_t ary[] = {16, 1};
    MPInt a(ary);

    size_t an = a.NTZ();
    TEST_ASSERT(a > 0);
    TEST_EQ(an, 4);
    a = -a;
    TEST_ASSERT(a < 0);
    an = a.NTZ();
    TEST_EQ(an, 4);
  }

  {
    const uint64_t ary[] = {1, 0};
    MPInt a(ary);

    size_t an = a.NTZ();
    TEST_ASSERT(a > 0);
    TEST_EQ(an, 0);
    a = -a;
    TEST_ASSERT(a < 0);
    an = a.NTZ();
    TEST_EQ(an, 0);
  }

  {
    const uint64_t ary[] = {32, 0};
    MPInt a(ary);

    size_t an = a.NTZ();
    TEST_ASSERT(a > 0);
    TEST_EQ(an, 5);
    a = -a;
    TEST_ASSERT(a < 0);
    an = a.NTZ();
    TEST_EQ(an, 5);
  }

  {
    const uint64_t ary[] = {0, 0, 1};
    MPInt a(ary);

    size_t an = a.NTZ();
    TEST_ASSERT(a > 0);
    TEST_EQ(an, 64*2);
    a = -a;
    TEST_ASSERT(a < 0);
    an = a.NTZ();
    TEST_EQ(an, 64*2);
  }

  {
    const uint64_t ary[] = {0, 0, 4};
    MPInt a(ary);

    size_t an = a.NTZ();
    TEST_ASSERT(a > 0);
    TEST_EQ(an, 64*2 + 2);
    a = -a;
    TEST_ASSERT(a < 0);
    an = a.NTZ();
    TEST_EQ(an, 64*2 + 2);
  }
}

void test_mpint_shr()
{
  PUTSERR(__func__);

  using namespace mpint;

  {
    MPInt a(123), b;

    MPInt::shr(b, a, 64);
    TEST_EQ(b.size(), 0);
    TEST_EQ(b.sign(), 0);
    TEST_ASSERT(b.isZero());

    MPInt::shr(b, a, 2);
    TEST_EQ(b.size(), 1);
    TEST_EQ(b.sign(), 1);
    TEST_ASSERT(! b.isZero());
    TEST_EQ(b, 30);

    MPInt::shr(b, a, 10);
    TEST_EQ(b.size(), 0);
    TEST_EQ(b.sign(), 0);
    TEST_ASSERT(b.isZero());
  }

  {
    MPInt a(-123), b;

    MPInt::shr(b, a, 64);
    TEST_EQ(b.size(), 0);
    TEST_EQ(b.sign(), 0);
    TEST_ASSERT(b.isZero());

    MPInt::shr(b, a, 2);
    TEST_EQ(b.size(), 1);
    TEST_EQ(b.sign(), -1);
    TEST_ASSERT(! b.isZero());
    TEST_EQ(b, -30);

    MPInt::shr(b, a, 10);
    TEST_EQ(b.size(), 0);
    TEST_EQ(b.sign(), 0);
    TEST_ASSERT(b.isZero());
  }

  {
    const uint64_t ary[] = {123, 123};
    MPInt a(ary), b;

    MPInt::shr(b, a, 64);
    TEST_EQ(b.size(), 1);
    TEST_EQ(b.sign(), 1);
    TEST_ASSERT(! b.isZero());

    MPInt::shr(b, a, 2);
    TEST_EQ(b.size(), 2);
    TEST_EQ(b.sign(), 1);
    TEST_ASSERT(! b.isZero());
    const uint64_t aryb[] = {0xc00000000000001eULL, 0x1eULL};
    TEST_EQ(b, MPInt(aryb));

    MPInt::shr(b, a, 10);
    TEST_EQ(b.size(), 1);
    TEST_EQ(b.sign(), 1);
    TEST_ASSERT(! b.isZero());
    TEST_EQ(b, 0x1ec0000000000000ULL);

    MPInt::shr(b, a, 90);
    TEST_EQ(b.size(), 0);
    TEST_EQ(b.sign(), 0);
    TEST_ASSERT(b.isZero());
  }

  {
    const uint64_t ary[] = {123, 123};
    MPInt a(ary, true), b;

    MPInt::shr(b, a, 64);
    TEST_EQ(b.size(), 1);
    TEST_EQ(b.sign(), -1);
    TEST_ASSERT(! b.isZero());

    MPInt::shr(b, a, 2);
    TEST_EQ(b.size(), 2);
    TEST_EQ(b.sign(), -1);
    TEST_ASSERT(! b.isZero());
    const uint64_t aryb[] = {0xc00000000000001eULL, 0x1eULL};
    TEST_EQ(b, MPInt(aryb, true));

    MPInt::shr(b, a, 10);
    TEST_EQ(b.size(), 1);
    TEST_EQ(b.sign(), -1);
    TEST_ASSERT(! b.isZero());
    TEST_EQ(b, -0x1ec0000000000000LL);

    MPInt::shr(b, a, 90);
    TEST_EQ(b.size(), 0);
    TEST_EQ(b.sign(), 0);
    TEST_ASSERT(b.isZero());
  }
}

void test_mpint_kronecker()
{
  PUTSERR(__func__);

  using namespace std;
  using namespace mpint;
  using namespace integer;

  {
    int64_t x = 0x6b;
    int64_t y = 0x119;
    MPInt mx(x);
    MPInt my(y);

    TEST_EQ(x, mx);
    TEST_EQ(y, my);

    int r = kronecker(x, y);
    int mr = impl::kronecker(mx, my);
  TEST_EQ(r, mr);
  }

  const unsigned long test_seed = 0;
  gmp_randclass rng(gmp_randinit_default);
  rng.seed(test_seed);

  const size_t numOfLoop = 100;
  const size_t multiLen = 100;
  const size_t offsetLen = 100;
  for (size_t i = 0; i < numOfLoop; ++i) {
    size_t l = multiLen*i + offsetLen;
    mpz_class gx = rng_odd(rng, l);
    mpz_class gy = rng_odd(rng, l);
    if (gx > gy) {
      mpz_swap(gx.get_mpz_t(), gy.get_mpz_t());
    }

    MPInt mx(gx), my(gy);

    int gr = mpz_kronecker(gx.get_mpz_t(), gy.get_mpz_t());
    TEST_LESSEQ(-1, gr);
    TEST_LESSEQ(gr, 1);

    int mr = impl::kronecker(mx, my);
    TEST_LESSEQ(-1, mr);
    TEST_LESSEQ(mr, 1);
    TEST_EQ(gr, mr);
  }
}

void bench_NTZ()
{
  printf("\n\n# %s\n", __func__);

  using namespace std;
  using namespace mpint;

  const unsigned long test_seed = 0;
  gmp_randclass rng(gmp_randinit_default);
  rng.seed(test_seed);

  const size_t numOfLoop = 100;
  const size_t multOfLen = 1000;
  const size_t offsetLen = 1000;
  {
    for (size_t i = 0; i < numOfLoop; ++i) {
      const size_t len = multOfLen * i + offsetLen;
#ifdef OUTPUT_GNUPLOT
      /*
        @note: Output is:
        length gmp_timing my_impl_timing
      */
      cout << len << " ";
#else
      PUT(len);
#endif

      mpz_class gx;
      mpz_ui_pow_ui(gx.get_mpz_t(), 2, len);
      MPInt mx;
      size_t gn, mn;
      mx.set(gx);

      {
        string gx_str, mx_str;
        gx_str = toString(gx);
        mx_str = mx.toString();
        TEST_EQ(gx_str, mx_str);
      }

      double mpz_time;
      {
        Xbyak::util::Clock clk;
        for (int j = 0; j < N; ++j) {
          clk.begin();
          gn = mpz_scan1(gx.get_mpz_t(), 0);
          clk.end();
          TEST_ASSERT(gn > 0);
        }
        mpz_time = (double)clk.getClock() / clk.getCount();
#ifdef OUTPUT_GNUPLOT
        printf(GNUPLOTF, mpz_time);
#else
        printf(BENCHF, "\tmpz_scan1", mpz_time);
#endif
      }

      {
        double mpint_time;
        Xbyak::util::Clock clk;
        for (int j = 0; j < N; ++j) {
          clk.begin();
          mn = MPInt::NumTrailZero(mx);
          clk.end();
          TEST_ASSERT(mn > 0);
        }
        mpint_time = (double)clk.getClock() / clk.getCount();
#ifdef OUTPUT_GNUPLOT
        printf(GNUPLOTF, mpint_time);
#else
        printf(BENCHF, "MPInt::NumTrailZero", mpint_time);
        printf("ratio:\t%f\n", mpint_time / mpz_time);
#endif
      }

      TEST_EQ(gn, mn);

#ifdef OUTPUT_GNUPLOT
      printf("\n");
#endif
    }
  }
}

void bench_shr()
{
  printf("\n\n# %s\n", __func__);

  using namespace std;
  using namespace mpint;

  const unsigned long test_seed = 0;
  gmp_randclass rng(gmp_randinit_default);
  rng.seed(test_seed);

  const size_t numOfLoop = 100;
  const size_t multOfLen = 1000;
  const size_t offsetLen = 1000;
  {
    for (size_t i = 0; i < numOfLoop; ++i) {
      const size_t len = multOfLen * i + offsetLen;
#ifdef OUTPUT_GNUPLOT
      /*
        @note: Output is:
        length gmp_timing my_impl_timing
      */
      cout << len << " ";
#else
      PUT(len);
#endif
      mpz_class gx = rng.get_z_bits(len);
      mpz_class gz;
      MPInt mx, mz;
      mx.set(gx);

      {
        string gx_str, mx_str;
        gx_str = toString(gx);
        mx_str = mx.toString();
        TEST_EQ(gx_str, mx_str);
      }

      const size_t shr_w = len >> 2;
      double mpz_time;
      {
        Xbyak::util::Clock clk;
        for (int j = 0; j < N; ++j) {
          clk.begin();
          mpz_tdiv_q_2exp(gz.get_mpz_t(), gx.get_mpz_t(), shr_w);
          /*
            mpn_rshift((gz.get_mpz_t())->_mp_d,
            (gx.get_mpz_t())->_mp_d,
            (gx.get_mpz_t())->_mp_size,
            len >> 2);
          */
          clk.end();
        }
        mpz_time = (double)clk.getClock() / clk.getCount();
#ifdef OUTPUT_GNUPLOT
        printf(GNUPLOTF, mpz_time);
#else
        printf(BENCHF, "mpz_tdiv_q_2exp", mpz_time);
#endif
      }

      {
        double mpint_time;
        Xbyak::util::Clock clk;
        for (int j = 0; j < N; ++j) {
          clk.begin();
          MPInt::shr(mz, mx, shr_w);
          clk.end();
        }
        mpint_time = (double)clk.getClock() / clk.getCount();
#ifdef OUTPUT_GNUPLOT
        printf(GNUPLOTF, mpint_time);
#else
        printf(BENCHF, "\tMPInt::shr", mpint_time);
        printf("ratio:\t%f\n", mpint_time / mpz_time);
#endif
      }

      {
        string gz_str,  mz_str;
        gz_str = toString(gz);
        mz_str = mz.toString();
        TEST_EQ(gz_str, mz_str);
      }

#ifdef OUTPUT_GNUPLOT
      puts("");
#endif
    }
  }
}

void bench_sub()
{
  printf("\n\n# %s\n", __func__);

  using namespace std;
  using namespace integer;
  using namespace mpint;

  const unsigned long test_seed = 0;
  gmp_randclass rng(gmp_randinit_default);
  rng.seed(test_seed);

  const size_t numOfLoop = 100;
  const size_t multOfLen = 1000;
  const size_t offsetLen = 1000;
  {
    for (size_t i = 0; i < numOfLoop; ++i) {
      const size_t len = multOfLen * i + offsetLen;
#ifdef OUTPUT_GNUPLOT
      /*
        @note: Output is:
        length gmp_timing my_impl_timing
      */
      cout << len << " ";
#else
      PUT(len);
#endif
      mpz_class gx = rng.get_z_bits(len);
      mpz_class gy = rng.get_z_bits(len);
      mpz_class gz;
      MPInt mx, my, mz;
      mx.set(gx);
      my.set(gy);

      {
        string gx_str, gy_str, mx_str, my_str;
        gx_str = toString(gx);
        gy_str = toString(gy);
        mx_str = mx.toString();
        my_str = my.toString();
        TEST_EQ(gx_str, mx_str);
        TEST_EQ(gy_str, my_str);
      }

      double mpz_time;
      {
        Xbyak::util::Clock clk;
        for (int j = 0; j < N; ++j) {
          clk.begin();
          mpz_sub(gz.get_mpz_t(), gx.get_mpz_t(), gy.get_mpz_t());
          clk.end();
        }
        mpz_time = (double)clk.getClock() / clk.getCount();
#ifdef OUTPUT_GNUPLOT
        printf(GNUPLOTF, mpz_time);
#else
        printf(BENCHF, "mpz_sub", mpz_time);
#endif
      }

      double mpint_time;
      {
        Xbyak::util::Clock clk;
        for (int j = 0; j < N; ++j) {
          clk.begin();
          MPInt::sub(mz, mx, my);
          clk.end();
        }
        mpint_time = (double)clk.getClock() / clk.getCount();
#ifdef OUTPUT_GNUPLOT
        printf(GNUPLOTF, mpint_time);
#else
        printf(BENCHF, "MPInt::sub", mpint_time);
        printf("ratio: %f\n", mpint_time / mpz_time);
#endif
      }

      {
        string gz_str,  mz_str;
        gz_str = toString(gz);
        mz_str = mz.toString();
        TEST_EQ(gz_str, mz_str);
      }

#ifdef OUTPUT_GNUPLOT
      puts("");
#endif
    }
  }
}

void bench_kronecker()
{
  printf("\n\n# %s\n", __func__);

  using namespace std;
  using namespace integer;
  using namespace mpint;

  const unsigned long test_seed = 0;
  gmp_randclass rng(gmp_randinit_default);
  rng.seed(test_seed);

  const size_t numOfLoop = 10;
  const size_t multOfLen = 100;
  const size_t offsetLen = 100;
  for (size_t i = 0; i < numOfLoop; ++i) {
    size_t len = multOfLen * i + offsetLen;
    mpz_class gx = rng_odd(rng, len);
    mpz_class gy = rng_odd(rng, len);
    if (gx > gy)
      mpz_swap(gx.get_mpz_t(), gy.get_mpz_t());
#ifdef OUTPUT_GNUPLOT
    /*
      @note: Output is:
      length gmp_timing my_impl_timing
    */
    cout << len << " ";
#else
    size_t x_len = mpz_sizeinbase(gx.get_mpz_t(), 2);
    size_t y_len = mpz_sizeinbase(gy.get_mpz_t(), 2);
    PUT(len);
    PUT(x_len);
    PUT(y_len);
#endif

    MPInt mx(gx), my(gy);

    double mpz_time, mpint_time;
    int gr;
    {
      Xbyak::util::Clock clk;
      for (int j = 0; j < N; ++j) {
        clk.begin();
        gr = mpz_kronecker(gx.get_mpz_t(), gy.get_mpz_t());
        clk.end();
      }
      TEST_LESSEQ(-1, gr);
      TEST_LESSEQ(gr, 1);
      mpz_time = (double)clk.getClock() / clk.getCount();
#ifdef OUTPUT_GNUPLOT
      printf(GNUPLOTF, mpz_time);
#else
      printf(BENCHF, "\tmpz_kronecker", mpz_time);
#endif
    }

    int mr;
    {
      Xbyak::util::Clock clk;
      for (int j = 0; j < N; ++j) {
        clk.begin();
        mr = impl::kronecker(mx, my);
        clk.end();
      }
      TEST_LESSEQ(-1, mr);
      TEST_LESSEQ(mr, 1);
      mpint_time = (double)clk.getClock() / clk.getCount();
#ifdef OUTPUT_GNUPLOT
      printf(GNUPLOTF, mpint_time);
#else
      printf(BENCHF, "impl::kronecker", mpint_time);
#endif
    }

#ifdef OUTPUT_GNUPLOT
    puts("");
#else
    TEST_EQ(gr, mr);
    printf("ratio: %f\n", mpint_time / mpz_time);
#endif
  }
}

void info_gmp()
{
  using namespace std;

  cerr << "GMP Version is " << gmp_version << endl
       << "number of bits in mp_limb is " << mp_bits_per_limb << endl;
}

void test_all()
{
  using namespace std;
  using namespace mpint;

  MPInt::codeGen(0);

  test_mpint();
  test_mpint_NTZ();
  test_mpint_shr();
  test_mpint_sub();
  test_mpint_kronecker();

  cout.flush();

  MPInt::codeGen(1);

  test_mpint();
  test_mpint_NTZ();
  test_mpint_shr();
  test_mpint_sub();
  test_mpint_kronecker();

  cout.flush();

  MPInt::codeGen();

  test_mpint();
  test_mpint_NTZ();
  test_mpint_shr();
  test_mpint_sub();
  test_mpint_kronecker();

  cout.flush();
}

void bench_for_gnuplot()
{
  using namespace std;
  using namespace mpint;

  MPInt::codeGen(0);
  bench_NTZ();
  bench_shr();
  bench_sub();
  bench_kronecker();

  MPInt::codeGen(1);
  bench_NTZ();
  bench_shr();
  bench_sub();
  bench_kronecker();

  MPInt::codeGen();
  bench_NTZ();
  bench_shr();
  bench_sub();
  bench_kronecker();
}

int main()
{
  using namespace std;
  using namespace mpint;

#ifndef NDEBUG
  cerr << "NDEBUG is undefined" << endl;
#endif
  cerr << "Number of sampling loop: " << N << endl;

  info_gmp();
  MPIntCodeGen();

  test_all();

  bench_for_gnuplot();

  return testsAreSucceeded() ? 0 : 1;
}
