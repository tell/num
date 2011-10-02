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

#include <climits>
#include <x86intrin.h>

#include <xbyak/xbyak.h>

#include "mpint.hpp"

namespace mpint {

static inline size_t emu_in_NumTrailZero1_bsfq(const MPInt::value_type x)
{
  return __bsfq(x);
}

static inline size_t emu_in_NumTrailZero1_popcnt(const MPInt::value_type x)
{
  return _mm_popcnt_u64((~x)&(x-1));
}

size_t MPInt::NumTrailZero(const MPInt& x)
{
  typedef MPInt::value_type value_type;
  static const size_t nbits = sizeof(value_type) * CHAR_BIT;

  assert(! x.isZero()); // @note: requirement

#if 1
  const value_type* ptr = x.d_ptr_.get();
  for (;;) {
    if (*ptr != 0) {
      break;
    }
    ++ptr;
  }
  return (ptr - x.d_ptr_.get())*nbits + in_NumTrailZero1(*ptr);
#else

  size_t i = 0;
  const size_t n = x.size();
  for (;;) {
    if (x.d_ptr_[i] != 0) {
      break;
    }
    ++i;
  }
  return i*nbits + in_NumTrailZero1(x.d_ptr_[i]);
#endif
}

static inline bool emu_in_shr_shift(MPInt::value_type* z, const MPInt::value_type* x, const size_t xn, const size_t sn)
{
  typedef MPInt::value_type value_type;

  if (sn == 0) {
    std::copy(x, x + xn, z);
  } else {
    const size_t bit_w = sizeof(value_type) * CHAR_BIT;
    for (size_t i = 1; i < xn; i++) {
      z[i - 1] = (x[i - 1] >> sn) | (x[i] << (bit_w - sn));
    }
    z[xn - 1] >>= sn;
  }
  return z[xn - 1] == 0;
}

template<int N>
static inline void shr_SIMD(MPInt::value_type* const z, const MPInt::value_type* const x, const size_t xn)
{
  assert(! "under construction");

  const size_t q = N >> 3;
  const size_t r = N & 7;
  __m128i* vz = z;
  __m128i* vx = x;

  /*
    shift right if N >= 8.
  */
  if (q > 0) {
    for (size_t i = 0; i < xn; ++i) {
      vz[i] = _mm_alignr_epi8(vx[i + 1], vx[i], q);
    }
  }

  for (size_t i = 0; i < (xn - 1); ++i) {
    __m128i t = _mm_alignr_epi8(vz[i + 1], vz[i], q + 1);
    vz[i] = _mm_srli_epi64(vz[i], r);
    t = _mm_slli_epi64(t, 8 - r);
    vz[i] = _mm_or_si128(vz[i], t);
  }

  /*
    under construction from here.
  */
}

void MPInt::shr(MPInt& z, const MPInt& x, const size_t n)
{
  typedef MPInt::value_type value_type;

  const size_t digit_w = 6;
  const size_t digit_mask = 0x3f;

  const size_t move_d = n >> digit_w;
  const size_t move_shift = n & digit_mask;
  const size_t x_size = x.size();

  if (move_d >= x_size) {
    // if x_size == 0, here is always true.
    assert(x_size >= 0);
    z.clear();
    return;

  } else {
    assert(x_size > 0);

    if (z.allocated_ < x.allocated_) {
      z.allocated_ = x.allocated_;
      z.d_ptr_.reset(new value_type[z.allocated_]);
    }

    std::copy(x.d_ptr_.get() + move_d, x.d_ptr_.get() + x_size, z.d_ptr_.get());
    z.sign_size_ = x_size - move_d;
    assert(z.sign_size_ > 0);

    // #define INSPECT
#ifdef INSPECT
    PUT(x);
    PUT(x.sign_size_);
    PUT(move_shift);
#endif
    bool lastIsZero = in_shr_shift(z.d_ptr_.get(), z.d_ptr_.get(), z.sign_size_, move_shift);
#ifdef INSPECT
    PUT(lastIsZero);
    if (z.sign_size_ > 0) {
      PUT(z); // may be invalid form.
    }
#undef INSPECT
#endif
    if (lastIsZero) {
      --z.sign_size_;
    }

    // zero clear.
    if (z.allocated_ > (size_t)x.sign_size_) {
      for (size_t i = x.allocated_; i < z.allocated_; ++i) {
        z.d_ptr_[i] = 0;
      }
    }

    if (x.sign_size_ < 0) {
      z.sign_size_ = -z.sign_size_;
    }
  }
}

/*
  @require:
  x.abs() >= y.abs().
  x.size() == xn.
  y.size() == yn.
  x.size() >= y.size().
  z.capacity() >= x.capacity() >= y.capacity() >= n.
*/
static inline bool emu_in_sub_nc(MPInt::value_type* z, const MPInt::value_type* x, const size_t xn, const MPInt::value_type* y, const size_t yn)
{
  typedef MPInt::value_type value_type;
  typedef MPInt::value_type value_type;

  value_type c = 0;
  for (size_t i = 0; i < yn; ++i) {
    value_type yc = y[i] + c;
    if (yc < c) {
      // y[i] = value_type(0xf...f) and c = 1
      z[i] = x[i];
    } else {
      c = x[i] < yc ? 1 : 0;
      z[i] = x[i] - yc;
    }
  }

  for (size_t i = yn; i < xn; ++i) {
    value_type c0 = c;
    c = x[i] < c ? 1 : 0;
    z[i] = x[i] - c0;
    if (c == 0 && (i == xn - 1)) {
      return z[i] == 0;
    }
  }
  return z[xn - 1] == 0;
}

/*
  z = x - y
*/
void MPInt::sub(MPInt& z, const MPInt& in_x, const MPInt& in_y)
{
  if ((in_x.sign_size_ ^ in_y.sign_size_) >= 0) {
    assert(in_x.size() >= 0);
    assert(in_x.size() >= 0);
    assert(in_x.sign() * in_y.sign() >= 0);
    bool lastIsZero;
    bool isNeg = in_x.sign() < 0;
    MPInt x, y;
    absolute(x, in_x);
    absolute(y, in_y);
    assert(x.sign() >= 0);
    assert(y.sign() >= 0);

    if ((x.sign_size_ == y.sign_size_
         && compare_SameSize(x, y) < 0)
        || x.sign_size_ < y.sign_size_) {
      x.swap(y);
      isNeg = ! isNeg;
    }
    assert(x.size() >= y.size());
    assert(x >= y);
    // @note: x - y, does not generate carry.
    if (z.allocated_ < x.allocated_) {
      z.allocated_ = x.allocated_;
      z.d_ptr_.reset(new value_type[z.allocated_]);
    }
    assert(x.sign_size_ >= 0);
    assert(y.sign_size_ >= 0);

    // #define INSPECT
#ifdef INSPECT
    PUT(x);
    PUT(x.sign_size_);
    PUT(y);
    PUT(y.sign_size_);
#endif
    lastIsZero = in_sub_nc(z.d_ptr_.get(), x.d_ptr_.get(), x.sign_size_, y.d_ptr_.get(), y.sign_size_);
    z.sign_size_ = x.sign_size_;
#ifdef INSPECT
    PUT(lastIsZero);
    if (z.sign_size_ > 1) {
      PUT(z.d_ptr_[0]);
    }
#undef INSPECT
#endif

    // zero clear.
    // @note: it does not guarantee that not used parts are always zero.
    if (z.allocated_ > (size_t)x.sign_size_) {
      for (size_t i = x.sign_size_; i < z.allocated_; ++i) {
        z.d_ptr_[i] = 0;
      }
    }

    // @note: Check size must be changed or not.
    if (z.sign_size_ > 0 && lastIsZero) {
      --z.sign_size_;
    }

    if (isNeg) {
      z.sign_size_ = -z.sign_size_;
    }
    return;

  } else {
    assert(in_x.sign() != in_y.sign());
    throw std::invalid_argument("not implemented for different sign sub");
    return;
  }
}

/*
  Assignment internal functions.
*/

MPInt::in_prop_op MPInt::in_NumTrailZero1 =
  emu_in_NumTrailZero1_bsfq;
  // emu_in_NumTrailZero1_popcnt;
MPInt::in_shift_op MPInt::in_shr_shift = emu_in_shr_shift;
MPInt::in_bin_op MPInt::in_sub_nc = emu_in_sub_nc;

class MPIntCode : public Xbyak::CodeGenerator {
public:
  typedef MPInt::value_type value_type;
  typedef Xbyak::Reg64 Reg64;

private:

  /*
    @require:
    xn > 0.
    0 <= rcx (4th operand, denotes shift width) < 64.
  */
  void genEntry_in_shr_shift()
  {
    const int bytes = sizeof(value_type);
    assert(bytes == 8);

    const Reg64& pz = rdi;
    const Reg64& px = rsi;
    const Reg64& xn = rdx;
    // @note: 4th operand is rcx;
    //const Reg64& sw = rcx;

    // working registers.
    const Reg64& t0 = r8;

inLocalLabel();

    cmp(xn, 1);
    je(".xn == 1");

    // xn > 1.
    dec(xn);

    align(16);
L(".loop");
    mov(rax, ptr [px]);
    mov(t0, ptr [px + bytes]);
    shrd(rax, t0, cl);
    mov(ptr [pz], rax);
    lea(px, ptr [px + bytes]);
    lea(pz, ptr [pz + bytes]);
    dec(xn);
    jnz(".loop");

    // Loop is over.

L(".xn == 1");

    mov(t0, ptr [px]);
    shr(t0, cl);
    mov(ptr [pz], t0);
    // t0 has last result.

    cmp(t0, 0);
    mov(rax, 0);
    sete(al);

outLocalLabel();

    ret();
  }

  /*
    @require:
    xn > 0.
  */
  void genEntry_in_shr_shift_4()
  {
    fprintf(stderr, "%s\n", __func__);

    const int bytes = sizeof(value_type);
    assert(bytes == 8);

    const Reg64& pz = rdi;
    const Reg64& px = rsi;
    const Reg64& xn = rdx;
    // @note: 4th operand is rcx;
    //const Reg64& sw = rcx;

    // working registers.
    const Reg64& t0 = r8;
    const Reg64& t1 = r9;
    const Reg64& t2 = r10;
    const Reg64& t3 = r11;

    push(r12);
    const Reg64& t4 = r12;

    auto ret_proc = [&]() -> void
      {
        pop(r12);
        ret();
      };

inLocalLabel();

    mov(rax, xn);
    and(rax, 3); // rax <- xn % 4.
    shr(xn, 2);  // xn  <- xn / 4.

    cmp(xn, 0);
    jne(".xn/4 > 0");
    // xn < 4.
    // @note: rax != 0, rax == 1, 2, or 3.
    // there exists: px -> [0,(unknown)].

    cmp(rax, 1);
    jne("@f");
    // rax == 1.

    mov(t2, ptr [px]);
    shr(t2, cl);
    mov(ptr [pz], t2);
    // t2 has last result.
    cmp(t2, 0);
    mov(rax, 0);
    sete(al);
    ret_proc();

L("@@");

    cmp(rax, 2);
    jne("@f");
    // rax == 2.

    mov(t1, ptr [px]);
    mov(t2, ptr [px + bytes]);
    shrd(t1, t2, cl);
    shr(t2, cl);
    mov(ptr [pz], t1);
    mov(ptr [pz + bytes], t2);
    // t2 has last result.
    cmp(t2, 0);
    mov(rax, 0);
    sete(al);
    ret_proc();

L("@@");

    // rax == 3.

    mov(t0, ptr [px]);
    mov(t1, ptr [px + bytes]);
    mov(t2, ptr [px + bytes*2]);
    shrd(t0, t1, cl);
    shrd(t1, t2, cl);
    shr(t2, cl);
    mov(ptr [pz], t0);
    mov(ptr [pz + bytes], t1);
    mov(ptr [pz + bytes*2], t2);
    // t2 has last result.
    cmp(t2, 0);
    mov(rax, 0);
    sete(al);
    ret_proc();

L(".xn/4 > 0");
    // @note: quite different from here.

    // xn >= 4.
    // @note: rax == 0, 1, 2, or 3.
    // there exists: px -> [0,1,2,3,(unknown)].

    cmp(rax, 0);
    je(".xn mod 4 == 0");
    // rax != 0.

    cmp(rax, 1);
    jne("@f");
    // rax == 1, and then px -> [0,1,(unknown),4x,4x+1].

    mov(t2, ptr [px]);
    mov(t3, ptr [px + bytes]);
    shrd(t2, t3, cl);
    mov(ptr [pz], t2);
    // t2 has last result.
    lea(px, ptr [px + bytes]);
    lea(pz, ptr [pz + bytes]);
    jmp(".xn mod 4 == 0");

L("@@");

    cmp(rax, 2);
    jne("@f");
    // rax == 2, and then px -> [0,1,(unknown),4x,4x+1,4x+2].

    mov(t1, ptr [px]);
    mov(t2, ptr [px + bytes]);
    mov(t3, ptr [px + bytes*2]);
    shrd(t1, t2, cl);
    shrd(t2, t3, cl);
    mov(ptr [pz], t1);
    mov(ptr [pz + bytes], t2);
    // t2 has last result.
    lea(px, ptr [px + bytes*2]);
    lea(pz, ptr [pz + bytes*2]);
    jmp(".xn mod 4 == 0");

L("@@");

    // rax == 3, and then px -> [0,1,(unknown),4x,4x+1,4x+2,4x+3].

    mov(t0, ptr [px]);
    mov(t1, ptr [px + bytes]);
    mov(t2, ptr [px + bytes*2]);
    mov(t3, ptr [px + bytes*3]);
    shrd(t0, t1, cl);
    shrd(t1, t2, cl);
    shrd(t2, t3, cl);
    mov(ptr [pz], t0);
    mov(ptr [pz + bytes], t1);
    mov(ptr [pz + bytes*2], t2);
    // t2 has last result.
    lea(px, ptr [px + bytes*3]);
    lea(pz, ptr [pz + bytes*3]);

L(".xn mod 4 == 0");

    cmp(xn, 1);
    je(".xn/4 == 1");
    // @note: xn/4 > 1.
    // There exists, px -> [0,1,2,3,4,(unknown),4x-2,4x-1].

    dec(xn); // @note: extra decrement.

    align(16);
L(".loop for xn/4 > 1");
    mov(t0, ptr [px]);
    mov(t1, ptr [px + bytes]);
    mov(t2, ptr [px + bytes*2]);
    mov(t3, ptr [px + bytes*3]);
    mov(t4, ptr [px + bytes*4]);
    shrd(t0, t1, cl);
    shrd(t1, t2, cl);
    shrd(t2, t3, cl);
    shrd(t3, t4, cl);
    mov(ptr [pz], t0);
    mov(ptr [pz + bytes], t1);
    mov(ptr [pz + bytes*2], t2);
    mov(ptr [pz + bytes*3], t3);
    lea(px, ptr [px + bytes*4]);
    lea(pz, ptr [pz + bytes*4]);
    // t3 has last result.
    dec(xn);
    jnz(".loop for xn/4 > 1");

    // @note: xn == 0,
    // But, there exists, px -> [0,1,2,3].
    // by extra decrement.

L(".xn/4 == 1");

    // @note: xn/4 == 1.
    // There exists, px -> [0,1,2,3].

    mov(t0, ptr [px]);
    mov(t1, ptr [px + bytes]);
    mov(t2, ptr [px + bytes*2]);
    mov(t3, ptr [px + bytes*3]);
    shrd(t0, t1, cl);
    shrd(t1, t2, cl);
    shrd(t2, t3, cl);
    shr(t3, cl);
    mov(ptr [pz], t0);
    mov(ptr [pz + bytes], t1);
    mov(ptr [pz + bytes*2], t2);
    mov(ptr [pz + bytes*3], t3);
    // t3 has last result.

    cmp(t3, 0);
    mov(rax, 0);
    sete(al);

outLocalLabel();

    ret_proc();
  }

  /*
    @require:
    xn >= yn.
    px[] - py[] never generate carry.

    @return: pz[xn - 1] == 0.

    @note: if xn == yn == 0, then never modify memory fields referred by pz.
  */
  void genEntry_in_sub_nc()
  {
    const int bytes = sizeof(value_type);
    assert(bytes == 8);

    const Reg64& pz = rdi;
    const Reg64& px = rsi;
    const Reg64& xn = rdx;
    const Reg64& py = rcx;
    const Reg64& yn = r8;

    // working space.
    const Reg64& t0 = r9;

inLocalLabel();

    sub(xn, yn);
    xor(rax, rax); // @note: clear CF.
    cmp(yn, 0);
    je(".EndSubXY");

    mov(t0, ptr [px]);
    sub(t0, ptr [py]);
    mov(ptr [pz], t0);
    lea(px, ptr [px + bytes]);
    lea(py, ptr [py + bytes]);
    lea(pz, ptr [pz + bytes]);
    dec(yn);
    jz(".EndSubXY");

    align(16);
L(".LoopSubXY");
    mov(t0, ptr [px]);
    sbb(t0, ptr [py]);
    mov(ptr [pz], t0);
    lea(px, ptr [px + bytes]);
    lea(py, ptr [py + bytes]);
    lea(pz, ptr [pz + bytes]);
    dec(yn);
    jnz(".LoopSubXY");

    align(16);
L(".EndSubXY"); // @note: run over of y.
    setc(al); // @note: save carry.
    movzx(rax, al);
    cmp(xn, 0);
    je(".End");

    align(16);
L(".StartSubX");
    mov(t0, ptr [px]);
    sub(t0, rax); // @note: saved carry propagation.
    mov(ptr [pz], t0);
    lea(pz, ptr [pz + bytes]);
    lea(px, ptr [px + bytes]);
    dec(xn);
    jz(".End");

    align(16);
L(".LoopSubX");
    mov(t0, ptr [px]);
    sbb(t0, 0);
    mov(ptr [pz], t0);
    lea(pz, ptr [pz + bytes]);
    lea(px, ptr [px + bytes]);
    dec(xn);
    jnz(".EndSubXY");

    align(16);
L(".End");
    cmp(t0, 0);
    sete(al);
    movzx(rax, al);

outLocalLabel();

    ret();
  }

  /*
    @require:
    xn >= yn > 0.
    px[] - py[] never generate carry.

    @return: pz[xn - 1] == 0.
  */
  void genEntry_in_sub_nc_4()
  {
    fprintf(stderr, "%s\n", __func__);

    const int bytes = sizeof(value_type);
    assert(bytes == 8);

    const Reg64& pz = rdi;
    const Reg64& px = rsi;
    const Reg64& xn = rdx;
    const Reg64& py = rcx;
    const Reg64& yn = r8;

    // working space.
    const Reg64& t0 = r9;
    const Reg64& t1 = r10;
    const Reg64& t2 = r11;

    push(r12);
    const Reg64& test = r12;
    mov(test, 0);

    auto ret_proc = [&]() -> void
      {
        pop(r12);
        ret();
      };

inLocalLabel();

    sub(xn, yn);

    cmp(yn, 4);
    jge(".yn >= 4");

    // yn < 4.

    xor(t2, t2);
    cmp(yn, 0);
    je(".yn == 0", T_NEAR);

    mov(t2, ptr [px]);
    sub(t2, ptr [py]);
    mov(ptr [pz], t2);
    lea(px, ptr [px + bytes]);
    lea(py, ptr [py + bytes]);
    lea(pz, ptr [pz + bytes]);
    dec(yn);
    // t2 has most significant value.
    jz(".yn == 0", T_NEAR);

    align(16);
L(".loop for yn < 4");
    mov(t2, ptr [px]);
    sbb(t2, ptr [py]);
    mov(ptr [pz], t2);
    lea(px, ptr [px + bytes]);
    lea(py, ptr [py + bytes]);
    lea(pz, ptr [pz + bytes]);
    dec(yn);
    // t2 has most significant value.
    jnz(".loop for yn < 4");
    jmp(".yn == 0", T_NEAR);

L(".yn >= 4");
    mov(rax, yn);
    shr(yn, 2); // yn <- yn / 4.
    and(rax, 3); // rax <- rax % 4, and CF <- 0.
    jz(".yn mod 4 = 0", T_NEAR);

    cmp(rax, 1);
    jne("@f");
    // rax == 1.
    mov(t2, ptr [px]);
    sub(t2, ptr [py]);
    mov(ptr [pz], t2);
    lea(px, ptr [px + bytes]);
    lea(py, ptr [py + bytes]);
    lea(pz, ptr [pz + bytes]);
    // t2 has most significant value.
    jmp(".yn mod 4 = 0");

L("@@");

    cmp(rax, 2);
    jne("@f");
    // rax == 2.
    mov(t1, ptr [px]);
    mov(t2,  ptr [px + bytes]);
    sub(t1, ptr [py]);
    sbb(t2,  ptr [py + bytes]);
    mov(ptr [pz],         t1);
    mov(ptr [pz + bytes], t2);
    lea(px, ptr [px + bytes*2]);
    lea(py, ptr [py + bytes*2]);
    lea(pz, ptr [pz + bytes*2]);
    // t2 has most significant value.
    jmp(".yn mod 4 = 0");

L("@@");

    // rax == 3.
    mov(t0, ptr [px]);
    mov(t1, ptr [px + bytes]);
    mov(t2, ptr [px + bytes*2]);
    sub(t0, ptr [py]);
    sbb(t1,  ptr [py + bytes]);
    sbb(t2,  ptr [py + bytes*2]);
    mov(ptr [pz],           t0);
    mov(ptr [pz + bytes],   t1);
    mov(ptr [pz + bytes*2], t2);
    lea(px, ptr [px + bytes*3]);
    lea(py, ptr [py + bytes*3]);
    lea(pz, ptr [pz + bytes*3]);
    // t2 has most significant value.

    align(16);
L(".yn mod 4 = 0");

    mov(rax, ptr [px]);
    mov(t0,  ptr [px + bytes]);
    mov(t1,  ptr [px + bytes*2]);
    mov(t2,  ptr [px + bytes*3]);
    sbb(rax, ptr [py]);
    sbb(t0,  ptr [py + bytes]);
    sbb(t1,  ptr [py + bytes*2]);
    sbb(t2,  ptr [py + bytes*3]);
    mov(ptr [pz],          rax);
    mov(ptr [pz + bytes],   t0);
    mov(ptr [pz + bytes*2], t1);
    mov(ptr [pz + bytes*3], t2);
    lea(px, ptr [px + bytes*4]);
    lea(py, ptr [py + bytes*4]);
    lea(pz, ptr [pz + bytes*4]);
    dec(yn);
    // t2 has most significant value.
    jnz(".yn mod 4 = 0");

L(".yn == 0");

    // py is not needed.
    // yn is not needed.

    setc(al);
    movzx(rax, al); // save carry to rax.

    cmp(xn, 4);
    jge(".xn >= 4");

    // t2 has mast significant value.
    cmp(xn, 0);
    je(".xn == 0", T_NEAR);

    mov(t2, ptr [px]);
    sub(t2, rax); // rax has saved carry propagation.
    mov(ptr [pz], t2);
    lea(px, ptr [px + bytes]);
    lea(pz, ptr [pz + bytes]);
    dec(xn);
    // t2 has most significant value.
    jz(".xn == 0", T_NEAR);

    align(16);
L(".loop on xn with xn < 4");
    mov(t2, ptr [px]);
    sbb(t2, 0);
    mov(ptr [pz], t2);
    lea(px, ptr [px + bytes]);
    lea(pz, ptr [pz + bytes]);
    dec(xn);
    // t2 has most significant value.
    jnz(".loop on xn with xn < 4");
    jmp(".xn == 0", T_NEAR);

L(".xn >= 4");
    // rax keep carry.
    mov(yn, xn); // yn <- xn.
    shr(xn, 2); // xn <- xn / 4.
    and(yn, 3); // yn <- yn % 4, and CF <- 0.
    jnz("@f");

    // yn == 0, propagate rax carry.
    mov(t1, 0);
    sub(t1, rax); // set CF from rax.
    jmp(".xn mod 4 == 0");

L("@@");

    cmp(yn, 1);
    jne("@f");
    // yn == 1.
    mov(t2, ptr [px]);
    sub(t2, rax); // rax has carry.
    mov(ptr [pz], t2);
    lea(px, ptr [px + bytes]);
    lea(pz, ptr [pz + bytes]);
    // t2 has most significant value.
    jmp(".xn mod 4 == 0");

L("@@");

    cmp(yn, 2);
    jne("@f");
    // yn == 2.
    mov(t1, ptr [px]);
    mov(t2, ptr [px + bytes]);
    sub(t1, rax); // rax has carry.
    sbb(t2, 0);
    mov(ptr [pz], t1);
    mov(ptr [pz + bytes], t2);
    lea(px, ptr [px + bytes*2]);
    lea(pz, ptr [pz + bytes*2]);
    // t2 has most significant value.
    jmp(".xn mod 4 == 0");

L("@@");

    // yn == 3.
    mov(t0, ptr [px]);
    mov(t1, ptr [px + bytes]);
    mov(t2, ptr [px + bytes*2]);
    sub(t0, rax); // rax has carry.
    sbb(t1, 0);
    sbb(t2, 0);
    mov(ptr [pz], t0);
    mov(ptr [pz + bytes], t1);
    mov(ptr [pz + bytes*2], t2);
    lea(px, ptr [px + bytes*3]);
    lea(pz, ptr [pz + bytes*3]);
    // t2 has most significant value.
    jmp(".xn mod 4 == 0");

    align(16);
L(".xn mod 4 == 0");

    mov(yn, ptr [px]);
    mov(t0, ptr [px + bytes]);
    mov(t1, ptr [px + bytes*2]);
    mov(t2, ptr [px + bytes*3]);
    sbb(yn, 0);
    sbb(t0, 0);
    sbb(t1, 0);
    sbb(t2, 0);
    mov(ptr [pz],           yn);
    mov(ptr [pz + bytes],   t0);
    mov(ptr [pz + bytes*2], t1);
    mov(ptr [pz + bytes*3], t2);
    lea(px, ptr [px + bytes*4]);
    lea(pz, ptr [pz + bytes*4]);
    dec(xn);
    // t2 has most significant value.
    jnz(".xn mod 4 == 0");

L(".xn == 0");
    // t2 must have most significant value.

    // evaluated loop counter,
    // so, flags do not correspond to last result is zero or not.
    cmp(t2, 0);
    mov(rax, 0);
    sete(al);

outLocalLabel();

    ret_proc();
  }

  void genDemo_andWithFlag()
  {
    xor(rax, rax);
    sub(rax, 1);
    mov(rax, 0);
    // and(rax, 1);

    mov(rax, 0);
    setc(al);
    ret();
  }

  void genDemo_pushPop()
  {
    const Reg64& op1 = rdi;
    push(op1);
    pop(rax);
    ret();
  }

  MPInt::in_shift_op code_shr_;
  MPInt::in_bin_op code_sub_;
  MPInt::in_shift_op code_shr_4_;
  MPInt::in_bin_op code_sub_4_;

public:

  void demo_andWithEflags()
  {
    auto func = (uint64_t (*)()) getCurr();
    genDemo_andWithFlag();
    align(16);
    std::cout << "demo " << __func__ << " : " << func() << std::endl;
  }

  void demo_pushPop()
  {
    auto func = (uint64_t (*)(const int)) getCurr();
    genDemo_pushPop();
    align(16);
    std::cout << "demo " << __func__ << " : " << func(10) << std::endl;
  }

  MPIntCode()
    : code_shr_(MPInt::in_shr_shift),
      code_sub_(MPInt::in_sub_nc)
  {
    assert((uintptr_t(getCurr()) & 0xf) == 0);

    code_shr_ = (MPInt::in_shift_op) getCurr();
    genEntry_in_shr_shift();
    align(16);
    assert((uintptr_t(getCurr()) & 0xf) == 0);

    code_shr_4_ = (MPInt::in_shift_op) getCurr();
    genEntry_in_shr_shift_4();
    align(16);
    assert((uintptr_t(getCurr()) & 0xf) == 0);

    MPInt::in_shr_shift = code_shr_4_;

    code_sub_ = (MPInt::in_bin_op) getCurr();
    genEntry_in_sub_nc();
    align(16);
    assert((uintptr_t(getCurr()) & 0xf) == 0);

    code_sub_4_ = (MPInt::in_bin_op) getCurr();
    genEntry_in_sub_nc_4();
    align(16);
    assert((uintptr_t(getCurr()) & 0xf) == 0);

    MPInt::in_sub_nc = code_sub_4_;

    MPIntCodeGen_();
  }

  void setGenedCode(const int version = -1) const
  {
    using namespace std;

    switch (version) {
    case -1:
    default:
      {
        MPInt::in_shr_shift = code_shr_4_;
        MPInt::in_sub_nc = code_sub_4_;
      }
      break;

    case 1:
      {
        MPInt::in_shr_shift = code_shr_;
        MPInt::in_sub_nc = code_sub_;
      }
      break;

    case 0:
      {
        MPInt::in_shr_shift = emu_in_shr_shift;
        MPInt::in_sub_nc = emu_in_sub_nc;
      }
    break;
    }

#if 1
    ostringstream oss;
    oss << hex;
    oss << "MPInt::in_shr_shift=" << uintptr_t(MPInt::in_shr_shift) << endl;
    oss << hex;
    oss << "MPInt::in_sub_nc=" << uintptr_t(MPInt::in_sub_nc) << endl;
    cerr << oss.str();
#endif
  }
};

static const MPIntCode& makeCodeGen()
{
  static const MPIntCode code;
  return code;
}

void MPInt::codeGen(const int version)
{
  const MPIntCode& code = makeCodeGen();
  code.setGenedCode(version);
}

void MPIntCodeGen()
{
  static bool isInited = false;
  if (isInited) return;
  isInited = true;

  try {
    fprintf(stderr, "MPIntCodeGen ");
#ifdef XBYAK32
#error "32bit is not supported"
#elif XBYAK64_WIN
#error "Windows is not supported"
#else
    fprintf(stderr, "64bit\n");

    makeCodeGen();
    return;
#endif
  } catch (Xbyak::Error err) {
    fprintf(stderr, "Xbyak ERROR: %s (%d)\n", Xbyak::ConvertErrorToString(err), err);
  } catch (...) {
    fprintf(stderr, "ERROR: unkown error\n");
  }
  ::exit(1);
}

void MPIntCodeGen_()
{
  fprintf(stderr, "MPIntCodeGen generated");
}

} // namespace mpint
