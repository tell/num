/* -*- mode: c++; mode: flymake; coding: utf-8-unix -*- */

#ifndef MPINT_HPP
#define MPINT_HPP

#include <cstdint>
#include <vector>
#include <iomanip>
#include <string>
#include <sstream>
#include <memory>

#include <boost/operators.hpp>

#define USE_GMP
#ifdef USE_GMP
#include <gmpxx.h>
#endif

#include "util.hpp"

namespace mpint {

namespace interface {

template<class T>
struct Empty {};

template<class T, class E = Empty<T> >
struct shiftable
  : E
{
  inline T operator<<(const size_t n)
  { T z; T::shl(z, static_cast<const T&>(*this), n); return z;}

  inline T operator>>(const size_t n)
  { T z; T::shr(z, static_cast<const T&>(*this), n); return z;}

  inline T& operator<<=(const size_t n)
  {
    T& ref = static_cast<T&>(*this);
    T::shl(ref, ref, n);
    return ref;
  }

  inline T& operator>>=(const size_t n)
  {
    T& ref = static_cast<T&>(*this);
    T::shr(ref, ref, n);
    return ref;
  }
};

template<class T, class E = Empty<T> >
struct addsubmul
  : E
{
  friend inline T operator+(const T& lhs, const T& rhs)
  { T z; T::add(z, lhs, rhs); return z; }

  friend inline T operator-(const T& lhs, const T& rhs)
  { T z; T::sub(z, lhs, rhs); return z; }

  friend inline T operator*(const T& lhs, const T& rhs)
  { T z; T::mul(z, lhs, rhs); return z; }

  inline T& operator+=(const T& rhs)
  {
    T& ref = static_cast<T&>(*this);
    T::add(ref, ref, static_cast<T&>(rhs));
    return ref;
  }

  inline T& operator-=(const T& rhs)
  {
    T& ref = static_cast<T&>(*this);
    T::sub(ref, ref, static_cast<T&>(rhs));
    return ref;
  }

  inline T& operator*=(const T& rhs)
  {
    T& ref = static_cast<T&>(*this);
    T::mul(ref, ref, static_cast<T&>(rhs));
    return ref;
  }
};

} // namespace interface

// Macro definitions.

#define MPINT_SIGN_(x) (((x) == 0 ? 0 : ((x) >= 0 ? 1 : -1)))
#define MPINT_ABS_(x) ((x) >= 0 ? (x) : -(x))

class MPInt
  : public interface::shiftable< MPInt,
                                 interface::addsubmul< MPInt > >,
    private boost::equality_comparable< MPInt >,
    private boost::equality_comparable< MPInt, int64_t >,
    private boost::less_than_comparable< MPInt >,
    private boost::less_than_comparable< MPInt, int64_t >
{
public:
  typedef size_t capacity_t;
  typedef int64_t sign_size_t;
  // typedef size_t size_t;
  typedef int sign_t;

  typedef uint64_t value_type;
  typedef std::unique_ptr<value_type[]> buffer_ptr;
  //typedef std::unique_ptr<ALIGN_(16) value_type[]> buffer_ptr;

  capacity_t allocated_;
  sign_size_t sign_size_;
  buffer_ptr d_ptr_;

  MPInt()
    : allocated_(0), sign_size_(0), d_ptr_(nullptr)
  {}

  MPInt(const MPInt& x)
    : allocated_(x.allocated_), sign_size_(x.sign_size_),
      d_ptr_(new value_type[x.allocated_])
  {
    std::copy(x.d_ptr_.get(), x.d_ptr_.get() + x.allocated_, d_ptr_.get());
  }

  explicit MPInt(const int x)
    : allocated_(1), sign_size_(MPINT_SIGN_(x)), d_ptr_(new value_type[1])
  {
    d_ptr_[0] = MPINT_ABS_(x);
  }

  explicit MPInt(const int64_t x)
    : allocated_(1), sign_size_(MPINT_SIGN_(x)), d_ptr_(new value_type[1])
  {
    d_ptr_[0] = MPINT_ABS_(x);
  }

  template<size_t N>
  explicit MPInt(const uint64_t(& digits)[N], bool setNegative = false)
    : d_ptr_(nullptr)
  {
    set(digits, setNegative);
  }

  explicit MPInt(const std::string& str)
    : d_ptr_(nullptr)
  { set(str); }

  void reserve(const size_t n)
  {
    allocated_ = n;
    d_ptr_.reset(new value_type[allocated_]);
    clear();
  }

  void clear()
  {
    for (size_t i = 0; i < allocated_; ++i) {
      d_ptr_[i] = 0;
    }
    sign_size_ = 0;
  }

  void release()
  {
    d_ptr_.release();
    sign_size_ = 0;
    allocated_ = 0;
  }

  MPInt& operator=(const MPInt& x)
  {
    /*
      @note: copy assignment.
    */

    if (this == &x) return *this;

    allocated_ = x.allocated_;
    sign_size_ = x.sign_size_;

    d_ptr_.reset(new value_type[allocated_]);
    std::copy(x.d_ptr_.get(), x.d_ptr_.get() + allocated_, d_ptr_.get());

    return *this;
  }

  void swap(MPInt& x)
  {
    std::swap(allocated_, x.allocated_);
    std::swap(sign_size_, x.sign_size_);
    d_ptr_.swap(x.d_ptr_);
  }

private:
  /*
    concat 32:32 -> 64 bit int.
  */
  static uint64_t concat_(const uint32_t low, const uint32_t high)
  {
    return (((uint64_t)high) << 32) | low;
  }

  static void convert_(std::vector<uint64_t>& dst, const std::vector<uint32_t>& src)
  {
    const size_t n = src.size();
    if (n > 1) {
      for (size_t i = 1; i < n; i += 2) {
        dst.push_back(concat_(src[i - 1], src[i]));
      }
    }
    if (n & 0x1) {
      dst.push_back(concat_(src[n - 1], 0));
    }
  }

  /*
    @return: position of most significant non-zero digit.
  */
  size_t scan_size_() const
  {
    assert(d_ptr_);
    for (size_t i = (size_t)allocated_; i > 0; --i) {
      if (d_ptr_[i - 1] != 0) {
        return i;
      }
    }
    return 0;
  }
public:

  void set(int x)
  {
    allocated_ = 1;
    sign_size_ = x == 0 ? 0 : (x > 0 ? 1 : -1);
    d_ptr_.reset(new value_type[allocated_]);
    d_ptr_[0] = MPINT_ABS_(x);
  }

  void set(unsigned int x)
  {
    allocated_ = 1;
    sign_size_ = 1;
    d_ptr_.reset(new value_type[allocated_]);
    d_ptr_[0] = x;
  }

  void set(const value_type* digits, size_t n, bool setNegative = false)
  {
    allocated_ = n;

    d_ptr_.reset(new value_type[allocated_]);
    std::copy(digits, digits + allocated_, d_ptr_.get());

    sign_size_ = scan_size_();
    if (setNegative) {
      sign_size_ = -sign_size_;
    }
  }

  template<size_t N>
  void set(const value_type (&digits)[N], bool setNegative = false)
  {
    allocated_ = N;

    d_ptr_.reset(new value_type[allocated_]);
    std::copy(digits, digits + allocated_, d_ptr_.get());

    sign_size_ = scan_size_();
    if (setNegative) {
      sign_size_ = -sign_size_;
    }
  }

  void set(const std::vector<value_type>& digits, bool setNegative = false)
  {
    set(&digits[0], digits.size(), setNegative);
  }

  void set(const std::vector<uint32_t>& digits, bool setNegative = false)
  {
    std::vector<value_type> x;
    convert_(x, digits);
    set(x, setNegative);
  }

  void set(const std::string& str, int base = 16)
  {
    std::string t = str;
    bool isNeg = false;

    if (t.size() >= 2 && t[0] == '-') {
      isNeg = true;
      t = t.substr(1);
    }

    if (t.size() >= 2 && t[0] == '0') {
      switch (t[1]) {
      case 'x':
        base = 16;
        t = t.substr(2);
        break;

      default:
        throw std::invalid_argument("not support base in set(str)");
        break;
      }
    }

    switch (base) {
    case 16:
      {
        std::vector<uint32_t> x;
        while (!t.empty()) {
          size_t remain = std::min((int)t.size(), 8);
          char* endp;
          uint32_t v = (uint32_t) strtoul(&t[t.size() - remain], &endp, 16);
          if (*endp) throw std::invalid_argument("bad hex str");
          x.push_back(v);
          t = t.substr(0, t.size() - remain);
        }
        set(x, isNeg);
      }
      break;

    default:
      break;
    }
  }

#ifdef USE_GMP
  MPInt(const mpz_class& x)
    : d_ptr_(nullptr)
  {
    set(x);
  }

  void set(const mpz_class& x)
  {
    std::ostringstream oss;
    mpz_class t = ::abs(x);
    oss << t.get_str(16);
    set(oss.str(), 16);
    sign_size_ = scan_size_();
    if (x < 0) {
      sign_size_ = -sign_size_;
    }
  }
#endif

  bool operator==(const MPInt& rhs) const
  {
    if (sign_size_ != rhs.sign_size_) {
      return false;
    } else {
      const size_t n = (size_t)MPINT_ABS_(sign_size_);
      for (size_t i = 0; i < n; ++i) {
        if (d_ptr_[i] != rhs.d_ptr_[i]) {
          return false;
        }
      }
      return true;
    }
  }

  bool operator==(const int64_t rhs) const
  {
    return *this == MPInt(rhs);
  }

  /*
    if x == y, then 0,
        x > y, then +,
        x < y, then -.
  */
  static sign_t compare_SameSize(const MPInt& lhs, const MPInt& rhs)
  {
    assert(lhs.sign_size_ == rhs.sign_size_);
    const size_t n = (size_t)MPINT_ABS_(lhs.sign_size_);
    for (size_t i = n; i > 0; --i) {
      if (lhs.d_ptr_[i - 1] != rhs.d_ptr_[i - 1]) {
        return lhs.d_ptr_[i - 1] < rhs.d_ptr_[i - 1] ? -1 : 1;
      }
    }
    return 0;
  }

  /*
    if x == y, then 0,
        x > y, then +,
        x < y, then -.
  */
  static sign_t compare(const MPInt& lhs, const MPInt& rhs)
  {
    const sign_size_t d_size = lhs.sign_size_ - rhs.sign_size_;
    if (d_size != 0) {
      return d_size;
    }
    return compare_SameSize(lhs, rhs);
  }

  /*
    if x == y, then 0,
        x > y, then +,
        x < y, then -.
  */
  static sign_t compare(const MPInt& lhs, const int64_t rhs)
  {
    if (lhs.sign_size_ == 0) {
      assert(lhs == 0);
      return -rhs;

    } else if (MPINT_ABS_(lhs.sign_size_) == 1) {
      if (lhs.d_ptr_[0] == rhs) {
        return 0;
      } else {
        if (lhs.sign_size_ > 0) {
          return lhs.d_ptr_[0] > rhs ? 1 : -1;
        } else {
          return lhs.d_ptr_[0] > rhs ? -1 : 1;
        }
      }
    } else {
      assert(MPINT_ABS_(lhs.sign_size_) > 1);
      return lhs.sign_size_;
    }
  }
  sign_t cmp(const MPInt& x) const { return compare(*this, x); }
  bool operator<(const MPInt& rhs) const { return compare(*this, rhs) < 0; }
  bool operator<(const int64_t rhs) const { return compare(*this, rhs) < 0; }
  bool operator>(const int64_t rhs) const { return compare(*this, rhs) > 0; }

  bool isZero() const { return sign_size_ == 0; }
  bool isPos() const { return sign_size_ > 0; }
  bool isNeg() const { return sign_size_ < 0; }
  bool isOdd() const { return (! isZero()) && d_ptr_[0] & 0x1; }
  bool isEven() const { return isZero() || (! d_ptr_[0] & 0x1); }

  capacity_t capacity() const { return allocated_; }
  size_t size() const { return MPINT_ABS_(sign_size_); }
  sign_t sign() const { return MPINT_SIGN_(sign_size_); }
  const value_type& operator[](size_t i) const { return d_ptr_[i]; }
  value_type& operator[](size_t i) { return d_ptr_[i]; }
  const value_type* get() const { return d_ptr_.get(); }
  value_type* get() { return d_ptr_.get(); }

  std::string toString(int base = 16) const
  {
    std::ostringstream oss;
    switch (base) {
    default:
    case 16:
      {
        if (sign_size_ < 0) {
          oss << "-";
        }
        oss << "0x" << std::hex;
        const size_t n = MPINT_ABS_(sign_size_);
        if (n == 0) {
          oss << "0";
        } else {
          oss << (*this)[n - 1];
          if (n > 1) {
            for (size_t i = n - 1; i > 0; --i) {
              oss << std::setfill('0')
                  << std::setw(sizeof(value_type)*2)
                  << (*this)[i - 1];
            }
          }
        }
      }
    }
    return oss.str();
  }

  friend std::ostream& operator<<(std::ostream& os, const MPInt& x)
  {
    return os << x.toString();
  }

  static void absolute(MPInt& z, const MPInt& x)
  {
    z = x;
    if (z.sign_size_ < 0) {
      z.sign_size_ = -z.sign_size_;
    }
  }
  MPInt abs() const { MPInt z; absolute(z, *this); return z; }

  static void negation(MPInt& z, const MPInt& x)
  {
    z = x;
    z.sign_size_ = -x.sign_size_;
  }
  MPInt operator-() const { MPInt z; negation(z, *this); return z; }

  typedef size_t (*in_prop_op)(const value_type);

  static in_prop_op in_NumTrailZero1;

  /*
    @return: number of trailing zero.
  */
  static size_t NumTrailZero(const MPInt& x);
  size_t NTZ() const { return NumTrailZero(*this); }

  typedef bool (*in_shift_op)(value_type*, const value_type*, const size_t, const size_t);

  static in_shift_op in_shr_shift;
  static void shr(MPInt& z, const MPInt& x, const size_t n);

  typedef bool (*in_bin_op)(value_type*, const value_type*, const size_t, const value_type*, const size_t);

  static in_bin_op in_sub_nc;

  /*
    z = x - y
  */
  static void sub(MPInt& z, const MPInt& in_x, const MPInt& in_y);

  /*
    Call code generator.
  */
  static void codeGen(const int version = -1);
};

void MPIntCodeGen();
void MPIntCodeGen_();

} // namespace mpint

#endif // MPINT
