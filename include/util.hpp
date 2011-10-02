/* mode: c++; mode: flymake; coding: utf-8-unix*/

#ifndef UTIL_HPP
#define UTIL_HPP

#include <iostream>
#include <sstream>

#if defined(__GNUC__)
#define ALIGN_(x) __attribute__((aligned(x)))
#endif

namespace {
int ctrTest = 0, ctrErr = 0;

bool testsAreSucceeded()
{
  std::cerr << "tests=" << ctrTest << std::endl
            << "errors=" << ctrErr << std::endl;
  return ctrErr == 0;
}
} // namespace

#define PUT(e) do {                             \
    std::cout << #e << "=" << (e) << std::endl; \
  } while (0)

#define PUTSERR(e) do {            \
    std::cerr << (e) << std::endl; \
  } while (0)

#define PUTHEX(e) do {                          \
    std::ostringstream oss;                     \
    oss << std::hex;                            \
    oss << #e << "=" << (e) << std::endl;       \
    std::cout << oss.str();                     \
  } while (0)

#define PUTHEXERR(e) do {                       \
    std::ostringstream oss;                     \
    oss << std::hex;                            \
    oss << #e << "=" << (e) << std::endl;       \
    std::cerr << oss.str();                     \
  } while (0)

#define TEST_ASSERT(p) do {                             \
    ++ctrTest;                                          \
    if (! (p)) {                                        \
      ++ctrErr;                                         \
      printf("%s:%d: error: Assertion failed: %s\n",    \
             __FILE__, __LINE__, #p);                   \
    }                                                   \
  } while (0)

#define TEST_EQ(x,y) do {                                   \
    ++ctrTest;                                              \
    if (!((x) == (y))) {                                    \
      ++ctrErr;                                             \
      printf("%s:%d: error: %s != %s\n",                    \
             __FILE__, __LINE__, #x, #y);                   \
      std::cout << "lhs=" << (x) << std::endl               \
                << "rhs=" << (y) << std::endl;              \
    }                                                       \
  } while (0)

#define TEST_LESS(x,y) do {                                 \
    ++ctrTest;                                              \
    if (!((x) < (y))) {                                     \
      ++ctrErr;                                             \
      printf("%s:%d: error: %s < %s\n",                    \
             __FILE__, __LINE__, #x, #y);                   \
      std::cout << "lhs=" << (x) << std::endl               \
                << "rhs=" << (y) << std::endl;              \
    }                                                       \
  } while (0)

#define TEST_LESSEQ(x,y) do {                                  \
    ++ctrTest;                                                 \
    if (!((x) <= (y))) {                                       \
      ++ctrErr;                                                \
      printf("%s:%d: error: %s <= %s\n",                       \
             __FILE__, __LINE__, #x, #y);                      \
      std::cout << "lhs=" << (x) << std::endl                  \
                << "rhs=" << (y) << std::endl;                 \
    }                                                          \
  } while (0)

namespace ff_util {

#ifdef USE_GMP

namespace {

std::string toString(const mpz_class& x)
{
  using namespace std;

  ostringstream oss;
  if (x < 0) {
    oss << "-0x";
  } else {
    oss << "0x";
  }
  oss << mpz_class(abs(x)).get_str(16);
  return oss.str();
}

mpz_class rng_odd(gmp_randclass& rng, const size_t l)
{
  using namespace std;

  mpz_class g = rng.get_z_bits(l);
  if (mpz_even_p(g.get_mpz_t())) {
    g += 1;
  }
  return g;
}

}

#endif // USE_GMP

} // namespace ff

#endif /* UTIL_HPP */
