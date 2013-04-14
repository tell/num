
WARNFLAGS := -Wall -Wextra -Wabi -Wformat=2 -Wcast-qual \
	-Wcast-align -Wwrite-strings -Wfloat-equal -Wpointer-arith \
	-Wconversion -Wfatal-errors
# WARNFLAGS += -Winline

OPTFLAGS := -O3 -fomit-frame-pointer
# OPTFLAGS += --param inline-unit-growth=1000
# OPTFLAGS += --param max-inline-insns-single=1000
# OPTFLAGS += --param large-function-growth=1000
# OPTFLAGS += -funroll-loops

DEBUGFLAGS := -DNDEBUG
# DEBUGFLAGS +=-g3

INCLUDEFLAGS := -I../xbyak -I../include

CXXFLAGS = -std=c++11 -stdlib=libc++ -m64 -msse4.2 -pedantic \
	-fno-operator-names \
	$(WARNFLAGS) $(OPTFLAGS) $(DEBUGFLAGS) $(INCLUDEFLAGS)

ifeq ($(sell uname),Darwin)
CC = $(CXX)
endif

LDFLAGS := -lstdc++ -lgmp -lgmpxx

MAKEDEPEND := makedepend
MAKEDEPENDFLAGS = -mv -Y. -- $(CFLAGS) $(CXXFLAGS) -- *.c *.cpp

.PHONY: all install clean
