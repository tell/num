include ../config.mk

MKDIR = mkdir
MKDIRFLAGS = -p

OBJFILES = kronecker-binary.o kronecker-binary_long.o \
	kronecker-jacobi.o mpint.o

TARGETs = ../lib/libint.a

all: $(TARGETs)

clean:
	$(RM) *.o $(TARGETs)

depend:
	$(MAKEDEPEND) $(MAKEDEPENDFLAGS)


../lib/libint.a: $(OBJFILES)
	$(MKDIR) $(MKDIRFLAGS) ../lib
	-$(RM) $@
	$(AR) $(ARFLAGS) $@ $(OBJFILES)

# DO NOT DELETE

kronecker-binary.o: ../include/util.hpp ../include/kronecker-jacobi.hpp
kronecker-binary.o: ../include/mpint.hpp
kronecker-binary_long.o: ../include/util.hpp ../include/mpint.hpp
kronecker-binary_long.o: ../include/kronecker-jacobi.hpp
kronecker-jacobi.o: ../include/util.hpp ../include/kronecker-jacobi.hpp
kronecker-jacobi.o: ../include/mpint.hpp
mpint.o: ../xbyak/xbyak/xbyak.h ../xbyak/xbyak/xbyak_bin2hex.h
mpint.o: ../xbyak/xbyak/xbyak_mnemonic.h ../include/mpint.hpp
mpint.o: ../include/util.hpp
