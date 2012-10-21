include ../config.mk

MKDIR = mkdir
MKDIRFLAGS = -p

OBJFILES = kronecker-binary.o kronecker-binary_long.o \
	kronecker-jacobi.o mpint.o

TARGETs = ../lib/libint.a

all: $(TARGETs)

clean:
	$(RM) *.o $(TARGETs)

../lib/libint.a: $(OBJFILES)
	$(MKDIR) $(MKDIRFLAGS) ../lib
	-$(RM) $@
	$(AR) $(ARFLAGS) $@ $(OBJFILES)
