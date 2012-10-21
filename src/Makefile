include ../config.mk

MKDIR = mkdir
MKDIRFLAGS = -p

OBJFILES = kronecker-binary.o kronecker-binary_long.o \
	kronecker-jacobi.o mpint.o

TARGET = ../lib/libint.a

all: $(TARGET)

../lib/libint.a: $(OBJFILES)
	$(MKDIR) $(MKDIRFLAGS) ../lib
	-$(RM) $@
	$(AR) $@ $(OBJFILES)

clean:
	$(RM) *.o $(TARGET)
