# make: MPI when mpicc is available; make MPI=0: serial; make MPI=1: require MPI.
MPI ?= auto
ifeq ($(MPI),auto)
MPI := $(if $(shell command -v mpicc 2>/dev/null),1,0)
endif
ifeq ($(origin CC),default)
CC := $(if $(filter 1,$(MPI)),mpicc,cc)
endif
ifeq ($(MPI),1)
CPPFLAGS += -DMPI_ON
endif
GSL_PREFIX ?= $(SIROCCO)
ifneq ($(strip $(GSL_PREFIX)),)
CPPFLAGS += -I$(GSL_PREFIX)/include
LDFLAGS += -L$(GSL_PREFIX)/lib
endif
CFLAGS ?= -O3 -g
CFLAGS += -std=gnu99 -Wall -Wextra
LDLIBS += -lgsl -lgslcblas -lm
.PHONY: all cak clean
all: cak_v4
cak: cak_v4
cak_v4: cak_v4.o
	$(CC) $(CFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@
cak_v4.o: cak_v4.c cak_v4.h
clean:
	rm -f cak_v4.o cak_v4
