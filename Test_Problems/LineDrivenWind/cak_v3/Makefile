INCLUDE = $(PYTHON)/include
LIB =  $(PYTHON)/lib
BIN =  $(PYTHON)/bin


CC = mpicc
CFLAGS = -O3 -m64 $(MPI_FLAG) -I$(INCLUDE)
LDFLAGS= -lm -m64  -L$(LIB) -lgsl -lgslcblas

cak_objects = cak_v3.o

MPI_FLAG = -DMPI_ON

cak: $(cak_objects)
	$(CC) $(CFLAGS) $(cak_objects) $(LDFLAGS) -o cak_v3
	
clean :
	rm -f *.o  *~

