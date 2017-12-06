CC = gcc
MPI_CC = mpicc
CFLAGS = 
CPOST_FLAGS  = -lm
NUMPROCS = 
seq: aquadsequential.o
	$(CC) $(CFLAGS) -o aquadsequential aquadsequential.c $(CPOST_FLAGS)

par:	aquadPartA.c stack.h stack.c
	$(MPI_CC) -o aquadPartA aquadPartA.c stack.h stack.c

exPar:
	mpirun -c $(NUMPROCS) ./aquadPartA

cleanSeq:
	rm aquadsequential.o aquadsequential

cleanPar:
	rm aquadPartA

clean:
	rm aquadsequential.o aquadsequential
	rm  aquadPartA.o stack.o aquadPartA
