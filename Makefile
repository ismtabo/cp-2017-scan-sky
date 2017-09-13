#
# Computación paralela (Grado Ingeniería Informática)
#
# Makefile: Contar cuerpos celestes
#

MPICC=mpicc
CFLAGS=-O3

EXES=ScanSky_mpiseq ScanSky_mpi

all: $(EXES)

clean: 
	rm -f $(EXES)

ScanSky_mpiseq: ScanSky_mpiseq.c cputils.h
	$(MPICC) $(CFLAGS) $< -o $@ 

ScanSky_mpi: ScanSky_mpi.c cputils.h
	$(MPICC) $(CFLAGS) $< -o $@ 
