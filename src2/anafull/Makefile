#////////////////////////////////////////////////////#
# * Makefile
#////////////////////////////////////////////////////#

## libfile ##
target=libanafull.a

## Complier and options ##
FC = ifort
FLAG = -O3 -fpp
#DBAG = -check all -std -gen_interfaces -fpe0 -ftrapuv -traceback
FLAGS = $(FLAG) $(DBAG)

## Directories ##
libroot = ../../
MyLib   = ${libroot}/
healpix = ${libroot}/pub/Healpix
cfitsio = ${libroot}/pub/cfitsio

MODD = -I$(healpix)/include -I$(MyLib)/mod
LIBD = -L$(cfitsio) -L$(healpix)/lib -L$(MyLib)/lib
LINK = -lutils -lhealpix -lcfitsio

OPTION = $(MODD) $(LIBD) $(LINK)

OBJ = anafull.o recfull.o delensfull.o

$(target): $(OBJ)
	ar cr $@ $^

%.o: %.f90
	$(FC) $(FLAGS) $(OPTION) -c $*.f90

install:
	mv *.mod ${MyLib}/mod; mv $(target) ${MyLib}/lib/

clean:
	rm -f *.o*

recfull.o : anafull.o

