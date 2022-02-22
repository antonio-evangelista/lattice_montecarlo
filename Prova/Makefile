TARGET = lattice

FC = gfortran

#il compilatore si ferma sempre dopo tre errori se per caso dovesse
#sbroccare e farne vedere dumila
#la lettera o maiusc = ottimizzazione, buon livello è 2, per far andare
#il codice più veloce
FFLAGS = -fmax-errors=3 -g -O0

#GCC = gcc
#qui commento con il cancelletto
#unico posto dove devo cambiare le cose è il nome dei sorgenti e nome TARGET
#che è il nome del file
#possiamo pure dirgli, prendi tutti i file .f90 della cartella
SOURCES_MAIN = precision.f90 function.f90 solver.f90 main.f90
SOURCES_TEST = precision.f90 test.f90

OBJECTS_MAIN = $(SOURCES_MAIN:.f90=.o)
OBJECTS_TEST = $(SOURCES_TEST:.f90=.o)

.SUFFIXES:
.SUFFIXES: .f90 .o #.c

#scriviamo la regola per passare da .f90 -> .o
.f90.o:
	$(FC) $(FFLAGS) -c $<

#.c.o:
#uguale a sopra ma con gcc

#valuta gfortran -o mach precision.o machine.o
main: $(OBJECTS_MAIN)
	$(FC) -o $(TARGET) $(OBJECTS_MAIN)

test: $(OBJECTS_TEST)
	$(FC) -o $(TARGET) $(OBJECTS_TEST)

#per pulire tutta la cartella di lavoro dai file inutili
clean:
	rm *.mod *.o $(TARGET)

#manca specificare l'albero delle dipendenze
main.o: function.o solver.o precision.o
solver.o: function.o precision.o
function.o: precision.o

test.o: precision.o
