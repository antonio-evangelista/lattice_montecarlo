TARGET = lattice

FC = gfortran

#il compilatore si ferma sempre dopo tre errori se per caso dovesse
#sbroccare e farne vedere dumila
#la lettera o maiusc = ottimizzazione, buon livello è 2, per far andare
#il codice più veloce
FFLAGS = -fmax-errors=3 # -O2

#GCC = gcc
#qui commento con il cancelletto
#unico posto dove devo cambiare le cose è il nome dei sorgenti e nome TARGET
#che è il nome del file
#possiamo pure dirgli, prendi tutti i file .f90 della cartella
SOURCES = $(wildcard *.f90) 

OBJECTS = $(SOURCES:.f90=.o)

.SUFFIXES:
.SUFFIXES: .f90 .o #.c

#scriviamo la regola per passare da .f90 -> .o
.f90.o:
	$(FC) -c $<

#.c.o:
#uguale a sopra ma con gcc


#valuta gfortran -o mach precision.o machine.o
all: $(OBJECTS)
	$(FC) -o $(TARGET) $(OBJECTS)

#per pulire tutta la cartella di lavoro dai file inutili
clean:
	rm *.mod *.o $(TARGET)

#manca specificare l'albero delle dipendenze
main.o: function.o solver.o
solver.o: functions.o
