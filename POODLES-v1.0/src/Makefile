CC  = g++
CXXFLAGS = -c -O2 -Wall -Wextra
LDFLAGS = -lm -O2
OBJ = BioTables.o Files.o PoodleS.o server.o SVMsmplStatus.o
EXE = POODLES-EX

all: $(EXE)

$(EXE): $(OBJ)            
	$(CC) -o $@ $(OBJ) $(LDFLAGS)
	cp $(EXE) -t ..

%.o: %.cpp stdafx.h
	$(CC) $(CXXFLAGS) $< -o $@

clean:
	rm -f $(OBJ) *~ ../$(EXE)
