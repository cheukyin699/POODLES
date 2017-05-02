CC  = g++
OBJ = BioTables.o Files.o PoodleS.o server.o SVMsmplStatus.o

POODLES-EX: $(OBJ)            
	$(CC) -o $@ $(OBJ) -lm
	cp POODLES-EX ../

BioTables.o: BioTables.cpp stdafx.h
	$(CC) -c BioTables.cpp 

Files.o: Files.cpp stdafx.h
	$(CC) -c Files.cpp 

PoodleS.o: PoodleS.cpp stdafx.h
	$(CC) -c PoodleS.cpp 

server.o: server.cpp stdafx.h
	$(CC) -c server.cpp 

SVMsmplStatus.o: SVMsmplStatus.cpp stdafx.h
	$(CC) -c SVMsmplStatus.cpp 

clean:
	rm *.o *~ ../POODLES-EX
