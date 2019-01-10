INST_DIR = ./bin
EXECUTABLE  = Seed2Cor

#DBG = -g
fflags = -Wall
FFLAGS = $(DBG) $(fflags)
cflags = -O2 -Wall ${DBG} -I./include
cflag2 = -Wall ${DBG} -I./include

LDLIBS = -L/usr/lib -lfftw3 -lstdc++ -lpthread -fopenmp

CC = g++

all : $(EXECUTABLE)

FOBJS = driver.o SacOps.o SysTools.o MakeSDB.o Pthread.o ProduceSac.o RemoveResp.o TempSpecNorm.o Whiten.o CrossCorr.o FreqOps.o

$(EXECUTABLE) : $(FOBJS)
	$(CC) $^ -o $@ $(LDLIBS) $(cflags)

Whiten.o : Whiten.c
	$(CC) $(cflag2) -c $<

FreqOps.o : FreqOps.c
	$(CC) $(cflag2) -c $<

%.o : %.c
	$(CC) $(cflags) -c $<

install : $(EXECUTABLE)
	install -s $(EXECUTABLE) $(INST_DIR)

install_ariadne : $(EXECUTABLE)
	cp $(EXECUTABLE) $(INST_DIR)/$(EXECUTABLE)_ariadne

clean :
	rm -f $(EXECUTABLE) core $(FOBJS)
