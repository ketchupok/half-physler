#
#       makefile for the csound opcode plugins of IWK
#				by Sebastian Schmutzhard & Alex Hofmann

CPP    = g++
UNAME := $(shell uname)

	OBJS    = tube.o const.o #single_reed.o
	FILE    = src/opcode_registers.cpp
	OUTFILE = opcodes_iwk.so

ifeq ($(UNAME), Linux)
	DYL 	= -shared
	FLAGS   = -fPIC -std=c++11
	LDLIBS  = -I /usr/local/include/csound
else ifeq ($(UNAME), Darwin)
	DYL	= -dynamiclib
	FLAGS   = -std=c++11
	LDLIBS  = -I ~/Library/Frameworks/CsoundLib64.framework/Versions/6.0/Headers/
endif


#       targets
default_target:
	make $(OUTFILE)

$(OUTFILE): $(OBJS) $(FILE) ./resonators/* makefile
	$(CPP) $(FLAGS) $(OBJS) $(FILE) $(LDLIBS) $(DYL) -o $(OUTFILE)

tube.o :: ./src/tube.cpp ./src/tube.h const.o
	$(CPP) $(FLAGS) ./src/tube.cpp $(LDLIBS) -c

const.o :: ./src/const.cpp ./src/const.h
	$(CPP) $(FLAGS) ./src/const.cpp $(LDLIBS) -c

run ::
	make
	csound Demo_PC_Mac/demo_halfphysler.csd

runBela:
	make
	make installBela
	belacsound --csd=/Bela/projects/Bela_HalfPhysler_DemoProject/_main.csd --period=32

installBela:
	cp $(OUTFILE) ./Bela/projects/Bela_HalfPhysler_DemoProject
  # TODO copy this project folder into Belas structure automatically..

lcean:
	make clean

clean: .
	@ rm -vf *~  *.o  *.so

clena: .
	make clean
#.SECONDARY:
#.SILENT:
