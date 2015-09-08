#########################################
#										#
#										#
# 	Makefile for WLTest					#
#	Language: C++                  		#
#                                      	#
#	(c) 2012 Jerry Shi             		#
#										#
#########################################

### Suffixes

.SUFFIXES:
.SUFFIXES: .cpp .o


### Macro ###

#-- Compiler 
CC = mpicxx

#-- Compiler flags
CFLAGS = -O3 

#-- Directory of include-files
MYINCDIR  = ./include
SYSINCDIR = /usr/local/gsl/1.15/include

#-- Directory of libraries
LIBDIR = /usr/local/gsl/1.15/lib

#-- Libraries
LIBS = -lgsl -lgslcblas -lm 

#-- source directory
SRCDIR = ./src/

VPATH=  .:$(MYINCDIR):$(SRCDIR):$(LIBDIR):$(SYSINCDIR)

#-- Object files
OBJS = 	Main.o 						\
		LipidModel.o				\
		Random.o					\
		Global.o					\
		Monomer.o					\
		HistogramND.o				\
		RandomAccessNeighborList.o
	
#-- Target
TARGET = PWLSIM 


#############

### Inference rule ###

%.o: %.cpp 
	${CC} ${CFLAGS} -I${MYINCDIR} -I${SYSINCDIR} -L${LIBDIR} ${LIBS} $< -c 

######################

### Target ###

${TARGET}: ${OBJS}
	${CC} ${CFLAGS} -I${MYINCDIR} -I${SYSINCDIR} $^ -o $@ -L${LIBDIR} ${LIBS} 

##############

### Indirect dependencies ###

Main.o:  					Main.cpp 			WLFrame.h 	LipidModel.o
Global.o: 					Global.cpp 			Global.h
Random.o: 					Random.cpp 			Random.h
LipidModel.o: 				LipidModel.cpp 		LipidModel.h Random.o MathVector.h RandomAccessNeighborList.h
HistogramND.o:				HistogramND.cpp 	HistogramND.h
Monomer.o:					Monomer.cpp			Monomer.h
RandomAccessNeighborList.o:	RandomAccessNeighborList.cpp RandomAccessNeighborList.h

#############################

### generate documentation ###
docs:
	doxygen ./doc/configure 
web:
	xdg-open ./doc/html/index.html

### test ###

run1D:
	mpirun -n 4 ${TARGET} wang_landau.input 1 4 0.75 
run2D: 
	mpirun -n 20 ${TARGET} wang_landau.input 2 10 2 0.75 0.75 

### Clean ###
clean:
	rm -f *.o $(TARGET) ./include/*.gch 
	
