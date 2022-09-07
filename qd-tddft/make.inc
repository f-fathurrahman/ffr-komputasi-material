MAKE = make
AR = ar

F90 = gfortran
#F90_OPTS = -cpp -Wall -O0 -fbounds-check #-O3
F90_OPTS = -cpp -Wall -O3
F90_OPTS_NOWARN = -cpp -O3

CC = gcc -cpp -I./
CC_OPTS = 

LIB_LINALG = -lblas -llapack
LIB_FFTW3 = -lfftw3
LIBS = $(LIB_LINALG) $(LIB_FFTW3) 

