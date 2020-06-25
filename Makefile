#
#INST_DIR = $(HOME)/bin
BIN  = spectral_snr_s2c

fflags = -O -ffixed-line-length-none
cflags = -O

LDLIBS = -lfftw3

FFLAGS = $(DBG) $(fflags)
CFLAGS = $(DBG) $(cflags)

FC = gfortran
CC = gcc

DBG =
FOBJS = spectral_snr_s2c.c filter4_cv.o swapn.o

$(BIN) : $(FOBJS)
	$(FC) $(FFLAGS) $(CFLAGS) $(FOBJS) -o $(BIN) $(LDLIBS)

install :: $(BIN)
	install -s $(BIN) $(INST_DIR)

clean ::
	rm -f $(BIN) core $(FOBJS)
