#-------------------------------------------
# Makefile for any C program with ROOT lib
# Zhihong Ye, 01/09/2009
#-------------------------------------------

PROGRAM = getsidis
CFLAGS = $(shell $(ROOTSYS)/bin/root-config --cflags)
LIBS   = $(shell $(ROOTSYS)/bin/root-config --libs)
GLIBS  = $(shell $(ROOTSYS)/bin/root-config --glibs)

CFLAGS += -Wall -g
INCLUDE := -I${ROOTSYS}/include
RTLIBS := -L${ROOTSYS}/lib -lCore /usr/lib64/libg2c.so.0 # /usr/lib/libg2c.so.0 #use /usr/lib64/ for ifarm
#RTLIBS := -L${ROOTSYS}/lib -lCore  /usr/lib/libg2c.so.0 #use /usr/lib64/ for ifarm
all: $(PROGRAM)

#FFLAGS = -C -g -m32 -fbounds-check -finit-local-zero -fno-automatic 
FFLAGS = -C -g -m64 -fbounds-check -finit-local-zero -fno-automatic # use -m64 on ifarm, which is a 64bit system
GFORTRAN_OBJS := ./fortran_libs/bkk.o ./fortran_libs/Cteq61Pdf.o ./fortran_libs/fbar.o ./fortran_libs/fd.o ./fortran_libs/fmfac.o \
         ./fortran_libs/fm.o ./fortran_libs/fphenom.o ./fortran_libs/fyukawa.o ./fortran_libs/LWW90_F2.o ./fortran_libs/LWW90_R.o \
		 ./fortran_libs/nform.o ./fortran_libs/qfsrad.o ./fortran_libs/radiate.o ./fortran_libs/recoil.o ./fortran_libs/rom.o \
         ./fortran_libs/SIDIS.o ./fortran_libs/sig2n.o ./fortran_libs/sigbar.o ./fortran_libs/sigdel.o ./fortran_libs/sigmot.o\
         ./fortran_libs/sigqfs.o ./fortran_libs/sigr1.o ./fortran_libs/sigr2.o ./fortran_libs/sigx.o ./fortran_libs/spence.o\
         ./fortran_libs/terp.o ./fortran_libs/valy.o ./fortran_libs/whitlow.o ./fortran_libs/wiser.o ./fortran_libs/xsecft.o

$(PROGRAM): $(PROGRAM).o dxs_dis.o ${GFORTRAN_OBJS} 
	g++ -o $(PROGRAM) $(PROGRAM).o dxs_dis.o ${GFORTRAN_OBJS} $(CFLAGS) $(LIBS) $(RTLIBS)

#$(PROGRAM).o: $(PROGRAM).C dxs_dis.o ${GFORTRAN_OBJS} 
#	g++ $(CFLAGS) -c -o $(PROGRAM).o dxs_dis.o ${GFORTRAN_OBJS} $(PROGRAM).C ${INCLUDE}

$dxs_dis.o: dxs_dis.f ${GFORTRAN_OBJS} 
	rm -f $@
	f95 ${GFORTRAN_OBJS} $(FFLAGS) -c dxs_dis.f;

clean:
	rm -f *.o

