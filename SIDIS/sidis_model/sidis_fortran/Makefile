#-------------------------------------------
# Makefile for any C program with ROOT lib
# Zhihong Ye, 01/09/2009
#-------------------------------------------

PROGRAM = GetSIDIS
#PROGRAM = getsidis
CFLAGS = $(shell $(ROOTSYS)/bin/root-config --cflags)
LIBS   = $(shell $(ROOTSYS)/bin/root-config --libs)
GLIBS  = $(shell $(ROOTSYS)/bin/root-config --glibs)

CFLAGS += -Wall -g -ggdb -O2
INCLUDE := -I${ROOTSYS}/include
#RTLIBS := -L${ROOTSYS}/lib -lCore -lm /usr/lib64/libgfortran.so.3.0.0 
#RTLIBS := -L${ROOTSYS}/lib -lCore -lm /usr/lib/libg2c.so.0
RTLIBS := -L${ROOTSYS}/lib -lCore -lm /home/yez/mylib/libg2c.so.0
#RTLIBS := -L${ROOTSYS}/lib -lCore -lm /usr/lib/libgfortran.so.1

F77  = g77
disk = /apps/cernlib/i386_rhel3/2004
FFLAGS = -u -C -O

#disk = /apps/cernlib/x86_64_rhel6/2005
#FFLAGS = -u -C -O -g -m64 -fbounds-check -finit-local-zero -fno-automatic # use -m64 on ifarm, which is a 64bit system

GFORTRAN_OBJS = whitlow.o LWW90_F2.o LWW90_R.o xsecft.o wiser.o SIDIS.o bkk.o Cteq61Pdf.o qfsrad.o radiate.o rom.o valy.o sig2n.o sigbar.o sigdel.o sigmot.o sigqfs.o sigr1.o sigr2.o sigx.o fbar.o fd.o fm.o fmfac.o fphenom.o fyukawa.o recoil.o spence.o terp.o nform.o
Clibs  = -lm -ldl
Xlibs  = -L/usr/X11R6/lib -lX11
CERN   = -L/$(disk)/lib -lmathlib -lkernlib -lpacklib 

all: $(PROGRAM)

dxs_dis.o: dxs_dis.f
	$(F77) $(fFlags) -c $< -o $@

$(PROGRAM).o:	$(PROGRAM).C 
	g++ $(CFLAGS) -c -o $(PROGRAM).o $(PROGRAM).C $(INCLUDE) 

$(PROGRAM):	$(PROGRAM).o ./dxs_dis.o $(GFORTRAN_OBJS)
	g++ $(CFLAGS) -o $(PROGRAM) $(PROGRAM).o ./dxs_dis.o $(GFORTRAN_OBJS) $(LIBS) $(Clibs) $(Xlibs) $(RTLIBS) $(CERN)

clean:
	rm -f dxs_dis.o $(PROGRAM).o

