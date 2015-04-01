# \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\  oct'98
# \\\          ---------------------------------------------------       \\\
# \\\          e907 / e931 simulation                                    \\\
# \\\                                                   geant 3.21       \\\
# \\\          uses specially prepared geant library                     \\\
# \\\          cernlib at:   /disk/cernlib/98a/lib                       \\\
# \\\                                                                    \\\
# \\\          support for   Linux - Digital Unix                        \\\
# \\\          ---------------------------------------------------       \\\
# \\\                              ae  -->  empl@homer.phys.uh.edu       \\\


#-------------- Linux
# f77  = f77
# disk = users
#-------------- Linux
  f77  = g77
  disk = /apps/cernlib/i386_rhel3/2004


  cFlags  = -g
  fFlags  = -u -C -O

  Clibs  = -lm -ldl
  Xlibs  = -L/usr/X11R6/lib -lX11
  CERN   = -L/$(disk)/lib -lmathlib -lkernlib -lpacklib 

  both   = $(libG) $(libg) $(libG) $(libg)
  OBJ = NT20.o getkinp.o qgetkin.o hitraxy.o whitlow.o LWW90_F2.o LWW90_R.o xsecft.o wiser.o SIDIS.o bkk.o Cteq61Pdf.o qfsrad.o radiate.o rom.o valy.o sig2n.o sigbar.o sigdel.o sigmot.o sigqfs.o sigr1.o sigr2.o sigx.o fbar.o fd.o fm.o fmfac.o fphenom.o fyukawa.o recoil.o spence.o terp.o nform.o
  source =  main.f 
%.o: %.f
	$(F77) $(fFlags) -c $< -o $@
all:  main.o $(OBJ)
	$(f77) $(fFlags) -o main main.o $(OBJ) $(CERN)
clean: 
	rm -f *.o *~



