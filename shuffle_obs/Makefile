F90 = ifort -convert big_endian -O2 -i4 -r8 -g -traceback -fno-alias -assume byterecl 
L90 = $(F90)
NCDI = ${NETCDF}/include 
#/users/home/opt/netcdf/netcdf-3.6.3
NCDL = ${NETCDF}/lib
BIN = shuffle_obs.x
BIN2 = shuffle_obs_pre.x
EXTRALIBS = 

LIB = abor1.o io_obs.o mod_obsdef.o obs_str.o set_dim.o setfgfiles.o \
check.o mod_allocobs.o mod_allocobs2.o mod_readfg.o readobs103.o set_knd.o setfgflayout.o

LIB2 = abor1.o io_obs.o mod_obsdef.o obs_str.o set_dim.o setfgfiles.o \
check.o mod_allocobs.o mod_allocobs2.o mod_readfg.o set_knd.o \
setfgflayout.o readobs103_pre.o mod_phinterp.o int_par_ins_2.o mod_mynetcdf.o index_sort.o

%.o:	%.F90
	$(F90) -c -o $@ $< -I$(NCDI)

all:	$(BIN) $(BIN2)

$(BIN):	$(LIB)
	$(L90) -o $(BIN) $(LIB) -L$(NCDL) -lnetcdf -lnetcdff $(EXTRALIBS)

$(BIN2):	$(LIB2)
	$(L90) -o $(BIN2) $(LIB2) -L$(NCDL) -lnetcdf -lnetcdff $(EXTRALIBS)

clean:	
	rm -f *.o *.mod $(BIN) $(BIN2)

include dependencies.inc
