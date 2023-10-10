#... Compiler F90
DEBUG0 = 
F_C = ifort -c -openmp
#... Compiler F77
F_77 = ifort -c -openmp
#... Compiler flags
F_F = -convert big_endian -align all 
#... Compiler flags F77
F_F77 = $(F_F) 
#... Pre-processor options
P_P = -Dopt_netcdf -DNOMPI -DLINUX -DDEEPLY_FGAT -DSHARED_MEMORY
# -Wp,-DUSE_POINTERS
#... Optimization flags
F_O = -O3 -xSSSE3 -fno-alias -ip
#... Debug flags
F_D = -Cdebug
#... C Compiler
C_C = gcc -c
#... C Flags
C_F = -DLINUX -DSTATIC_LINKING -DXPRIVATE=PRIVATE \
  -UINTERCEPT_ALLOC -DINTEGER_IS_INT -D_ABI64 -UHAS_XMOTIF \
#... Linker
F_L = ifort -openmp -threads
#... Linking options
L_O = 
#... Debug flag (Debug on: 1 ; Debug off: whatelse)
DEBUG = 0
#... AR add
AR_I = ar -rv
#... AR delete
AR_E = ar -d
#.. Remove
RM   = rm -f

#... Include dirs
EXTINC = -I/S/home00/G5032/c0341/local/netcdf-4.1.3-ifort/include \
         -I../include
#... Libraries
DRHOOK_LIB = /S/home00/G5032/c0341/ODA/aux/drhook_CY33R1.007.ifort
ARPACK = /S/home00/G5032/c0341/ODA/aux/arpack/ARPACK
EXTLIB = -L/S/home00/G5032/c0341/local/netcdf-4.1.3-ifort/lib -lnetcdf -lnetcdff \
         -L$(DRHOOK_LIB) -ldrhook -lmpi_serial \
         -L$(ARPACK) -larpack_IFORT
#
