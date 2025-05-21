#-------------------------------------------------------------------------------
# compile the UMFPACK demos (for GNU make and original make)
#-------------------------------------------------------------------------------

# Modified from Timothy Davis's Makefile in VISCO2.5D/SuiteSparse/UMFPACK/Demo

default: libs

include ../SuiteSparse/SuiteSparse_config/SuiteSparse_config.mk

#-----------

#2pt5d = /Users/fpollitz/fred4/CIDER2015/VISCO2.5D/SOURCE
2pt5d = ../SOURCE

#-------------------------------------------------------------------------------
# the optional Partition module requires METIS, CAMD, and CCOLAMD
I_WITH_PARTITION = 
LIB_WITH_PARTITION =
CONFIG1 = 
CONFIG2 = -DNPARTITION
ifeq (,$(findstring -DNCHOLMOD, $(UMFPACK_CONFIG)))
    # CHOLMOD is requested.  See if it is available
    ifeq (../SuiteSparse/CHOLMOD, $(wildcard ../SuiteSparse/CHOLMOD))
        ifeq (../SuiteSparse/COLAMD, $(wildcard ../SuiteSparse/COLAMD))
            # CHOLMOD and COLAMD are available
            CONFIG1 =
            LIB_WITH_CHOLMOD = ../SuiteSparse/CHOLMOD/Lib/libcholmod.a \
                ../SuiteSparse/COLAMD/Lib/libcolamd.a
            # check if METIS is requested and available
            ifeq (,$(findstring -DNPARTITION, $(CHOLMOD_CONFIG)))
                # METIS is requested.  See if it is available
                ifeq ($(METIS_PATH), $(wildcard $(METIS_PATH)))
                    ifeq (../SuiteSparse/CAMD, $(wildcard ../SuiteSparse/CAMD))
                        ifeq (../SuiteSparse/CCOLAMD, $(wildcard ../SuiteSparse/CCOLAMD))
                            # METIS, CAMD, and CCOLAMD are available
                            LIB_WITH_PARTITION = $(METIS) \
                                ../SuiteSparse/CCOLAMD/Lib/libccolamd.a \
                                ../SuiteSparse/CAMD/Lib/libcamd.a
                            I_WITH_PARTITION = -I$(METIS_PATH)/Lib \
                                -I../SuiteSparse/CCOLAMD/Include -I../SuiteSparse/CAMD/Include
                            CONFIG2 =
                        endif
                    endif
                endif
            endif
        endif
    endif
endif
#-------------------------------------------------------------------------------

C = $(CC) $(CF) $(UMFPACK_CONFIG) $(CONFIG1) $(CONFIG2) \
    -I../SuiteSparse/UMFPACK/Include -I../SuiteSparse/AMD/Include -I../SuiteSparse/SuiteSparse_config

INC = ../SuiteSparse/UMFPACK/Include/umfpack.h ../SuiteSparse/AMD/Include/amd.h ../SuiteSparse/SuiteSparse_config/SuiteSparse_config.h

LIBS = $(LAPACK) $(BLAS) $(XERBLA) $(LIB) $(LIB_WITH_CHOLMOD) $(LIB_WITH_PARTITION) $(CUBLAS_LIB) $(CUDART_LIB)

../SuiteSparse/UMFPACK/Lib/libumfpack.a:
	( cd ../SuiteSparse/UMFPACK/Lib ; $(MAKE) )

../SuiteSparse/AMD/Lib/libamd.a:
	( cd ../SuiteSparse/AMD ; $(MAKE) library )

../SuiteSparse/SuiteSparse_config/libsuitesparseconfig.a:
	( cd ../SuiteSparse/SuiteSparse_config ; $(MAKE) library )

../SuiteSparse/CHOLMOD/Lib/libcholmod.a:
	- ( cd ../SuiteSparse/CHOLMOD && $(MAKE) library )

../SuiteSparse/COLAMD/Lib/libcolamd.a:
	- ( cd ../SuiteSparse/COLAMD && $(MAKE) library )

../SuiteSparse/CCOLAMD/Lib/libccolamd.a:
	- ( cd ../SuiteSparse/CCOLAMD && $(MAKE) library )

../SuiteSparse/CAMD/Lib/libcamd.a:
	- ( cd ../SuiteSparse/CAMD && $(MAKE) library )

${BLAS}:
	- ( cd ../OpenBLAS-0.3.13 && $(MAKE) )

$(METIS):
	( cd $(METIS_PATH) && $(MAKE) )

UMFPACK = ../SuiteSparse/UMFPACK/Lib/libumfpack.a ../SuiteSparse/AMD/Lib/libamd.a \
    ../SuiteSparse/SuiteSparse_config/libsuitesparseconfig.a \
    $(LIB_WITH_CHOLMOD) $(LIB_WITH_PARTITION) 

libs: $(UMFPACK) ${BLAS}

#-------------------------------------------------------------------------------

visco2pt5d: ${2pt5d}/visco2pt5d.f ${2pt5d}/visco2pt5dsubs.f ${2pt5d}/glsubs.f ${2pt5d}/svdcmp8.f 
	$(F77) $(F77FLAGS) -ffixed-form -ffixed-line-length-132 -fno-align-commons -o visco2pt5d ${2pt5d}/visco2pt5d.f \
 ${2pt5d}/visco2pt5dsubs.f ${2pt5d}/glsubs.f ${2pt5d}/svdcmp8.f \
 umf4_f77zwrapper64.o $(UMFPACK) $(LIBS)

visco3d: ${2pt5d}/visco3d.f ${2pt5d}/visco3dsubs.f ${2pt5d}/glsubs3.f ${2pt5d}/svdcmp8.f
	$(F77) $(F77FLAGS) -ffixed-form -ffixed-line-length-132 -fno-align-commons -o visco3d ${2pt5d}/visco3d.f \
 ${2pt5d}/visco3dsubs.f ${2pt5d}/glsubs3.f ${2pt5d}/svdcmp8.f \
 umf4_f77zwrapper64.o $(UMFPACK) $(LIBS)

umf4_f77wrapper.o: umf4_f77wrapper.c $(INC)
	$(C) -c umf4_f77wrapper.c -o umf4_f77wrapper.o

umf4_f77zwrapper.o: umf4_f77zwrapper.c $(INC)
	$(C) -c umf4_f77zwrapper.c -o umf4_f77zwrapper.o

umf4_f77wrapper64.o: umf4_f77wrapper.c $(INC)
	$(C) -DDLONG -c umf4_f77wrapper.c -o umf4_f77wrapper64.o

umf4_f77zwrapper64.o: umf4_f77zwrapper.c $(INC)
	$(C) -DZLONG -c umf4_f77zwrapper.c -o umf4_f77zwrapper64.o

#-------------------------------------------------------------------------------
# Remove all but the files in the original distribution
#-------------------------------------------------------------------------------

CLEAN = ../SuiteSparse/*/Lib/*.o ../SuiteSparse/*/Lib/*.a \
    ../SuiteSparse/SuiteSparse_config/libsuitesparseconfig.a \
   ../OpenBLAS-0.3.13/*.o ../OpenBLAS-0.3.13/*/*.o ../OpenBLAS-0.3.13/*/*/*.o \
   ../OpenBLAS-0.3.13/*/*/*/*.o ../OpenBLAS-0.3.13/*.a \
   ../MAINPROG/umf4_f77zwrapper64.o ../MAINPROG/visco2pt5d ../MAINPROG/visco3d 

clean:
	- $(RM) $(CLEAN)
