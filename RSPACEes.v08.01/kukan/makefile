# makefile for ver 8f 04/23/2023
###############################################################################

system = xeon
systm= $(strip $(system))

# [blank]   | [ analyzing hostname (not applicable to all systems) ]
# xeon      | Intel oneAPI
# ohtaka    | System B in ISSP (Intel compiler)
# ohtaka_gcc  | System B in ISSP (GCC compiler)
# ohtaka_aocc | System B in ISSP (AOCC compiler)
# ubuntu    | gfortran for Linux PC (serial)
# fugaku    | mpifrtpx for Fugaku
# octopus   | octopus in CMC
# ito       | mpifrt for ito
# ito-intel | mpiifort for ito


openmp    = yes
sclpck    = yes

# yes / no  (only relevant, if corresponding compiler flag is specified for that system)

###############################################################################

###############################################################################
ifeq ($(systm),xeon)
F90         = mpiifort
OPTION =
ifeq ($(strip $(openmp)),yes)
OPENMP      = -qopenmp
ifeq ($(strip $(sclpck)),yes)
LOPTION     = -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -pthread -lm
else
LOPTION     = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm
endif
else
OPENMP      =
ifeq ($(strip $(sclpck)),yes)
LOPTION     =  -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -pthread -lm
else
LOPTION     =  -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm
endif
endif
endif
###############################################################################
ifeq ($(systm),ohtaka)
F90         = mpiifort
OPTION = -O3 -march=core-avx2
ifeq ($(strip $(openmp)),yes)
OPENMP      = -qopenmp
ifeq ($(strip $(sclpck)),yes)
LOPTION     = -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -pthread -lm
else
LOPTION     = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm
endif
else
OPENMP      =
ifeq ($(strip $(sclpck)),yes)
LOPTION     =  -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -pthread -lm
else
LOPTION     =  -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm
endif
endif
endif
###############################################################################
ifeq ($(systm),ohtaka_gcc)
F90         = mpif90
F90OPT      = -ffree-line-length-none -w -fallow-argument-mismatch -fdollar-ok
FIXOPT      = -ffixed-form
OPTDOUBLE   = -fdefault-real-8
OPTION      =  -O3 -march=znver2
ifeq ($(strip $(openmp)),yes)
OPENMP      = -fopenmp
ifeq ($(strip $(sclpck)),yes)
LOPTION     = -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lpthread -lm
else
LOPTION     = -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lpthread -lm
endif
else
OPENMP      =
ifeq ($(strip $(sclpck)),yes)
LOPTION     = -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread -lm
else
LOPTION     = -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread -lm
endif
endif
endif
###############################################################################
ifeq ($(systm),ohtaka_aocc)
F90         = mpif90
OPTION = -I ${AOCL_ROOT}/include -O3 -march=znver2
ifeq ($(strip $(openmp)),yes)
OPENMP      = -fopenmp
ifeq ($(strip $(sclpck)),yes)
LOPTION     = -lblis-mt -lflame -lscalapack
else
LOPTION     = -lblis-mt -lflame
endif
else
OPENMP      =
ifeq ($(strip $(sclpck)),yes)
LOPTION     =  -lblis -lflame -lscalapack
else
LOPTION     =  -lblis -lflame
endif
endif
endif
###############################################################################
ifeq ($(systm),fugaku)
F90         = mpifrtpx
OPTION      = -Kfast -Kocl -Kilfunc -Kpreex
ifeq ($(strip $(openmp)),yes)
OPENMP      = -Kopenmp
ifeq ($(strip $(sclpck)),yes)
LOPTION     = -SSL2BLAMP -SCALAPACK
else
LOPTION     = -SSL2BLAMP
endif
else
OPENMP      =
ifeq ($(strip $(sclpck)),yes)
LOPTION     = -SSL2 -SCALAPACK
else
LOPTION     = -SSL2
endif
endif
endif
###############################################################################
ifeq ($(systm),ubuntu)
F90         = mpif90
F90OPT      = -ffree-line-length-none
FIXOPT      = -ffixed-form
OPTDOUBLE   = -fdefault-real-8
OPTION      = 
ifeq ($(strip $(openmp)),yes)
OPENMP      = -fopenmp
ifeq ($(strip $(sclpck)),yes)
LOPTION     = -llapack -lblas -lscalapack-openmpi -lblacs-openmpi
else
LOPTION     = -llapack -lblas
endif
else
OPENMP      =
ifeq ($(strip $(sclpck)),yes)
LOPTION     = -llapack -lblas -lscalapack-openmpi -lblacs-openmpi
else
LOPTION     = -llapack -lblas
endif
endif
endif
###############################################################################
ifeq ($(systm),octopus)
F90         = mpiifort -O3
OPTDOUBLE   = -r8
OPTION      = -I$(MKLROOT)/include
ifeq ($(strip $(openmp)),yes)
OPENMP      = -qopenmp
ifeq ($(strip $(sclpck)),yes)
LOPTION     = -I${MKLROOT}/include ${MKLROOT}/lib/intel64/libmkl_scalapack_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_lp64.a -Wl,--end-group -liomp5 -lpthread -lm -ldl
else
LOPTION     = -I${MKLROOT}/include ${MKLROOT}/lib/intel64/libmkl_scalapack_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_lp64.a -Wl,--end-group -liomp5 -lpthread -lm -ldl
endif
else
OPENMP      =
ifeq ($(strip $(sclpck)),yes)
LOPTION     = -I${MKLROOT}/include ${MKLROOT}/lib/intel64/libmkl_scalapack_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_lp64.a -Wl,--end-group -liomp5 -lpthread -lm -ldl
else
LOPTION     = -I${MKLROOT}/include ${MKLROOT}/lib/intel64/libmkl_scalapack_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_lp64.a -Wl,--end-group -liomp5 -lpthread -lm -ldl
endif
endif
endif
###############################################################################
ifeq ($(systm),ito)
F90         = mpifrt
OPTION      = -Kfast
ifeq ($(strip $(openmp)),yes)
OPENMP      = -Kopenmp
ifeq ($(strip $(sclpck)),yes)
LOPTION     = -SSL2BLAMP -SCALAPACK
else
LOPTION     = -SSL2BLAMP
endif
else
OPENMP      =
ifeq ($(strip $(sclpck)),yes)
LOPTION     = -SSL2 -SCALAPACK
else
LOPTION     = -SSL2
endif
endif
endif
###############################################################################
ifeq ($(systm),ito-intel)
F90         = mpiifort
OPTION      = -I${MKLROOT}/include
ifeq ($(strip $(openmp)),yes)
OPENMP      = -qopenmp
ifeq ($(strip $(sclpck)),yes)
LOPTION     = ${MKLROOT}/lib/intel64/libmkl_scalapack_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_lp64.a -Wl,--end-group -liomp5 -lpthread -lm -ldl
else
LOPTION     = -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl
endif
else
OPENMP      =
ifeq ($(strip $(sclpck)),yes)
LOPTION     = ${MKLROOT}/lib/intel64/libmkl_scalapack_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_lp64.a -Wl,--end-group -lpthread -lm -ldl
else
LOPTION     = -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl
endif
endif
endif
###############################################################################

###############################################################################
OBJECT-00= 
TEST-01=
OBJECT-01= keys8f.o mathfunctions8f.o modules8f.o radialhartree8f.o scf_diffcharge8f.o var_global8f.o 
TEST-02=
OBJECT-02= dim_global8f.o tools8f.o fuzzycell8f.o interpolation8f.o output8f.o overlap_finitedifference8f.o overlap_interpolation8f.o vxpot8f.o ggaxyz8f.o ggartp8f.o \
 ksprecondition8f.o kmesh8f.o
TEST-03=
OBJECT-03= kslaplacian8f.o read_ppfile8f.o
TEST-04=
OBJECT-04= read_input_kukan8f.o
TEST-05=
OBJECT-05= pseudocalc8f.o
TEST-06=
OBJECT-06= nonlocaloperation8f.o orthogonalization8f.o
TEST-07=
OBJECT-07= kscg8f.o ksdiis8f.o reorderstates8f.o
OBJECT-07sclp= rayleigh-ritz8f_sclpck.o
OBJECT-07lp= rayleigh-ritz8f.o
TEST-08=
OBJECT-08=
TEST-09=
OBJECT-09=
TEST-10=
OBJECT-10= scf_augcharge8f.o scf_fuzzycellmoment8f.o scf_hartree8f.o scf_onecentervh8f.o \
 scf_rhoxcalc8f.o scf_radialmoment8f.o scf_rhoaugdense8f.o scf_vhbound8f.o scf_vxcalc8f.o \
 trans8f.o
TEST-11=
OBJECT-11= filterpseudopotentials8f.o force8f.o rw_wfpot8f.o scf_charge8f.o \
 scf_chargemixing8f.o scf_diag8f.o scf_dij8f.o scf_occupation8f.o scf_polcon8f.o scf_potentials8f.o \
 setup_initialwave8f.o setup_pseudopot8f.o totalenergy8f.o potout8f.o atomicocc_pdos8f.o
TEST-12=
OBJECT-12= 
TEST-13=
OBJECT-13=
TEST-14=
OBJECT-14= kukan8f.o
TEST-15=
OBJECT-15= ele8f.o
TEST-16=
OBJECT-16= ldos8f.o
TEST-16=
OBJECT-17= pdos8f.o
TEST-17=
OBJECT-18= dos8f.o
TEST-18=
OBJECT-19= orbcharge8f.o
TEST-19=

ifeq ($(strip $(sclpck)),yes)
OBJECT=   $(OBJECT-00) $(TEST-01) $(OBJECT-01) $(TEST-02) $(OBJECT-02) $(TEST-03) $(OBJECT-03) \
 $(TEST-04) $(OBJECT-04) $(TEST-05) $(OBJECT-05) $(TEST-06) $(OBJECT-06) $(TEST-07) $(OBJECT-07) $(OBJECT-07sclp) \
 $(TEST-08) $(OBJECT-08) $(TEST-09) $(OBJECT-09) $(TEST-10) $(OBJECT-10) $(TEST-11) $(OBJECT-11) \
 $(TEST-12) $(OBJECT-12) $(TEST-13) $(OBJECT-13) $(TEST-14) $(OBJECT-14)
OBJECTele=  $(OBJECT-00) $(TEST-01) $(OBJECT-01) $(TEST-02) $(OBJECT-02) $(TEST-03) $(OBJECT-03) \
 $(TEST-04) $(OBJECT-04) $(TEST-05) $(OBJECT-05) $(TEST-06) $(OBJECT-06) $(TEST-07) $(OBJECT-07) $(OBJECT-07sclp) \
 $(TEST-08) $(OBJECT-08) $(TEST-09) $(OBJECT-09) $(TEST-10) $(OBJECT-10) $(TEST-11) $(OBJECT-11) \
 $(TEST-12) $(OBJECT-12) $(TEST-13) $(OBJECT-13) $(TEST-15) $(OBJECT-15)
OBJECTldos= $(OBJECT-00) $(TEST-01) $(OBJECT-01) $(TEST-02) $(OBJECT-02) $(TEST-03) $(OBJECT-03) \
 $(TEST-04) $(OBJECT-04) $(TEST-05) $(OBJECT-05) $(TEST-06) $(OBJECT-06) $(TEST-07) $(OBJECT-07) $(OBJECT-07sclp) \
 $(TEST-08) $(OBJECT-08) $(TEST-09) $(OBJECT-09) $(TEST-10) $(OBJECT-10) $(TEST-11) $(OBJECT-11) \
 $(TEST-12) $(OBJECT-12) $(TEST-13) $(OBJECT-13) $(TEST-16) $(OBJECT-16)
OBJECTpdos= $(OBJECT-00) $(TEST-01) $(OBJECT-01) $(TEST-02) $(OBJECT-02) $(TEST-03) $(OBJECT-03) \
 $(TEST-04) $(OBJECT-04) $(TEST-05) $(OBJECT-05) $(TEST-06) $(OBJECT-06) $(TEST-07) $(OBJECT-07) $(OBJECT-07sclp) \
 $(TEST-08) $(OBJECT-08) $(TEST-09) $(OBJECT-09) $(TEST-10) $(OBJECT-10) $(TEST-11) $(OBJECT-11) \
 $(TEST-12) $(OBJECT-12) $(TEST-13) $(OBJECT-13) $(TEST-17) $(OBJECT-17)
OBJECTdos=  $(OBJECT-00) $(TEST-01) $(OBJECT-01) $(TEST-02) $(OBJECT-02) $(TEST-03) $(OBJECT-03) \
 $(TEST-04) $(OBJECT-04) $(TEST-05) $(OBJECT-05) $(TEST-06) $(OBJECT-06) $(TEST-07) $(OBJECT-07) $(OBJECT-07sclp) \
 $(TEST-08) $(OBJECT-08) $(TEST-09) $(OBJECT-09) $(TEST-10) $(OBJECT-10) $(TEST-11) $(OBJECT-11) \
 $(TEST-12) $(OBJECT-12) $(TEST-13) $(OBJECT-13) $(TEST-18) $(OBJECT-18)
OBJECTorb=  $(OBJECT-00) $(TEST-01) $(OBJECT-01) $(TEST-02) $(OBJECT-02) $(TEST-03) $(OBJECT-03) \
 $(TEST-04) $(OBJECT-04) $(TEST-05) $(OBJECT-05) $(TEST-06) $(OBJECT-06) $(TEST-07) $(OBJECT-07) $(OBJECT-07sclp) \
 $(TEST-08) $(OBJECT-08) $(TEST-09) $(OBJECT-09) $(TEST-10) $(OBJECT-10) $(TEST-11) $(OBJECT-11) \
 $(TEST-12) $(OBJECT-12) $(TEST-13) $(OBJECT-13) $(TEST-19) $(OBJECT-19)
else
OBJECT=   $(OBJECT-00) $(TEST-01) $(OBJECT-01) $(TEST-02) $(OBJECT-02) $(TEST-03) $(OBJECT-03) \
 $(TEST-04) $(OBJECT-04) $(TEST-05) $(OBJECT-05) $(TEST-06) $(OBJECT-06) $(TEST-07) $(OBJECT-07) $(OBJECT-07lp) \
 $(TEST-08) $(OBJECT-08) $(TEST-09) $(OBJECT-09) $(TEST-10) $(OBJECT-10) $(TEST-11) $(OBJECT-11) \
 $(TEST-12) $(OBJECT-12) $(TEST-13) $(OBJECT-13) $(TEST-14) $(OBJECT-14)
OBJECTele=  $(OBJECT-00) $(TEST-01) $(OBJECT-01) $(TEST-02) $(OBJECT-02) $(TEST-03) $(OBJECT-03) \
 $(TEST-04) $(OBJECT-04) $(TEST-05) $(OBJECT-05) $(TEST-06) $(OBJECT-06) $(TEST-07) $(OBJECT-07) $(OBJECT-07lp) \
 $(TEST-08) $(OBJECT-08) $(TEST-09) $(OBJECT-09) $(TEST-10) $(OBJECT-10) $(TEST-11) $(OBJECT-11) \
 $(TEST-12) $(OBJECT-12) $(TEST-13) $(OBJECT-13) $(TEST-15) $(OBJECT-15)
OBJECTldos= $(OBJECT-00) $(TEST-01) $(OBJECT-01) $(TEST-02) $(OBJECT-02) $(TEST-03) $(OBJECT-03) \
 $(TEST-04) $(OBJECT-04) $(TEST-05) $(OBJECT-05) $(TEST-06) $(OBJECT-06) $(TEST-07) $(OBJECT-07) $(OBJECT-07lp) \
 $(TEST-08) $(OBJECT-08) $(TEST-09) $(OBJECT-09) $(TEST-10) $(OBJECT-10) $(TEST-11) $(OBJECT-11) \
 $(TEST-12) $(OBJECT-12) $(TEST-13) $(OBJECT-13) $(TEST-16) $(OBJECT-16)
OBJECTpdos= $(OBJECT-00) $(TEST-01) $(OBJECT-01) $(TEST-02) $(OBJECT-02) $(TEST-03) $(OBJECT-03) \
 $(TEST-04) $(OBJECT-04) $(TEST-05) $(OBJECT-05) $(TEST-06) $(OBJECT-06) $(TEST-07) $(OBJECT-07) $(OBJECT-07lp) \
 $(TEST-08) $(OBJECT-08) $(TEST-09) $(OBJECT-09) $(TEST-10) $(OBJECT-10) $(TEST-11) $(OBJECT-11) \
 $(TEST-12) $(OBJECT-12) $(TEST-13) $(OBJECT-13) $(TEST-17) $(OBJECT-17)
OBJECTdos=  $(OBJECT-00) $(TEST-01) $(OBJECT-01) $(TEST-02) $(OBJECT-02) $(TEST-03) $(OBJECT-03) \
 $(TEST-04) $(OBJECT-04) $(TEST-05) $(OBJECT-05) $(TEST-06) $(OBJECT-06) $(TEST-07) $(OBJECT-07) $(OBJECT-07lp) \
 $(TEST-08) $(OBJECT-08) $(TEST-09) $(OBJECT-09) $(TEST-10) $(OBJECT-10) $(TEST-11) $(OBJECT-11) \
 $(TEST-12) $(OBJECT-12) $(TEST-13) $(OBJECT-13) $(TEST-18) $(OBJECT-18)
OBJECTorb=  $(OBJECT-00) $(TEST-01) $(OBJECT-01) $(TEST-02) $(OBJECT-02) $(TEST-03) $(OBJECT-03) \
 $(TEST-04) $(OBJECT-04) $(TEST-05) $(OBJECT-05) $(TEST-06) $(OBJECT-06) $(TEST-07) $(OBJECT-07) $(OBJECT-07lp) \
 $(TEST-08) $(OBJECT-08) $(TEST-09) $(OBJECT-09) $(TEST-10) $(OBJECT-10) $(TEST-11) $(OBJECT-11) \
 $(TEST-12) $(OBJECT-12) $(TEST-13) $(OBJECT-13) $(TEST-19) $(OBJECT-19)
endif
###############################################################################

###############################################################################

%.o: %.F
	$(F90) $(OPENMP) $(FIXOPT) $(OPTION) -c $*.F

%.o: %.F90
	$(F90) $(OPENMP) $(F90OPT) $(OPTION) -c $*.F90

all: $(OBJECT) $(OBJECTele) $(OBJECTldos) $(OBJECTpdos) $(OBJECTdos) $(OBJECTorb)
	$(F90) $(OPENMP) $(OPTION) -o kukan8 $(OBJECT)     $(LOPTION)
	$(F90) $(OPENMP) $(OPTION) -o ele    $(OBJECTele)  $(LOPTION)
	$(F90) $(OPENMP) $(OPTION) -o ldos   $(OBJECTldos) $(LOPTION)
	$(F90) $(OPENMP) $(OPTION) -o pdos   $(OBJECTpdos) $(LOPTION)
	$(F90) $(OPENMP) $(OPTION) -o dos   $(OBJECTdos) $(LOPTION)
	$(F90) $(OPENMP) $(OPTION) -o orbcharge   $(OBJECTorb) $(LOPTION)

clean:
	rm *~ kukan8 ele ldos pdos dos orbcharge *.o *.mod

###############################################################################
