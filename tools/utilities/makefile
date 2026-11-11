LIBDIR = .
FC = ifx
FFLAGS1 = 
ifeq ($(FC),ifx)
FFLAGS2 = -extend-source
else ifeq ($(FC),gfortran)
# debug segmentation fault:
#FFLAGS2 = -ffixed-line-length-none -g -fcheck=all -Wall -O0
# debug NAN or Infinity:
#FFLAGS2 = -ffixed-line-length-none -g -fbacktrace -ffpe-trap=zero,overflow,underflow,invalid
# default:
FFLAGS2 = -ffixed-line-length-none
else
$(error Unknown FC '$(FC)'; use ifx or gfortran)
endif
AR = ar
RANLIB = ranlib
#CC = gcc

LIBBASES = libutilities.a

SRCSF1 = readinput.F90 global_vars.F90 global_xsection.F90 global_mom_xpdf.F90 global_cuts.F90 global_dists.F90 global_virt.F90\
         ps_VV_2to4.F90 ps_VV_2to5.F90 KinDistSetup.F90 Do_VecSet.F90 Do_VecSet_tilde.F90 ps_VV_2to5_decayV1.F90 ps_VV_2to5_decayV2.F90\
         ps_WW_2to5_singleT.F90 GetScales.F90 ps_WW_2to6_WWbb.F90 ps_VVV_2to6.F90
# Do_VecSet_FC9.F90
SRCSF2 = xsection.F Integrand.F print_info.F model_sm.F ps_tools.F GetPdfs.F OSmap_VV_2to4.F OSmap_VV_2to5.F Cut_tools.F FillHIS.F\
         RenConst.F softphoton_sub.F FuncDPcoll_had.F FuncDPsing_fin_rad.F FuncDPsing_fin_ind.F GetPolarVec.F softphoton_sub_Z.F\
         softphoton_sub_W.F OSmap_VV_2to5_decayV2.F OSmap_VV_2to5_decayV1.F FuncDPsing_decayZ.F FuncDPsing_decayW.F Hist_init_WZ.F\
         Hist_init_WW.F Hist_init_ZZ.F Hist_fill_WZ.F Hist_fill_WW.F Hist_fill_ZZ.F OSmap_TW_2to5.F OSmap_VV_2to6.F OSmap_VVV_2to6.F\
         Hist_fill_WWW.F Hist_init_WWW.F
#ps_2to4_WZlep.F
INC2 = model.h

OBJS1	= ${SRCSF1:.F90=.o}
OBJS2	= ${SRCSF2:.F=.o}

CFLAGS	= ${COPT} -m64

all:	${LIBBASES}

%.o : %.F90
	$(FC) $(FFLAGS1) -c $<
	
%.o : %.F $(INC2)
	$(FC) $(FFLAGS2) -c $<

${LIBBASES}: ${OBJS1} ${OBJS2}
	${AR} rc ${LIBBASES} $?
	${RANLIB} ${LIBBASES}

tags: ${OBJS1} ${OBJS2}

install: all
	cp ${LIBBASES} ${LIBDIR}

clean:
	rm -f ${LIBBASES} ${OBJS1} ${OBJS2} *.mod
