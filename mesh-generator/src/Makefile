#************************************************************
#                                                           *
#           Makefile                                        *
#                                                           *
#************************************************************
#*************  makefile def variables                      A
#************************************************************
# Define shell                                             A1
SHELL = /bin/sh
# END SECTION shell
#************************************************************
#************************************************************

CPPFLAGS =  -g 
#LOCDIR   = ../../petscexample/ex1/

include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules


# END SECTION compiler
#************************************************************
#*************  makefile core                               B
#************************************************************
#  Define obj source                                       B1
SOURCES = hex.cpp hexh.cpp wedge.cpp wedgeh.cpp tet.cpp teth.cpp quad.cpp quadh.cpp tri.cpp trih.cpp line.cpp lineh.cpp elem.cpp elem_type.cpp mesh.cpp lsysPDE.cpp multigrid.cpp main.cpp hFun.cpp trihFun.cpp triFun.cpp tetraFun.cpp linFun.cpp linehFun.cpp 

HEAD= main.hpp basis.hpp basish.hpp mesh.hpp elem.hpp elem_type.hpp lsysPDE.hpp multigrid.hpp hFun.hpp hGaussPoints.hpp

OBJS = hex.o hexh.o wedge.o wedgeh.o tet.o teth.o quad.o quadh.o tri.o trih.o line.o lineh.o elem.o mesh.o elem_type.o lsysPDE.o multigrid.o  main.o hFun.o trihFun.o triFun.o tetraFun.o linFun.o linehFun.o 
#************************************************************
# linker rules for OBJS from  OBJS a.out                   B2

zero:  $(OBJS) $(HEAD) chkopts
	$(CLINKER) $(CPPFLAGS)  -o mesh.out  $(OBJS) ${PETSC_KSP_LIB} -L$(LD_LIBRARY_PATH)
	rm -f *~ *.l *.j *.m #*#

