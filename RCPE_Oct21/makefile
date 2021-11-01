#
# fortran 
#

#FC 	= gfortran -fdefault-real-8 -O3

TARGET  = RCPEtest

#RCAC_obj 	= RCAC_functions.o \
#                  RCAC_variables.o
#RCACtest	= $(RCAC_obj)

#ALL: $(TARGET)

#RCACtest: ${RCACtest}
#	${FC} -o RCACtest ${RCACtest}
#RCAC_functions.o:

ALL: 
	#module load mkl
	gfortran -c RCPE_variables.F95
	gfortran -c RCPE_functions.F95
	gfortran -o RCPEtest test_Oct2021_RCPE.F95 RCPE_variables.o RCPE_functions.o 

clean:
	${RM} *.o $(TARGET)


#gfortran -o RCACtest test.F95 RCAC_variables.o RCAC_functions.o -lmkl_gf_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm
#gfortran -std=f95 -o RCACtest test.F95 RCAC_variables.o RCAC_functions.o -framework accelerate
#gfortran -o RCACtest test.F95 RCAC_variables.o RCAC_functions.o liblapack.so.3
#gfortran -o RCACtest test_May2018_V1.F95 RCAC_variables.o RCAC_functions.o liblapack.so.3
#gfortran -o RCACtest test_May2018_V1.F95 RCAC_variables.o RCAC_functions.o 