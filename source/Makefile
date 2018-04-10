#=======================================
# CPD cluster (ifort)
#=======================================
#FC = gfortran
FC = mpif90
#FLAGS   =  -DUSE_MPI -O3 -unroll -I/usr/local/amd64/magny-cours/gcc/fftw-3.3.3/include -I/usr/local/amd64/magny-cours/gcc/mvapich2-1.9-ib/include

FLAGS   = -DUSE_MPI -O3 -unroll -cpp -DPOINTER_SIZE=8 -I/home/evitali/sprng2.0/include -DUSE_MKL -fopenmp -m64 -I${MKLINCLUDE} -I/opt/fftw/include
LIBS=-L/usr/local/lib -lgmp -L/home/evitali/sprng2.0/lib -lsprng -Wl,--no-as-needed -L${MKLPATH} -lmkl_gf_lp64 -lmkl_core -lmkl_gnu_thread -ldl -lpthread -lm -L/opt/fftw/lib -lfftw3

#LIBS = -L/usr/local/acml-5.3.1/gfortran64/lib -lacml -L/usr/local/amd64/magny-cours/gcc/fftw-3.3.3/lib -lfftw3 -L/usr/local/amd64/magny-cours/gcc/mvapich2-1.9-ib/lib -lmpich -lmpichf90
##FLAGS   = -DUSE_MPI -warn all -pg -g -check all -O0 -debug-parameters all -traceback -ftrapuv -g -fpe0

#=======================================
# linux box (gfortran)
#=======================================
#FC = mpif90
#FLAGS   = -O3 
#FLAGS   = -Wall -fcheck=all -O0 -fbacktrace -g -ffpe-trap=zero,overflow

#all : bcslat_SOC shell_SOC

#all : bcslat_SOC_comp_w_MFpaper shell_SOC

#all : bcslat_SOC_w_zeeman shell_SOC

#all : bcslat_SOC_w_zeeman_test shell_SOC

all : bcs_SOC_openBCs shell_SOC

#bcslat_SOC : bcslat_SOC.F90
#	$(FC) $(FLAGS) -o bcslat_SOC.exe mpiwrap.F90 bcslat_SOC.F90 $(LIBS) #-llapack -lfftw3

#bcslat_SOC_comp_w_MFpaper : bcslat_SOC_comp_w_MFpaper.F90
#	$(FC) $(FLAGS) -o bcslat_SOC_comp_w_MFpaper.exe mpiwrap.F90 bcslat_SOC_comp_w_MFpaper.F90 $(LIBS) #-llapack -lfftw3

#bcslat_SOC_w_zeeman : bcslat_SOC_w_zeeman.F90
#	$(FC) $(FLAGS) -o bcslat_SOC_w_zeeman.exe caldet.f90 mpiwrap.F90 bcslat_SOC_w_zeeman.F90 $(LIBS)

#bcslat_SOC_w_zeeman_test : bcslat_SOC_w_zeeman_test.F90
#	$(FC) $(FLAGS) -o bcslat_SOC_w_zeeman_test.exe caldet.f90 mpiwrap.F90 bcslat_SOC_w_zeeman_test.F90 $(LIBS)

bcs_SOC_openBCs : bcs_SOC_openBCs.F90
	$(FC) $(FLAGS) -o bcs_SOC_openBCs.exe caldet.f90 mpiwrap.F90 current_current.F90 bcs_SOC_openBCs.F90 $(LIBS) #-llapack -lfftw3

shell_SOC : shell_SOC.F90
	$(FC) $(FLAGS) -o shell_SOC.exe shell_SOC.F90 $(LIBS)
clean :	
	rm -f *.o *.mod *.exe
