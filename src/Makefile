#
# Make file
#
#******************************************************
# BASH code
#******************************************************
define folder_set =
#
k=0
for i in $(ls ../); do
	if [ $i = "bin" ]; then
		k=1
	fi
done

if [ $k = 0 ]; then
	mkdir ../bin
	echo ../bin folder created
else
	echo ../bin folder exists
fi
#
endef
#
folder_set: ; @$(value folder_set)

.ONESHELL:

#******************************************************
bin_path=../bin

QP_pgr_src = QP_module.f90 QP_program.f90
QP_pgr_obj = QP_module.o QP_program.o

.SUFFIXES: 
.SUFFIXES: .o .f90 .oOMP

LIBS = -L/usr/local/lib/lapack95 -llapack95 -L/usr/local/lib -llapack95 -llapack
OPT  = -cpp -I/usr/local/lib/lapack95/lapack95_modules -I/usr/local/include/lapack95_modules
.f90.o:
	gfortran -c $(OPT)  "$<"

#
all: folder_set QP_pgr
#
.phony : all
#
#
QP_pgr: $(QP_pgr_obj) Makefile
	$(info)
	gfortran -o $(bin_path)/$@.out $(QP_pgr_obj) $(OPT) $(LIBS)

clean:
	@(rm -f *.o *.mod)
	@rm -f ../*out

