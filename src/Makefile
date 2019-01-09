#
# Make file
#
#******************************************************
# Required libraries and local paths
#******************************************************
FC      =  gfortran
mod_lib = -I/usr/local/lib/lapack95/lapack95_modules
optlib  = -lblas -llapack -llapack95
#
lapack_path   = /usr/local/lib
lapack95_path = $(lapack_path)
#
LIBS = -L$(lapack95_path) $(optlib)
OPT  = -cpp $(mod_lib)
#
#******************************************************
#

#******************************************************
# BASH code
#******************************************************
define folder_set =
#
k=0
q=0
for i in $(ls ../); do
	case $i in
	"bin")
		k=1
	;;
	"obj")
		q=1
	;;
	esac
done
#
if [ $k = 0 ]; then
	mkdir ../bin
	echo ../bin folder created
else
	echo ../bin folder exists
fi
#
if [ $q = 0 ]; then
	mkdir ../obj
	echo ../obj folder created
else
	echo ../obj folder exists
fi
#
endef
#
folder_set: ; @$(value folder_set)

.ONESHELL:
#
#*******************************************************
#
bin_path = ../bin
obj_path = ../obj
#
QP_pgr_src = QP_module.f90 QP_program.f90
QP_pgr_obj = QP_module.o   QP_program.o #$(obj_path)/$@QP_module.o $(obj_path)/$@QP_program.o

.SUFFIXES: 
.SUFFIXES: .o .f90 .oOMP 
#
.f90.o:
	$(FC)  -c -cpp $(mod_lib)  "$<" -o $(obj_path)/$@


#
all: folder_set QP_pgr
#
.phony : all
#
#
QP_pgr:	$(QP_pgr_obj) Makefile
	$(info)
	$(FC) -o $(bin_path)/$@.out $(obj_path)/*.o  -cpp $(mod_lib) $(LIBS)


clean:
	@(rm -f *.o *.mod)
	@(rm -f $(bin_path)/*out)
	@(rm -f $(obj_path)/*o)

