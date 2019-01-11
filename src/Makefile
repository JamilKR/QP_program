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
	echo
else
	echo ../bin folder exists
	echo
fi
#
if [ $q = 0 ]; then
	mkdir ../obj
	echo ../obj folder created
	echo
else
	echo ../obj folder exists
	echo
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
QP_pgr_obj = $(QP_pgr_src:.f90=.o) 
##

vpath %.f ./quadpack # tells make to look for any dependency whose name ends in `.f' in the directory ./quadpack.

Quadpack_src = dqc25s.f  dqcheb.f  dqk15w.f  dqwgts.f
Quadpack_obj = $(Quadpack_src:.f=.o)


.SUFFIXES: 
.SUFFIXES: .o .f90 .oOMP .f
#
.f90.o:
	$(info )
	$(FC)  -c -cpp $(mod_lib)  "$<" -o $(obj_path)/$@
	$(info )



.f.o:
	$(info )
	$(FC) -c -cpp  "$<" $(LIBS) -o $(obj_path)/$@
	$(info )


#
all: folder_set Quadpack QP_pgr 
#
.phony : all
#
#
#*******************************************************************************
QP_pgr:	$(QP_pgr_obj) Makefile
	$(info )
	$(FC) -o $(bin_path)/$@.out  $(addprefix $(obj_path)/,$(QP_pgr_obj))  -cpp $(mod_lib) $(LIBS)
	$(info )


#*******************************************************************************
Quadpack:	$(Quadpack_obj) Makefile
		$(info )
		ar cr $(obj_path)/quadpack.a $(addprefix $(obj_path)/,$(Quadpack_obj))
		$(info )


#*******************************************************************************
clean:
	@(rm -f *.o *.mod)
	@(rm -f $(bin_path)/*out)
	@(rm -f $(obj_path)/*o $(obj_path)/*a)

