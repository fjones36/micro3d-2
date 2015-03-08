.SUFFIXES:

SUFFIXlinux     = linux
SUFFIXlinuxMyrtle     = linuxMyrtle
SUFFIXlinuxMPI  = linuxMPI
# ------------  ------------
SUFFIX = ${SUFFIXlinux}


CPPlinux    = g++
CPPlinuxMyrtle    = g++
CPPlinuxMPI = mpiCC
# ------------  ------------
CPP = ${CPP${SUFFIX}}

MPIRUN      = mpirun
MPINUMCPU   = 4
MACHINEFILE = machines

OBJDIR = objects

DEFINESlinux    = 
DEFINESlinuxMyrtle    = -DMYRTLE
DEFINESlinuxMPI = -D__MPI__
# ------------  ------------
DEFINES = ${DEFINES${SUFFIX}}

INCLUDES = -I. -Ibase

LIBS = -lm 

OBJECTS = ${OBJDIR}/ran01.o ${OBJDIR}/data_array.o \
		${OBJDIR}/polygon.o ${OBJDIR}/generate.o \
		${OBJDIR}/neuron.o ${OBJDIR}/move.o ${OBJDIR}/vec.o  \
		${OBJDIR}/stringfunc.o ${OBJDIR}/streamfunc.o \
		${OBJDIR}/data.o ${OBJDIR}/funcs.o ${OBJDIR}/collisions.o \
		${OBJDIR}/grid.o ${OBJDIR}/integer.o ${OBJDIR}/pair_corr.o \
		${OBJDIR}/glia.o ${OBJDIR}/Custom.o

DO_OBJ = ${OBJDIR}/stringfunc.o 

OPT = -O3
#OPT = -O3 -fno-stack-protector
#OPT =

EXEPRE   = slicegen_
EXEPOS   = _0.13
EXEC     = ${EXEPRE}${SUFFIX}${EXEPOS}
BINDIR   = ~/bin
INPUT    = input_v0.13.txt
NUMCONFS = 1
SILENT   = > /dev/null
#SILENT   = 

# TO TRANSFER TO THE CLUSTER
RAND3DSRC = /home/ccruz/micro3D/v0.13/ 
RAND3DDEST= ccruz@cpsbu:~/micro3D/v0.13/
RAND3DEXCL=  --exclude="archive" --exclude="objects" --exclude="poster" \
	--exclude="steps2brain_1" --exclude="steps2brain_2" --exclude="test" 

first_make:
	${MAKE} -s -f makefile usage

usage:
	echo ""
	echo "TARGETS FOR THE SLICEGEN PROGRAM:"
	echo ""
	echo "   MAIN PROGRAM:"
	echo "		 1.  make build"
	echo "		 1a. make build_myrtle"
	echo "		 2.  make build_mpi"
	echo ""
	echo "		10. make run_mpi"
	echo ""
	echo "   UTILS:"
	echo "		 make clean"
	echo "		 make transfer_cpsbu_pretend"
	echo "		 make transfer_cpsbu"
	echo ""

1:build
1a:build_myrtle
2:build_mpi
10:run_mpi

build:
	make -f makefile all SUFFIX='linux'
build_myrtle:
	make -f makefile all SUFFIX='linuxMyrtle'
build_mpi:
	make -f makefile all SUFFIX='linuxMPI'

# -----------------------------------------------------------

all: ${OBJDIR} ${OBJECTS}
	${CPP} ${OPT} ${OBJECTS} -o ${EXEC} ${LIBS} 
	cp ${EXEC} ${BINDIR}

run_mpi:
	${MPIRUN} -np ${MPINUMCPU} -machinefile ${MACHINEFILE} ./${EXEPRE}linuxMPI${EXEPOS}  ${INPUT} ${NUMCONFS}  ${SILENT}

clean:
	rm -f ${OBJDIR}/*.o ${EXEC} *~ *.orig core* 
	find . -name '*rsyncBak*' -exec rm {} \;

transfer_cpsbu:
	rsync -e ssh -avz --delete ${RAND3DSRC} ${RAND3DDEST} ${RAND3DEXCL}

transfer_cpsbu_pretend:
	rsync -e ssh -avz -n --delete ${RAND3DSRC} ${RAND3DDEST} ${RAND3DEXCL}


${OBJDIR}:
	mkdir ${OBJDIR}

# -------------------------------------------------------
${OBJDIR}/generate.o:    generate.cc generate.h polygon.h data_array.h move.h defs.h data.h funcs.h collisions.h grid.h base/infobjarray.h base/macros.h neuron.h pair_corr.h glia.h Custom.h
	${CPP} -c ${OPT} ${DEFINES} ${INCLUDES} ${@F:.o=.cc} -o $@
${OBJDIR}/polygon.o:    polygon.cc polygon.h base/macros.h
	${CPP} -c ${OPT} ${DEFINES} ${INCLUDES} ${@F:.o=.cc} -o $@
${OBJDIR}/neuron.o:    neuron.cc neuron.h base/vec.h base/macros.h
	${CPP} -c ${OPT} ${DEFINES} ${INCLUDES} ${@F:.o=.cc} -o $@
${OBJDIR}/glia.o:    glia.cc glia.h base/vec.h base/macros.h
	${CPP} -c ${OPT} ${DEFINES} ${INCLUDES} ${@F:.o=.cc} -o $@
${OBJDIR}/data_array.o:    data_array.cc data_array.h base/macros.h
	${CPP} -c ${OPT} ${DEFINES} ${INCLUDES} ${@F:.o=.cc} -o $@
${OBJDIR}/move.o:    move.cc move.h neuron.h defs.h data.h collisions.h grid.h base/infobjarray.h polygon.h base/macros.h glia.h
	${CPP} -c ${OPT} ${DEFINES} ${INCLUDES} ${@F:.o=.cc} -o $@
${OBJDIR}/ran01.o:    ran01.cc ran01.h
	${CPP} -c ${OPT} ${DEFINES} ${INCLUDES} ${@F:.o=.cc} -o $@
${OBJDIR}/grid.o:    grid.cc grid.h defs.h base/integer.h base/infobjarray.h neuron.h data.h glia.h
	${CPP} -c ${OPT} ${DEFINES} ${INCLUDES} ${@F:.o=.cc} -o $@
${OBJDIR}/data.o:    data.cc data.h
	${CPP} -c ${OPT} ${DEFINES} ${INCLUDES} ${@F:.o=.cc} -o $@
${OBJDIR}/funcs.o:    funcs.cpp funcs.h data.h neuron.h data.h base/infobjarray.h defs.h base/macros.h glia.h
	${CPP} -c ${OPT} ${DEFINES} ${INCLUDES} ${@F:.o=.cpp} -o $@
${OBJDIR}/vec.o:    base/vec.cxx base/vec.h
	${CPP} -c ${OPT} ${DEFINES} ${INCLUDES} base/${@F:.o=.cxx} -o $@
${OBJDIR}/integer.o:    base/integer.cc base/integer.h
	${CPP} -c ${OPT} ${DEFINES} ${INCLUDES} base/${@F:.o=.cc} -o $@
${OBJDIR}/stringfunc.o:    base/stringfunc.cxx base/stringfunc.h
	${CPP} -c ${OPT} ${DEFINES} ${INCLUDES} base/${@F:.o=.cxx} -o $@
${OBJDIR}/streamfunc.o:    base/streamfunc.cxx base/streamfunc.h
	${CPP} -c ${OPT} ${DEFINES} ${INCLUDES} base/${@F:.o=.cxx} -o $@
${OBJDIR}/collisions.o:    collisions.cc collisions.h
	${CPP} -c ${OPT} ${DEFINES} ${INCLUDES} ${@F:.o=.cc} -o $@
${OBJDIR}/pair_corr.o:    pair_corr.cc pair_corr.h defs.h data.h base/vec.h base/macros.h
	${CPP} -c ${OPT} ${DEFINES} ${INCLUDES} ${@F:.o=.cc} -o $@
#${OBJDIR}/Custom.o:    Custom.cpp Custom.h
#	${CPP} -c ${OPT} ${DEFINES} ${INCLUDES} ${@F:.o=.cc} -o $@


# -------------------------------------------------------
