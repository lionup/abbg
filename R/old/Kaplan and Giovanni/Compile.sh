#! /bin/sh                                                                                                             
#remember: cant have spaces between equals signs in assignments
#input (1) Debug/Release (2) run
# Set up environment variables for Mac and Intel Fortan 11.0

source /opt/intel/Compiler/11.1/058/bin/ifortvars.sh intel64

#Debug version
if [ "$1" == "Debug" ]; then
	echo "Compiling in " $1 "mode..."
	
	# Delete existing files
	rm *.o *.mod
	
	# Compile and link debug version: produces all debugging output, no optimizations, for use with Fx3
	ifort -m64 -g -debug all -implicitnone -save-temps -warn all -fp-stack-check -ftrapuv -check all -traceback -check noarg_temp_created \
		-gen-interfaces -warn interfaces -openmp -O3 \
		Parameters.f90 Globals.f90 Procedures.f90  OtherCode/random.f90 \
		SetParameters.f90 LoadData.f90 Grids.f90 	\
		Decisions.f90 Simulate.f90 SaveOutput.f90 \
		Functions/FnTax.f90 Functions/FnTaxParam.f90 Functions/FnTaxParamNet.f90 \
		Functions/FnSSParam.f90 Functions/FnGrossInc.f90 Functions/FnGridTrans.f90 \
		Functions/FnEqumR.f90 Functions/FnBetaKY.f90 Functions/FnGridPerm.f90 \
		OtherCode/golden.f90 OtherCode/mnbrak.f90 OtherCode/cumnor.f90 \
		OtherCode/rtsecBetaKY.f90 OtherCode/rtsecEqumR.f90 OtherCode/rtsecSSParam.f90 \
		OtherCode/zbrentTaxParam.f90 OtherCode/zbrentTaxParamNet.f90 OtherCode/rtnewtGrossInc.f90 \
		Main.f90 -o MainDebug.out

	# Create debug file with all script information
	dsymutil MainDebug.out


	# -g tells compiler to produce symbol info
	# -save-temps tells compiler not to delete intermediate files and puts them in working directory

	if [ "$2" == "run" ]; then
		source /opt/intel/Compiler/11.1/058/bin/ifortvars.sh intel64
		./MainDebug.out
	fi

# Release version with openmp
elif [ "$1" == "Release" ]; then
		echo "Compiling in " $1 "mode..."
		
		ifort -m64 -O3 -openmp -traceback \
			Parameters.f90 Globals.f90 Procedures.f90  OtherCode/random.f90 \
			SetParameters.f90 LoadData.f90 Grids.f90 	\
			Decisions.f90 Simulate.f90 SaveOutput.f90 \
			Functions/FnTax.f90 Functions/FnTaxParam.f90 Functions/FnTaxParamNet.f90 \
			Functions/FnSSParam.f90 Functions/FnGrossInc.f90 Functions/FnGridTrans.f90 \
			Functions/FnEqumR.f90 Functions/FnBetaKY.f90 Functions/FnGridPerm.f90 \
			OtherCode/golden.f90 OtherCode/mnbrak.f90 OtherCode/cumnor.f90 \
			OtherCode/rtsecBetaKY.f90 OtherCode/rtsecEqumR.f90 OtherCode/rtsecSSParam.f90 \
			OtherCode/zbrentTaxParam.f90 OtherCode/zbrentTaxParamNet.f90 OtherCode/rtnewtGrossInc.f90 \
			Main.f90 -o MainRelease.out

			
	# -fast inclues-O3 plus some other optimization things
	# -openmp-link static forces the openmp libraries to be linked staticaly.

	
	if [ "$2" == "run" ]; then
		source /opt/intel/Compiler/11.1/058/bin/ifortvars.sh intel64
		./MainRelease.out
	fi


	fi

