#--------------------------------------
# Makefile
#--------------------------------------

#Use gfortran
F90 = gfortran
FFLAG := -O3 -fopenmp

#Include path - for modules created
IMOD=-Ibin 

#Objects
OBJ=bin/fv_arrays.obj \
bin/fv_grid_tools.obj \
bin/fv_duogrid.obj \
bin/fv_control.obj \
bin/test_cases.obj \
bin/tp_core.obj \
bin/sw_core.obj \
bin/dyn_core.obj \
bin/atmosphere.obj \

#Compile and build all
all: header config bin/main ending

#Make all and run executable
run: all
	./main

#Print heading
header:
	@echo --------------------------------------------
	@echo Compiling and building the software   
	@echo --------------------------------------------
	@echo 
	@echo 

#Configure Enviroment (directories)
config:
	chmod +x sh/*.sh
	. sh/dirs.sh

#fv_arrays
bin/fv_arrays.obj: model/fv_arrays.f90
	$(F90) $(FFLAG) -c  $^ -o $@ $(IMOD)
	mv fv_arrays.mod bin/.

#fv_grid_tools
bin/fv_grid_tools.obj: tools/fv_grid_tools.f90
	$(F90) $(FFLAG) -c  $^ -o $@ $(IMOD)
	mv fv_grid_tools.mod bin/.

#fv_duogrid
bin/fv_duogrid.obj: tools/fv_duogrid.f90
	$(F90) $(FFLAG) -c  $^ -o $@ $(IMOD)
	mv fv_duogrid.mod bin/.

#fv_control
bin/fv_control.obj: tools/fv_control.f90
	$(F90) $(FFLAG) -c  $^ -o $@ $(IMOD)
	mv fv_control.mod bin/.

#test_cases
bin/test_cases.obj: tools/test_cases.f90
	$(F90) $(FFLAG) -c  $^ -o $@ $(IMOD)
	mv test_cases.mod bin/.

#tp_core
bin/tp_core.obj: model/tp_core.f90
	$(F90) $(FFLAG) -c  $^ -o $@ $(IMOD)
	mv tp_core.mod bin/.

#sw_core
bin/sw_core.obj: model/sw_core.f90
	$(F90) $(FFLAG) -c  $^ -o $@ $(IMOD)
	mv sw_core.mod bin/.

#dyn_core
bin/dyn_core.obj: model/dyn_core.f90
	$(F90) $(FFLAG) -c  $^ -o $@ $(IMOD)
	mv dyn_core.mod bin/.
	
#atmosphere
bin/atmosphere.obj: driver/atmosphere.f90
	$(F90) $(FFLAG) -c  $^ -o $@ $(IMOD)
	mv atmosphere.mod bin/.
	
#Main executable
bin/main: main.f90 $(OBJ)
	$(F90) $(FFLAG)  main.f90 $(OBJ) -o $@ $(IMOD)

#Creates a link for executable and prints ending
ending: 
	chmod +x sh/link.sh
	sh/link.sh
	@echo End of compilation
	@echo
	@echo "Set parameter files (pars / *.par )" 
	@echo "   and then run 'main'"
	@echo "------------------------------------------------------------------"


# Create tarball and backup
archive: 
        #Backup all important files
	chmod +x sh/backup.sh
	./sh/backup.sh

#Clean targets
clean: 
	rm -rf bin/*.obj bin/*.o bin/*.mod	
	rm -rf bin/main*
	rm -rf *~

cleandata: clean
	rm -rf data/
	rm -rf graphs/
	rm -rf bin/
	rm main

cleangrids: clean
	rm -rf graphs/
	rm -rf grids/

cleanall: clean cleandata cleangrids
