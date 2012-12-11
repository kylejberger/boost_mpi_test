all: DEM.exe

DEM.exe: main.o 
  	mpic++   -Wl,-rpath,/curc/tools/x_86_64/rh5/boost/1.50/openmpi/1.6/intel/12.1.4/torque/2.5.11/lib, main.o  -o DEM.exe -L/curc/tools/x_86_64/rh5/boost/1.50/openmpi/1.6/intel/12.1.4/torque/2.5.11/lib -lboost_serialization -lboost_mpi -L/curc/tools/x_86_64/rh5/openmpi/1.6/intel/12.1.4/torque/2.5.11/ib/lib -lmpi
		
main.o: main.cpp
		mpic++ -I/curc/tools/x_86_64/rh5/boost/1.50/openmpi/1.6/intel/12.1.4/torque/2.5.11/include -c -g -Wall main.cpp
		
clean:
		rm -rf *o DEM.exe