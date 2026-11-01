CXX = g++
FC = gfortran
CXXFLAGS = -O3 -march=native -fopenmp -I /usr/include/eigen3
FCFLAGS = -ffree-line-length-none -fopenmp

all: test_fortran_trajectory

easywfo_cisd_bridge.o: easywfo_cisd_bridge.cpp
	$(CXX) $(CXXFLAGS) -c easywfo_cisd_bridge.cpp -o easywfo_cisd_bridge.o

test_fortran_trajectory: test_fortran_trajectory.f90 easywfo_cisd_bridge.o
	$(FC) $(FCFLAGS) test_fortran_trajectory.f90 easywfo_cisd_bridge.o -lstdc++ -o test_fortran_trajectory

clean:
	rm -f *.o test_fortran_trajectory st_adi_*.txt
