# Compiler settings
CXX = g++
NVCC = nvcc
CXXFLAGS = -O3 -std=c++17 -fopenmp
NVCCFLAGS = -O3 -std=c++17 -Xcompiler -fopenmp
LDFLAGS = -fopenmp

# CGAL flags (if needed)
CGAL_FLAGS = -lCGAL -lgmp -lmpfr

# Target executables
PT_EXEC = pt
CGAL_EXEC = cgal
DC_EXEC = dc
TESTING_EXEC = testing

# All executables
ALL_EXECS = $(PT_EXEC) $(CGAL_EXEC) $(DC_EXEC) $(TESTING_EXEC)

# Default target: compile pt first, run it, then compile and run others
all: run-pt run-others

# Compile and run pt.cpp first
run-pt: $(PT_EXEC)
	@echo "Executing pt..."
	./$(PT_EXEC)

# Compile pt.cpp
$(PT_EXEC): pt.cpp
	@echo "Compiling pt.cpp..."
	$(CXX) $(CXXFLAGS) -o $(PT_EXEC) pt.cpp $(LDFLAGS)

# Compile and run the rest after pt
run-others: run-cgal run-dc run-testing

# CGAL executable
$(CGAL_EXEC): cgal.cpp
	@echo "Compiling cgal.cpp..."
	$(CXX) $(CXXFLAGS) -o $(CGAL_EXEC) cgal.cpp $(CGAL_FLAGS) $(LDFLAGS)

# DC executable
$(DC_EXEC): dc.cpp
	@echo "Compiling dc.cpp..."
	$(CXX) $(CXXFLAGS) -o $(DC_EXEC) dc.cpp $(LDFLAGS)

# CUDA testing executable with jemalloc
$(TESTING_EXEC): testing.cu
	@echo "Compiling testing.cu with jemalloc..."
	$(NVCC) $(NVCCFLAGS) -o $(TESTING_EXEC) testing.cu -ljemalloc

# Run individual executables
run-cgal: $(CGAL_EXEC)
	@echo "Executing cgal..."
	./$(CGAL_EXEC)

run-dc: $(DC_EXEC)
	@echo "Executing dc..."
	./$(DC_EXEC)

run-testing: $(TESTING_EXEC)
	@echo "Executing testing..."
	./$(TESTING_EXEC)

clean:
	rm -f $(ALL_EXECS)
	rm -f *.csv

rebuild: clean all

.PHONY: all run-pt run-others run-cgal run-dc run-testing clean rebuild