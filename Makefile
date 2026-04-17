CXX      = g++
NVCC     = nvcc
CXXFLAGS  = -O3 -std=c++17 -fopenmp
NVCCFLAGS = -O3 -std=c++17 -arch=sm_70 -Xcompiler -fopenmp
LDFLAGS   = -fopenmp
CGAL_FLAGS = -lCGAL -lgmp -lmpfr

PT_EXEC      = pt
CGAL_EXEC    = cgal
DC_EXEC      = dc
TESTING_EXEC = testing
ALL_EXECS    = $(PT_EXEC) $(CGAL_EXEC) $(DC_EXEC) $(TESTING_EXEC)

all: run-others

run-others: run-pt
	@echo "Starting others..."
	$(MAKE) run-testing

run-pt: $(PT_EXEC)
	@echo "Executing pt..."
	./$(PT_EXEC)

$(PT_EXEC): pt.cpp
	@echo "Compiling pt.cpp..."
	$(CXX) $(CXXFLAGS) -o $(PT_EXEC) pt.cpp $(LDFLAGS)

# Chain: run-testing -> run-dc -> run-cgal (strict serial order)
run-testing: run-dc $(TESTING_EXEC)
	@echo "Executing testing..."
	./$(TESTING_EXEC)

run-dc: run-cgal $(DC_EXEC)
	@echo "Executing dc..."
	./$(DC_EXEC)

run-cgal: $(CGAL_EXEC)
	@echo "Executing cgal..."
	./$(CGAL_EXEC)

$(CGAL_EXEC): cgal.cpp
	@echo "Compiling cgal.cpp..."
	$(CXX) $(CXXFLAGS) -o $(CGAL_EXEC) cgal.cpp $(CGAL_FLAGS) $(LDFLAGS)

$(DC_EXEC): dc.cpp
	@echo "Compiling dc.cpp..."
	$(CXX) $(CXXFLAGS) -o $(DC_EXEC) dc.cpp $(LDFLAGS)

$(TESTING_EXEC): testing.cu
	@echo "Compiling testing.cu..."
	$(NVCC) $(NVCCFLAGS) -o $(TESTING_EXEC) testing.cu -ljemalloc

clean:
	rm -f $(ALL_EXECS) *.csv

rebuild: clean all

.PHONY: all run-pt run-others run-cgal run-dc run-testing clean rebuild
