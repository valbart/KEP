CCALL = g++ -DIL_STD -g -fopenmp
CC = g++ -DIL_STD -g
CFLAGS  = -I/opt/ibm/ILOG/CPLEX_Studio128/opl/include -I/opt/ibm/ILOG/CPLEX_Studio128/concert/include
LDFLAGS = -L/opt/ibm/ILOG/CPLEX_Studio128/concert/lib/x86-64_linux/static_pic/ -lconcert -L/opt/ibm/ILOG/CPLEX_Studio128/cplex/lib/x86-64_linux/static_pic/ -lilocplex -lcplex -lcplexdistmip -lm -lpthread -ldl -lboost_program_options

all: RobustKEPSolver

RobustKEPSolver: RobustKEPSolver.o CycleFormulation.o Graph_Generation.o HPIEF.o LimitedWorld.o MathHeuristics.o PICEF.o PreTestOpt.o SubsetRecourse.o TestSolutions.o UnlimitedCycle.o
		$(CCALL) -o $@ $^ $(LDFLAGS)

RobustKEPSolver.o: RobustKEPSolver.cpp Structures.h Functions.h
		$(CC) -c $(CFLAGS) $<

CycleFormulation.o: CycleFormulation.cpp Structures.h Functions.h
		$(CC) -c $(CFLAGS) $<

Graph_Generation.o: Graph_Generation.cpp Structures.h Functions.h
		$(CC) -c $(CFLAGS) $<

HPIEF.o: HPIEF.cpp Structures.h Functions.h
		$(CC) -c $(CFLAGS) $<

LimitedWorld.o: LimitedWorld.cpp Structures.h Functions.h
		$(CC) -c $(CFLAGS) $<

MathHeuristics.o: MathHeuristics.cpp Structures.h Functions.h
		$(CC) -c $(CFLAGS) $<

PICEF.o: PICEF.cpp Structures.h Functions.h
		$(CC) -c $(CFLAGS) $<

PreTestOpt.o: PreTestOpt.cpp Structures.h Functions.h
		$(CC) -c $(CFLAGS) $<

SubsetRecourse.o: SubsetRecourse.cpp Structures.h Functions.h
		$(CC) -c $(CFLAGS) $<

TestSolutions.o: TestSolutions.cpp Structures.h Functions.h
		$(CC) -c $(CFLAGS) $<

UnlimitedCycle.o: UnlimitedCycle.cpp Structures.h Functions.h
		$(CC) -c $(CFLAGS) $<

.PHONY: clean cleanest

clean:
		rm RobustKEPSolver

cleanall:
		rm *.o
