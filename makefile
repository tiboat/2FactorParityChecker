compiler=gcc
flags= -lm -O4 -march=native -fopenmp

all: 2FactorParityChecker

2FactorParityChecker: 2FactorParityChecker.c
	$(compiler) -g -o 2FactorParityChecker 2FactorParityChecker.c perfMatchings/enumPerfMatchings.c matchmaker/matching.c matchmaker/cheap.c nauty/nauty.a $(flags)

clean:
	rm 2FactorParityChecker
