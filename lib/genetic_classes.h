#ifndef LIB_GENETIC_CLASSSES_H_
#define LIB_GENETIC_CLASSES_H_

class Gene
{
    public:
	    Symbolic elem; //one symbolic gene
	    int oper; //operator of the gene: 1 - plus, 2 - minus, 3 - mult, 4 - div, 5 - pow
};

struct data
{
	double y;
	double x;
};

class Individ
{
    public:
	    Symbolic ind; //one symbolic individual
	    vector<Gene> genes; //genes of the individual
	    double fit; //Value of the fitness function for the individual
};

class Population
{
    public:
        std::vector<Individ> inds; //all of the individuals in the population
};


#endif // LIB_GENETIC_H_
