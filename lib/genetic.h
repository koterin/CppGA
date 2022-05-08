#ifndef GENETIC_H_
#define GENETIC_H_

#define PI 3.14159265
#define NGA_LOOPS_LIMIT 100 // number of loops limit (limit * number of elements)
#define GAsize 10 // number of inds in the numerical GA

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

	    //Function for calculation the fitness func for the individual
	    void CalcFit(vector<vector<struct data>> ExpData, Symbolic t, vector<Symbolic> Variables,
																vector<vector<double>> VarValues)
	    {
		    struct data buf;
		    double dev = 0.0; //Devitation ind current v from expdata v
		    double devSUM = 0.0; //Sum devitation
		    Symbolic y;

		    for (int f = 0; f < ExpData.size(); f++) //loops for files
		    {
			    for (int i = 0; i < ExpData[f].size(); i++) //loops for expdata lines
			    {
				    buf = ExpData[f][i];
				    y.auto_expand = 0;
				    y.simplified = 0;
				    y = ind[t == buf.x];
				    for (int j = 2; j < Variables.size(); j++) //always from the 2: 0 - v(t), 1 - t.
				    {
					    y = y[Variables[j] == VarValues[j-1][f]];
				    }
				    y.upr();

				    auto yClone = y->clone();
				    if ((typeid(*(yClone)) != typeid(Number<double>)) && (typeid(*(yClone)) != typeid(Number<int>)))
				    {
					    dev = 1e+30;
					    devSUM += dev;
				    } else
				    {
					    dev = (pow((double(y) - buf.y), 2) / (pow(buf.y, 2)));
					    devSUM += dev;
				    }

                    yClone->unreference(yClone);
				    y.auto_expand = 1;
				    y.simplified = 1;
			    }
		    }

		    //Checking if devSUM is too small so we won't get zeros in fits
		    if (devSUM > 1e9)
		    {
			    fit = 0;
		    } else
		    {	
			    fit = (1 / (1 + sqrt(devSUM))) * 100;
		    }
		    return;
	    }
};


#endif // GENETIC_H_
