#ifndef LIB_GENETIC_FUNCTIONS_H_
#define LIB_GENETIC_FUNCTIONS_H_

#include "genetic_classes.h"

//Checking if input y is already a complex formula
vector<Gene> InputGeneDecomposition(Symbolic y)
{
    vector<Gene> bufGenes, newGenes;
    bufGenes.clear();
    newGenes.clear();

    Gene outputGene;

    list <Symbolic> bufList;
    bufList.clear();

    auto bufY = y->clone();
    //SystemDiaDebug.Write(bufY);

    if (typeid(*bufY) == typeid(Sum))
    {
        bufList = (*((Sum*)(bufY))).summands;
        auto iter = bufList.begin();
        int oper = 1;

        for (int i = 0; i < bufList.size(); i++)
        {
            auto bufNum = Symbolic(*iter)->clone();
            if ((typeid(*bufNum) == typeid(Number<double>)) || (typeid(*bufNum) == typeid(Number<int>)))
            {
                outputGene.elem = *iter;
                outputGene.oper = oper;
                bufGenes.push_back(outputGene);
            }
            if (typeid(*bufNum) == typeid(Symbol))
            {
                outputGene.elem = *iter;
                outputGene.oper = oper;
                bufGenes.push_back(outputGene);
            }
            else if ((typeid(*bufNum) != typeid(Number<double>)) && (typeid(*bufNum) != typeid(Number<int>))
                && (typeid(*bufNum) != typeid(Symbol)))
            {
                newGenes.clear();
                newGenes = InputGeneDecomposition(*iter);
                newGenes[0].oper = oper;
                for (int j = 0; j < newGenes.size(); j++)
                {
                    bufGenes.push_back(newGenes[j]);
                }
            }

            iter++;
            bufNum->unreference(bufNum);
        }    
    }

    if (typeid(*bufY) == typeid(Product))
    {
        bufList = (*((Product*)(bufY))).factors;
        auto iter = bufList.begin();
        int oper = 3;

        for (int i = 0; i < bufList.size(); i++)
        {
            auto bufNum = Symbolic(*iter)->clone();
            if ((typeid(*bufNum) == typeid(Number<double>)) || (typeid(*bufNum) == typeid(Number<int>)))
            {
                outputGene.elem = *iter;
                outputGene.oper = oper;
                bufGenes.push_back(outputGene);
            }
            if (typeid(*bufNum) == typeid(Symbol))
            {
                outputGene.elem = *iter;
                outputGene.oper = oper;
                bufGenes.push_back(outputGene);
            }
            else if ((typeid(*bufNum) != typeid(Number<double>)) && (typeid(*bufNum) != typeid(Number<int>))
                && (typeid(*bufNum) != typeid(Symbol)))
            {
                newGenes.clear();
                newGenes = InputGeneDecomposition(*iter);
                newGenes[0].oper = oper;
                for (int j = 0; j < newGenes.size(); j++)
                {
                    bufGenes.push_back(newGenes[j]);
                }
            }

            iter++;
            bufNum->unreference(bufNum);
        }
    }

    if (typeid(*bufY) == typeid(Power))
    {
        bufList = (*((Symbol*)(&(*((Power*)(bufY)))))).parameters;
        auto iter = bufList.begin();
        int oper = 5;

        for (int i = 0; i < bufList.size(); i++)
        {
            auto bufNum = Symbolic(*iter)->clone();
            if ((typeid(*bufNum) == typeid(Number<double>)) || (typeid(*bufNum) == typeid(Number<int>)))
            {
                outputGene.elem = *iter;
                outputGene.oper = oper;
                bufGenes.push_back(outputGene);
            }
            if (typeid(*bufNum) == typeid(Symbol))
            {
                outputGene.elem = *iter;
                outputGene.oper = oper;
                bufGenes.push_back(outputGene);
            }
            else if ((typeid(*bufNum) != typeid(Number<double>)) && (typeid(*bufNum) != typeid(Number<int>))
                && (typeid(*bufNum) != typeid(Symbol)))
            {
                newGenes.clear();
                newGenes = InputGeneDecomposition(*iter);
                newGenes[0].oper = oper;
                for (int j = 0; j < newGenes.size(); j++)
                {
                    bufGenes.push_back(newGenes[j]);
                }
            }

            iter++;
            bufNum->unreference(bufNum);
        }
    }

    bufY->unreference(bufY);
    return (bufGenes);
}

//Function for creating new individual from the genes given
//WARNING! Input genes vector MUST BE SORTED and organized the way it should be in the output Individual
Individ IndFromGenes(vector<Gene> genes)
{
    Individ outputInd;
    int genNum = genes.size();
    Gene bufGene = genes[0];
    outputInd.ind = bufGene.elem;

    //outputInd.genes.reserve(genNum);
    //For the 1st element
    outputInd.genes.push_back(bufGene);
    //outputInd.ind = genes[0].elem;

    for (int i = 1; i < genNum; i++)
    {
        auto buf = genes[i].elem->clone();
        if (((typeid(*buf) == typeid(Number<double>)) || (typeid(*buf) == typeid(Number<int>)))
            && (double(genes[i].elem) > 100))
        {
            genes[i].elem = 99.9;
        }

        //So power won't be simplified
        if ((typeid(*buf) == typeid(Number<int>)))
        {
            genes[i].elem = genes[i].elem - 0.01;
        }

        if (genes[i].oper == 1)
        {
            outputInd.ind = (outputInd.ind) + genes[i].elem;
            outputInd.genes.push_back(genes[i]);
        }

        if (genes[i].oper == 2)
        {
            outputInd.ind = (outputInd.ind) - genes[i].elem;
            outputInd.genes.push_back(genes[i]);
        }

        if (genes[i].oper == 3)
        {
            outputInd.ind = (outputInd.ind) * genes[i].elem;
            outputInd.genes.push_back(genes[i]);
        }

        if (genes[i].oper == 4)
        {
            outputInd.ind = (outputInd.ind) / genes[i].elem;
            outputInd.genes.push_back(genes[i]);
        }

        if (genes[i].oper == 5)
        {
            outputInd.ind = pow((outputInd.ind), genes[i].elem);
            outputInd.genes.push_back(genes[i]);
        }
        buf->unreference(buf);
    }
    return(outputInd);
}

Individ ZeroIndTermination(Individ IndZero, Symbolic vars)
{
    Gene plusGene;
    plusGene.elem = vars;
    plusGene.oper = 1;

    std::cout << "\nZeroIndTermination, bufInd is " << IndZero.ind << " Missing var is " << vars << std::endl;
    std::cout << "previous genes length is " << IndZero.genes.size() << std::endl;

    if (IndZero.genes.size() == 0)
    {
        IndZero.genes.push_back(plusGene);
    }

    IndZero.genes.push_back(plusGene);
    std::cout << "Now genes size must be plus one: " << IndZero.genes.size() << std::endl;

    IndZero = IndFromGenes(IndZero.genes);
    std::cout << "ind after IndFromGenes " << IndZero.ind << std::endl;
    IndZero.genes = InputGeneDecomposition(IndZero.ind);

    return(IndZero);
}

Individ SearchForAbscentVars(Individ IndZero, vector<Symbolic> Variables)
{
    vector<int> varCount(Variables.size(), 0);

    //Looking for abscent variables
    for (int j = 1; j < Variables.size(); j++)
    {
        for (int i = 0; i < IndZero.genes.size(); i++)
        {
            if (IndZero.genes[i].elem == Variables[j])
            {
                varCount[j]++;
            }
        }

        if (varCount[j] == 0)
        {
            IndZero = ZeroIndTermination(IndZero, Variables[j]);
        }

        //Checking if there're still abscent variables
        for (int i = 0; i < IndZero.genes.size(); i++)
        {
            if (IndZero.genes[i].elem == Variables[j])
            {
                varCount[j]++;
            }
        }

        if (varCount[j] == 0)
        {
            IndZero = SearchForAbscentVars(IndZero, Variables);
        }
    } 

    ////Checking if there're still abscent variables
    //for (int j = 1; j < Variables.size(); j++)
    //{
    //    for (int i = 0; i < IndZero.genes.size(); i++)
    //    {
    //        if (IndZero.genes[i].elem == Variables[j])
    //        {
    //            varCount[j]++;
    //        }
    //    }

    //    if (varCount[j] == 0)
    //    {
    //        IndZero = SearchForAbscentVars(IndZero, Variables);
    //    }
    //}

    return(IndZero);
}

double CalcFit(Symbolic ind, vector<vector<struct data>> ExpData, Symbolic x, vector<Symbolic> Variables,
																vector<vector<double>> VarValues)
{
    struct data buf;
    double dev = 0.0; //Devitation ind current v from expdata v
    double devSUM = 0.0; //Sum devitation
    Symbolic y;
    double fit;

    // loops for files
    for (int f = 0; f < ExpData.size(); f++)
    {
        // loops for expdata file lines
		for (int i = 0; i < ExpData[f].size(); i++)
        {
		    buf = ExpData[f][i];
		    y.auto_expand = 0;
		    y.simplified = 0;
		    y = ind[x == buf.x];
		    for (int j = 2; j < Variables.size(); j++) //always from the 2: 0 - y, 1 - x
			    y = y[Variables[j] == VarValues[j-1][f]];
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
        fit = 0;
    else	
	    fit = (1 / (1 + sqrt(devSUM))) * 100;
	return (fit);
}

Population CreatePop(Population pop, vector<Symbolic> Variables, int numInd, std::string foutname)
{
    double coeff = 0.0; //random coefficient for 1st pop creation
    Individ indZero;
    Gene genZero;
    Individ outputInd;
    Gene outputGene;

    pop.inds.clear();
    pop.inds.reserve(numInd);

    std::ofstream fout;
    fout.open(foutname, std::ios_base::app);

    vector<Gene> startGenes;
    startGenes = InputGeneDecomposition(Variables[0]);
    std::cout << "Final input is [ ";
    for (int i = 0; i < startGenes.size(); i++)
    {
        std::cout << startGenes[i].elem << " ";
    }
    std::cout << " ]" << Variables[0] << std::endl;

    //1st element of the Population should be the input individual
    indZero.genes = startGenes;
    pop.inds.push_back(indZero);
    //Cycle for sequential creating individuals
    for (int i = 1; i < numInd; i++)
    {
        int dist = rand() % 7 + 1;
        indZero.genes.clear();
        indZero.genes = startGenes;
        int varChoice = rand() % (Variables.size() - 1) + 1; //variable for addition

        //Creating individuals with the operands decided by the probability distribution
        //(weights can be adjusted according to the problem)

        if (dist == 1)
        {
            genZero.elem = Variables[varChoice]; //setting the 2nd gene as the operation inititated - "+ x"
            genZero.oper = 1;
            indZero.genes.push_back(genZero);
            pop.inds.push_back(indZero); //Putting new individual in the vector of the Population
        }

        if (dist == 2)
        {
            genZero.elem = Variables[varChoice];
            genZero.oper = 3;
            indZero.genes.push_back(genZero);
            pop.inds.push_back(indZero);
        }

        if (dist == 3)
        {
            genZero.elem = Variables[varChoice];
            genZero.oper = 4;
            indZero.genes.push_back(genZero);
            pop.inds.push_back(indZero);
        }

        if (dist == 4)
        {
            genZero.elem = Variables[varChoice];
            genZero.oper = 5;
            indZero.genes.push_back(genZero);
            pop.inds.push_back(indZero);
        }

        if (dist == 5)
        {
            coeff = ((double(rand() % 12) + 1.0) / (double(rand() % 10 + 1.0))) - 0.001;
            genZero.elem = coeff;
            genZero.oper = 5;
            indZero.genes.push_back(genZero);
            pop.inds.push_back(indZero);
        }
            
        if (dist == 6)
        {
            coeff = ((double(rand() % 12) + 1.0) / (double(rand() % 10 + 1.0))) - 0.001;
            genZero.elem = coeff;
            genZero.oper = 3;
            indZero.genes.push_back(genZero);
            pop.inds.push_back(indZero);
        }

        if (dist == 7)
        {
            coeff = ((double(rand() % 12) + 1.0) / (double(rand() % 10 + 1.0))) - 0.001;
            int coeffBool = rand() % 2;
            if (coeffBool == 0)
                coeff = -coeff;
            genZero.elem = coeff;
            genZero.oper = 1;
            indZero.genes.push_back(genZero);
            pop.inds.push_back(indZero);
        }
    }

    //Outputting all of the population, checking if there's zeros
    std::cout << "INITIAL POPULATION " << 1 << " = { ";
    fout << "\nINITIAL POPULATION " << 1 << " = { ";

    for (int i = 0; i < numInd; i++)
    {
        pop.inds[i] = IndFromGenes(pop.inds[i].genes);
        pop.inds[i].genes = InputGeneDecomposition(pop.inds[i].ind);
        pop.inds[i] = SearchForAbscentVars(pop.inds[i], Variables);
   
        fout << pop.inds[i].ind;
        std::cout << pop.inds[i].ind;

        if (i != (numInd - 1))
        {
            fout << " | ";
            std::cout << " | ";
        }
    }
    fout << " }; \n";
    std::cout << " }; \n";

    fout.close();
    return pop;
}

#endif // LIB_GENETIC_FUNCTIONS_H_
