#define _SCL_SECURE_NO_WARNINGS
#include <iostream>
#include <string>
#include <fstream> //for files
#include "symbolicc++.h" //for symbolic
#include <ctime> //for runtime calculations
#include <list>
#include <iterator>

class Gene {

public:
	Symbolic elem; //one symbolic gene
	int oper; //operator of the gene: 1 - plus, 2 - minus, 3 - mult, 4 - div, 5 - pow
};

struct data
{
	double v;
	double t;
};

class Individ {

public:
	Symbolic ind; //one symbolic individual
	int pop; //popluation where the individual was last modified or created
	vector<Gene> genes; //genes of the individual
	double fit; //Value of the fitness function for the individual

	//Function for calculation the fitness func for the individual
	void CalcFit(vector<struct data> ExpData, Symbolic t) {

		struct data buf;
		//double velocity = 0.0; //ind current v
		double dev = 0.0; //Devitation ind current v from expdata v
		double devSUM = 0.0; //Sum devitation
		Symbolic velocity;
		
		for (int i = 0; i < ExpData.size(); i++)
		{
			buf = ExpData[i];

			velocity.auto_expand = 0;
			velocity.simplified = 0;
			velocity = ind[t == buf.t];
			velocity.upr();

			auto yy = velocity->clone();
			//Если текущий yy - не число (например, комплексное число), то считаем efr как бесконечно большое число
			if ((typeid(*(yy)) != typeid(Number<double>)) && (typeid(*(yy)) != typeid(Number<int>)))
			{
				dev = 1e+20; //тоже можно заменить на efr +=
				devSUM += dev;
			}
			else
			{
				dev = pow((double(velocity) - buf.v), 2);
				devSUM += dev;
			}
			yy->unreference(yy);
			velocity.auto_expand = 1;
			velocity.simplified = 1;

		}
		
		//Checking if devSUM is too small so we won't get zeros in fits
		if (devSUM < 1e-10)
		{
			fit = 1e-20;
		}

		else
		{	
			fit = (1 / (1 + sqrt(devSUM))) * 100;
		}
		return;
	}
};

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
	outputInd.ind = "";
	outputInd.genes.resize(genes.size());

	//For the 1st element
	outputInd.genes[0] = genes[0];
	outputInd.ind = genes[0].elem;

	for (int i = 1; i < genes.size(); i++)
	{
		if (genes[i].oper == 1)
		{
			outputInd.ind = (outputInd.ind) + genes[i].elem;
			outputInd.genes[i] = genes[i];
		}

		if (genes[i].oper == 2)
		{
			outputInd.ind = (outputInd.ind) - genes[i].elem;
			outputInd.genes[i] = genes[i];
		}

		if (genes[i].oper == 3)
		{
			outputInd.ind = (outputInd.ind) * genes[i].elem;
			outputInd.genes[i] = genes[i];
		}

		if (genes[i].oper == 4)
		{
			outputInd.ind = (outputInd.ind) / genes[i].elem;
			outputInd.genes[i] = genes[i];
		}

		if (genes[i].oper == 5)
		{
			outputInd.ind = pow((outputInd.ind), genes[i].elem);
			outputInd.genes[i] = genes[i];
		}
	}
	return(outputInd);
}

Individ ZeroIndTermination(Individ IndZero, Symbolic t)
{
	Gene plusGene;
	plusGene.elem = t;
	plusGene.oper = 1;
	Individ bufInd;

	bufInd = IndZero;
	std::cout << "\nZeroIndTermination, bufInd is " << bufInd.ind << std::endl;
	std::cout << "previous genes length is " << bufInd.genes.size() << std::endl;

	if (bufInd.genes.size() == 0)
	{
		bufInd.genes.push_back(plusGene);
	}

	bufInd.genes.push_back(plusGene);
	std::cout << "Now genes size is " << bufInd.genes.size() << std::endl;
	
	for (int i = 0; i < bufInd.genes.size(); i++)
	{
		std::cout << "Gene " << i << " is " << bufInd.genes[i].elem << std::endl;
	}

	bufInd = IndFromGenes(bufInd.genes);
	std::cout << "ind after IndFromGenes " << bufInd.ind << std::endl;
	bufInd.genes = InputGeneDecomposition(bufInd.ind);

	return(bufInd);
}

Individ ZeroMultIndCheck(Individ IndZero, vector<Symbolic> Variables)
{
	vector<int> varCount(Variables.size() - 1, 0);
	
	//Looking for abscent variables
	for (int i = 0; i < IndZero.genes.size(); i++)
	{
		auto bufG = IndZero.genes[i].elem->clone();
		if ((typeid(bufG) != typeid(Number<int>)) || (typeid(bufG) != typeid(Number<double>)))
		{
			for (int j = 1; j < Variables.size(); j++)
			{
				if (IndZero.genes[i].elem == Variables[j])
				{
					varCount[j]++;
				}
			}
		}
		bufG->unreference(bufG);
	}

	for (int i = 1; i < varCount.size(); i++)
	{
		if (varCount[i] == 0)
		{
			IndZero = ZeroIndTermination(IndZero, Variables[i]);
		}
	}

	return(IndZero);
}

class Population {
private:
	Individ indZero;
	Gene genZero;
	Individ outputInd;
	Gene outputGene;
public:
	std::vector<Individ> inds; //all of the individuals in the population
	int pop; //number of the population

	//Function for creating the first population
	void CreatePop(vector<Symbolic> Variables, int numInd, std::string foutname)
	{
		//srand(time(0)); //turning on the random distribution
		double coeff = 0.0; //random coefficient for 1st pop dreation

		inds.clear();
		inds.reserve(numInd);

		std::ofstream fout;
		fout.open(foutname, std::ios_base::app);

		vector<Gene> startGenes;
		startGenes = InputGeneDecomposition(Variables[0]);
		std::cout << "Final input is [ ";
		for (int i = 0; i < startGenes.size(); i++)
		{
			std::cout << startGenes[i].elem << " ";
		}
		std::cout << " ]" << std::endl;

		//1st element of the Population should be the input individual
		indZero.genes = startGenes;
		inds.push_back(indZero);
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
				inds.push_back(indZero); //Putting new individual in the vector of the Population
			}

			if (dist == 2)
			{
				genZero.elem = Variables[varChoice];
				genZero.oper = 3;
				indZero.genes.push_back(genZero);
				inds.push_back(indZero);
			}

			if (dist == 3)
			{
				genZero.elem = Variables[varChoice];
				genZero.oper = 4;
				indZero.genes.push_back(genZero);
				inds.push_back(indZero);
			}

			if (dist == 4)
			{
				genZero.elem = Variables[varChoice];
				genZero.oper = 5;
				indZero.genes.push_back(genZero);
				inds.push_back(indZero);
			}

			if (dist == 5)
			{
				coeff = (double(rand() % 12) + 1.0) / (double(rand() % 10 + 1.0));
				genZero.elem = coeff;
				genZero.oper = 5;
				indZero.genes.push_back(genZero);
				inds.push_back(indZero);
			}
			
			if (dist == 6)
			{
				coeff = (double(rand() % 12) + 1.0) / (double(rand() % 10 + 1.0));
				genZero.elem = coeff;
				genZero.oper = 3;
				indZero.genes.push_back(genZero);
				inds.push_back(indZero);
			}

			if (dist == 7)
			{
				coeff = (double(rand() % 12) + 1.0) / (double(rand() % 10 + 1.0));
				int coeffBool = rand() % 2;
				if (coeffBool == 0)
				{
					coeff = -coeff;
				}
				genZero.elem = coeff;
				genZero.oper = 1;
				indZero.genes.push_back(genZero);
				inds.push_back(indZero);
			}
		}

		//Outputting all of the population, checking if there's zeros
		std::cout << "INITIAL POPULATION " << 1 << " = { ";
		fout << "\nINITIAL POPULATION " << 1 << " = { ";

		for (int i = 0; i < numInd; i++)
		{
			for (int j = 0; j < inds[i].genes.size(); j++)
			{
				std::cout << "\nCurrent gene is " << inds[i].genes[j].elem << " oper "
					<< inds[i].genes[j].oper << std::endl;

			}

			inds[i] = IndFromGenes(inds[i].genes);
			inds[i].genes = InputGeneDecomposition(inds[i].ind);
			inds[i] = ZeroMultIndCheck(inds[i], Variables);

			//If one of the inds in the 1st population is a number
			while (inds[i].genes.size() <= (Variables.size() - 1))
			{
				inds[i] = ZeroMultIndCheck(inds[i], Variables);
			}
			
			fout << inds[i].ind;
			std::cout << inds[i].ind;

			if (i != (numInd - 1))
			{
				fout << " | ";
				std::cout << " | ";
			}
		}
		fout << " }; \n";
		std::cout << " }; \n";

		fout.close();
		return;
	}

};

//Function for writing dataset from the file to the array
vector<struct data> SetData(std::string fileroute, int len)
{
	std::ifstream datafile(fileroute);
	vector<struct data> ExpData;
	ExpData.clear();
	ExpData.reserve(len);
	struct data curdata;
	int i = 0;
	double num1, num2;

	std::ofstream datain;
	datain.open("Data\\expdata.txt");
	std::string currentLine;
	int totalLen = 0;
	
	//Counting the number of lines in the input file
	while (!datafile.eof())
	{
		getline(datafile, currentLine);
		totalLen++;
	}

	int freq = (totalLen / len); //number of lines which will be repeatedly skipped
	if (freq < 1)
	{
		datain << "\nDesired length is bigger than the input file=" << std::endl;
		return(ExpData);
	}

	std::cout << "\nDatafile is " << totalLen << " lines long, input frequiency will be "
		<< freq << " lines" << std::endl;

	datafile.seekg(std::ios_base::beg);

	while (i < totalLen)
	{
		if (i % freq == 0)
		{
			datafile >> curdata.t >> curdata.v;
			datain << curdata.t << "	" << curdata.v << std::endl;
			ExpData.push_back(curdata);
		}

		else
		{
			datafile >> num1 >> num2;
		}

		i += 1;
	}

	ExpData.shrink_to_fit();

	double tMax, vMax;
	tMax = ExpData[0].t;
	vMax = ExpData[0].v;
	for (int i = 1; i < ExpData.size(); i++)
	{
		if (ExpData[i].t > tMax)
		{
			tMax = ExpData[i].t;
		}
		if (ExpData[i].v > vMax)
		{
			vMax = ExpData[i].v;
		}
	}

	std::ofstream expdataNORM;
	expdataNORM.open("Data\\expdataNORM.txt");

	for (int i = 0; i < ExpData.size(); i++)
	{
		ExpData[i].t = ExpData[i].t / tMax;
		ExpData[i].v = ExpData[i].v / vMax;
		expdataNORM << ExpData[i].t << "	" << ExpData[i].v << std::endl;
	}

	datafile.close();
	datain.close();
	expdataNORM.close();
	return(ExpData);
};

Individ NumMutation(Individ KID)
{
	Individ newKID;
	newKID = KID;

	int m = 10;
	double arg_st = 0.5;

	for (int i = 0; i < newKID.genes.size(); i++)
	{
		auto buf = newKID.genes[i].elem->clone();

		if ((typeid(*buf) == typeid(Number<double>)) || (typeid(*buf) == typeid(Number<int>)))
		{
			double nmax = newKID.genes[i].elem - 0.5 * newKID.genes[i].elem; //lower limit
			double nmin = newKID.genes[i].elem + 0.5 * newKID.genes[i].elem; //upper limit
			double dx = fabs(nmax - nmin); //limit line
			double m_ver = 0.0;

			for (int j = 0; j < m; j++)
			{
				double randNum = (double(rand() % 100) + 1) / double(100);

				double a = 0.0;
				if (randNum < (1 / double(m)))
				{
					a = 1.0;
					m_ver = m_ver + a * pow(arg_st, j + 1);
				}
			}

			int boolNum;
			boolNum = rand() % 2;
			if (boolNum == 0)
			{
				newKID.genes[i].elem += 0.5 * dx * m_ver;
			}
			else
			{
				newKID.genes[i].elem += - 0.5 * dx * m_ver;
			}
		}

		buf->unreference(buf);
	}

	return(newKID);
}

Individ SymbMutation(Individ KID, string foutname, Symbolic t)
{
	std::ofstream fout;
	fout.open(foutname, std::ios_base::app);
	int boolOper = rand() % 4 + 1;
	Gene mutGene;

	fout << "Mutation" << std::endl;
	fout << "\nTHE KID WAS " << KID.ind;

	if (boolOper == 1)
	{
		mutGene.elem = t;
		mutGene.oper = 3;
		fout << "\nmult mutation";
	}
	if (boolOper == 2)
	{
		mutGene.elem = t;
		mutGene.oper = 4;
		fout << "\ndiv mutation";
	}
	if (boolOper == 3)
	{
		mutGene.elem = t;
		mutGene.oper = 5;
		fout << "\npow mutation";
	}
	if (boolOper == 4)
	{
		mutGene.elem = double((rand() % 12 + 1) / (rand() % 10 + 1));
		mutGene.oper = 5;
		fout << "\nNum pow mutatuion ";
	}

	KID.genes.push_back(mutGene);
	KID = IndFromGenes(KID.genes);
	KID.genes = InputGeneDecomposition(KID.ind);

	fout.close();
	return(KID);
}

//Function for GA coefficient optimization
Individ numGA(Individ inputInd, vector<struct data> ExpData, Symbolic t, std::string foutname)
{
	Individ outputInd;
	vector<Individ> numPop; //Population for numeric GA
	int GAsize = 40; //number of inds in GA
	double resCoef = 0.0;
	int bol = 0;
	int dec = 0; //decision - is there any genes to optimize?

	std::ofstream fout;
	fout.open(foutname, std::ios_base::app);

	std::ofstream fitfile;
	fitfile.open("Data\\FitfileNUM.txt");

	numPop.clear();
	numPop.resize(GAsize);
	std::cout << "Numeric GA for V(t) = " << inputInd.ind << std::endl;
	fout << "------------NUMERIC GA STARTED-------------";
	fout << "\n Ind " << inputInd.ind << std::endl;

	outputInd.ind = inputInd.ind;
	outputInd.genes = inputInd.genes;
	outputInd.pop = inputInd.pop;

	for (int i = 0; i < numPop.size(); i++)
	{
		numPop[i] = outputInd;
	}

	//Looking for numeric genes in the Individual
	for (int i = 0; i < inputInd.genes.size(); i++)
	{
		auto buf = inputInd.genes[i].elem->clone();
		
		if ((typeid(*buf) == typeid(Number<double>)) || (typeid(*buf) == typeid(Number<int>)))
		{
			dec += 1;
			//creating new individual
			for (int j = 1; j < numPop.size(); j++)
			{
				bol = rand() % 2; //boolean imitation

				if (bol == 0)
				{
					resCoef = double(rand() % 10);
				}
				else if (bol == 1)
				{
					resCoef = double(-(rand() % 10));
				};

				numPop[j].genes[i].elem = resCoef;
			}
		}

		buf->unreference(buf);
	}

	//No numeric genes were found
	if (dec == 0)
	{
		std::cout << "No numeric GA optimization needed for the " << outputInd.ind << std::endl;
		fout << "No numeric GA optimization needed for the " << outputInd.ind << std::endl;
		fout.close();
		return(inputInd);
	}

	//Found numeric genes
	else if (dec > 0)
	{
		//Preparing the Population
		for (int i = 0; i < numPop.size(); i++)
		{
			numPop[i] = IndFromGenes(numPop[i].genes);
			fout << "New ind is " << numPop[i].ind << std::endl;
			numPop[i].CalcFit(ExpData, t);
			auto buf = numPop[i].ind->clone();
			if ((typeid(*buf) == typeid(Number<double>)) || (typeid(*buf) == typeid(Number<int>)))
			{
				numPop[i].fit = 0;
			}
			buf->unreference(buf);
		}

		fout << "The initial ind fit is " << numPop[0].fit << std::endl;
		double fitAVG = 0.0;

		//MAIN GA LOOP

		int limit = dec*1000; //manual limit for GA loops
		std::cout << "NumGA limits = " << limit << std::endl;
		int mind, maxd;
		double coef1, coef2;

		for (int f = 0; f < limit; f++) {

			mind = 0;
			maxd = 0;
			fitAVG = 0.0;

			//Displaying current population
			fout << "\nNumeric GA Population " << f + 1 << std::endl;
			for (int g = 0; g < numPop.size(); g++)
			{
				fout << numPop[g].ind << " and fit " << numPop[g].fit << std::endl;
				//Founding minimum
				if (numPop[g].fit < numPop[mind].fit)
				{
					mind = g;
				}

				//Founding maximum
				if (numPop[g].fit > numPop[maxd].fit)
				{
					maxd = g;
				}

				fitAVG += numPop[g].fit;
			}

			fitAVG = fitAVG / numPop.size();
			fitfile << f + 1 << " " << numPop[mind].fit << " " << fitAVG << " " << numPop[maxd].fit << std::endl;

			fout << "\nThe minimum fit in NumPopulation " << f + 1 << " is " << numPop[mind].fit << std::endl;
			fout << "The maximum fit in NumPopulation " << f + 1 << " is " << numPop[maxd].fit << std::endl;

			//Checking if the current population is converged
			if (abs(1 - (numPop[maxd].fit / numPop[mind].fit)) < 0.005)
			{
				fout << "\nNumPop was converged" <<
					"Optimum coefficients found after " << f + 1 << " loops, final ind is "
					<< numPop[maxd].ind << " and fit = " << numPop[maxd].fit << std::endl;
				outputInd = numPop[maxd];
				fout.close();
				fitfile.close();
				return(outputInd);
			}

			//Setting Mom and Dad as 2 random elements of the population
			//MOM
			int numMOM = rand() % (GAsize - 1);
			Individ MOM = numPop[numMOM];
			fout << "num MOM is " << MOM.ind << " and fit is " << MOM.fit << std::endl;

			//DAD
			int numDAD = rand() % (GAsize - 1);
			while (numDAD == numMOM)
			{
				numDAD = rand() % (GAsize - 1);
			}
			Individ DAD = numPop[numDAD];
			fout << "num DAD is " << DAD.ind << " and fit is " << DAD.fit << std::endl;

			//KID
			coef1 = (rand() % 100) / double(100);
			Individ KID = MOM;
			
			//Crossing over
			for (int q = 0; q < KID.genes.size(); q++)
			{
				auto buf = KID.genes[q].elem->clone();
				if ((typeid(*buf) == typeid(Number<double>)) || (typeid(*buf) == typeid(Number<int>)))
				{
					KID.genes[q].elem = coef1 * MOM.genes[q].elem + (1 - coef1) * DAD.genes[q].elem;
				}
				buf->unreference(buf);
			}
			
			//Numeric mutation
			double randProb = (rand() % 100) / double(100);
			if (randProb < 0.4)
			{
				KID = NumMutation(KID);
				fout << "NumMutation" << std::endl;
			}

			KID = IndFromGenes(KID.genes);

			//Checking if ind is converged to a number
			auto buf = KID.ind->clone();
			if ((typeid(*buf) == typeid(Number<double>)) || (typeid(*buf) == typeid(Number<int>)))
			{
				KID.fit = 0;
			}
			else
			{
				KID.CalcFit(ExpData, t);
			}
			buf->unreference(buf);

			fout << "num KID is " << KID.ind << " and fit is " << KID.fit << std::endl;

			//Replacing the worst element of the population with the KID
			numPop[mind] = KID;
			outputInd = numPop[maxd];
		};
	}

	std::cout << "numGA optimization failed, ending the loop" << std::endl;
	fout << "numGA optimization failed, ending the loop" << std::endl;
	fout.close();
	fitfile.close();
	return(outputInd);
}

vector<int> MomDadChoice(Population popul)
{
	vector<int> nums;
	nums.clear();
	nums.reserve(3);
	
	int numMOM = 0;

	//MOM
	while (popul.inds[numMOM].genes.size() == 1)
	{
		numMOM = rand() % popul.inds.size();
	}
	Individ MOM = popul.inds[numMOM];
	
	int nodeCross = (rand() % (MOM.genes.size() - 1)) + 1; //always excluding the 1st gene
	
	//DAD
	int numDAD = rand() % popul.inds.size();

	while (numDAD == numMOM) //Checking if DAD is the same as MOM
	{
		numDAD = rand() % popul.inds.size();
	}	
	Individ DAD = popul.inds[numDAD];
	
	nums.push_back(numMOM);
	nums.push_back(numDAD);
	nums.push_back(nodeCross);
	return(nums);

}

Individ StakingCrossover(Individ MOM, Individ DAD, int nodeCross)
{
	Individ KID;
	vector<Gene> newKID;
	newKID.clear();
	Gene bufGene;

	for (int i = 0; i < nodeCross; i++)
	{
		newKID.push_back(MOM.genes[i]);
	}

	bufGene = DAD.genes[0];
	int boolR = rand() % 5 + 1;

	switch (boolR)
	{
	case 1:
		bufGene.oper = 1;
	case 2:
		bufGene.oper = 2;
	case 3:
		bufGene.oper = 3;
	case 4:
		bufGene.oper = 4;
	case 5:
		bufGene.oper = 5;
	}

	newKID.push_back(bufGene);

	for (int i = 1; i < DAD.genes.size(); i++)
	{
		newKID.push_back(DAD.genes[i]);
	}

	KID = IndFromGenes(newKID);
	return(KID);
}

Individ RandPlusCrossover(Individ MOM, Individ DAD)
{
	Individ KID;
	vector<Gene> newKID;
	newKID.clear();
	int countMOM = MOM.genes.size();
	int countDAD = DAD.genes.size();

	int countLimit = countMOM;
	if (countMOM < countDAD)
	{
		countLimit = countDAD;
	}

	for (int i = 0; i < countLimit; i++)
	{
		if (i < countMOM)
		{
			newKID.push_back(MOM.genes[i]);
		}

		if (i < countDAD)
		{
			newKID.push_back(DAD.genes[(double(countDAD) - 1) - i]);
		}
	}

	KID = IndFromGenes(newKID);
	return(KID);
}

Individ MixCrossover(Individ MOM, Individ DAD)
{
	Individ KID;
	vector<Gene> newKID;
	newKID.clear();
	int countMOM = MOM.genes.size();
	int countDAD = DAD.genes.size();

	int countLimit = countMOM;
	if (countMOM < countDAD)
	{
		countLimit = countDAD;
	}

	for (int i = 0; i < countLimit; i++)
	{
		if ((i < countMOM) && (i % 2 == 0))
		{
			newKID.push_back(MOM.genes[i]);
		}

		if ((i >= countMOM) && (i % 2 == 0))
		{
			newKID.push_back(DAD.genes[i]);
		}

		if ((i < countDAD) && (i % 2 == 1))
		{
			newKID.push_back(DAD.genes[i]);
		}
		if ((i >= countDAD) && (i % 2 == 1))
		{
			newKID.push_back(MOM.genes[i]);
		}

	}

	KID = IndFromGenes(newKID);
	return(KID);
}

Individ ClassicCrossover(Individ MOM, Individ DAD, Symbolic t)
{
	Individ KID;
	vector<Gene> newKID;
	newKID.clear();

	int MOMpart = (rand() % (MOM.genes.size() - 1)) + 1; //size of MOM individual to be swapped
	int DADpart = (rand() % (DAD.genes.size() - 1)) + 1; //size of DAD individual to be swapped

	int MOMdiv = rand() % (MOM.genes.size() - MOMpart + 1); //Start gene for swap in MOM
	int DADdiv = rand() % (DAD.genes.size() - DADpart + 1); //Start gene for swap in DAD

	for (int i = 0; i < MOMdiv; i++)
	{
		newKID.push_back(MOM.genes[i]);
	}

	for (int i = 0; i < DADpart; i++)
	{
		newKID.push_back(DAD.genes[int(DADdiv + i)]);
	}

	if ((int(MOMdiv + MOMpart)) < MOM.genes.size())
	{
		for (int i = 0; i < (MOM.genes.size() - MOMdiv - MOMpart); i++)
		{
			newKID.push_back(MOM.genes[int(MOMdiv + MOMpart + i)]);
		}
	}

	KID = IndFromGenes(newKID);
	KID.genes = InputGeneDecomposition(KID.ind);

	while (KID.genes.size() <= 1)
	{
		KID = ZeroIndTermination(KID, t);
	}

	std::cout << "Final KID is " << KID.ind << std::endl;
	return(KID);
}

//Function for GA symbolic optimization (main)
Individ symbGA(Population popul, vector<struct data> ExpData, Symbolic t, unsigned int startime,
	std::string foutname, std::string fitfilename)
{
	Individ outputInd;
	double stopPoint = 99.9;
	
	//MAIN GA LOOP
	int mind, maxd;
	double coef1, coef2;
	double fitAVG = 0.0;
	Gene bufGene;
	bufGene.elem = t;
	bufGene.oper = 3;

	std::ofstream fout;
	fout.open(foutname, std::ios_base::app);

	std::ofstream fitfile;
	fitfile.open(fitfilename, std::ios_base::app);

	std::cout << "\n--------------------------THE SYMBOLIC GA STARTED--------------------------\n" << std::endl;
	fout << "\n--------------------------THE SYMBOLIC GA STARTED--------------------------\n" << std::endl;
	fout << "Stop point is " << stopPoint << "% accuracy" << std::endl;
	
	int f = 0;
	maxd = 0;
	while (popul.inds[maxd].fit < 100) {

		unsigned int nowtime = clock();
		double curtime = (nowtime - startime) / (double)CLOCKS_PER_SEC;
		std::cout << "\ncurrent runtime is " << curtime << " seconds" << std::endl;
		fout << "\ncurrent runtime is " << curtime << " seconds" << std::endl;
		
		mind = 0;
		maxd = 0;
		fitAVG = 0.0;

		//Displaying current population
		std::cout << "\nsymbGA Population " << f + 1 << std::endl;
		fout << "\nsymbGA Population " << f + 1 << std::endl;
		for (int g = 0; g < popul.inds.size(); g++)
		{
			fout << popul.inds[g].ind << " and fit " << popul.inds[g].fit << std::endl;
			//Founding minimum
			if (popul.inds[g].fit < popul.inds[mind].fit)
			{
				mind = g;
			}

			//Founding maximum
			if (popul.inds[g].fit > popul.inds[maxd].fit)
			{
				maxd = g;
			}

			fitAVG += popul.inds[g].fit;
		}

		fitAVG = fitAVG / popul.inds.size();
		fitfile << f + 1 << " " << popul.inds[mind].fit << " " << fitAVG << " " << popul.inds[maxd].fit << std::endl;

		std::cout << "\nThe minimum fit in Population " << f + 1 << " is " <<
			popul.inds[mind].ind << " with fit " << popul.inds[mind].fit << std::endl;
		std::cout << "The maximum fit in Population " << f + 1 << " is " << 
			popul.inds[maxd].ind << " with fit " << popul.inds[maxd].fit << std::endl;

		fout << "\nThe minimum fit in Population " << f + 1 << " is " <<
			popul.inds[mind].ind << " with fit " << popul.inds[mind].fit << std::endl;
		fout << "The maximum fit in Population " << f + 1 << " is " <<
			popul.inds[maxd].ind << " with fit " << popul.inds[maxd].fit << std::endl;

		//Checking if the current minimum fit is the desired one
		if (popul.inds[maxd].fit > stopPoint)
		{
			std::cout << "\n---------------FINAL---------------" <<
				"Optimum coefficients found after " << f + 1 << " loops, final ind is "
				<< popul.inds[maxd].ind << " and fit = " << popul.inds[maxd].fit << std::endl;
			fout << "\n---------------FINAL---------------" <<
				"Optimum coefficients found after " << f + 1 << " loops, final ind is "
				<< popul.inds[maxd].ind << " and fit = " << popul.inds[maxd].fit << std::endl;
			fout.close();
			fitfile.close();
			outputInd = popul.inds[maxd];
			return(outputInd);
		}

		//Setting Mom and Dad as 2 random elements of the population
		int numMOM = 0;
		int numDAD = 0;
		int nodeCross = 0;
		vector<int> nums;
		nums.resize(3);
		nums = MomDadChoice(popul);
		numMOM = nums[0];
		numDAD = nums[1];
		nodeCross = nums[2];

		Individ MOM = popul.inds[numMOM];
		Individ DAD = popul.inds[numDAD];
		fout << "MOM is " << MOM.ind << " and fit is " << MOM.fit << std::endl;
		fout << "DAD is " << DAD.ind << " and fit is " << DAD.fit << std::endl;
		//fout << "cross node is " << nodeCross << std::endl;

		//KID
		Individ KID = MOM; //Подумать над простой инициализацией класса

		fout << "ClassicCrossover" << std::endl;
		int numAmount = 11;
		while (numAmount > 10) //Limiting the amount of numerical coefficients up to 10
		{
			numAmount = 0;
			KID = ClassicCrossover(MOM, DAD, t);

			double boolCross = (1 / (double(rand() % 10) + 1)); //Mutation probability
			if (boolCross < 0.18)
			{
				KID = SymbMutation(KID, foutname, t);
			}

			for (int i = 0; i < KID.genes.size(); i++) //Numerical genes search
			{
				auto buf = KID.genes[i].elem->clone();
				if ((typeid(*buf) == typeid(Number<double>)) || (typeid(*buf) == typeid(Number<int>)))
				{
					numAmount += 1;
				}
				buf->unreference(buf);
			}
		}

		fout << "KID is " << KID.ind << std::endl;
		KID = numGA(KID, ExpData, t, foutname);
		KID.genes = InputGeneDecomposition(KID.ind);

		if (KID.genes.size() <= 1)
		{
			//KID.fit = 0;
			fout << "KID is 1 gene long, fit = 0" << std::endl;
			KID = ZeroIndTermination(KID, t);
			fout << "now ind is " << KID.ind << std::endl;
		}
		KID.CalcFit(ExpData, t);
		fout << "Optimimzed KID is " << KID.ind << " and fit is " << KID.fit << std::endl;
		//Replacing the worst element of the population with the KID

		popul.inds[mind] = KID;

		//Different logic: adding new ind to population only if it's better than the worst
		//if (KID.fit > popul.inds[mind].fit)
		//{
		//	popul.inds[mind] = KID;
		//}
		//else
		//{
		//	fout << "This KID's fit is worse than minimum, skipping that individual" << std::endl;
		//}

		f++;
		outputInd = popul.inds[maxd];
	};

	std::cout << "symbGA optimization failed, ending the loop" << std::endl;
	fout << "symbGA optimization failed, ending the loop" << std::endl;
	fout.close();
	fitfile.close();

	return(outputInd);
}

void main(void) {

	std::clock_t start;
	unsigned int startime = clock();
	std::fixed;
	
	Symbolic v("v"); //V - velocity
	Symbolic t("t"); //t - time
	Symbolic m("m"); //m - bullet mass
	double mass = 0.02; //Mass of the bullet
	//Area of the bullet; area of the fabric; impedance of the fabric; width of the layer; amount of layers

	v = t + m + 1; //starting ind must contain all of the variables
	int numInd = 20; //number of individuals in the population
	int len = 100; //number of lines in ExpData to read

	vector<Symbolic> Variables; //First variable must be the wanted one
	Variables.push_back(v);
	Variables.push_back(t);
	Variables.push_back(m);

	Individ outputInd; //Buffer for Individ class
	vector<struct data> ExpData; //Vector of experimental data
	ExpData.clear();
	ExpData.reserve(len);

	//Output file with all the logs
	std::ofstream fout;
	std::string foutname = "Data\\logs.txt";
	fout.open(foutname);
	fout << "Program started\n";
	fout.close();

	//Output file for Fitness function
	std::ofstream fitfile;
	std::string fitfilename = "Data\\fitfile.txt";
	fitfile.open(fitfilename);
	fitfile.close();

	Population popul;
	popul.CreatePop(Variables, numInd, foutname); //Creating 1st population

	//Insert here the path to the input data file
	//WARNING! All the phrases must be deleted from the file
	std::string inputroute = "Data\\1.txt";
	std::ifstream datafile(inputroute);

	fout.open(foutname, std::ios_base::app);

	if (datafile.is_open())
	{
		ExpData = SetData(inputroute, len);
		datafile.close();
	}
	else
	{
		std::cout << "\nInput data file not found" << std::endl;
		fout << "\nInput data file not found" << std::endl;
		std::cout << "The Program will be terminated\n";
		fout << "The Program will be terminated\n";
		return;
	}

	fout << "\nExperimental data file size is " << ExpData.size() << " lines" << std::endl;
	fout.close();

	//Fitness function calculations for the 1st gen
	for (int i = 0; i < numInd; i++)
	{
		popul.inds[i] = numGA(popul.inds[i], ExpData, t, foutname);
		while (popul.inds[i].genes.size() <= 1)
		{
			popul.inds[i] = ZeroIndTermination(popul.inds[i], t);
		}
		popul.inds[i].CalcFit(ExpData, t); //t - the agrument which should be substituted
	}

	outputInd = symbGA(popul, ExpData, t, startime, foutname, fitfilename);

	fout.open(foutname, std::ios_base::app);
	std::cout << "PROGRAM RESULT IS " << outputInd.ind << " with fit = " << outputInd.fit << std::endl;
	fout << "PROGRAM RESULT IS " << outputInd.ind << " with fit = " << outputInd.fit << std::endl;

	//Program runtime calculation
	unsigned int endtime = clock();
	double runtime = (endtime - startime) / (double)CLOCKS_PER_SEC;
	std::cout << "\nTotal runtime is " << runtime << " seconds" << std::endl;
	fout << "\n Total runtime is " << runtime << " seconds" << std::endl;

	fout.close();

	return;
}