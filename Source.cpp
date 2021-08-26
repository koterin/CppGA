#define _SCL_SECURE_NO_WARNINGS
#include <iostream>
#include <string>
#include <fstream> //for files
#include "symbolicc++.h" //for symbolic
#include <ctime> //for runtime calculations

class Gene {

public:
	Symbolic elem; //one symbolic gene
	int depth; //depth of the gene (number of operations past the intial
	int oper; //operator of the gene: 1 - plus, 2 - minus, 3 - mult, 4 - div, 5 - pow
	int type; //type of the gene: 1 - numeric, 2 - symbolic
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
	int numCoef; //Number of numeric coefficients in the individual

	//Function for calculation the fitness func for the individual
	void CalcFit(vector<struct data> ExpData, Symbolic t) {

		struct data buf;
		double res = 0.0; //ind current v
		double dev = 0.0; //Devitation ind current v from expdata v
		double devf = 0.0; //Sum devitation
		auto iter = ExpData.begin();

		while (iter != ExpData.end())
		{
			buf = *iter;
			res = ind[t == buf.t];
			dev = pow((res - buf.v),2);
			devf += dev;
			++iter;
		}

		fit = (1 / (1 + sqrt(devf))) * 100;
		return;
	}
};

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
	void CreatePop(Symbolic y, Symbolic x, int n, int nc)
	{
		//srand(time(0)); //turning on the random distribution
		Symbolic z = y; //temporary variable for creating individuals
		double coeff = 0.0; //random coefficient for 1st pop dreation
		inds.clear();
		inds.reserve(n);
		indZero.genes.clear();
		indZero.genes.reserve(2);

		//Cycle for sequential creating individuals
		for (int i = 0; i < n; i++)
		{
			int dist = rand() % 10 +1;

			//Creating individuals with the operands decided by the probability distribution
			//(weights can be adjusted according to the problem)

			if (dist == 1)
			{
				z = y + x;

				indZero.ind = z; //initialization of y as an individual
				indZero.pop = 1; //Ind population - 1st
				indZero.numCoef = 0; //Number of numeric coefficients in the individual

				genZero.elem = y; //setting the 1st gene as the imput y
				genZero.depth = 1;
				genZero.oper = 1;
				genZero.type = 2;
				indZero.genes.clear();
				indZero.genes.push_back(genZero);

				genZero.elem = x; //setting the 2nd gene as the operation inititated - "+ x"
				genZero.depth = 2;
				genZero.oper = 1;
				genZero.type = 2;
				indZero.genes.push_back(genZero);
				inds.push_back(indZero); //Putting new individual in the vector of the Population

			}

			if (dist == 2)
			{
				z = y - x;

				indZero.ind = z;
				indZero.pop = 1;
				indZero.numCoef = 0;

				genZero.elem = y;
				genZero.depth = 1;
				genZero.oper = 1;
				genZero.type = 2;
				indZero.genes.clear();
				indZero.genes.push_back(genZero);

				genZero.elem = x;
				genZero.depth = 2;
				genZero.oper = 2;
				genZero.type = 2;
				indZero.genes.push_back(genZero);
				inds.push_back(indZero);

			}

			if (dist == 3)
			{
				z = y * x;

				indZero.ind = z;
				indZero.pop = 1;
				indZero.numCoef = 0;

				genZero.elem = y;
				genZero.depth = 1;
				genZero.oper = 1;
				genZero.type = 2;
				indZero.genes.clear();
				indZero.genes.push_back(genZero);

				genZero.elem = x;
				genZero.depth = 2;
				genZero.oper = 3;
				genZero.type = 2;
				indZero.genes.push_back(genZero);
				inds.push_back(indZero);

			}

			if (dist == 4)
			{
				z = y / x;

				indZero.ind = z;
				indZero.pop = 1;
				indZero.numCoef = 0;

				genZero.elem = y;
				genZero.depth = 1;
				genZero.oper = 1;
				genZero.type = 2;
				indZero.genes.clear();
				indZero.genes.push_back(genZero);

				genZero.elem = x;
				genZero.depth = 2;
				genZero.oper = 4;
				genZero.type = 2;
				indZero.genes.push_back(genZero);
				inds.push_back(indZero);

			}

			if (dist == 5)
			{
				z = pow(y, x);

				indZero.ind = z;
				indZero.pop = 1;
				indZero.numCoef = 0;

				genZero.elem = y;
				genZero.depth = 1;
				genZero.oper = 1;
				genZero.type = 2;
				indZero.genes.clear();
				indZero.genes.push_back(genZero);

				genZero.elem = x;
				genZero.depth = 2;
				genZero.oper = 5;
				genZero.type = 2;
				indZero.genes.push_back(genZero);
				inds.push_back(indZero);

			}

			if (dist == 6)
			{
				coeff = double(rand() % 10) + 1.0;
				z = y ^ coeff;

				indZero.ind = z;
				indZero.pop = 1;
				indZero.numCoef = 1;

				genZero.elem = y;
				genZero.depth = 1;
				genZero.oper = 1;
				genZero.type = 2;
				indZero.genes.clear();
				indZero.genes.push_back(genZero);

				genZero.elem = coeff;
				genZero.depth = 2;
				genZero.oper = 5;
				genZero.type = 1;
				indZero.genes.push_back(genZero);
				inds.push_back(indZero);

			}
			
			if (dist == 7)
			{
				coeff = double(rand() % 10) + 1.0;
				z = y * coeff;

				indZero.ind = z;
				indZero.pop = 1;
				indZero.numCoef = 1;

				genZero.elem = y;
				genZero.depth = 1;
				genZero.oper = 1;
				genZero.type = 2;
				indZero.genes.clear();
				indZero.genes.push_back(genZero);

				genZero.elem = coeff;
				genZero.depth = 2;
				genZero.oper = 3;
				genZero.type = 1;
				indZero.genes.push_back(genZero);
				inds.push_back(indZero);

			}

			if (dist == 8)
			{
				coeff = double(rand() % 10) + 1.0;
				z = y + coeff;

				indZero.ind = z;
				indZero.pop = 1;
				indZero.numCoef = 1;

				genZero.elem = y;
				genZero.depth = 1;
				genZero.oper = 1;
				genZero.type = 2;
				indZero.genes.clear();
				indZero.genes.push_back(genZero);

				genZero.elem = coeff;
				genZero.depth = 2;
				genZero.oper = 1;
				genZero.type = 1;
				indZero.genes.push_back(genZero);
				inds.push_back(indZero);

			}

			if (dist == 9)
			{
				coeff = double(rand() % 10) + 1.0;
				z = y - coeff;

				indZero.ind = z;
				indZero.pop = 1;
				indZero.numCoef = 1;

				genZero.elem = y;
				genZero.depth = 1;
				genZero.oper = 1;
				genZero.type = 2;
				indZero.genes.clear();
				indZero.genes.push_back(genZero);

				genZero.elem = coeff;
				genZero.depth = 2;
				genZero.oper = 2;
				genZero.type = 1;
				indZero.genes.push_back(genZero);
				inds.push_back(indZero);

			}

			if (dist == 10)
			{
				coeff = double(rand() % 10) + 1.0;
				z = y / coeff;

				indZero.ind = z;
				indZero.pop = 1;
				indZero.numCoef = 1;

				genZero.elem = y;
				genZero.depth = 1;
				genZero.oper = 1;
				genZero.type = 2;
				indZero.genes.clear();
				indZero.genes.push_back(genZero);

				genZero.elem = coeff;
				genZero.depth = 2;
				genZero.oper = 4;
				genZero.type = 1;
				indZero.genes.push_back(genZero);
				inds.push_back(indZero);;

			}
		}

		//Outputting all of the population, checking if there's zeros
		std::cout << "Population " << 1 << " = { ";

		for (int i = 0; i < n; i++)
		{
			outputInd = inds.at(i);
			//If one of the inds in the 1st population is zero
			if (outputInd.ind == 0)
			{
				outputInd.ind = y;
				outputInd.numCoef = 0;
				outputGene.elem = y;
				outputGene.depth = 1;
				outputGene.oper = 1;
				genZero.type = 2;
				outputInd.genes.clear();
				outputInd.genes.push_back(outputGene);
				inds.at(i) = outputInd;
				//std::cout << outputInd.genes.at(0).elem << " Gene at the i ";
			}

			std::cout << outputInd.ind;

			if (i != (n - 1))
			{
				std::cout << " | ";
			}
		}
		std::cout << " }; \n";
		//std::cout << inds.size();

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
	
	while (i < len)
	{
		if (datafile.eof())
		{
			std::cout << "\nEnd of file on the line " << i << std::endl;
			break;
		}

		datafile >> curdata.t >> curdata.v;
		ExpData.push_back(curdata);
		//std::cout << "curent data is " << curdata.t << " | " << curdata.v << std::endl;
		i += 1;
	}

	ExpData.shrink_to_fit();
	datafile.close();
	return(ExpData);
};

//NOT WRITTEN YET
void OutputPlot()
{
	//creates file with current ExpData
	//and best function results from the calculated individual
	//MathCAD makes one plot with both of them to compare
}

//Function for creating new individual from the genes given
//WARNING! Input genes vector MUST BE SORTED and organized the way it should be in the output Individual
Individ IndFromGenes(vector<Gene> genes)
{
	Individ outputInd;
	outputInd.ind = "";
	int vsize = genes.size();
	outputInd.genes.resize(vsize);
	int depth = 1;

	//For the 1st element
	outputInd.genes.at(0) = genes.at(0);
	outputInd.genes.at(0).depth = depth;
	depth += 1;
	outputInd.ind = genes.at(0).elem;

	for (int i = 1; i < vsize; i++)
	{
		if (genes.at(i).oper == 1)
		{
			outputInd.ind = (outputInd.ind) + genes.at(i).elem;
			outputInd.genes.at(i) = genes.at(i);
			outputInd.genes.at(i).depth = depth;
			depth += 1;
		}

		if (genes.at(i).oper == 2)
		{
			outputInd.ind = (outputInd.ind) - genes.at(i).elem;
			outputInd.genes.at(i) = genes.at(i);
			outputInd.genes.at(i).depth = depth;
			depth += 1;
		}

		if (genes.at(i).oper == 3)
		{
			outputInd.ind = (outputInd.ind) * genes.at(i).elem;
			outputInd.genes.at(i) = genes.at(i);
			outputInd.genes.at(i).depth = depth;
			depth += 1;
		}

		if (genes.at(i).oper == 4)
		{
			outputInd.ind = (outputInd.ind) / genes.at(i).elem;
			outputInd.genes.at(i) = genes.at(i);
			outputInd.genes.at(i).depth = depth;
			depth += 1;
		}

		if (genes.at(i).oper == 5)
		{
			outputInd.ind = pow((outputInd.ind), genes.at(i).elem);
			outputInd.genes.at(i) = genes.at(i);
			outputInd.genes.at(i).depth = depth;
			depth += 1;
		}
	}

	return(outputInd);
}

//Function for GA coefficient optimization
Individ numGA(Individ inputInd, vector<struct data> ExpData, Symbolic t)
{
	Individ outputInd;
	vector<Individ> numPop; //Population for numeric GA
	int GAsize = 15; //number of inds in GA
	double resCoef = 0.0;
	int bol = 0;
	int dec = 0; //decision - is there any genes to optimize?

	numPop.clear();
	numPop.resize(GAsize);
	std::cout << "\n Ind " << inputInd.ind << std::endl;

	outputInd.ind = inputInd.ind;
	outputInd.genes = inputInd.genes;
	outputInd.numCoef = inputInd.numCoef;
	outputInd.pop = inputInd.pop;

	for (int i = 0; i < numPop.size(); i++)
	{
		numPop[i] = outputInd;
	}

	for (int i = 0; i < inputInd.genes.size(); i++)
	{
		if (inputInd.genes.at(i).type == 1)
		{
			dec += 1;
			//creating new individual
			for (int j = 0; j < numPop.size(); j++)
			{
				resCoef = numPop[j].genes[i].elem;
				bol = rand() % 2; //boolean imitation

				if (bol == 0)
				{
					//resCoef = resCoef / double(rand() % 40);
					resCoef = double(rand() % 100);
				}
				else if (bol == 1)
				{
					//resCoef = resCoef / double(-(rand() % 40));
					resCoef = double(-(rand() % 100));
				};

				numPop[j].genes[i].elem = resCoef;

				std::cout << "new gene is " << numPop[j].genes.at(i).elem << std::endl;
			}
		}
	}

	if (dec == 0)
	{
		std::cout << "No numeric GA optimization needed" << std::endl;
		return(outputInd);
	}

	else if (dec > 0)
	{
		std::cout << "The initial ind fit is " << inputInd.fit << std::endl;
		double multCoef = 2; //Condition: if multCoef * current_fit < initial_fit then GA stops

		for (int i = 0; i < numPop.size(); i++)
		{
			numPop[i] = IndFromGenes(numPop[i].genes);
			std::cout << "New ind is " << numPop[i].ind << std::endl;
		}

		//MAIN GA LOOP

		int limit = 2; //manual limit for GA loops
		int mind, maxd;
		double coef1, coef2;

		for (int f = 0; f < limit; f++) {

			mind = 0;
			maxd = 0;
			numPop[0].CalcFit(ExpData, t);

			//Displaying current population
			std::cout << "\nNumeric GA Population " << f + 1 << std::endl;
			for (int g = 0; g < numPop.size(); g++)
			{
				numPop[g].CalcFit(ExpData, t);
				std::cout << numPop[g].ind << " and fit " << numPop[g].fit << std::endl;
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
			}

			std::cout << "\nThe minimum fit in Population " << f << " is " << numPop[mind].fit << std::endl;
			std::cout << "The maximum fit in Population " << f << " is " << numPop[maxd].fit << std::endl;

			//Checking if the current minimum fit is twice smaller than the initial
			if (((numPop[maxd].fit) / multCoef) > inputInd.fit)
			{
				std::cout << "Optimum coefficients found after " << f << " loops, final ind is "
					<< numPop[maxd].ind << " and fit = " << numPop[maxd].fit << std::endl;
				outputInd = numPop[maxd];
				return(outputInd);
			}

			//Setting Mom and Dad as 2 random elements of the population
			//MOM
			int numMOM = rand() % (GAsize - 1);
			std::cout << "numMOM is " << numMOM << std::endl;
			Individ MOM = numPop[numMOM];
			std::cout << "MOM is " << MOM.ind << " and fit is " << MOM.fit << std::endl;

			//DAD
			int numDAD = rand() % (GAsize - 1);
			if (numDAD == numMOM)
			{
				while (numDAD == numMOM)
				{
					numDAD = rand() % (GAsize - 1);
				}
			}

			std::cout << "numDAD is " << numDAD << std::endl;
			Individ DAD = numPop[numDAD];
			std::cout << "DAD is " << DAD.ind << " and fit is " << DAD.fit << std::endl;

			//KID
			coef1 = (rand() % 100) / double(100);
			coef2 = (rand() % 100) / double(100);
			Individ KID = MOM;
			
			for (int q = 0; q < MOM.genes.size(); q++)
			{
				if (KID.genes[q].type == 1)
				{
					KID.genes[q].elem = coef1 * MOM.genes[q].elem + coef2 * DAD.genes[q].elem;
				}
			}
			
			KID = IndFromGenes(KID.genes);
			KID.CalcFit(ExpData, t);
			
			std::cout << "KID is " << KID.ind << " and fit is " << KID.fit << std::endl;

			//Replacing the worst element of the population with the KID
			numPop[mind] = KID;
			
		};
	}

	std::cout << "GA optimization failed, ending the loop" << std::endl;
	return(outputInd);
}

vector<int> MomDadChoice(Population popul)
{
	vector<int> nums;
	nums.clear();
	nums.reserve(3);
	
	//MOM
	int numMOM = rand() % (popul.inds.size() - 1);
	Individ MOM = popul.inds[numMOM];
	std::cout << "\nMOM is " << MOM.ind << " numMom is " << numMOM << std::endl;

	int nodeCross = (rand() % (MOM.genes.size() - 1)) + 1; //always excluding the 1st gene
	std::cout << "nodeMom is " << nodeCross << std::endl;

	//DAD
	int numDAD = rand() % (popul.inds.size() - 1);

	if (numDAD == numMOM)
	{
		while (numDAD == numMOM)
		{
			numDAD = rand() % (popul.inds.size() - 1);
		}
	}

	Individ DAD = popul.inds[numDAD];
	std::cout << "numDAD is " << numDAD << " DAD is " << DAD.ind << " size is " << DAD.genes.size() << std::endl;

	if (DAD.genes.size() < (double(nodeCross) + 1))
	{
		std::cout << "nodeDAD not found, again" << std::endl;
		nums = MomDadChoice(popul);

		if (DAD.genes.size() < (double(nodeCross) + 1))
		{
			return(nums);
		}
	}
	else
	{
		nums.push_back(numMOM);
		nums.push_back(numDAD);
		nums.push_back(nodeCross);
		return(nums);
	}
}

//Function for GA symbolic optimization (main)
Individ symbGA(Population popul, vector<struct data> ExpData, Symbolic t, unsigned int startime)
{
	Individ outputInd;
	double stopPoint = 99.9;
	
	//MAIN GA LOOP
	int limit = 2; //manual limit for GA loops
	int mind, maxd;
	double coef1, coef2;

	std::cout << "\n-------------THE SYMBOLIC GA STARTS-------------\n" << std::endl;
	
	for (int f = 0; f < limit; f++) {

		unsigned int nowtime = clock();
		double curtime = (nowtime - startime) / (double)CLOCKS_PER_SEC;
		std::cout << "\ncurrent runtime is " << curtime << " seconds" << std::endl;
		
		mind = 0;
		maxd = 0;
		popul.inds[0].CalcFit(ExpData, t);

		//Displaying current population
		std::cout << "\nGA Population " << f + 1 << std::endl;
		for (int g = 0; g < popul.inds.size(); g++)
		{
			popul.inds[g].CalcFit(ExpData, t);
			std::cout << popul.inds[g].ind << " and fit " << popul.inds[g].fit << std::endl;
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
		}

		std::cout << "\nThe minimum fit in Population " << f << " is " << popul.inds[mind].fit << std::endl;
		std::cout << "The maximum fit in Population " << f << " is " << popul.inds[maxd].fit << std::endl;

		//Checking if the current minimum fit is the desired one
		if (popul.inds[maxd].fit > stopPoint)
		{
			std::cout << "\n-----FINAL-----" <<
				"Optimum coefficients found after " << f << " loops, final ind is "
				<< popul.inds[maxd].ind << " and fit = " << popul.inds[maxd].fit << std::endl;
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

		std::cout << "numMom is " << numMOM << std::endl;
		Individ MOM = popul.inds[numMOM];
		std::cout << "nodeMom is " << nodeCross << std::endl;
		
		Individ DAD = popul.inds[numDAD];
		std::cout << "numDAD is " << numDAD << std::endl;

		std::cout << "MOM is " << MOM.ind << " and fit is " << MOM.fit << std::endl;
		std::cout << "DAD is " << DAD.ind << " and fit is " << DAD.fit << std::endl;

		////KID
		//coef1 = (rand() % 100) / double(100);
		//coef2 = (rand() % 100) / double(100);
		//Individ KID = MOM;

		//for (int q = 0; q < MOM.genes.size(); q++)
		//{
		//	if (KID.genes.at(q).type == 1)
		//	{
		//		KID.genes.at(q).elem = coef1 * MOM.genes.at(q).elem + coef2 * DAD.genes.at(q).elem;
		//	}
		//}

		//KID = IndFromGenes(KID.genes);
		//KID.CalcFit(ExpData, t);

		//std::cout << "KID is " << KID.ind << " and fit is " << KID.fit << std::endl;

		//Replacing the worst element of the population with the KID
		//numPop.at(mind) = KID;

	};

	//std::cout << "GA optimization failed, ending the loop" << std::endl;
	return(outputInd);
}

void main(void) {

	std::clock_t start;
	unsigned int startime = clock();
	
	Symbolic v("v"); //V - velocity
	Symbolic t("t"); //t - time
	v = t; //v depending on t
	int k = 1; //Individual serial number
	int numInd = 6; //number of individuals in the population
	int numCoef = 0; //number of coefficients in the origin individual
	int len = 1000; //number of lines in ExpData to read
	Individ outputInd; //Buffer for Individ class
	vector<struct data> ExpData; //Vector of experimental data
	ExpData.clear();
	ExpData.reserve(len);

	Population popul;
	popul.CreatePop(v, t, numInd, numCoef); //Creating 1st population

	//Insert here the path to the input data file
	//WARNING! All the phrases must be deleted from the file
	std::string fileroute = "InputData\\1.txt";
	std::ifstream datafile(fileroute);

	if (datafile.is_open())
	{
		ExpData = SetData(fileroute, len);
		datafile.close();
	}
	else
	{
		std::cout << "\nInput data file not found" << std::endl;
		std::cout << "The Program will be terminated\n";
		return;
	}

	std::cout << "\nExperimental data file size is " << ExpData.size() << std::endl;
	std::cout << "Number of individuals in the population is " << numInd << "\n" << std::endl;
	
	//Fitness function calculations for the 1st gen
	for (int i = 0; i < numInd; i++)
	{
		popul.inds.at(i).CalcFit(ExpData, t); //t - the agrument which should be substituted
		popul.inds.at(i) = numGA(popul.inds.at(i), ExpData, t);
	}

	outputInd = symbGA(popul, ExpData, t, startime);

	//Program runtime calculation
	unsigned int endtime = clock();
	double runtime = (endtime - startime) / (double)CLOCKS_PER_SEC;
	std::cout << "\nRuntime is " << runtime << " seconds" << std::endl;

	return;
}