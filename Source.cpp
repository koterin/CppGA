#define _SCL_SECURE_NO_WARNINGS
#include <iostream>
#include <string>
#include <fstream> //for files
#include "symbolicc++.h" //for symbolic
#include <ctime> //for runtime calculations

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
	void CreatePop(Symbolic y, Symbolic x, int numInd, std::string foutname)
	{
		//srand(time(0)); //turning on the random distribution
		Symbolic z = y; //temporary variable for creating individuals
		double coeff = 0.0; //random coefficient for 1st pop dreation
		inds.clear();
		inds.reserve(numInd);
		indZero.genes.clear();
		indZero.genes.reserve(2);
		std::ofstream fout;
		fout.open(foutname, std::ios_base::app);

		//Cycle for sequential creating individuals
		for (int i = 0; i < numInd; i++)
		{
			int dist = rand() % 10 +1;

			//Creating individuals with the operands decided by the probability distribution
			//(weights can be adjusted according to the problem)

			if (dist == 1)
			{
				z = y + x;

				indZero.ind = z; //initialization of y as an individual
				indZero.pop = 1; //Ind population - 1st

				genZero.elem = y; //setting the 1st gene as the imput y
				genZero.oper = 1;
				indZero.genes.clear();
				indZero.genes.push_back(genZero);

				genZero.elem = x; //setting the 2nd gene as the operation inititated - "+ x"
				genZero.oper = 1;
				indZero.genes.push_back(genZero);
				inds.push_back(indZero); //Putting new individual in the vector of the Population

			}

			if (dist == 2)
			{
				z = y - x;

				indZero.ind = z;
				indZero.pop = 1;

				genZero.elem = y;
				genZero.oper = 1;
				indZero.genes.clear();
				indZero.genes.push_back(genZero);

				genZero.elem = x;
				genZero.oper = 2;
				indZero.genes.push_back(genZero);
				inds.push_back(indZero);

			}

			if (dist == 3)
			{
				z = y * x;

				indZero.ind = z;
				indZero.pop = 1;

				genZero.elem = y;
				genZero.oper = 1;
				indZero.genes.clear();
				indZero.genes.push_back(genZero);

				genZero.elem = x;
				genZero.oper = 3;
				indZero.genes.push_back(genZero);
				inds.push_back(indZero);

			}

			if (dist == 4)
			{
				z = y / x;

				indZero.ind = z;
				indZero.pop = 1;

				genZero.elem = y;
				genZero.oper = 1;
				indZero.genes.clear();
				indZero.genes.push_back(genZero);

				genZero.elem = x;
				genZero.oper = 4;
				indZero.genes.push_back(genZero);
				inds.push_back(indZero);

			}

			if (dist == 5)
			{
				z = pow(y, x);

				indZero.ind = z;
				indZero.pop = 1;

				genZero.elem = y;
				genZero.oper = 1;
				indZero.genes.clear();
				indZero.genes.push_back(genZero);

				genZero.elem = x;
				genZero.oper = 5;
				indZero.genes.push_back(genZero);
				inds.push_back(indZero);

			}

			if (dist == 6)
			{
				coeff = double(rand() % 100) + 1.0;
				z = y ^ coeff;

				indZero.ind = z;
				indZero.pop = 1;

				genZero.elem = y;
				genZero.oper = 1;
				indZero.genes.clear();
				indZero.genes.push_back(genZero);

				genZero.elem = coeff;
				genZero.oper = 5;
				indZero.genes.push_back(genZero);
				inds.push_back(indZero);

			}
			
			if (dist == 7)
			{
				coeff = double(rand() % 100) + 1.0;
				z = y * coeff;

				indZero.ind = z;
				indZero.pop = 1;

				genZero.elem = y;
				genZero.oper = 1;
				indZero.genes.clear();
				indZero.genes.push_back(genZero);

				genZero.elem = coeff;
				genZero.oper = 3;
				indZero.genes.push_back(genZero);
				inds.push_back(indZero);

			}

			if (dist == 8)
			{
				coeff = double(rand() % 100) + 1.0;
				z = y + coeff;

				indZero.ind = z;
				indZero.pop = 1;

				genZero.elem = y;
				genZero.oper = 1;
				indZero.genes.clear();
				indZero.genes.push_back(genZero);

				genZero.elem = coeff;
				genZero.oper = 1;
				indZero.genes.push_back(genZero);
				inds.push_back(indZero);

			}

			if (dist == 9)
			{
				coeff = double(rand() % 100) + 1.0;
				z = y - coeff;

				indZero.ind = z;
				indZero.pop = 1;

				genZero.elem = y;
				genZero.oper = 1;
				indZero.genes.clear();
				indZero.genes.push_back(genZero);

				genZero.elem = coeff;
				genZero.oper = 2;
				indZero.genes.push_back(genZero);
				inds.push_back(indZero);

			}

			if (dist == 10)
			{
				coeff = double(rand() % 100) + 1.0;
				z = y / coeff;

				indZero.ind = z;
				indZero.pop = 1;

				genZero.elem = y;
				genZero.oper = 1;
				indZero.genes.clear();
				indZero.genes.push_back(genZero);

				genZero.elem = coeff;
				genZero.oper = 4;
				indZero.genes.push_back(genZero);
				inds.push_back(indZero);;

			}
		}

		//Outputting all of the population, checking if there's zeros
		std::cout << "INITIAL POPULATION " << 1 << " = { ";
		fout << "\nINITIAL POPULATION " << 1 << " = { ";

		for (int i = 0; i < numInd; i++)
		{
			outputInd = inds.at(i);
			//If one of the inds in the 1st population is zero
			auto buf = outputInd.ind->clone();
			if ((typeid(*buf) == typeid(Number<double>)) || (typeid(*buf) == typeid(Number<int>)))
			{
				outputInd.ind = y;
				outputGene.elem = y;
				outputGene.oper = 1;
				outputInd.genes.clear();
				outputInd.genes.push_back(outputGene);
				inds.at(i) = outputInd;
			}
			buf->unreference(buf);

			fout << outputInd.ind;
			std::cout << outputInd.ind;

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

	std::ofstream datain;
	datain.open("Data\\expdata.txt");
	
	while (i < len)
	{
		if (datafile.eof())
		{
			datain << "\nDesired length if bigger than the input file, end of file on the line " << i << std::endl;
			break;
		}

		datafile >> curdata.t >> curdata.v;
		datain << curdata.t << "	" << curdata.v << std::endl;
		ExpData.push_back(curdata);
		i += 1;
	}

	ExpData.shrink_to_fit();
	datafile.close();
	datain.close();
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

	//For the 1st element
	outputInd.genes[0] = genes[0];
	outputInd.ind = genes[0].elem;

	for (int i = 1; i < vsize; i++)
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

//Function for GA coefficient optimization
Individ numGA(Individ inputInd, vector<struct data> ExpData, Symbolic t, std::string foutname)
{
	Individ outputInd;
	vector<Individ> numPop; //Population for numeric GA
	int GAsize = 10; //number of inds in GA
	double resCoef = 0.0;
	int bol = 0;
	int dec = 0; //decision - is there any genes to optimize?

	std::ofstream fout;
	fout.open(foutname, std::ios_base::app);

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

	for (int i = 0; i < inputInd.genes.size(); i++)
	{
		auto buf = inputInd.genes[i].elem->clone();
		
		if ((typeid(*buf) == typeid(Number<double>)) || (typeid(*buf) == typeid(Number<int>)))
		{
			dec += 1;
			//creating new individual
			for (int j = 0; j < numPop.size(); j++)
			{
				resCoef = numPop[j].genes[i].elem;
				bol = rand() % 2; //boolean imitation

				if (bol == 0)
				{
					resCoef = double(rand() % 100);
				}
				else if (bol == 1)
				{
					resCoef = double(-(rand() % 100));
				};

				numPop[j].genes[i].elem = resCoef;
			}
		}

		buf->unreference(buf);
	}

	if (dec == 0)
	{
		std::cout << "No numeric GA optimization needed for the " << outputInd.ind << std::endl;
		fout << "No numeric GA optimization needed for the " << outputInd.ind << std::endl;
		fout.close();
		return(inputInd);
	}

	else if (dec > 0)
	{
		fout << "The initial ind fit is " << inputInd.fit << std::endl;
		double multCoef = 2; //Condition: if multCoef * current_fit < initial_fit then GA stops

		for (int i = 0; i < numPop.size(); i++)
		{
			numPop[i] = IndFromGenes(numPop[i].genes);
			fout << "New ind is " << numPop[i].ind << std::endl;
		}

		//MAIN GA LOOP

		int limit = 50; //manual limit for GA loops
		int mind, maxd;
		double coef1, coef2;

		for (int f = 0; f < limit; f++) {

			mind = 0;
			maxd = 0;
			numPop[0].CalcFit(ExpData, t);

			//Displaying current population
			fout << "\nNumeric GA Population " << f + 1 << std::endl;
			for (int g = 0; g < numPop.size(); g++)
			{
				numPop[g].CalcFit(ExpData, t);
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
			}

			fout << "\nThe minimum fit in NumPopulation " << f + 1 << " is " << numPop[mind].fit << std::endl;
			fout << "The maximum fit in NumPopulation " << f + 1<< " is " << numPop[maxd].fit << std::endl;

			//Checking if the current minimum fit is twice smaller than the initial
			if (((numPop[maxd].fit) / multCoef) > inputInd.fit)
			{
				fout << "\nOptimum coefficients found after " << f + 1 << " loops, final ind is "
					<< numPop[maxd].ind << " and fit = " << numPop[maxd].fit << std::endl;
				outputInd = numPop[maxd];
				fout.close();
				return(outputInd);
			}

			//Setting Mom and Dad as 2 random elements of the population
			//MOM
			int numMOM = rand() % (GAsize - 1);
			Individ MOM = numPop[numMOM];
			fout << "num MOM is " << MOM.ind << " and fit is " << MOM.fit << std::endl;

			//DAD
			int numDAD = rand() % (GAsize - 1);
			if (numDAD == numMOM)
			{
				while (numDAD == numMOM)
				{
					numDAD = rand() % (GAsize - 1);
				}
			}

			Individ DAD = numPop[numDAD];
			fout << "num DAD is " << DAD.ind << " and fit is " << DAD.fit << std::endl;

			//KID
			coef1 = (rand() % 100) / double(100);
			coef2 = (rand() % 100) / double(100);
			Individ KID = MOM;
			
			for (int q = 0; q < MOM.genes.size(); q++)
			{
				auto buf = KID.genes[q].elem->clone();
				if (typeid(*buf) == typeid(Number<double>))
				{
					KID.genes[q].elem = coef1 * MOM.genes[q].elem + coef2 * DAD.genes[q].elem;
				}
				buf->unreference(buf);
			}
			
			KID = IndFromGenes(KID.genes);
			KID.CalcFit(ExpData, t);
			
			fout << "num KID is " << KID.ind << " and fit is " << KID.fit << std::endl;

			//Replacing the worst element of the population with the KID
			numPop[mind] = KID;
			outputInd = numPop[maxd];
			
		};
	}

	std::cout << "numGA optimization failed, ending the loop" << std::endl;
	fout << "numGA optimization failed, ending the loop" << std::endl;
	fout.close();
	return(outputInd);
}

vector<int> MomDadChoice(Population popul)
{
	vector<int> nums;
	nums.clear();
	nums.reserve(3);
	
	int numMOM = 0;
	//MOM
	for (int i = 0; i < 1000; i++) //1000 - random number of loops, just so it would be enough
	{
		numMOM = rand() % (popul.inds.size() - 1);
		if (popul.inds[numMOM].genes.size() != 1)
		{
			break;
		};
	}
	Individ MOM = popul.inds[numMOM];
	
	int nodeCross = (rand() % (MOM.genes.size() - 1)) + 1; //always excluding the 1st gene
	
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
	
	if (DAD.genes.size() < (double(nodeCross) + 1))
	{
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

Individ OnePointCrossover(Individ MOM, Individ DAD, int nodeCross)
{
	Individ KID;
	vector<Gene> newKID;
	newKID.clear();

	for (int i = 0; i < DAD.genes.size(); i++)
	{
		if (i < nodeCross)
		{
			newKID.push_back(MOM.genes[i]);
		}
		else if (i >= nodeCross)
		{
			newKID.push_back(DAD.genes[i]);
		}
	}

	KID = IndFromGenes(newKID);
	return(KID);
}

Individ StakingCrossover(Individ MOM, Individ DAD, int nodeCross)
{
	Individ KID;
	vector<Gene> newKID;
	newKID.clear();

	for (int i = 0; i < nodeCross; i++)
	{
		newKID.push_back(MOM.genes[i]);
	}

	for (int i = 0; i < DAD.genes.size(); i++)
	{
		newKID.push_back(DAD.genes[i]);
	}

	KID = IndFromGenes(newKID);
	return(KID);
}

//Function for GA symbolic optimization (main)
Individ symbGA(Population popul, vector<struct data> ExpData, Symbolic t, unsigned int startime, std::string foutname)
{
	Individ outputInd;
	double stopPoint = 99.9;
	
	//MAIN GA LOOP
	int limit = 1000; //manual limit for GA loops
	int mind, maxd;
	double coef1, coef2;

	std::ofstream fout;
	fout.open(foutname, std::ios_base::app);

	std::cout << "\n--------------------------THE SYMBOLIC GA STARTED--------------------------\n" << std::endl;
	fout << "\n--------------------------THE SYMBOLIC GA STARTED--------------------------\n" << std::endl;
	fout << "Stop point is " << stopPoint << "% accuracy, loop limit is " << limit << std::endl;
	
	for (int f = 0; f < limit; f++) {

		unsigned int nowtime = clock();
		double curtime = (nowtime - startime) / (double)CLOCKS_PER_SEC;
		std::cout << "\ncurrent runtime is " << curtime << " seconds" << std::endl;
		fout << "\ncurrent runtime is " << curtime << " seconds" << std::endl;
		
		mind = 0;
		maxd = 0;
		popul.inds[0].CalcFit(ExpData, t);

		//Displaying current population
		std::cout << "\nsymbGA Population " << f + 1 << std::endl;
		for (int g = 0; g < popul.inds.size(); g++)
		{
			popul.inds[g].CalcFit(ExpData, t);
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
		}

		std::cout << "\nThe minimum fit in Population " << f << " is " <<
			popul.inds[mind].ind << " with fit " << popul.inds[mind].fit << std::endl;
		std::cout << "The maximum fit in Population " << f << " is " << 
			popul.inds[maxd].ind << " with fit " << popul.inds[maxd].fit << std::endl;

		fout << "\nThe minimum fit in Population " << f << " is " <<
			popul.inds[mind].ind << " with fit " << popul.inds[mind].fit << std::endl;
		fout << "The maximum fit in Population " << f << " is " <<
			popul.inds[maxd].ind << " with fit " << popul.inds[maxd].fit << std::endl;

		//Checking if the current minimum fit is the desired one
		if (popul.inds[maxd].fit > stopPoint)
		{
			std::cout << "\n---------------FINAL---------------" <<
				"Optimum coefficients found after " << f << " loops, final ind is "
				<< popul.inds[maxd].ind << " and fit = " << popul.inds[maxd].fit << std::endl;
			fout << "\n---------------FINAL---------------" <<
				"Optimum coefficients found after " << f << " loops, final ind is "
				<< popul.inds[maxd].ind << " and fit = " << popul.inds[maxd].fit << std::endl;
			fout.close();
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
		fout << "cross node is " << nodeCross << std::endl;

		//KID
		Individ KID = MOM;
		//KID = OnePointCrossover(MOM, DAD, nodeCross);
		KID = StakingCrossover(MOM, DAD, nodeCross);
		std::cout << "KID is " << KID.ind << std::endl;
		KID.CalcFit(ExpData, t);
		KID = numGA(KID, ExpData, t, foutname);
		KID.pop = f + 2;
		fout << "Optimimzed KID is " << KID.ind << " and fit is " << KID.fit << std::endl;
		//Replacing the worst element of the population with the KID
		popul.inds[mind] = KID;
		outputInd = popul.inds[maxd];
	};

	std::cout << "GA optimization failed, ending the loop" << std::endl;
	fout << "GA optimization failed, ending the loop" << std::endl;
	fout.close();

	return(outputInd);
}

void main(void) {

	std::clock_t start;
	unsigned int startime = clock();
	
	Symbolic v("v"); //V - velocity
	Symbolic t("t"); //t - time
	v = t; //v depending on t
	int k = 1; //Individual serial number
	int numInd = 20; //number of individuals in the population
	int numCoef = 0; //number of coefficients in the origin individual
	int len = 1000; //number of lines in ExpData to read
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

	Population popul;
	popul.CreatePop(v, t, numInd, foutname); //Creating 1st population

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
		popul.inds.at(i).CalcFit(ExpData, t); //t - the agrument which should be substituted
		popul.inds.at(i) = numGA(popul.inds.at(i), ExpData, t, foutname);
	}

	outputInd = symbGA(popul, ExpData, t, startime, foutname);

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