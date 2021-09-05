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
	void CalcFit(vector<struct data> ExpData, Symbolic t, std::string foutname) {

		struct data buf;
		double velocity = 0.0; //ind current v
		double dev = 0.0; //Devitation ind current v from expdata v
		double devSUM = 0.0; //Sum devitation

		std::ofstream fout;
		fout.open(foutname, std::ios_base::app);
		
		for (int i = 0; i < ExpData.size(); i++)
		{
			buf = ExpData[i];
			velocity = ind[t == buf.t];
			dev = pow((velocity - buf.v), 2);
			devSUM += dev;
		}

		fit = (1 / (1 + sqrt(devSUM))) * 100;
		return;
	}
};

//Checking if input y is already a complex formula
vector<Gene> GeneDecomposition(Symbolic y)
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
		std::cout << "It's a sum " << (*((Sum*)(bufY))).summands << std::endl;
		bufList = (*((Sum*)(bufY))).summands;
		auto iter = bufList.begin();
		int oper = 1;

		for (int i = 0; i < bufList.size(); i++)
		{
			auto bufNum = Symbolic(*iter)->clone();
			if ((typeid(*bufNum) == typeid(Number<double>)) || (typeid(*bufNum) == typeid(Number<int>)))
			{
				std::cout << *iter << " is a num" << std::endl;
				outputGene.elem = *iter;
				outputGene.oper = oper;
				bufGenes.push_back(outputGene);
			}
			if (typeid(*bufNum) == typeid(Symbol))
			{
				std::cout << *iter << " is a symb" << std::endl;
				outputGene.elem = *iter;
				outputGene.oper = oper;
				bufGenes.push_back(outputGene);
			}
			else if ((typeid(*bufNum) != typeid(Number<double>)) && (typeid(*bufNum) != typeid(Number<int>))
				&& (typeid(*bufNum) != typeid(Symbol)))
			{
				std::cout << *iter << " is a complex function" << std::endl;
				newGenes.clear();
				newGenes = GeneDecomposition(*iter);
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
		std::cout << "It's a profuct " << (*((Product*)(bufY))).factors << std::endl;
		bufList = (*((Product*)(bufY))).factors;
		auto iter = bufList.begin();
		int oper = 3;

		for (int i = 0; i < bufList.size(); i++)
		{
			auto bufNum = Symbolic(*iter)->clone();
			if ((typeid(*bufNum) == typeid(Number<double>)) || (typeid(*bufNum) == typeid(Number<int>)))
			{
				std::cout << *iter << " is a num" << std::endl;
				outputGene.elem = *iter;
				outputGene.oper = oper;
				bufGenes.push_back(outputGene);
			}
			if (typeid(*bufNum) == typeid(Symbol))
			{
				std::cout << *iter << " is a symb" << std::endl;
				outputGene.elem = *iter;
				outputGene.oper = oper;
				bufGenes.push_back(outputGene);
			}
			else if ((typeid(*bufNum) != typeid(Number<double>)) && (typeid(*bufNum) != typeid(Number<int>))
				&& (typeid(*bufNum) != typeid(Symbol)))
			{
				std::cout << *iter << " is a complex function" << std::endl;
				newGenes.clear();
				newGenes = GeneDecomposition(*iter);
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
		std::cout << "It's a power " << (*((Symbol*)(&(*((Power*)(bufY)))))).parameters << std::endl;
		bufList = (*((Symbol*)(&(*((Power*)(bufY)))))).parameters;
		auto iter = bufList.begin();
		int oper = 5;

		for (int i = 0; i < bufList.size(); i++)
		{
			auto bufNum = Symbolic(*iter)->clone();
			if ((typeid(*bufNum) == typeid(Number<double>)) || (typeid(*bufNum) == typeid(Number<int>)))
			{
				std::cout << *iter << " is a num" << std::endl;
				outputGene.elem = *iter;
				outputGene.oper = oper;
				bufGenes.push_back(outputGene);
			}
			if (typeid(*bufNum) == typeid(Symbol))
			{
				std::cout << *iter << " is a symb" << std::endl;
				outputGene.elem = *iter;
				outputGene.oper = oper;
				bufGenes.push_back(outputGene);
			}
			else if ((typeid(*bufNum) != typeid(Number<double>)) && (typeid(*bufNum) != typeid(Number<int>))
				&& (typeid(*bufNum) != typeid(Symbol)))
			{
				std::cout << *iter << " is a complex function" << std::endl;
				newGenes.clear();
				newGenes = GeneDecomposition(*iter);
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
		std::ofstream fout;
		fout.open(foutname, std::ios_base::app);
		vector<Gene> startGenes;
		startGenes = GeneDecomposition(y);
		std::cout << "Start genes are [ ";
		for (int k = 0; k < startGenes.size(); k++)
		{
			std::cout << startGenes[k].elem << " ";
		}
		std::cout << "]" << std::endl;
		Individ sampleind;
		sampleind = IndFromGenes(startGenes);
		std::cout << "An ind from them is " << sampleind.ind << "\n";


		//Cycle for sequential creating individuals
		for (int i = 0; i < numInd; i++)
		{
			int dist = rand() % 9 + 1;
			indZero.genes.clear();
			indZero.genes = startGenes;

			//Creating individuals with the operands decided by the probability distribution
			//(weights can be adjusted according to the problem)

			if (dist == 1)
			{
				z = y + x;

				indZero.ind = z; //initialization of y as an individual
				indZero.pop = 1; //Ind population - 1st

				genZero.elem = x; //setting the 2nd gene as the operation inititated - "+ x"
				genZero.oper = 1;
				indZero.genes.push_back(genZero);
				inds.push_back(indZero); //Putting new individual in the vector of the Population

			}

			if (dist == 2)
			{
				z = y * x;

				indZero.ind = z;
				indZero.pop = 1;

				genZero.elem = x;
				genZero.oper = 3;
				indZero.genes.push_back(genZero);
				inds.push_back(indZero);

			}

			if (dist == 3)
			{
				z = y / x;

				indZero.ind = z;
				indZero.pop = 1;

				genZero.elem = x;
				genZero.oper = 4;
				indZero.genes.push_back(genZero);
				inds.push_back(indZero);

			}

			if (dist == 4)
			{
				z = pow(y, x);

				indZero.ind = z;
				indZero.pop = 1;

				genZero.elem = x;
				genZero.oper = 5;
				indZero.genes.push_back(genZero);
				inds.push_back(indZero);

			}

			if (dist == 5)
			{
				coeff = double(rand() % 100) + 1.0;
				z = y ^ coeff;

				indZero.ind = z;
				indZero.pop = 1;

				genZero.elem = coeff;
				genZero.oper = 5;
				indZero.genes.push_back(genZero);
				inds.push_back(indZero);

			}
			
			if (dist == 6)
			{
				coeff = double(rand() % 100) + 1.0;
				z = y * coeff;

				indZero.ind = z;
				indZero.pop = 1;

				genZero.elem = coeff;
				genZero.oper = 3;
				indZero.genes.push_back(genZero);
				inds.push_back(indZero);

			}

			if (dist == 7)
			{
				coeff = double(rand() % 100) + 1.0;
				z = y + coeff;

				indZero.ind = z;
				indZero.pop = 1;

				genZero.elem = coeff;
				genZero.oper = 1;
				indZero.genes.push_back(genZero);
				inds.push_back(indZero);

			}

			if (dist == 8)
			{
				coeff = double(rand() % 100) + 1.0;
				z = y - coeff;

				indZero.ind = z;
				indZero.pop = 1;

				genZero.elem = coeff;
				genZero.oper = 2;
				indZero.genes.push_back(genZero);
				inds.push_back(indZero);

			}

			if (dist == 9)
			{
				coeff = double(rand() % 100) + 1.0;
				z = y / coeff;

				indZero.ind = z;
				indZero.pop = 1;

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
			outputInd = inds[i];
			//If one of the inds in the 1st population is zero
			auto buf = outputInd.ind->clone();
			if ((typeid(*buf) == typeid(Number<double>)) || (typeid(*buf) == typeid(Number<int>)))
			{
				outputInd.ind = y;
				outputGene.elem = y;
				outputGene.oper = 1;
				outputInd.genes.clear();
				outputInd.genes.push_back(outputGene);
				inds[i] = outputInd;
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

	std::cout << "Datafile is " << totalLen << " lines long, input frequiency will be "
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
	datafile.close();
	datain.close();
	return(ExpData);
};

Individ Mutation(Individ KID)
{
	Individ newKID;
	newKID = KID;

	int m = 100;
	double arg_st = 0.5;

	for (int i = 0; i < newKID.genes.size(); i++)
	{
		auto buf = newKID.genes[i].elem->clone();

		if ((typeid(*buf) == typeid(Number<double>)) || (typeid(*buf) == typeid(Number<int>)))
		{
			double nmax = newKID.genes[i].elem - 0.5 * newKID.genes[i].elem;
			double nmin = newKID.genes[i].elem + 0.5 * newKID.genes[i].elem;
			double dx = fabs(nmax - nmin);
			double m_ver = 0.0;

			for (int j = 0; j < m; j++)
			{
				double randNum = (rand() % 100 + 1) / double(100);

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

//Function for GA coefficient optimization
Individ numGA(Individ inputInd, vector<struct data> ExpData, Symbolic t, std::string foutname)
{
	Individ outputInd;
	vector<Individ> numPop; //Population for numeric GA
	int GAsize = 20; //number of inds in GA
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
		double fitAVG = 0.0;
		
		for (int i = 0; i < numPop.size(); i++)
		{
			numPop[i] = IndFromGenes(numPop[i].genes);
			fout << "New ind is " << numPop[i].ind << std::endl;
		}

		//MAIN GA LOOP

		int limit = 100; //manual limit for GA loops
		int mind, maxd;
		double coef1, coef2;

		for (int f = 0; f < limit; f++) {

			mind = 0;
			maxd = 0;
			numPop[mind].CalcFit(ExpData, t, foutname); //for the min calc later
			fitAVG = 0.0;

			//Displaying current population
			fout << "\nNumeric GA Population " << f + 1 << std::endl;
			for (int g = 0; g < numPop.size(); g++)
			{
				numPop[g].CalcFit(ExpData, t, foutname);
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
			fout << "The maximum fit in NumPopulation " << f + 1<< " is " << numPop[maxd].fit << std::endl;

			//Checking if the current population is converged
			if (abs(1 - (numPop[maxd].fit / numPop[mind].fit)) < 0.01)
			{
				fout << "\nOptimum coefficients found after " << f + 1 << " loops, final ind is "
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
			//coef2 = (rand() % 100) / double(100);
			Individ KID = MOM;
			
			for (int q = 0; q < MOM.genes.size(); q++)
			{
				auto buf = KID.genes[q].elem->clone();
				if (typeid(*buf) == typeid(Number<double>))
				{
					KID.genes[q].elem = coef1 * MOM.genes[q].elem + (1 - coef1) * DAD.genes[q].elem;
				}
				buf->unreference(buf);
			}
			
			KID = IndFromGenes(KID.genes);
			KID.CalcFit(ExpData, t, foutname);
			
			fout << "num KID is " << KID.ind << " and fit is " << KID.fit << std::endl;

			double randProb = 0.0;
			randProb = (rand() % 100) / double(100);
			if (randProb < 0.9)
			{
				KID = Mutation(KID);
			}

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

//Individ OnePointCrossover(Population popul, Individ MOM, Individ DAD, int nodeCross)
//{
//	Individ KID;
//	vector<Gene> newKID;
//	newKID.clear();
//
//	vector<int> nums;
//	nums.clear();
//	nums.reserve(3);
//	if ((MOM.genes.size() - nodeCross) > DAD.genes.size())
//	{
//		nums = MomDadChoice(popul);
//
//		if (DAD.genes.size() < nodeCross)
//		{
//			break;
//		}
//	}
//	else
//	{
//		nums.push_back(numMOM);
//		nums.push_back(numDAD);
//		nums.push_back(nodeCross);
//		return(nums);
//	}
//
//	for (int i = 0; i < MOM.genes.size(); i++)
//	{
//		if (i < nodeCross)
//		{
//			newKID.push_back(MOM.genes[i]);
//		}
//		
//		if (i >= nodeCross)
//		{
//			newKID.push_back(DAD.genes[i]);
//		}
//	}
//
//	KID = IndFromGenes(newKID);
//	return(KID);
//}

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

//Function for GA symbolic optimization (main)
Individ symbGA(Population popul, vector<struct data> ExpData, Symbolic t, unsigned int startime,
	std::string foutname, std::string fitfilename)
{
	Individ outputInd;
	double stopPoint = 99.9;
	
	//MAIN GA LOOP
	int limit = 1000; //manual limit for GA loops
	int mind, maxd;
	double coef1, coef2;
	double fitAVG = 0.0;

	std::ofstream fout;
	fout.open(foutname, std::ios_base::app);

	std::ofstream fitfile;
	fitfile.open(fitfilename, std::ios_base::app);

	std::setprecision(9);

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
		fitAVG = 0.0;
		popul.inds[0].CalcFit(ExpData, t, foutname);

		//Displaying current population
		std::cout << "\nsymbGA Population " << f + 1 << std::endl;
		for (int g = 0; g < popul.inds.size(); g++)
		{
			popul.inds[g].CalcFit(ExpData, t, foutname);
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

		auto bufMOM = MOM.ind->clone();
		auto bufDAD = DAD.ind->clone();



		if (abs(1 - (MOM.fit / DAD.fit)) < 0.05)
		{
			double boolCross = (1 / (double(rand() % 10) + 1)); //Particular crossover probability

			if ((boolCross >= 0.2) && (boolCross <= 0.45))
			{
				fout << "RandPlusCrossover" << std::endl;
				KID = RandPlusCrossover(MOM, DAD);
			}
			if ((boolCross > 0.15) && (boolCross < 0.2))
			{
				fout << "StakingCrossover" << std::endl;
				KID = StakingCrossover(MOM, DAD, nodeCross);
			}
			else
			{
				fout << "StakingCrossover" << std::endl;
				KID = MixCrossover(MOM, DAD);
			}
		}
		else
		{
			fout << "OnePointCrossover" << std::endl;
			KID = MixCrossover(MOM, DAD);
		}

		bufMOM->unreference(bufMOM);
		bufDAD->unreference(bufDAD);

		fout << "KID is " << KID.ind << std::endl;
		KID.CalcFit(ExpData, t, foutname);
		KID = numGA(KID, ExpData, t, foutname);
		KID.pop = f + 2;
		fout << "Optimimzed KID is " << KID.ind << " and fit is " << KID.fit << std::endl;
		//Replacing the worst element of the population with the KID
		popul.inds[mind] = KID;
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
	v = t+50; //v depending on t
	int k = 1; //Individual serial number
	int numInd = 5; //number of individuals in the population
	int numCoef = 0; //number of coefficients in the origin individual
	int len = 200; //number of lines in ExpData to read
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
		popul.inds.at(i).CalcFit(ExpData, t, foutname); //t - the agrument which should be substituted
		popul.inds.at(i) = numGA(popul.inds.at(i), ExpData, t, foutname);
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