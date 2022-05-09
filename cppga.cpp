#define _SCL_SECURE_NO_WARNINGS
#include <iostream>
#include <string>
#include <fstream> //for files
#include <ctime> //for runtime calculations
#include <list>
#include <cmath>
#include <iterator>
#include "lib/symbolicc++.h" //for symbolic
#include "lib/cppga.h"

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

Individ SymbMutation(Individ KID, string foutname, vector<Symbolic> Variables)
{
    std::ofstream fout;
    fout.open(foutname, std::ios_base::app);
    int boolOper = rand() % 4 + 1;
    Gene mutGene;

    fout << "Mutation" << std::endl;
    fout << "\nTHE KID WAS " << KID.ind;

    if (boolOper == 1)
    {
        int boolGene = rand() % (Variables.size() - 1) + 1;
        mutGene.elem = Variables[boolGene];
        mutGene.oper = 3;
        fout << "\nmult mutation";
    }
    if (boolOper == 2)
    {
        int boolGene = rand() % (Variables.size() - 1) + 1;
        mutGene.elem = Variables[boolGene];
        mutGene.oper = 4;
        fout << "\ndiv mutation";
    }
    if (boolOper == 3)
    {
        int boolGene = rand() % (Variables.size() - 1) + 1;
        mutGene.elem = Variables[boolGene];
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
Individ numGA(Individ inputInd, vector<vector<struct data>> ExpData, Symbolic x, vector<Symbolic> Variables,
                                                        vector<vector<double>> VarValues, std::string foutname)
{
    Individ outputInd;
    vector<Individ> numPop; //Population for numeric GA
    double resCoef = 0.0;
    int bol = 0;
    int dec = 0; //decision - "is there any genes to optimize?"

    std::ofstream fout;
    fout.open(foutname, std::ios_base::app);

    std::ofstream fitfile;
    fitfile.open("Data/Diploma/FitfileNUM.txt");

    numPop.clear();
    numPop.resize(NGAsize);
    std::cout << "Numeric GA for V(t) = " << inputInd.ind << std::endl;
    fout << "------------NUMERIC GA STARTED-------------";
    fout << "\n Ind " << inputInd.ind << std::endl;

    outputInd.ind = inputInd.ind;
    outputInd.genes = inputInd.genes;

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
                    resCoef = double(rand() % 9) + 0.999;
                }
                else if (bol == 1)
                {
                    resCoef = double(-(rand() % 9)) + 0.999;
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
            fout << "New ind for NGA is " << numPop[i].ind << std::endl;
            numPop[i].fit = CalcFit(numPop[i].ind, ExpData, x, Variables, VarValues);
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

        int limit = dec*NGA_LOOPS_LIMIT; //manual limit for GA loops
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
            if ((numPop[maxd].fit == 0) || (abs(1 - (numPop[maxd].fit / numPop[mind].fit)) < 0.00001))
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
            int numMOM = rand() % (NGAsize - 1);
            Individ MOM = numPop[numMOM];
            fout << "num MOM is " << MOM.ind << " and fit is " << MOM.fit << std::endl;

            //DAD
            int numDAD = rand() % (NGAsize - 1);
            while (numDAD == numMOM)
            {
                numDAD = rand() % (NGAsize - 1);
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
            if (randProb < 0.6)
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
                KID.fit = CalcFit(KID.ind, ExpData, x, Variables, VarValues);
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

vector<int> MomDadChoice(Population popul, vector<Symbolic> Variables)
{
    vector<int> nums;
    nums.clear();
    nums.reserve(2);
    
    int numMOM = 0;

    //MOM
    while (popul.inds[numMOM].genes.size() < (Variables.size()-1))
    {
        numMOM = rand() % popul.inds.size();
    }
    Individ MOM = popul.inds[numMOM];
    
    //DAD
    int numDAD = rand() % popul.inds.size();

    while (numDAD == numMOM) //Checking if DAD is the same as MOM
    {
        numDAD = rand() % popul.inds.size();
    }    
    Individ DAD = popul.inds[numDAD];
    
    nums.push_back(numMOM);
    nums.push_back(numDAD);
    return(nums);

}

Individ ClassicCrossover(Individ MOM, Individ DAD, vector<Symbolic> Variables)
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

    std::cout << "Final KID is " << KID.ind << std::endl;
    return(KID);
}

//Function for GA symbolic optimization (main)
Individ symbGA(Population popul, vector<vector<struct data>> ExpData, Symbolic x, vector<Symbolic> Variables,
    vector<vector<double>> VarValues, unsigned int startime, std::string foutname, std::string fitfilename)
{
    Individ outputInd;
    double stopPoint = 95; //Another stopping mechanism - desired accuracy
    
    //MAIN GA LOOP
    int mind, maxd;
    double coef1, coef2;
    double fitAVG = 0.0;
    Gene bufGene;
    bufGene.elem = Variables[1];
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
        nums = MomDadChoice(popul, Variables);
        numMOM = nums[0];
        numDAD = nums[1];

        Individ MOM = popul.inds[numMOM];
        Individ DAD = popul.inds[numDAD];
        fout << "MOM is " << MOM.ind << " and fit is " << MOM.fit << std::endl;
        fout << "DAD is " << DAD.ind << " and fit is " << DAD.fit << std::endl;

        //KID
        Individ KID;

        fout << "ClassicCrossover" << std::endl;
        int numAmount = ((Variables.size() - 1) * 5) + 1; //Checking if there're too many genes
        while (numAmount > ((Variables.size() - 1) * 7))
        {
            KID = ClassicCrossover(MOM, DAD, Variables);
            numAmount = KID.genes.size();
        }

        double boolCross = (1 / (double(rand() % 10) + 1)); //Mutation probability
        if (boolCross < 0.5)
        {
            KID = SymbMutation(KID, foutname, Variables);
        }

        KID = SearchForAbscentVars(KID, Variables);
        fout << "KID is " << KID.ind << std::endl;
        KID = numGA(KID, ExpData, x, Variables, VarValues, foutname);
        KID.genes = InputGeneDecomposition(KID.ind);
        KID = SearchForAbscentVars(KID, Variables);        
        KID.fit = CalcFit(KID.ind, ExpData, x, Variables, VarValues);
        fout << "Optimimzed KID is " << KID.ind << " and fit is " << KID.fit << std::endl;
        //Replacing the worst element of the population with the KID

        popul.inds[mind] = KID;

        f++;
        outputInd = popul.inds[maxd];
    };

    std::cout << "symbGA optimization failed, ending the loop" << std::endl;
    fout << "symbGA optimization failed, ending the loop" << std::endl;
    fout.close();
    fitfile.close();

    return(outputInd);
}

int main() {

    std::clock_t start;
    unsigned int startime = clock();
    std::fixed;

    Symbolic y("y"); // angle
    Symbolic x("x"); // h/d
    Symbolic ss("ss"); // sigma ss / sr
    
    
    // y = x * 0.5 + ss;
    // y = atan((1/(2*x))*(-(1+0.5*x)+((1+0.5*x)^2+4*x*(0.58*ss))^0.5)) * 180.0 / PI;
    //y = (1 / (2 * x)) * (-(1 + 0.5 * x) + ((1 + 0.5 * x) ^ 2 + 4 * x * (0.58 * ss)) ^ 0.5);
    y = (0.5 * x) * (-(1 + 0.5 * x) + ((1 + 0.5 * x) ^ 2 + 4 * x * (0.58 * ss)) ^ 0.5);
    int numInd = 10; // number of inds in te population

    vector<Symbolic> Variables; // First variable must be the wanted one
    Variables.push_back(y);
    Variables.push_back(x);
    Variables.push_back(ss);

    vector<vector<double>> VarValues; // Values of the parameters
    vector<double> hdValues; // from the text files
    vector<double> ssValues;
    ssValues.push_back(0.848); // B4C 6.154
    ssValues.push_back(1); // Al2O3 7.252
    VarValues.push_back(hdValues);
    VarValues.push_back(ssValues);

    Individ outputInd; // Buffer for Individ class
    vector <vector<struct data>> ExpData; // vector of expdata vectors
    ExpData.clear();
    
   // outputInd.ind = y;

    // TEST
//    char test[100];
  //  test << outputInd.ind;
   // std::cout << "\nsprintf is " << test << std::endl;

    // Logger initilization
    std::ofstream fout;
    fout.open(LOGFILE);
    logger("Program started\n", term_and_file);
    fout.close();

    //Output file for Fitness function
    std::ofstream fitfile;
    std::string fitfilename = "Data/Diploma/fitfile.txt";
    fitfile.open(fitfilename);
    fitfile.close();

    Population popul;
    popul = CreatePop(popul, Variables, numInd, LOGFILE); //Creating 1st population

    //Insert here the path to the input data file
    //WARNING! All the phrases must be deleted from the file
    vector<std::string> datafiles;
    datafiles.push_back("Data/Diploma/1.txt");
    datafiles.push_back("Data/Diploma/2.txt");
    
    vector<std::string> expNormFiles;
    expNormFiles.push_back("Data/Diploma/1expNorm.txt");
    expNormFiles.push_back("Data/Diploma/2expNorm.txt");

    vector<struct data> bufData;
    for (int f = 0; f < datafiles.size(); f++)
    {
        std::ifstream datafile(datafiles[f]);
        if (!datafile.is_open())
        {
            logger("\nInput data file not found.\
                    The Program will be terimenated\n", term_and_file);
            return (1);
        }
    }

    ExpData = SetData(datafiles, expNormFiles);

    //Fitness function calculations for the 1st gen
    for (int i = 0; i < numInd; i++)
    {
        popul.inds[i] = numGA(popul.inds[i], ExpData, x, Variables, VarValues, LOGFILE);
        popul.inds[i] = SearchForAbscentVars(popul.inds[i], Variables);
        popul.inds[i].fit = CalcFit(popul.inds[i].ind, ExpData, x, Variables, VarValues);
    }

    // TODO: delete logfile from input parameters
    outputInd = symbGA(popul, ExpData, x, Variables, VarValues, startime, LOGFILE, fitfilename);

    std::cout << "PROGRAM RESULT IS " << outputInd.ind << " with fit = " << outputInd.fit << std::endl;
    fout << "PROGRAM RESULT IS " << outputInd.ind << " with fit = " << outputInd.fit << std::endl;

    //Program runtime calculation
    unsigned int endtime = clock();
    double runtime = (endtime - startime) / (double)CLOCKS_PER_SEC;
    std::cout << "\nTotal runtime is " << runtime << " seconds" << std::endl;
    fout << "\n Total runtime is " << runtime << " seconds" << std::endl;

    fout.close();

    return (0);
}
