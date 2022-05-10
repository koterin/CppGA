#ifndef LIB_SYMBGA_H_
#define LIB_SYMBGA_H_

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

    logger(terminal, (char *)("\nFinal KID is "));
    symb_logger(terminal, KID.ind);
    return(KID);
}


Individ SymbMutation(Individ KID, vector<Symbolic> Variables)
{
    std::ofstream fout;
    int boolOper = rand() % 4 + 1;
    Gene mutGene;

    logger(file, (char *)("Mutation\nThe KID was "));
    symb_logger(file, KID.ind);

    if (boolOper == 1)
    {
        int boolGene = rand() % (Variables.size() - 1) + 1;
        mutGene.elem = Variables[boolGene];
        mutGene.oper = 3;
        logger(file, (char *)("mult mutation"));
    }
    if (boolOper == 2)
    {
        int boolGene = rand() % (Variables.size() - 1) + 1;
        mutGene.elem = Variables[boolGene];
        mutGene.oper = 4;
        logger(file, (char *)("div mutation"));
    }
    if (boolOper == 3)
    {
        int boolGene = rand() % (Variables.size() - 1) + 1;
        mutGene.elem = Variables[boolGene];
        mutGene.oper = 5;
        logger(file, (char *)("pow mutation"));
    }
    if (boolOper == 4)
    {
        mutGene.elem = double((rand() % 12 + 1) / (rand() % 10 + 1));
        mutGene.oper = 5;
        logger(file, (char *)("num pow mutation"));
    }

    KID.genes.push_back(mutGene);
    KID = IndFromGenes(KID.genes);
    KID.genes = InputGeneDecomposition(KID.ind);

    return(KID);
}

// Function for GA symbolic optimization (main)
Individ symbGA(Population popul, vector<vector<struct data>> ExpData, Symbolic x, vector<Symbolic> Variables,
    vector<vector<double>> VarValues, unsigned int startime, std::string fitfilename)
{
    Individ outputInd;
    double stopPoint = 95;  // Another stopping mechanism - desired accuracy
    
    // MAIN GA LOOP
    int mind, maxd;
    double coef1, coef2;
    double fitAVG = 0.0;
    Gene bufGene;
    bufGene.elem = Variables[1];
    bufGene.oper = 3;

    std::ofstream fitfile;
    fitfile.open(fitfilename, std::ios_base::app);

    logger(term_and_file, (char *)("Stop point is %.f%% accuracy\n"), stopPoint);
    logger(term_and_file, (char *)("\n--------------------------THE SYMBOLIC GA STARTED"));
    logger(term_and_file, (char *)("--------------------------\n"));
 
    int f = 0;
    maxd = 0;
    while (popul.inds[maxd].fit < 100) {

        unsigned int nowtime = clock();
        double curtime = (nowtime - startime) / (double)CLOCKS_PER_SEC;
        logger(term_and_file, (char *)("\ncurrent time is %.f seconds\n"), curtime);

        mind = 0;
        maxd = 0;
        fitAVG = 0.0;

        // Displaying current population
        logger(term_and_file, (char *)("\nsymbGA Population %d\n"), f + 1);
        for (int g = 0; g < popul.inds.size(); g++)
        {
            symb_logger(file, popul.inds[g].ind);
            logger(term_and_file, (char *)(" and fit %.2f\n"), popul.inds[g].fit);

            // Finding minimum
            if (popul.inds[g].fit < popul.inds[mind].fit)
            {
                mind = g;
            }

            // Finding maximum
            if (popul.inds[g].fit > popul.inds[maxd].fit)
            {
                maxd = g;
            }

            fitAVG += popul.inds[g].fit;
        }

        fitAVG = fitAVG / popul.inds.size();
        fitfile << f + 1 << " " << popul.inds[mind].fit << " " << fitAVG << " " << popul.inds[maxd].fit << std::endl;

        logger(term_and_file, (char *)("\nMIN in Population %d: %.2f%% "), f + 1, popul.inds[mind].fit);
        symb_logger(term_and_file, popul.inds[mind].ind);
        logger(term_and_file, (char *)("\nMAX in Population %d: %.2f%% "), f + 1, popul.inds[maxd].fit);
        symb_logger(term_and_file, popul.inds[maxd].ind);
        
        //Checking if the current minimum fit is the desired one
        if (popul.inds[maxd].fit > stopPoint)
        {
            logger(term_and_file, (char *)("\n---------------FINAL---------------\n"));
            logger(term_and_file, (char *)("Optimum coefficients found after %d loops, "), f + 1);
            logger(term_and_file, (char *)("final fit is %.2f%% from ind \n"), popul.inds[maxd].fit);
            symb_logger(term_and_file, popul.inds[maxd].ind);

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

        logger(term_and_file, (char *)("\nMOM is %.2f "), MOM.fit);
        symb_logger(file, MOM.ind);
        logger(term_and_file, (char *)("\nDAD is %.2f "), DAD.fit);
        symb_logger(file, DAD.ind);

        //KID
        Individ KID;

        logger(file, (char *)("ClassicCrossover"));
        int numAmount = ((Variables.size() - 1) * 5) + 1; //Checking if there're too many genes
        while (numAmount > ((Variables.size() - 1) * 7))
        {
            KID = ClassicCrossover(MOM, DAD, Variables);
            numAmount = KID.genes.size();
        }

        double boolCross = (1 / (double(rand() % 10) + 1)); //Mutation probability
        if (boolCross < 0.5)
        {
            KID = SymbMutation(KID, Variables);
        }

        KID = SearchForAbscentVars(KID, Variables);
        logger(file, (char *)("\nKID is "));
        symb_logger(file, KID.ind);
        KID = numGA(KID, ExpData, x, Variables, VarValues);
        KID.genes = InputGeneDecomposition(KID.ind);
        KID = SearchForAbscentVars(KID, Variables);        
        KID.fit = CalcFit(KID.ind, ExpData, x, Variables, VarValues);

        logger(term_and_file, (char *)("\nOptimimzed KID is %.2f%% "), KID.fit);
        symb_logger(file, KID.ind);

        //Replacing the worst element of the population with the KID
        popul.inds[mind] = KID;

        f++;
        outputInd = popul.inds[maxd];
    };

    logger(term_and_file, (char *)("\n symbGA optimization failed, ending the loop"));
    fitfile.close();

    return(outputInd);
}

#endif // LIB_SYMBGA_H_
