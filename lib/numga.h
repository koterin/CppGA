#ifndef LIB_NUMGA_H_
#define LIB_NUMGA_H_

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

//Function for GA coefficient optimization
Individ numGA(Individ inputInd, vector<vector<struct data>> ExpData, Symbolic x,
                vector<Symbolic> Variables,vector<vector<double>> VarValues)
{
    Individ outputInd;
    vector<Individ> numPop; //Population for numeric GA
    double resCoef = 0.0;
    int bol = 0;
    int dec = 0; //decision - "is there any genes to optimize?"

    std::ofstream fitfile;
    fitfile.open("Data/Diploma/FitfileNUM.txt");

    numPop.clear();
    numPop.resize(NGAsize);
    logger(file, (char *)("\n\n------------NUMERIC GA STARTED-------------"));
    logger(file, (char *)("\nNumeric GA for V(t) = "));
    symb_logger(file, inputInd.ind);

    outputInd.ind = inputInd.ind;
    outputInd.genes = inputInd.genes;

    for (int i = 0; i < numPop.size(); i++)
    {
        numPop[i] = outputInd;
    }

    // Looking for numeric genes in the Individual
    for (int i = 0; i < inputInd.genes.size(); i++)
    {
        auto buf = inputInd.genes[i].elem->clone();
        
        if ((typeid(*buf) == typeid(Number<double>)) || (typeid(*buf) == typeid(Number<int>)))
        {
            dec += 1;
            // Creating new individual
            for (int j = 1; j < numPop.size(); j++)
            {
                bol = rand() % 2;  // boolean imitation

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

    // No numeric genes were found
    if (dec == 0)
    {
        logger(term_and_file, (char *)("\nNo numeric GA optimization needed for the "));
        symb_logger(term_and_file, outputInd.ind);
        return(inputInd);
    }

    // Found numeric genes
    else if (dec > 0)
    {
        // Preparing the Population
        for (int i = 0; i < numPop.size(); i++)
        {
            numPop[i] = IndFromGenes(numPop[i].genes);
            logger(file, (char *)("\nNew ind for NGA is "));
            symb_logger(file, numPop[i].ind);;
            numPop[i].fit = CalcFit(numPop[i].ind, ExpData, x, Variables, VarValues);
            auto buf = numPop[i].ind->clone();
            if ((typeid(*buf) == typeid(Number<double>)) || (typeid(*buf) == typeid(Number<int>)))
            {
                numPop[i].fit = 0;
            }
            buf->unreference(buf);
        }

        logger(file, (char *)("\nThe initial ind fit is %.2f%%"), numPop[0].fit);
        double fitAVG = 0.0;

        // MAIN GA LOOP

        int limit = dec * NGA_LOOPS_LIMIT;  // manual limit for GA loops
        logger(term_and_file, (char *)("\nNumGA limits = %f"), limit);
        int mind, maxd;
        double coef1, coef2;

        for (int f = 0; f < limit; f++) {

            mind = 0;
            maxd = 0;
            fitAVG = 0.0;

            // Displaying current population
            logger(file, (char *)("\nNumeric GA Population %d\n"), (f + 1));
            for (int g = 0; g < numPop.size(); g++)
            {
                symb_logger(file, numPop[g].ind);
                logger(file, (char *)(" and fit %.2f%%\n"), numPop[g].fit);

                // Fiding minimum
                if (numPop[g].fit < numPop[mind].fit)
                {
                    mind = g;
                }
                // Finding maximum
                if (numPop[g].fit > numPop[maxd].fit)
                {
                    maxd = g;
                }
                fitAVG += numPop[g].fit;
            }

            fitAVG = fitAVG / numPop.size();
            fitfile << f + 1 << " " << numPop[mind].fit << " " << fitAVG << " " << numPop[maxd].fit << std::endl;

            logger(file, (char *)("\nThe minimum fit in NumPopulation %d is %.2f%%"),
                                                            f + 1, numPop[mind].fit);
            logger(file, (char *)("\nThe maximum fit in NumPopulation %d is %.2f%%"),
                                                            f + 1, numPop[maxd].fit);

            // Checking if the current population is converged
            if ((numPop[maxd].fit == 0) || (abs(1 - (numPop[maxd].fit / numPop[mind].fit)) < 0.00001))
            {
                logger(file, (char *)("\nNumPop was converged"));
                logger(file, (char *)("\nOptimum coefficients found after %d loops, "), f + 1);
                logger(file, (char *)("final ind is %.2f%% - "), numPop[maxd].fit);
                symb_logger(file, numPop[maxd].ind);

                outputInd = numPop[maxd];
                fitfile.close();
                return(outputInd);
            }

            // Setting Mom and Dad as 2 random elements of the population
            // MOM
            int numMOM = rand() % (NGAsize - 1);
            Individ MOM = numPop[numMOM];
            logger(file, (char *)("\nnumMOM is %.2f%% - "), MOM.fit);
            symb_logger(file, MOM.ind);

            //DAD
            int numDAD = rand() % (NGAsize - 1);
            while (numDAD == numMOM)
            {
                numDAD = rand() % (NGAsize - 1);
            }
            Individ DAD = numPop[numDAD];
            logger(file, (char *)("\nnumDAD is %.2f%% - "), DAD.fit);
            symb_logger(file, DAD.ind);

            // KID
            coef1 = (rand() % 100) / double(100);
            Individ KID = MOM;
            
            // Crossing over
            for (int q = 0; q < KID.genes.size(); q++)
            {
                auto buf = KID.genes[q].elem->clone();
                if ((typeid(*buf) == typeid(Number<double>)) || (typeid(*buf) == typeid(Number<int>)))
                {
                    KID.genes[q].elem = coef1 * MOM.genes[q].elem + (1 - coef1) * DAD.genes[q].elem;
                }
                buf->unreference(buf);
            }
            
            // Numeric mutation
            double randProb = (rand() % 100) / double(100);
            if (randProb < 0.6)
            {
                KID = NumMutation(KID);
                logger(file, (char *)("\nNumMutation\n"));
            }

            KID = IndFromGenes(KID.genes);

            // Checking if ind is converged to a number
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

            logger(file, (char *)("\nnumKID is %.2f%% - "), KID.fit);
            symb_logger(file, KID.ind);
            
            // Replacing the worst element of the population with the KID
            numPop[mind] = KID;
            outputInd = numPop[maxd];
        };
    }

    logger(term_and_file, (char *)("\nnumGA optimization failed, ending the loop"));
    fitfile.close();
    return(outputInd);
}

#endif // LIB_NUMGA_H_
