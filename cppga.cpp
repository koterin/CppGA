#define _SCL_SECURE_NO_WARNINGS

#include "lib/symbolicc++.h" //for symbolic
#include "lib/cppga.h"

int main() {

    std::clock_t start;
    unsigned int startime = clock();
    Population popul;
    vector<struct data> bufData;
    Individ outputInd;  // Buffer for Individ class
    vector <vector<struct data>> ExpData;  // vector of expdata vectors
    vector<std::string> datafiles;
    vector<std::string> expNormFiles;

    Symbolic y("y"); // angle
    Symbolic x("x"); // h/d
    Symbolic ss("ss"); // sigma ss / sr
     
    // y = x * 0.5 + ss;
    // y = atan((1/(2*x))*(-(1+0.5*x)+((1+0.5*x)^2+4*x*(0.58*ss))^0.5)) * 180.0 / PI;
    //y = (1 / (2 * x)) * (-(1 + 0.5 * x) + ((1 + 0.5 * x) ^ 2 + 4 * x * (0.58 * ss)) ^ 0.5);
    y = (0.5 * x) * (-(1 + 0.5 * x) + ((1 + 0.5 * x) ^ 2 + 4 * x * (0.58 * ss)) ^ 0.5);

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

    ExpData.clear();
    
    // Logger initilization
    std::ofstream fout;
    fout.open(LOGFILE);
    fout.close();

    // Output file for Fitness function
    std::ofstream fitfile;
    std::string fitfilename = "Data/Diploma/fitfile.txt";
    fitfile.open(fitfilename);
    fitfile.close();

    // Insert here the path to the input data file
    // WARNING! All the phrases must be deleted from the file
    datafiles.push_back("Data/Diploma/1.txt");
    datafiles.push_back("Data/Diploma/2.txt");

    // Path to the normalized input data (will be created by the SetData func)
    expNormFiles.push_back("Data/Diploma/1expNorm.txt");
    expNormFiles.push_back("Data/Diploma/2expNorm.txt");

    if (!checkDatafiles(datafiles))
        return (1);

    ExpData = SetData(datafiles, expNormFiles);
    popul = CreatePop(popul, Variables, SGAsize);  // Creating 1st population

    // Fitness function calculations for the 1st gen
    for (int i = 0; i < SGAsize; i++)
    {
        popul.inds[i] = numGA(popul.inds[i], ExpData, x, Variables, VarValues);
        popul.inds[i] = SearchForAbscentVars(popul.inds[i], Variables);
        popul.inds[i].fit = CalcFit(popul.inds[i].ind, ExpData, x, Variables, VarValues);
    }

    outputInd = symbGA(popul, ExpData, x, Variables, VarValues, startime, fitfilename);

    logger(term_and_file, (char *)("PROGRAM RESULT IS %.2f "), outputInd.fit);
    symb_logger(file, outputInd.ind);

    // Program runtime calculation
    unsigned int endtime = clock();
    double runtime = (endtime - startime) / (double)CLOCKS_PER_SEC;
    logger(term_and_file, (char *)("Total runtime is %.f seconds"), runtime);

    return (0);
}
