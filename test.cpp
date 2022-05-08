#define _SCL_SECURE_NO_WARNINGS
#include <iostream>
#include <string>
#include <fstream> //for files
#include "symbolicc++.h" //for symbolic
#include <ctime> //for runtime calculations
#include <list>
#include <cmath>
#include <iterator>
//#include <vector>

struct data
{
	double y;
	double x;
};

int main(void) {
    int len = 10;
	vector<std::string> fileroutes;
	fileroutes.push_back("Data/Diploma/1.txt");
	fileroutes.push_back("Data/Diploma/2.txt");
	
	vector<std::string> exproutes;
	exproutes.push_back("Data/Diploma/1expNorm.txt");
	exproutes.push_back("Data/Diploma/2expNorm.txt");
	
	vector<vector<struct data>> ExpData;
	ExpData.clear();
	struct data curdata;
	double num1, num2;

	for (int f = 0; f < fileroutes.size(); f++)
	{
		std::ifstream datafile(fileroutes[f]);
		std::string currentLine;
		int totalLen = 0;
		vector<struct data> bufData;

		//Counting the number of lines in the input file
		while (!datafile.eof())
		{
			getline(datafile, currentLine);
			std::cout << currentLine << std::endl;
            totalLen++;
		}

		int freq = (totalLen / len); //number of lines which will be repeatedly skipped
		if (freq < 1)
		{
			std::cout << "\nDesired length is bigger than the input file=" << std::endl;
			return(1);
		}

		std::cout << "\nDatafile " << f << " is " << totalLen << " lines long, input will be every "
			<< freq << " lines" << std::endl;

		datafile.clear();
        datafile.seekg(std::ios_base::beg);

		int i = 0;
		while (i < totalLen)
		{
			if (i % freq == 0)
			{
				datafile >> curdata.x >> curdata.y;
				std::cout << "WRTTING " << curdata.x << " " << curdata.y << std::endl;
                bufData.push_back(curdata);
			}

			else
			{
				datafile >> num1 >> num2;
			}

			i += 1;
		}

		ExpData.push_back(bufData);
		bufData.clear();
	}

	double tMax, vMax;
	tMax = ExpData[0][0].x;
	vMax = ExpData[0][0].y;
	for (int f = 0; f < fileroutes.size(); f++)
	{
		for (int i = 1; i < ExpData[f].size(); i++)
		{
			if (ExpData[f][i].x > tMax)
			{
				tMax = ExpData[f][i].x;
			}
			if (ExpData[f][i].y > vMax)
			{
				vMax = ExpData[f][i].y;
			}
		}
	}
	
    if (tMax == 0)
       tMax = 1;
    if (vMax == 0)
        vMax = 1;

	std::ofstream expfile;
	for (int f = 0; f < fileroutes.size(); f++)
	{
		expfile.open(exproutes[f]);
		for (int i = 0; i < ExpData[f].size(); i++)
		{
			ExpData[f][i].x = ExpData[f][i].x / tMax;
			ExpData[f][i].y = ExpData[f][i].y / vMax;
			expfile << ExpData[f][i].x << "	" << ExpData[f][i].y << std::endl;
		}
		expfile.close();
	}

    return (0);
}
