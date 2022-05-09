#ifndef LIB_LOGGER_H_
#define LIB_LOGGER_H_

#include <iostream>
#include <string>
#include <fstream> //for files
#include <ctime> //for runtime calculations


enum outputType
{
    terminal = 0, 
    file = 1,
    term_and_file = 2
};

void logger(std::string message, outputType flag)
{
    if ((flag == 0) || (flag == 2))
        std::cout << message << std::endl;
    if ((flag == 1) || (flag == 2))
    {
        std::ofstream fout;
        fout.open(LOGFILE, std::ios::app);
        fout << message;
        fout.close();
    }   
}

std::ofstream logger_initialization(void)
{
    std::ofstream fout;
    fout.open(LOGFILE);
    fout.close();
    return (fout);
}

#endif // LIB_LOGGER_H_
