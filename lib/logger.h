#ifndef LIB_LOGGER_H_
#define LIB_LOGGER_H_

#include <iostream>
#include <string>
#include <fstream> //for files
#include <ctime> //for runtime calculations

#define FOUTNAME "Data/Diploma/logs.txt"

void logger(std::string message, unsigned int terminal_flag, unsigned int file_flag)
{
    if (terminal_flag == 1)
        std::cout << message << std::endl;
    if (file_flag == 1)
    {
        std::ofstream fout;
        fout.open(FOUTNAME, std::ios::app);
        fout << message;
        fout.close();
    }   
}

#endif // LIB_LOGGER_H_
