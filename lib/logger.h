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

void print_logger(outputType flag, char *message)
{
    if ((flag == 0) || (flag == 2))
    {
        std::cout << message;
    }
    if ((flag == 1) || (flag == 2))
    {
        std::ofstream fout;
        fout.open(LOGFILE, std::ios::app);
        fout << message;
        fout.close();
    }   
}

void symb_logger(outputType flag, Symbolic symb)
{
    if ((flag == 0) || (flag == 2))
    {
        std::cout << symb;
    }
    if ((flag == 1) || (flag == 2))
    {
        std::ofstream fout;
        fout.open(LOGFILE, std::ios::app);
        fout << symb;
        fout.close();
    }   
}

void logger(outputType flag, char *format, ...)
{
    va_list args;
    char str[100];

    memset(str, 0, strlen(str));
    va_start(args, format);
    sprintf(str, format, args);
    print_logger(flag, str);
    va_end(args);
}

#endif // LIB_LOGGER_H_
