#ifndef LIB_CPPGA_H_
#define LIB_CPPGA_H_

#define PI 3.14159265
#define NGA_LOOPS_LIMIT 100 // Number of loops limit (limit * number of elements)
#define NGAsize 10 // Number of inds in the numerical GA
#define SGAsize 10 // Number of inds in the symbolic GA
#define FILEDIR "Data/Diploma/"
#define LOGFILE "logs.txt"
#define DATA_LINES_NUM 10 // Number of data lines to read from file

#include <iostream>
#include <string>
#include <fstream> //for files
#include <ctime> //for runtime calculations
#include <list>
#include <cmath>
#include <iterator>
#include <cstring>
#include <cstdarg>

#include "logger.h"
#include "processing.h"
#include "genetic_classes.h"
#include "genetic_functions.h"
#include "expdata.h"
#include "numga.h"
#include "symbga.h"

// TODO logger level (INFO, DEBUG);
// TODO logger fitfiles
// TODO decomposition
// TODO file names + filedir

#endif // LIB_CPPGA_H_
