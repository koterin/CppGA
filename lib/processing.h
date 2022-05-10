#ifndef LIB_PROCESSING_H_
#define LIB_PROCESSING_H_

int checkDatafiles(vector<std::string> datafiles)
{
    for (int f = 0; f < datafiles.size(); f++)
    {
        std::ifstream datafile(datafiles[f]);
        if (!datafile.is_open())
        {
            logger(term_and_file, (char *)("\nInput data file not found.\
                    The Program will be terimenated\n"));
            return (0);
        }
    }
    return (1);
}

#endif // LIB_PROCESSING_H_
