#include "routine.h"

std::vector<std::string>
StringSplit(std::string line, char delim)
{
    stringstream sline(line);
    string word;
    vector<string> linelist;
    while (getline(sline, word, delim)) {
        if (word.size() > 0) {
            linelist.push_back(word);
        }
    }
    return linelist;
}




