#ifndef BASIC_H
#define BASIC_H

#include <iostream>
#include <string>
#include <sstream>
#include <vector>

using namespace std;

vector<string>
StringSplit(string line, char delim);

template <class Type>
Type Str2Num(const string &str)
{
    stringstream is;
    is << str;
    Type Num = -1;
    is >> Num;
    return Num;
}

#endif // _COMMON_H
