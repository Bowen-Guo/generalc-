/*
   Read parameters from parameter file
*/
#include "parser.h"
#include <fstream>
#include <string>
#include <sstream>
#include <stdexcept>
#include <iostream>
using namespace std;

int readParameters(const char *filename, int &n3, int &n2, int &n1, int &n, char *fin_refl,
                   char * fin_seed, char * fout_t, char * fout_p, char * fout_q, char * flog)
{
    string line;
    string paramName;
    int paramValue = 0;
    string paramStr;

    ifstream fin(filename);
    if (!fin.good())
    {
        string msg("parameters file not found");
        msg.append(filename);
        throw runtime_error(msg);
    }
    while (fin.good())
    {
        getline(fin,line);

        
        if(line[0] != '#')
        {
          istringstream buffer(line);
          buffer >> paramName;
          if(paramName.compare("n") == 0)
          {
            buffer >> paramValue;
            n = paramValue;
          }
          else if(paramName.compare("n3") == 0)
          {
            buffer >> paramValue;
            n3 = paramValue;
          }
          else if(paramName.compare("n2") == 0)
          {
            buffer >> paramValue;
            n2 = paramValue;
          }
          else if(paramName.compare("n1") == 0)
          {
            buffer >> paramValue;
            n1 = paramValue;
          }
          else if(paramName.compare("fin_refl") == 0)
          {
            buffer >> paramStr;
            paramStr.copy(fin_refl,paramStr.length(), 0);
          }
          else if(paramName.compare("fin_seed") == 0)
          {
            buffer >> paramStr;
            paramStr.copy(fin_seed,paramStr.length(), 0);
          }
          else if(paramName.compare("fout_t") == 0)
          {
            buffer >> paramStr;
            paramStr.copy(fout_t,paramStr.length(), 0);
          }
          else if(paramName.compare("fout_p") == 0)
          {
            buffer >> paramStr;
            paramStr.copy(fout_p,paramStr.length(), 0);
          }
          else if(paramName.compare("fout_q") == 0)
          {
            buffer >> paramStr;
            paramStr.copy(fout_q,paramStr.length(), 0);
          }
          else if(paramName.compare("flog") == 0)
          {
            buffer >> paramStr;
            paramStr.copy(flog,paramStr.length(), 0);
          }
          else 
            cout << "WARNING: Key parameter is missing!!!" << endl;
        } 
      }
    fin.close();

}
