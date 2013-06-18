#include <string>
#include <iostream>
#include <fstream>

#include "SchematicEye.hh"


void _parse_args(int argc, const char * argv[], int * option,
                float * object_distance, float * offaxis, std::string * model)
{    
    for (int i = 1; i < argc; i++) 
    {   
        std::string word = argv[i];
        // transform to lower case if not already
        std::transform(word.begin(), word.end(), word.begin(), ::tolower);

        if (word == "plot") { *option = 0; }
        if (word == "lsa") { *option = 1; } 
        if (word == "series") { *option = 2; }       

        if (word.substr(0, 8) == "distance") 
            { *object_distance = ::atof(word.substr(9).c_str()); }
        if (word.substr(0, 7) == "offaxis") 
            { *offaxis = ::atof(word.substr(8).c_str()); }

        // add in model change param (i.e. Navarro vs Dubbelman)
        // and age option for loop.
        if (word.substr(0, 5) == "model")
        {
            if (word.substr(6) == "navarro" or word.substr(6) == "dubbelman") 
                { *model = word.substr(6); }
            else 
                {std::cout << "sorry model option not understood, \
                 using dubbelman" << std::endl;}
        }

    }
}

int main(int argc, const char * argv[])
{
    // set defaults:
    int * option = new int (0);
    float * object_distance = new float (100000);
    float * offaxis = new float (0);
    std::string * model = new std::string ("dubbelman");

    _parse_args(argc, argv, option, object_distance, offaxis, model);

    Eye::Eye eye;

    if (*option == 0)
    {
        eye.SimplePlot(*object_distance, *offaxis, *model);
    }
    
    if (*option == 1)
    {
        eye.LSAanalysis(*object_distance, *offaxis, *model);
    }

    if (*option == 2)
    {
        eye.AccommodationAnalysis(*object_distance, *offaxis, *model);
    }
    /*
    if (*option == 3)
    {
        eye.ImageSeriesAnalysis(*offaxis, *model);
    }
    */
    
    return 0;
}
