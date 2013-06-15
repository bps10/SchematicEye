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

        if (word == "plot") { *option = 0; }
        if (word == "loop") { *option = 1; } 
        if (word == "series") { *option = 2; }       

        if (word.substr(0, 8) == "distance") 
            { *object_distance = ::atof(word.substr(9).c_str()); }
        if (word.substr(0, 7) == "offaxis") 
            { *offaxis = ::atof(word.substr(8).c_str()); }

        // add in model change param (i.e. Navarro vs Dubbelman) and age option for loop.
        if (i == 2 && option == 0)
        {
            std::string word = argv[2];
            if (argv[i] == "navarro" or argv[i] == "dubbelman") 
                { *model = argv[i]; }
            else 
                {std::cout << "sorry model option not understood" << std::endl;}
        }

    }
}

int main(int argc, const char * argv[])
{
    // set defaults:
    int * option;
    float * object_distance;
    float * offaxis;
    std::string * model;

    option = new int (0);
    object_distance = new float (10000);
    offaxis = new float (0);
    model = new std::string ("dubbelman");

    std::cout << *object_distance << std::endl;

    _parse_args(argc, argv, option, object_distance, offaxis, model);
    
    std::cout << *object_distance << std::endl;

    Eye::Eye eye;

    if (argc < 1 || *option == 0)
    {
        
        eye.SimplePlot(*object_distance, *offaxis);
    }
    
    if (*option == 1)
    {

        eye.LSAanalysis(*option, *object_distance, *offaxis, *model);
    }
    if (*option == 2)
    {
    std::cout << "starting values... " << std::endl;
    std::cout << "  "  << std::endl;
    
    Eye::Eye        eye;
    eye.set_params("dubbelman");
    eye.EyePlots(2, *object_distance, *offaxis);
    eye.Intensity(1, *object_distance, *offaxis);
    }
    
    
    return 0;
}
