#include <string>
#include <iostream>
#include <fstream>

#include "eye_analysis.hh"


void _parse_args(int argc, const char * argv[], int * option,
                float * object_distance, float * offaxis, std::string * model,
                float * age, float * pupil, float * diopters, std::string * param,
                int * iters, float * wavelength)
{    
    for (int i = 1; i < argc; i++) 
    {   
        std::string word = argv[i];
        // transform to lower case if not already
        std::transform(word.begin(), word.end(), word.begin(), ::tolower);

        if (word == "plot") { *option = 0; }
        if (word == "lsa") { *option = 1; } 
        if (word == "series") { *option = 2; }  
        if (word == "spot") { *option = 3; }   

        if (word.substr(0, 5) == "--age")
            { *age = atof(word.substr(6).c_str()); }
        if  (word.substr(0, 2) == "-a")
            { *age = atof(word.substr(3).c_str()); }

        if (word.substr(0, 7) == "--pupil")
            { *pupil = atof(word.substr(8).c_str()); }
        if (word.substr(0, 2) == "-t")
            { *pupil = atof(word.substr(3).c_str()); }

        if (word.substr(0, 10) == "--diopters") 
            { *diopters = atof(word.substr(11).c_str()); } 
        if (word.substr(0,2) == "-d") 
            { *diopters = atof(word.substr(3).c_str()); } 

        if (word.substr(0, 7) == "--param")
            { *param = word.substr(8); }
        if (word.substr(0,2) == "-p")
            { *param = word.substr(3); }

        if (word.substr(0, 7) == "--iters")
            { *iters = atof(word.substr(8).c_str()); }
        if (word.substr(0,2) == "-i")
            { *iters = atof(word.substr(3).c_str()); }

        if (word.substr(0, 7) == "--wvlen")
            { *wavelength = atof(word.substr(8).c_str()); }
        if (word.substr(0,2) == "-w")
            { *wavelength = atof(word.substr(3).c_str()); }

        if (word.substr(0, 10) == "--distance") 
            { *object_distance = atof(word.substr(11).c_str()); }
        if (word.substr(0,2) == "-l")
            { *object_distance = atof(word.substr(3).c_str()); }

        if (word.substr(0, 9) == "--offaxis") 
            { *offaxis = atof(word.substr(10).c_str()); }
        if (word.substr(0,2) == "-o")
            { *offaxis = atof(word.substr(3).c_str()); }

        // and age option for loop.
        if (word.substr(0, 6) == "-model")
        {
            if (word.substr(7) == "navarro" or word.substr(7) == "dubbelman") 
                { *model = word.substr(7); }
            else 
                {std::cout << "sorry model option not understood, \
                 using dubbelman" << std::endl;}
        }

    }
    if (*param != "age" && *param != "pupil" && *param != "angle" && *param != "focus" && 
        *param != "distance")
    {
        std::cout << "sorry option not understood, using object angle." << std::endl;
        *param = "angle";
    }
}

void _help()
{
    std::cout << " " << std::endl;
    std::cout << "Schematic Eye" << std::endl;
    std::cout << " " << std::endl;
    std::cout << "options:" << std::endl;
    std::cout << "============================================================" << std::endl;
    std::cout << "plot \t\t\t choose plot option" << std::endl;
    std::cout << "lsa \t\t\t choose LSA option, not working" << std::endl;
    std::cout << "series \t\t\t chose series option, add additional flags" << std::endl;
    std::cout << "spot \t\t\t choose a spot plot" << std::endl;
    std::cout << " " << std::endl;
    std::cout << "-i\t--iters=ITERS \t\t set number of iterations for a series plot" << std::endl;
    std::cout << "-a\t--age=AGE \t\t set age parameter - Dubbelman only" << std::endl;
    std::cout << "-t\t--pupil=PUPIL \t\t set pupil width" << std::endl;
    std::cout << "-d\t--diopters=DIOPTERS \t set lens accommodation in diopters" << std::endl;
    std::cout << "-l\t--distance=DISTANCE \t set object distance (mm)" << std::endl;
    std::cout << "-o\t--offaxis=OFFAXIS \t set object offaxis in degrees" << std::endl;
    std::cout << "-m\t--model=MODEL \t\t set model [dubbelman or navarro]" << std::endl;
    std::cout << "-w\t--wvlen=WVLEN \t\t set wavelength traced (nm)" << std::endl;
    std::cout << "-p\t--param=PARAM \t\t choose a parameter to iterate - series only" << std::endl;
    std::cout << "\t\t\t [age, pupil, angle, focus, distance, wavelength]" << std::endl;
}

int main(int argc, const char * argv[])
{
    // set defaults:
    int * option = new int (0);
    int * iters = new int (4);
    float * object_distance = new float (1000000);
    float * offaxis = new float (0);
    float * age = new float (20);
    float * pupil = new float (3);
    float * diopters = new float (0);
    float * wavelength = new float (555.0);
    std::string * param = new std::string ("angle");
    std::string * model = new std::string ("navarro");

    std::string word = "";
    if (argc == 2) { word = argv[1]; }

    if (word == "--help" or word == "help" or word == "-h") 
    {
        _help();
        return 0;
    }
    else
    {
        _parse_args(argc, argv, option, object_distance, offaxis, model, 
                age, pupil, diopters, param, iters, wavelength);
        
        SchematicEye::Analysis analysis;
        
        if (*option == 0)
        {
            analysis.EyePlots(1, *object_distance, *offaxis, *model, 
                        *age, *pupil, *diopters, *wavelength);
        }
        
        if (*option == 1)
        {
            analysis.LSAanalysis(*object_distance, *offaxis, *model);
        }

        if (*option == 2)
        {
            analysis.IntensityAnalysis(*param, *object_distance, *offaxis, *model, 
                        *age, *pupil, *diopters, *iters, *wavelength);
        }

        if (*option == 3)
        {
            analysis.SpotPlot(*object_distance, *offaxis, *model, 
                        *age, *pupil, *diopters, *wavelength);
        }

        return 0;
    }
}
