#include <string>
#include <iostream>
#include <fstream>

#include "SchematicEye.hh"


int main(int argc, const char * argv[])
{
    int option (0);
    float object_distance = 10000;
    float offaxis = 5;
    
    for (int i = 1; i < argc; i++) {
        
        std::string word = argv[1];

        if (word == "plot") {option = 0;}
        if (word == "loop") {option = 1;} 
        if (word == "series") {option = 2;}       

        if (word.substr(0, 8) == "distance") 
            { object_distance = ::atof(word.substr(9).c_str()); }
        if (word.substr(0, 7) == "offaxis") 
            { offaxis = ::atof(word.substr(8).c_str()); }

        // add in model change param (i.e. Navarro vs Dubbelman) and age option for loop.
        if (i == 2 && option == 0)
        {
            std::string word = argv[2];
            if (argv[i] == "navarro" or argv[i] == "dubbelman") 
                {std::string  mod;
                mod = argv[i];}
            else 
                {std::cout << "sorry model option not understood" << std::endl;}
        }

        if (i > 2) {std::cout << "sorry additional options not understood" << std::endl;}
    }
    
    if (argc < 1 || option == 0)
    {
    Eye::Eye     eye;
    eye.set_params(option, "dubbelman");
    eye.SchematicEye();
    eye.EyeTracer(object_distance, offaxis);
    eye.EyePlots(1, object_distance, offaxis); 
    eye.Diopters(0);
    }
    
    if (option == 1)
    {

    std::ofstream outputfile;
    outputfile.open ("dat/EyeLSA.csv", std::ios::trunc);
    float AGE, LensAccomm, PupilSize, AccommOptPower, RelaxedOptPower, Defocus_diopters;
    std::string mod;
    
    for (int i = 0; i < 20; i++)
    {
        
        for (int j = 1; j < 21; j++)
        {
            for (int k = 10; k <25; k++)
            {
                LensAccomm = i / 2.0;
                PupilSize = j / 4.0;
                AGE = k;
                
                Eye::Eye     eye;
                eye.set_params( LensAccomm, PupilSize, mod = "dubbelman", AGE );
                eye.SchematicEye();
                eye.EyeTracer(object_distance, offaxis);
                
                AccommOptPower = eye.FindOpticalPower(1);
                RelaxedOptPower = eye.FindOpticalPower(2);
                Defocus_diopters = eye.Diopters(option);
                
                outputfile << LensAccomm << "," << PupilSize << "," << AGE << "," << 
                            AccommOptPower << "," << RelaxedOptPower << "," << 
                            Defocus_diopters << std::endl;
            }
        }
    }
    
    outputfile.close();
     }
    
    if (option == 2)
    {
    std::cout << "starting values... " << std::endl;
    std::cout << "  "  << std::endl;
    
    Eye::Eye        eye;
    eye.set_params(option, "dubbelman");
    eye.EyePlots(2, object_distance, offaxis);
    eye.Intensity(1, object_distance, offaxis);
    }
    
    
    return 0;
}