#include <Goptical/Data/Set>
#include <Goptical/Data/PlotData>
#include <Goptical/Data/DiscreteSet>
#include <Goptical/Data/Plot>

#include <Goptical/Analysis/RayFan>
#include <Goptical/Analysis/Spot>
#include <Goptical/Analysis/Focus>

#include <Goptical/Trace/Result>
#include <Goptical/Trace/Ray>

#include <Goptical/Io/RendererSvg>
#include <Goptical/Io/Rgb>

#include <armadillo>
#include <math.h>

#include "eye_eye.hh"


namespace SchematicEye 
{
    class Analysis
    {

    private:

        Goptical::Analysis::RayFan * fan;

    public:
        
        Analysis();
        ~Analysis();

        double _max_intensity;

        void LSAanalysis(float object_distance, float off_axis,
            std::string model);
        void EyePlots(float object_distance, 
                    float off_axis, std::string model, float age,
                    float pupil_size, float diopters, float wavelength);
        void EyePlots(int best_focus, float object_distance, 
                    float off_axis, std::string model, float age,
                    float pupil_size, float diopters, float wavelength);
        float *IntensityAnalysis(float object_distance, float off_axis,
                    float pupil_size, float diopters, float wavelength);
        void IntensityAnalysis(std::string param, float object_distance, 
                    float off_axis, std::string model, float age,
                    float pupil_size, float diopters, int iter,
                    float wavelength);
        void SpotPlot(float object_distance, 
                    float off_axis, std::string model, float age,
                    float pupil_size, float diopters, float wavelength);


        arma::vec PSF(arma::vec intensity, arma::vec xvals, bool symmetric);
        arma::vec MTF(arma::vec psf);

        inline float _get_input(std::string message)
        {
            float out;
            std::cout << "enter " << message << ": " << std::endl;
            std::cin >> out;
            return out;
        }

    };
}