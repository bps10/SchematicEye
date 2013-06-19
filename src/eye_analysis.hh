#include <Goptical/Data/Set>
#include <Goptical/Data/PlotData>
#include <Goptical/Data/DiscreteSet>
#include <Goptical/Data/Plot>

#include <Goptical/Analysis/RayFan>
#include <Goptical/Analysis/Spot>
#include <Goptical/Analysis/Focus>


#include <Goptical/Io/RendererSvg>
#include <Goptical/Io/Rgb>

#include "eye_eye.hh"

namespace SchematicEye 
{
    class Analysis
    {

    private:

        Goptical::Analysis::RayFan * fan;

        void Intensity(int best_focus, float object_distance, float off_axis, std::string model);
        void EyePlots(int best_focus, float object_distance, float off_axis, std::string model);
        void SpotPlot(int option, float object_distance, float off_axis, std::string model);

    public:
        
        Analysis();
        ~Analysis();


        void SimplePlot(float object_distance, float off_axis, 
            std::string model);
        void AccommodationAnalysis(float object_distance, float off_axis, 
            std::string model);
        void LSAanalysis(float object_distance, float off_axis,
                std::string model);

    };
}