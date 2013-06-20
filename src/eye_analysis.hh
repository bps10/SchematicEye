#include <Goptical/Data/Set>
#include <Goptical/Data/PlotData>
#include <Goptical/Data/DiscreteSet>
#include <Goptical/Data/Plot>

#include <Goptical/Analysis/RayFan>
#include <Goptical/Analysis/Spot>
#include <Goptical/Analysis/Focus>

#include <Goptical/Math/Vector>
#include <Goptical/Trace/Result>
#include <Goptical/Trace/Ray>

#include <Goptical/Sys/Surface>
#include <Goptical/Sys/Image>

#include <Goptical/Io/RendererSvg>
#include <Goptical/Io/Rgb>

#include "eye_eye.hh"


//#define _FOREACH(var, cont)             \
//    for (typeof((cont).begin()) var = (cont).begin(); var !=(cont).end(); ++var)

namespace SchematicEye 
{
    class Analysis
    {

    private:

        Goptical::Analysis::RayFan * fan;

        void Intensity(int best_focus, float object_distance, float off_axis, std::string model);
        void EyePlots(int best_focus, float object_distance, float off_axis, std::string model);
        void ReturnIntercepts(const Trace::Result &result, const Sys::Surface &s);

    public:
        
        Analysis();
        ~Analysis();

        double _max_intensity;

        void SimplePlot(float object_distance, float off_axis, 
            std::string model);
        void AccommodationAnalysis(float object_distance, float off_axis, 
            std::string model);
        void LSAanalysis(float object_distance, float off_axis,
            std::string model);
        void SpotPlot(int option, float object_distance, float off_axis, 
            std::string model);

        inline float _get_input(std::string message)
        {
            float out;
            std::cout << "enter " << message << ": " << std::endl;
            std::cin >> out;
            return out;
        }

    };
}