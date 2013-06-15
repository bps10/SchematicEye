#include <math.h>

#include <algorithm>
#include <string>
#include <iostream>
#include <fstream>

#include <Goptical/Math/Vector>
#include <Goptical/Math/VectorPair>

#include <Goptical/Material/Abbe>
#include <Goptical/Material/Base>

#include <Goptical/Data/Set>
#include <Goptical/Data/PlotData>
#include <Goptical/Data/DiscreteSet>
#include <Goptical/Data/Plot>

#include <Goptical/Sys/System>
#include <Goptical/Sys/OpticalSurface>
#include <Goptical/Sys/Source>
#include <Goptical/Sys/SourceRays>
#include <Goptical/Sys/SourcePoint>
#include <Goptical/Sys/Image>
#include <Goptical/Sys/Stop>

#include <Goptical/Curve/Sphere>
#include <Goptical/Curve/Conic>
#include <Goptical/Shape/Disk>

#include <Goptical/Trace/Tracer>
#include <Goptical/Trace/Result>
#include <Goptical/Trace/Distribution>
#include <Goptical/Trace/Sequence>
#include <Goptical/Trace/Params>

#include <Goptical/Light/SpectralLine>

#include <Goptical/Analysis/RayFan>
#include <Goptical/Analysis/Spot>
#include <Goptical/Analysis/Focus>


#include <Goptical/Io/RendererSvg>
#include <Goptical/Io/Rgb>

#define PI 4.0 * atan(1.0)

using namespace Goptical;

class Eye  {
    
private:    
    float diopters, pupil_size, age, init_diop;
    int print;

    std::string  model;
    
    Sys::System * sys;
    Trace::Tracer * tracer;
    
    Sys::SourceRays  * source_rays;
    Sys::SourcePoint  *  source_point; 
    
    Sys::OpticalSurface * anterior_cornea;
    Sys::OpticalSurface * posterior_cornea;
    Sys::Stop * pupil;
    Sys::OpticalSurface * anterior_lens;
    Sys::OpticalSurface * posterior_lens;
    Sys::Image * image;
    
    Curve::Conic * anterior_cornea_curve;
    Curve::Conic * posterior_cornea_curve;
    Curve::Conic * anterior_lens_curve;
    Curve::Conic * posterior_lens_curve;
    Shape::Disk * ant_cornea_shape;
    Shape::Disk * post_cornea_shape;
    Shape::Disk * lens_shape;
    Shape::Disk * EyeShape;
    Curve::Sphere * EyeCurve;
    
    Analysis::RayFan * fan;
    
public:

    Eye();
    ~Eye();
    void set_params(int option, std::string mod);
    void set_params(float diop, float pup, std::string mod, float a);
    void SchematicEye( );
    void EyeTracer(float object_distance, float offaxis);
    void Intensity(int option, float object_distance, float off_axis);
    void EyePlots(int option, float object_distance, float off_axis);
    void SpotPlot(int option, float object_distance, float off_axis);
    
    float GetCornealThickness(std::string model);
    float GetAnteriorChamber(std::string model, float age, float diopters);
    float GetLensThickness(std::string model, float age, float diopters);
    float GetAxialLength(std::string model, float age, float diopters);
    float GetVitreousLen(std::string model);
    float DegreesToMM(float object_distance, float degrees);
    float FindOpticalPower(int opt);
    float Diopters(int option);
}; 
