//
//  SchematicEye.cpp
//
//  Created by Brian Schmidt on 10/8/12.
//

#include "SchematicEye.hh"

/*
1
- MTF, PSF.
- Create object series option.
2
- Add spectacle lens option - introduce chromatic analysis.
    + can define wavelengths in Goptical: add_spectral_line.
    + need to figure out how to create a lens with specific wavelength pass.
3
- Allow user to change model parameters (dubbelman vs navarro): navarro
    needs work, see python script.
4
- Work out GRIN model.
5
- Off axis - change location of point source for marginal rays.
- GUI.
*/

Eye::Eye() {}

Eye::~Eye() {}


void Eye::set_params(std::string mod = "dubbelman")
{

    std::cout << "enter lens diopters (mm): " << std::endl;
    std::cin >> diopters;
    std::cout << "enter pupil size (mm): " << std::endl;
    std::cin >> pupil_size;
    std::cout << "enter age (years): " << std::endl;
    std::cin >> age;
    
    //age = a;
    model = mod;
    print = 1.0;
    init_diop = diopters;
}

void Eye::set_params(float diop = 0.0, float pup = 5.0, 
    std::string mod = "dubbelman", float a = 20.0)
{

    diopters = diop;
    pupil_size = pup;
    init_diop = diopters;
    age = a;
    model = mod;
    print = 0.0;
    
}


void Eye::SchematicEye()
{
    // create new system
    sys = new Sys::System();
    
    //**********************************************************************
    // set material refraction indexes. Eventually include GRIN for lens:
    //**********************************************************************
    
    ref<Material::AbbeVd> cornea_refract =      ref<Material::AbbeVd>::create(1.367, 56.50);
    
    ref<Material::AbbeVd> extraOcular_refract = ref<Material::AbbeVd>::create(1.3374, 49.61);
    
    ref<Material::AbbeVd> lens_refract =        ref<Material::AbbeVd>::create(1.42, 48.00);
    
    ref<Material::AbbeVd> intraOcular_refract = ref<Material::AbbeVd>::create(1.336, 50.90);
    
    //**********************************************************************
    // Define system parameters using Dubbelman 2005 or Navarro 1985:
    //**********************************************************************
    std::transform(model.begin(), model.end(), model.begin(), ::tolower);
    
    float corneal_thickness, anterior_chamber, lens_thickness, axial_length, 
    pupil_rad, cornea_ant_k, cornea_post_k, cornea_radius_ant, cornea_radius_post, 
    lens_ant_k, lens_post_k, lens_ant_radius, lens_post_radius, vitreous_length;
    
    corneal_thickness = GetCornealThickness(model);
    anterior_chamber = GetAnteriorChamber(model, age, diopters);
    lens_thickness = GetLensThickness(model, age, diopters);
    vitreous_length = GetVitreousLen(model);
    axial_length = GetAxialLength(model, age, diopters);
    pupil_rad = pupil_size / 2.0;
        
    if (model == "dubbelman") 
    {
        cornea_ant_k = 0.82;             // Schwarzschild constant (k).
        cornea_post_k = 0.66;            // Schwarzschild constant (k).
        cornea_radius_ant = 7.87;        // radius of curvature (c).
        cornea_radius_post = 6.40;    // radius of curvature (c).
        
        lens_ant_k = -(4.0 - ( 0.5 * diopters ) );                                
        lens_post_k = -3.0;                                                        
        lens_ant_radius =   1.0 / ( 1.0 / (12.7 - 0.058 * age ) + (0.0077 * diopters) ); 
        lens_post_radius = -1.0 / ( 1.0 / (5.9 -  0.013 * age ) + (0.0043 * diopters) ); 
    }
    
    if (model == "navarro")
    { 

        pupil_rad = pupil_size / 2.0;
        
        cornea_ant_k = -0.26;             // Schwarzschild constant (k).
        cornea_post_k = 0;                // Schwarzschild constant (k).
        cornea_radius_ant = 7.72;        // radius of curvature (c).
        cornea_radius_post = 6.50;        // radius of curvature (c).
        
        lens_ant_k = -3.1316;                                
        lens_post_k = -1.0;                                                        
        lens_ant_radius =   10.2; 
        lens_post_radius = -6.0; 
    }
    
    if (print == 1)
    {
        std::cout << "  " << std::endl;
        
        std::cout << "pupil size (mm): "              << pupil_size          << std::endl;
        std::cout << "corneal thickness (mm): "       << corneal_thickness   << std::endl;
        std::cout << "anterior chamber depth (mm): "  << anterior_chamber    << std::endl;
        std::cout << "lens thickness (mm): "          << lens_thickness      << std::endl;
        std::cout << "vitreous length (mm): "         << vitreous_length     << std::endl;
        std::cout << "axial length (mm): "            << axial_length        << std::endl;
        
        std::cout << "cornea: " << std::endl;
        std::cout << "anterior surface, k (mm): "     << cornea_ant_k        << std::endl;
        std::cout << "posterior surface, k (mm): "    << cornea_post_k       << std::endl;
        std::cout << "anterior radius, c (mm): "      << cornea_radius_ant   << std::endl;
        std::cout << "posterior radius, c (mm): "     << cornea_radius_post  << std::endl;
        
        std::cout << "lens: " << std::endl;
        std::cout << "anterior surface, k (mm): "     << lens_ant_k          << std::endl;
        std::cout << "posterior surface, k (mm): "    << lens_post_k         << std::endl;
        std::cout << "anterior radius, c (mm): "      << lens_ant_radius     << std::endl;
        std::cout << "posterior radius, c (mm): "     << lens_post_radius    << std::endl;
    }
    
    //**********************************************************************
    // Cornea: 
    //**********************************************************************
    
    // update cornea parameters
    // radius of curvature (c), Schwarzschild constant (k).
    anterior_cornea_curve = new Curve::Conic( cornea_radius_ant,  cornea_ant_k);  
    posterior_cornea_curve = new Curve::Conic(cornea_radius_post, cornea_post_k); 
    ant_cornea_shape = new Shape::Disk(5.8);
    post_cornea_shape = new Shape::Disk(4.9); // cornea diameter in mm
    
    anterior_cornea = new Sys::OpticalSurface(Math::Vector3(0, 0, 0),   // position.
                                              *anterior_cornea_curve,   // curve.
                                              *ant_cornea_shape,        // aperture shape.
                                              Material::none,           // material to left.
                                              cornea_refract);          // material to right.
    
    posterior_cornea = new Sys::OpticalSurface(Math::Vector3(0, 0, 
                                               corneal_thickness), // position.
                                               *posterior_cornea_curve, // curve.
                                               *post_cornea_shape,      // aperature shape.
                                               cornea_refract,          // material to left.
                                               extraOcular_refract);    // material to right.
    
    
    //**********************************************************************
    // add pupil (set location and radius (mm) of pupil (default = 1.5mm).
    //**********************************************************************
    
    pupil = new Sys::Stop(Math::Vector3(0, 0, 
                        corneal_thickness + anterior_chamber),  pupil_rad);      
    pupil->set_external_radius( 6.0); // make sure pupil radius is at least as large as cornea.
    
    //**********************************************************************
    // Crystalline lens
    //**********************************************************************
    
    // update lens parameters
    anterior_lens_curve = new Curve::Conic(    lens_ant_radius,     lens_ant_k);     
    posterior_lens_curve = new Curve::Conic(    lens_post_radius,     lens_post_k);             
    lens_shape = new Shape::Disk(6.0); // lens diameter in mm
    
    anterior_lens = new Sys::OpticalSurface(Math::Vector3(0, 0, 
                        corneal_thickness + anterior_chamber), 
                                            *anterior_lens_curve,   // curve.
                                            *lens_shape,            // aperture shape.
                                            extraOcular_refract,    // material to left.
                                            lens_refract);          // material to right.
    // position.
    posterior_lens = new Sys::OpticalSurface(Math::Vector3(0, 0, 
        corneal_thickness + anterior_chamber + lens_thickness),     
                                             *posterior_lens_curve, // curve.
                                             *lens_shape,           // aperture shape.
                                             lens_refract,          // material to left.
                                             intraOcular_refract);  // material to right.
    
    //**********************************************************************
    // add all of the optical components.
    //**********************************************************************
    
    sys->add(*anterior_cornea);
    sys->add(*posterior_cornea);
    sys->add(*pupil);
    sys->add(*anterior_lens);
    sys->add(*posterior_lens);
    
    //**********************************************************************
    // set the eye shape:
    //**********************************************************************
    
    EyeCurve = new Curve::Sphere(-12.0); 
    EyeShape = new Shape::Disk(10.0); // eye radius (mm)
    
    image = new Sys::Image(Math::Vector3(0, 0, axial_length), *EyeCurve, *EyeShape); 
    
    sys->add(*image);
    sys->set_entrance_pupil(*anterior_cornea);
    
}


void Eye::EyeTracer(float object_distance=10000, float offaxis=0)
{
    //**********************************************************************
    // Setup light sources and ray tracer
    //**********************************************************************
    
    float object_dist, object_offaxis;
    object_dist = -1 * abs(object_distance); // make sure negative
    object_offaxis = DegreesToMM(object_distance, offaxis);    

    source_rays  = new Sys::SourceRays(
                Math::Vector3(0, object_offaxis, object_dist));
    //source_point = new Sys::SourcePoint(Sys::SourceAtInfinity, 
    //                                    Math::vector3_001);  
    source_point = new Sys::SourcePoint(Sys::SourceAtFiniteDistance, 
                Math::Vector3(0, object_offaxis, object_dist));
    
    // add sources to system
    sys->add(*source_rays);
    sys->add(*source_point);
    
    // configure sources
    source_rays->add_chief_rays(*sys);
    source_rays->add_marginal_rays(*sys, pupil_size / 2.0);
    source_rays->add_marginal_rays(*sys, -pupil_size / 2.0);
    
    source_point->clear_spectrum();

    // add wavelengths of light
    source_point->add_spectral_line(Light::SpectralLine::e);
    
    // ray tracer
    tracer = new Trace::Tracer(*sys);
    
}


float Eye::GetCornealThickness(std::string model)
{    
    float corneal_thickness;
    if (model == "dubbelman")    {corneal_thickness = 0.574;}
    if (model == "navarro")     {corneal_thickness = 0.55;}
    return corneal_thickness; 
}


float Eye::GetAnteriorChamber(std::string model, float age, float diopters)
{
    float anterior_chamber;
    if (model == "dubbelman") // in one paper is written as 3.87 !PLUS! (0.010 ...
    {anterior_chamber = 3.87 - ( 0.010  * age ) - ( diopters * ( 0.048 - 0.0004 * age) );}
    if (model == "navarro") {anterior_chamber = 3.05;}
    
    return anterior_chamber; 
}


float Eye::GetLensThickness(std::string model, float age, float diopters)
{
    float lens_thickness;
    if (model == "dubbelman")
    {lens_thickness = 2.93 + ( 0.0236 * age ) + ( diopters * ( 0.058 - 0.0005 * age));}
    if (model == "navarro") {lens_thickness = 4.0;}
    return lens_thickness;
}


float Eye::GetAxialLength(std::string model, float age, float diopters)
{    
    float a, b, c;
    a = GetAnteriorChamber(model, age, diopters);
    b = GetLensThickness(model, age, diopters);
    c = GetVitreousLen(model);
    return a + b + c;
}


float Eye::GetVitreousLen(std::string model)
{ 
    float vitreous_length;
    if (model == "dubbelman")    {vitreous_length = 16.9935;}
    if (model == "navarro")     {vitreous_length = 16.6;}
    return vitreous_length;
}

float Eye::ReturnPupilSize()
{
    return pupil_size;
}

float Eye::DegreesToMM(float object_distance, float degrees)
{
    float radians;
    radians = degrees * PI / 180;
    return abs(object_distance * tan(radians));
}


float Eye::FindOpticalPower(int opt = 1)
{
    Analysis::Focus        focus(*sys);
    float power, focal_len;
    if (opt == 1) { focal_len = focus.get_best_focus()[0][2]; }
    if (opt == 2) { focal_len = GetAxialLength(model, age, diopters); }
    
    power = 1.0 / (focal_len / 1000.0);
    return power;
}


float Eye::Diopters(bool print = true)
{
    Analysis::Focus     focus(*sys);
    float relaxed_power, accomm_power, defocus_diopter;
    relaxed_power = FindOpticalPower(2);
    accomm_power = FindOpticalPower(1);
    defocus_diopter = accomm_power - relaxed_power;
    
    if (print == true)
    {
    std::cout << " " << std::endl;
    std::cout << "focal plane: " << std::endl;
    std::cout << focus.get_best_focus()[0][2] << std::endl;
    std::cout << "defocus diopters: " << std::endl;
    std::cout << defocus_diopter << std::endl;
    return 0.0;
    }
    
    if (print == false)
    {
    return defocus_diopter;
    }
    
}


