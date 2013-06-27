//
//  SchematicEye.cpp
//
//  Created by Brian Schmidt on 10/8/12.
//

#include "eye_eye.hh"

/*
1
- Add additional options: # of iterations, wavelength
2
- Iterative solve for best diopter accommodation of lens.
- Add spectacle lens option - introduce chromatic analysis.
    + can define wavelengths in Goptical: add_spectral_line.
    + need to figure out how to create a lens with specific wavelength pass.
    + take mean of wavelengths during analysis.
3
- double check navarro: refractive indices need updating.
- create a comprehensive analysis option - schematic, abberations, psf, etc.
4
- Work out GRIN model.
- Allow sagittal and tang MTFs by specifiying geometry of ray tracer dist.
5
- Off axis - change location of point source for marginal rays.
- GUI - web based. Big issue is getting everything installed on AWS.
- add verbose option and single letter options (e.g. -h = help)
*/

namespace SchematicEye {

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
    // Define system parameters using Dubbelman 2005 or Navarro 1985/1999:
    //**********************************************************************
    
    // make sure model is lower case:
    std::transform(model.begin(), model.end(), model.begin(), ::tolower);
    
    float corneal_thickness, anterior_chamber, lens_thickness, 
    axial_length, pupil_rad, cornea_ant_k, cornea_post_k, 
    cornea_radius_ant, cornea_radius_post, lens_ant_k, lens_post_k,
    lens_ant_radius, lens_post_radius, vitreous_length, lens_refractive_index;
    
    ref<Material::AbbeVd> cornea_refract, extraOcular_refract,
    lens_refract, intraOcular_refract;

    corneal_thickness = GetCornealThickness(model);
    anterior_chamber = GetAnteriorChamber(model, age, diopters);
    lens_thickness = GetLensThickness(model, age, diopters);
    vitreous_length = GetVitreousLen(model);
    axial_length = GetAxialLength(model, age, diopters);
    pupil_rad = pupil_size / 2.0;
    lens_refractive_index = GetLensRefractiveIndex(model, age, diopters);

    if (model == "dubbelman") 
    {
        cornea_refract = ref<Material::AbbeVd>::create(1.376, 56.50);
        extraOcular_refract = ref<Material::AbbeVd>::create(1.336, 49.61);
        lens_refract = ref<Material::AbbeVd>::create(lens_refractive_index, 48.00);
        intraOcular_refract = ref<Material::AbbeVd>::create(1.336, 50.90);

        cornea_ant_k = 0.82;          // Schwarzschild constant (k).
        cornea_post_k = 0.66;         // Schwarzschild constant (k).
        cornea_radius_ant = 7.87;     // radius of curvature (c).
        cornea_radius_post = 6.40;    // radius of curvature (c).
        
        lens_ant_k = -(4.0 - ( 0.5 * diopters ) );
        lens_post_k = -3.0;
        lens_ant_radius =   1.0 / ( 1.0 / (12.7 - 0.058 * age ) + 
                (0.0077 * diopters) ); 
        lens_post_radius = -1.0 / ( 1.0 / (5.9 -  0.013 * age ) + 
                (0.0043 * diopters) ); 
    }
    
    if (model == "navarro")
    { 
        // refract indices from Navarro 1985 for lambda = 656.3nm
        cornea_refract = ref<Material::AbbeVd>::create(1.37405, 56.50);
        extraOcular_refract = ref<Material::AbbeVd>::create(1.3354, 49.61);
        lens_refract = ref<Material::AbbeVd>::create(lens_refractive_index,
                                                     48.00);
        intraOcular_refract = ref<Material::AbbeVd>::create(1.3407, 50.90);
    
        cornea_ant_k = -0.26;            // Schwarzschild constant (k).
        cornea_post_k = 0;               // Schwarzschild constant (k).
        cornea_radius_ant = 7.72;        // radius of curvature (c).
        cornea_radius_post = 6.50;       // radius of curvature (c).
        
        lens_ant_k = -3.1316 - ( 0.34 * log(diopters + 1.0) ); 
        lens_post_k = -1.0 - ( 0.125 * log(diopters + 1.0) );            
        lens_ant_radius =   10.2 - ( 1.75 * log(diopters + 1.0) ); 
        lens_post_radius = -6.0 + ( 0.2294 * log(diopters + 1.0) ); 
    }
    
    if (print == 1)
    {
        _simplePrint("  ");

        _simplePrint("pupil size (mm): ", pupil_size);
        _simplePrint("corneal thickness (mm): ", corneal_thickness);
        _simplePrint("anterior chamber depth (mm): ", anterior_chamber);
        _simplePrint("lens thickness (mm): ", lens_thickness);
        _simplePrint("vitreous length (mm): ", vitreous_length);
        _simplePrint("axial length (mm): ", axial_length);

        _simplePrint("cornea");
        _simplePrint("anterior surface, k (mm): ", cornea_ant_k);
        _simplePrint("posterior surface, k (mm): ", cornea_post_k);
        _simplePrint("anterior radiuc, c (mm): ", cornea_radius_ant);
        _simplePrint("posterior radius, c (mm): " , cornea_radius_post);

        _simplePrint("lens");
        _simplePrint("anterior surface, k (mm): ", lens_ant_k);
        _simplePrint("posterior surface, k (mm): ", lens_post_k);
        _simplePrint("anterior radius, c (mm): ", lens_ant_radius);
        _simplePrint("posterior radius, c (mm): ", lens_post_radius);
    }
    
    //**********************************************************************
    // Cornea: 
    //**********************************************************************
    
    // update cornea parameters
    // radius of curvature (c), Schwarzschild constant (k).
    anterior_cornea_curve = new Curve::Conic(cornea_radius_ant,  
        cornea_ant_k);
    posterior_cornea_curve = new Curve::Conic(cornea_radius_post, 
        cornea_post_k);
    ant_cornea_shape = new Shape::Disk(5.8);
    post_cornea_shape = new Shape::Disk(4.9); // cornea diameter in mm
    
                                                        // position.
    anterior_cornea = new Sys::OpticalSurface(Math::Vector3(0, 0, 0),   
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
    // make sure pupil radius is at least as large as cornea.
    pupil->set_external_radius( 6.0); 
    
    //**********************************************************************
    // Crystalline lens
    //**********************************************************************
    
    // update lens parameters
    anterior_lens_curve = new Curve::Conic(lens_ant_radius, lens_ant_k);     
    posterior_lens_curve = new Curve::Conic(lens_post_radius, lens_post_k);
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
    
    image = new Sys::Image(Math::Vector3(0, 0, axial_length), 
        *EyeCurve, *EyeShape); 
    
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

    // add wavelengths of light
    source_point->clear_spectrum();
    //source_point->add_spectral_line(Light::SpectralLine::e);
    source_point->add_spectral_line(Light::SpectralLine::C);
    //source_point->add_spectral_line(Light::SpectralLine::F);
    
    // ray tracer
    tracer = new Trace::Tracer(*sys);
    
}

float Eye::GetLensRefractiveIndex(std::string model, float age, float diopters)
{

    float lens_ref_ind;
    if (model == "dubbelman") 
    { 
            lens_ref_ind = 1.441 - 0.00039 * age + 0.0013 * diopters; 
    }
    if (model == "navarro") 
    { 
        lens_ref_ind = 1.42 + ( 9.0 * pow(10, -5) * (10.0 * diopters + pow(diopters, 2)) ); 
    }
    return lens_ref_ind;
}

float Eye::GetCornealThickness(std::string model)
{    
    float corneal_thickness;
    if (model == "dubbelman") { corneal_thickness = 0.574; }
    if (model == "navarro") { corneal_thickness = 0.55; }
    return corneal_thickness; 
}


float Eye::GetAnteriorChamber(std::string model, float age, float diopters)
{
    float anterior_chamber;
    if (model == "dubbelman") 
    // in one paper is written as 3.87 !PLUS! (0.010 ...
    {
        anterior_chamber = 3.87 - ( 0.010  * age ) - ( 
            diopters * ( 0.048 - 0.0004 * age) );
    }
    if (model == "navarro") 
    {
        anterior_chamber = 3.05 - ( 0.05 * log(diopters + 1.0) );
    }
    
    return anterior_chamber; 
}


float Eye::GetLensThickness(std::string model, float age, float diopters)
{
    float lens_thickness;
    if (model == "dubbelman")
    {
        lens_thickness = 2.93 + ( 0.0236 * age ) + ( diopters * 
            ( 0.058 - 0.0005 * age));
    }
    if (model == "navarro")
    {
        lens_thickness = 4.0 + ( 0.1 * log(diopters + 1.0) );
    }
    return lens_thickness;
}


float Eye::GetAxialLength(std::string model, float age, float diopters)
{    
    float a, b, c, d;
    a = GetCornealThickness(model);
    b = GetAnteriorChamber(model, age, diopters);
    c = GetLensThickness(model, age, diopters);
    d = GetVitreousLen(model);
    return a + b + c + d;
}

float Eye::ReturnAxialLength()
{
    return GetAxialLength(model, age, diopters);
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

}