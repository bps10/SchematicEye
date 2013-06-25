#include "eye_analysis.hh" 

namespace SchematicEye {

Analysis::Analysis() {};

Analysis::~Analysis() {};


float *Analysis::IntensityAnalysis(float object_distance=1000000, float off_axis=5, 
        float pupil_size=3, float diopters=0)
{
    std::string model="navarro";
    std::string param="foo";
    int _diop = diopters;

    Eye eye = Eye();  
    eye.set_params(_diop, pupil_size, model, 20);
    eye.SchematicEye();
    eye.EyeTracer(object_distance, off_axis);
                
    eye.sys->get_tracer_params().set_default_distribution(
                                Trace::Distribution(Trace::HexaPolarDist, 250)); 
                    
    Goptical::Analysis::Spot spot(*eye.sys);  

    double maxRadius = 0.2;
    double radiusStep = 0.0005;
    int len = (maxRadius / radiusStep) - 1;
    float *outfile = new float [len];
    int i(0);
    double radius = 0.0005; // in mm.
    while ( radius <= maxRadius)
    {
        // use .get_ray_wavelen_set() to process spot for each wavelength? 
        // then wighted average?
        outfile[i] = spot.get_encircled_intensity( radius );
        i++;
        radius += radiusStep;
    }
    return outfile;
}

void Analysis::IntensityAnalysis(std::string param, float object_distance=1000000, 
                    float off_axis=5, std::string model="dubbleman", 
                    float age=20, float pupil_size=3, float diopters=0, int iter=4)
{
    int _diop = diopters;
    std::ofstream outputfile;

    // add ability to change name of Intensity
    // setup output files and headers:
    outputfile.open ("dat/Intensity.csv", std::ios::trunc);
    outputfile << "encircled_intensity" << "," 
                << "radius_mm" << "," 
                << "lens_focus_D" << "," 
                << "pupil_size_mm" << ","
                << "offaxis_deg" << ","
                << "age_y" << ","
                << "obj_distance_mm" << ","
                << "model" << ","
                << "axial_len_mm" << ","
                << "iterations" << std::endl;

    for (int i = 0; i < iter; i++)
    {    
        Eye eye = Eye();  
        eye.set_params(_diop, pupil_size, model, age);
        eye.SchematicEye();
        eye.EyeTracer(object_distance, off_axis);

        std::cout << " " <<std::endl;
        std::cout << "lens (diopters): " << _diop << "  Age (years): " 
                    << eye.age << "  pupil size (mm):" << eye.pupil_size << std::endl;

        std::cout << "off axis (deg):" << off_axis
                    << "  obj dist (mm): " << object_distance << std::endl;
                    
        std::cout << "unaccommodated defocus: " << eye.FindOpticalPower(2)
                    << "  best focus (diopter): " << eye.FindOpticalPower(1) << std::endl;
                    
        eye.sys->get_tracer_params().set_default_distribution(
                                    Trace::Distribution(Trace::HexaPolarDist, 250)); 
                        
        Goptical::Analysis::Spot spot(*eye.sys);  

        float eye_len, pup_s;
        eye_len = eye.ReturnAxialLength();
        pup_s = eye.ReturnPupilSize();
        double radius = 0.0005; // in mm.
        while ( radius <= 0.2)
        {
            // use .get_ray_wavelen_set() to process spot for each wavelength? 
            // then wighted average?
            outputfile << spot.get_encircled_intensity( radius ) 
                << "," << radius
                << "," << _diop 
                << "," << pup_s
                << "," << off_axis 
                << "," << age 
                << "," << object_distance 
                << "," << model 
                << "," << eye_len
                << "," << iter << std::endl;
            
            radius += 0.0005;
        }
        // increment the chosen variable.
        if (param == "age") { age += 2; }
        if (param == "pupil") { pupil_size++; }
        if (param == "focus") { _diop += 2; }
        if (param == "angle") { off_axis += 5; }
        if (param == "distance") { object_distance = pow(10, 4 - i); }
    }

}


void Analysis::EyePlots(int best_focus = 1, float object_distance=100000, 
                    float off_axis=5, std::string model="dubbleman", 
                    float age=20, float pupil_size=3, float diopters=0)
{    

    switch (best_focus) 
    {
        case 1:
        {
            Eye eye = Eye(); 
            eye.set_params(diopters, pupil_size, model, age);
            eye.SchematicEye();
            eye.EyeTracer(object_distance, off_axis);

            Goptical::Io::RendererSvg renderer("img/eye.svg", 1200,1200);
            renderer.set_margin_ratio(0.35, 0.25, 0.1, 0.1);

            // layout plot
            renderer.set_page_layout(1,2);
            
            // draw 2d system layout: 
            eye.sys->draw_2d_fit(renderer);
            eye.sys->draw_2d(renderer);

            // trace and draw rays from rays source
            eye.sys->enable_single<Sys::Source>(*eye.source_rays);
            eye.tracer->get_trace_result().set_generated_save_state(*eye.source_rays);
            
            renderer.set_page(0);
            eye.tracer->trace();
            eye.tracer->get_trace_result().draw_2d(renderer, true); //, eye.image);
            
            // longitudinal aberration
            fan = new Goptical::Analysis::RayFan(*eye.sys);
            
            eye.sys->enable_single<Sys::Source>(*eye.source_point);
            ref<Data::Plot> abber_plot = fan->get_plot(Goptical::Analysis::RayFan::EntranceHeight,
                                                Goptical::Analysis::RayFan::LongitudinalDistance);
            
            renderer.set_page(1);
            abber_plot->draw(renderer);

            eye.Diopters(true);
            break;
        }
        case 2:
        {
            int _diop = 0;
            Goptical::Io::RendererSvg renderer("img/eye.svg", 1200,800);
            renderer.set_margin_ratio(0.1, 0.1, 0.1, 0.1);
            renderer.set_page_layout(1,4);
            for (int i = 0; i < 4; i++)
            {    
                Eye eye = Eye(); 
                eye.set_params(model);
                eye.SchematicEye();
                eye.EyeTracer(object_distance, off_axis);
                //EyePlots();
                
                renderer.set_page(i);
                
                eye.sys->enable_single<Sys::Source>(*eye.source_rays);
                eye.tracer->get_trace_result().set_generated_save_state(*eye.source_rays);
                eye.sys->draw_2d_fit(renderer);
                eye.sys->draw_2d(renderer);                
                eye.tracer->trace();
                eye.tracer->get_trace_result().draw_2d(renderer);

                //delete &eye;
                _diop += 2;
                
            }
        }    
        
    }
}


void Analysis::LSAanalysis(float object_distance, float off_axis,
    std::string model)
{
    std::ofstream outputfile;
    outputfile.open ("dat/EyeLSA.csv", std::ios::trunc);
    float AGE, LensAccomm, PupilSize, AccommOptPower, 
            RelaxedOptPower, Defocus_diopters;

    for (int i = 0; i < 20; i++)
    {   
        for (int j = 1; j < 21; j++)
        {
            for (int k = 10; k <25; k++)
            {
                LensAccomm = i / 2.0;
                PupilSize = j / 4.0;
                AGE = k;
                
                Eye eye = Eye(); 
                eye.set_params( LensAccomm, PupilSize, model, AGE );
                eye.SchematicEye();
                eye.EyeTracer(object_distance, off_axis);
                
                AccommOptPower = eye.FindOpticalPower(1);
                RelaxedOptPower = eye.FindOpticalPower(2);
                Defocus_diopters = eye.Diopters(true);
                
                outputfile << LensAccomm << "," << PupilSize << "," << AGE << "," << 
                            AccommOptPower << "," << RelaxedOptPower << "," << 
                            Defocus_diopters << std::endl;

            }
        }
    }
    
    outputfile.close();
}


void Analysis::SpotPlot(float object_distance=100000, 
                    float off_axis=5, std::string model="dubbleman", 
                    float age=20, float pupil_size=3, float diopters=0)
{
    Eye eye = Eye(); 
    //float diopters, pupil_size, age;
    //diopters = _get_input("lens diopters (mm)");
    //pupil_size = _get_input("pupil size (mm)");
    //age = _get_input("age (years)");

    eye.set_params(diopters, pupil_size, model, age);
    eye.SchematicEye();
    eye.EyeTracer(object_distance, off_axis);

    eye.sys->enable_single<Sys::Source>(*eye.source_point);
    Goptical::Io::RendererSvg     renderer("img/spot.svg", 600,600, 
        Io::rgb_black);
    //renderer.set_margin_ratio(0.1, 0.1, 0.1, 0.1);
    
    eye.sys->get_tracer_params().set_default_distribution(
                            Trace::Distribution(Trace::HexaPolarDist, 100)); 
    //renderer.set_page_layout(1, 1);

    Goptical::Analysis::Spot spot(*eye.sys);
    
    //renderer.set_page(0);
    spot.draw_diagram(renderer);

}

void Analysis::EyePlots(float object_distance=100000, float off_axis=5, 
                    std::string model="navarro", float age=20, float pupil_size=3, 
                    float diopters=0)
{ 
    Eye eye = Eye(); 
    Goptical::Io::RendererSvg renderer("img/eye_analysis.svg", 300,1200, Io::rgb_black);

    renderer.set_margin_ratio(0.1, 0.1, 0.1, 0.1);
    renderer.set_page_layout(1,4);

    int _diop = 0;
    for (int i = 0; i < 4; i++)
    {    
        
        eye.set_params(_diop, eye.pupil_size, eye.model, eye.age);
        eye.SchematicEye();
        eye.EyeTracer(object_distance, off_axis);
        
        renderer.set_page(0);
        eye.sys->get_tracer_params().set_default_distribution(
                            Trace::Distribution(Trace::HexaPolarDist, 2550)); 
        
        
        Goptical::Analysis::Spot spot(*eye.sys);
        
        spot.draw_diagram(renderer);    

        ref<Data::Plot> plot = spot.get_encircled_intensity_plot(100);

        plot->draw(renderer);                

        renderer.set_page(1);
        renderer.set_page(2);
        
        //Goptical::Analysis::Spot spot(*eye.sys);
        Goptical::Analysis::Focus focus(*eye.sys);
        eye.image->set_plane(focus.get_best_focus());

        spot.draw_diagram(renderer);    

        //ref<Data::Plot> plot = spot.get_encircled_intensity_plot(100);

        renderer.set_page(3);
        plot->draw(renderer);                
    }    
}


arma::vec Analysis::PSF(arma::vec intensity, arma::vec xvals, bool symmetric=false)
{
    float samples;
    samples = intensity.n_elem;
    arma::vec psf = arma::zeros<arma::vec>(1, samples);
    arma::vec psftotal = arma::zeros<arma::vec>(1, samples * 2);

    // we have an integral, therefore take the deriv to get rays / bin
    arma::vec deriv = arma::zeros<arma::vec>(samples);

    deriv(0) = intensity(0);
    deriv.rows(1, samples) = intensity.rows(1, samples) - intensity.rows(0, samples - 1);

    float radius0, radius1, area;
    for (int i = 0; i < samples - 1; i++)
    {
        // account for increasing size of area
        radius0 = xvals(i);
        radius1 = xvals(i + 1);
        
        // subtract inner and outer circle area to get sliver of interest
        area = (PI * pow(radius1, 2)) - (PI * pow(radius0, 2));

        // deriv = amount in each circle; then divide by area
        psf(i) = deriv(i) / area;
    } 

    // normalize so that PSF has integral of 1.
    psf = psf / arma::sum(psf);

    psftotal.rows(1, samples + 1) = arma::fliplr(psf);
    psftotal.rows(0, samples) = psf.rows(1, samples - 1);

    if (symmetric) { return psftotal; }
    else { return psf; }
}

arma::vec Analysis::MTF(arma::vec psf)
{
    arma::vec temp;

    arma::vec mtf = arma::zeros<arma::vec>(psf.n_elem);

    temp = arma::abs(arma::fft(psf));
    temp = arma::real(temp);

    // normalize MTF
    mtf = temp / psf.max();
    return mtf;
}
}