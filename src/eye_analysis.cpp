//#include "eye_eye.hh"
#include "eye_analysis.hh" 

namespace SchematicEye {

Analysis::Analysis() {};

Analysis::~Analysis() {};


void Analysis::Intensity(int best_focus = 1, float object_distance=10000, 
                    float off_axis=5, std::string model="dubbleman")
{
    int _diop;
    float diopters, pupil_size, age, init_diop;
    std::ofstream outputfile;

    switch (best_focus)
    { 
        case 1: // not best focus
        {
            diopters = _get_input("lens diopters (mm)");
            pupil_size = _get_input("pupil size (mm)");
            age = _get_input("age (years)");
            _diop = diopters;

            // setup output files and headers:
            outputfile.open ("dat/Intensity.csv", std::ios::trunc);
            outputfile << "encircled intensity" << "," 
                        << "radius (mm)" << "," 
                        << "lens focus (D)" << "," 
                        << "pupil size (mm)" << ","
                        << "offaxis (deg)" << std::endl;
            for (int i = 0; i < 4; i++)
            {    
                Eye eye = Eye();  
                eye.set_params(_diop, pupil_size, model, age);
                eye.SchematicEye();
                eye.EyeTracer(object_distance, off_axis);

                std::cout << "lens (diopters): " << _diop << "  Age (years): " 
                            << eye.age << "  pupil size (mm):" << eye.pupil_size << std::endl;
                            
                std::cout << "unaccommodated defocus: " << eye.FindOpticalPower(2)
                            << "  best focus (diopter): " << eye.FindOpticalPower(1) << std::endl;
                            
                eye.sys->get_tracer_params().set_default_distribution(
                                            Trace::Distribution(Trace::HexaPolarDist, 200)); 
                                
                Goptical::Analysis::Spot spot(*eye.sys);  

                double    radius = 0.005; // in mm.
                while ( radius < 1.01)
                {
                    // use .get_ray_wavelen_set() to process spot for each wavelength, then average
                    outputfile << spot.get_encircled_intensity( radius ) 
                        << "," << radius
                        << "," << _diop 
                        << "," << eye.ReturnPupilSize() 
                        << "," << off_axis << std::endl;
                    
                    radius += 0.005;
                }
                //delete &eye;
                _diop += 2;
            }
            break;
        }
        case 2: // best focus
        {
            outputfile.open ("dat/IntensityBestFocus.csv", std::ios::trunc);

            for (int i = 0; i < 4; i++)
            {    
                Eye eye = Eye(); 
                eye.set_params(diopters, pupil_size, model, age);
                eye.SchematicEye();
                eye.EyeTracer(object_distance, off_axis);

                Goptical::Analysis::Spot spot(*eye.sys);
                Goptical::Analysis::Focus focus(*eye.sys);
                eye.image->set_plane(focus.get_best_focus());            

                double    radius = 0.0; // in mm.
                while ( radius < 1.0)
                {
                    outputfile << spot.get_encircled_intensity( radius ) 
                        << "," << radius
                        << "," << _diop 
                        << "," << eye.ReturnPupilSize() << std::endl;
                    
                    radius += 0.005;
                }
                //delete &eye;
                _diop += 2;
            }
        }
        outputfile.close();
           
    }
}

void Analysis::EyePlots(int best_focus = 1, float object_distance=10000, 
                    float off_axis=5, std::string model="dubbleman")
{    

    switch (best_focus) 
    {
        case 1:
        {
            Eye eye = Eye(); 
            eye.set_params(model);
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
            
            /* use below to get out actual spot data - need to create new Result() */
            //Trace::Result &result = eye.tracer->get_trace_result();
            //std::cout << result.get_max_ray_intensity() << std::endl;
            //std::cout << *eye.image << std::endl;
            //Goptical::Sys::Image &s = *eye.image;
            //std::cout << result.get_intercepted(s).get_intercept_point().project_xy() << std::endl;
            //ReturnIntercepts(eye.tracer->get_trace_result(), s);

            renderer.set_page(0);
            eye.tracer->trace();
            eye.tracer->get_trace_result().draw_2d(renderer, false); //, eye.image);
            
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

void Analysis::ReturnIntercepts(const Goptical::Trace::Result &result, const Goptical::Sys::Surface &s)
{
    _max_intensity = result.get_max_ray_intensity();
    //std::cout <<  << std::endl;
    
    //_FOREACH(i, result.get_intercepted(s))
    //{

    //}

    //for(typeof(Sys::Surface i = (result.get_intercepted(s)).begin(); i != (result.get_intercepted(s)).end(); i++)
     //   {
            
            //const Trace::Ray &ray = **i;
            // dont need global transform here, draw ray intercept points in
            // surface local coordinates.
            //ray.get_intercept_point().project_xy();
            //std::cout << ray.get_intercept_point().project_xy();
           
     //   }
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


void Analysis::SpotPlot(int option = 1, float object_distance=10000, 
                    float off_axis=5, std::string model="dubbleman")
{
    Eye eye = Eye(); 
    float diopters, pupil_size, age;
    diopters = _get_input("lens diopters (mm)");
    pupil_size = _get_input("pupil size (mm)");
    age = _get_input("age (years)");

    switch (option)
    { 
        case 1:
        {
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

            break;
        }
        case 2:
        {
            Goptical::Io::RendererSvg renderer1("img/spot.svg", 300,1200, Io::rgb_black);
            Goptical::Io::RendererSvg renderer2("img/spot_intensity.svg", 640, 480*4);
            Goptical::Io::RendererSvg renderer3("img/spotBestFocus.svg", 300,1200, Io::rgb_black);
            Goptical::Io::RendererSvg renderer4("img/spot_intensityBestFocus.svg", 640, 480*4);

            renderer1.set_margin_ratio(0.1, 0.1, 0.1, 0.1);
            renderer2.set_margin_ratio(0.1, 0.1, 0.1, 0.1);
            renderer3.set_margin_ratio(0.1, 0.1, 0.1, 0.1);
            renderer4.set_margin_ratio(0.1, 0.1, 0.1, 0.1);

            renderer1.set_page_layout(1,4);
            renderer2.set_page_layout(1,4);
            renderer3.set_page_layout(1,4);
            renderer4.set_page_layout(1,4);

            int _diop = 0;
            for (int i = 0; i < 4; i++)
            {    
                
                eye.set_params(_diop, eye.pupil_size, eye.model, eye.age);
                eye.SchematicEye();
                eye.EyeTracer(object_distance, off_axis);
                EyePlots(1);
                
                renderer1.set_page(i);
                renderer2.set_page(i);
                eye.sys->get_tracer_params().set_default_distribution(
                                    Trace::Distribution(Trace::HexaPolarDist, 100)); 
                
                
                Goptical::Analysis::Spot spot1(*eye.sys);
                
                spot1.draw_diagram(renderer1);    

                ref<Data::Plot> plot1 = spot1.get_encircled_intensity_plot(100);

                plot1->draw(renderer2);                

                // now repeat measurements at point of best focus.
                
                renderer3.set_page(i);
                renderer4.set_page(i);
                
                Goptical::Analysis::Spot spot2(*eye.sys);
                Goptical::Analysis::Focus focus(*eye.sys);
                eye.image->set_plane(focus.get_best_focus());

                spot2.draw_diagram(renderer3);    

                ref<Data::Plot> plot2 = spot2.get_encircled_intensity_plot(100);

                plot2->draw(renderer4);                
                _diop += 2;

            }
        }    
    }
}


void Analysis::SimplePlot(float object_distance, float off_axis,
            std::string model)
{
    //Eye::Eye eye;
    //set_params(model);
    //SchematicEye();
    //EyeTracer(object_distance, off_axis);
    EyePlots(1, object_distance, off_axis, model); 
    //Diopters(true);

}

void Analysis::AccommodationAnalysis(float object_distance, float off_axis, 
            std::string model)
{
    std::cout << "starting values... " << std::endl;
    std::cout << "  "  << std::endl;
    
    //set_params(model);
    //EyePlots(1, object_distance, off_axis, model);
    Intensity(1, object_distance, off_axis, model);
}


}