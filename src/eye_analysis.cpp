#include "SchematicEye.hh" 


void Eye::Intensity(int option = 1, float object_distance=10000, 
                    float off_axis=5)
{
    int init_diop = 0;
    std::ofstream outputfile;
    switch (option)
    { 
        case 1: // not best focus
        {
            outputfile.open ("dat/Intensity.csv", std::ios::trunc);
            // headers:
            outputfile << "encircled intensity" << "," << "radius (mm)"
                        << "," << "lens focus (D)" << "pupil size (mm)"<< std::endl;
            for (int i = 0; i < 4; i++)
            {    
                
                set_params(init_diop + (i * 2.0), pupil_size, model, age);
                SchematicEye();
                EyeTracer(object_distance, off_axis);

                sys->get_tracer_params().set_default_distribution(
                                            Trace::Distribution(Trace::HexaPolarDist, 100)); 
                                
                Analysis::Spot spot(*sys);  

                double    radius = 0.0; // in mm.
                while ( radius < 1.0)
                {
                    outputfile << spot.get_encircled_intensity( radius ) << "," << radius
                        << "," << init_diop + (i * 2.0) << pupil_size << std::endl;
                    
                    radius += 0.005;
                }
            }
        }
        case 2: // best focus
        {
            outputfile.open ("dat/IntensityBestFocus.csv", std::ios::trunc);

            for (int i = 0; i < 4; i++)
            {    

                Analysis::Spot spot(*sys);
                Analysis::Focus focus(*sys);
                image->set_plane(focus.get_best_focus());            

                double    radius = 0.0; // in mm.
                while ( radius < 1.0)
                {
                    outputfile << spot.get_encircled_intensity( radius ) << "," << radius
                        << "," << init_diop + (i * 2.0) <<std::endl;
                    
                    radius += 0.005;
                }
            }
        }
        outputfile.close();
           
    }
}

void Eye::EyePlots(int option = 1, float object_distance=10000, 
                    float off_axis=5)
{    

            
    switch (option) 
    {
        case 1:
        {
            Io::RendererSvg renderer("img/eye.svg", 1200,1200);
            renderer.set_margin_ratio(0.35, 0.25, 0.1, 0.1);
    
            const int number (1);
            // layout plot
            

            renderer.set_page_layout(1,2);
            
            // draw 2d system layout: 
            sys->draw_2d_fit(renderer);
            sys->draw_2d(renderer);
            for (int i = 0; i < number; i++)
            {
                // trace and draw rays from rays source
                sys->enable_single<Sys::Source>(*source_rays);
                tracer->get_trace_result().set_generated_save_state(*source_rays);
                
                renderer.set_page(0);
                tracer->trace();
                tracer->get_trace_result().draw_2d(renderer);
                
                // longitudinal aberration
                fan = new Analysis::RayFan(*sys);
                
                sys->enable_single<Sys::Source>(*source_point);
                ref<Data::Plot> abber_plot = fan->get_plot(Analysis::RayFan::EntranceHeight,
                                                    Analysis::RayFan::LongitudinalDistance);
                
                renderer.set_page(1);
                abber_plot->draw(renderer);

            }
                break;
        }
        case 2:
        {
            Io::RendererSvg renderer("img/eye.svg", 1200,800);
            renderer.set_margin_ratio(0.1, 0.1, 0.1, 0.1);
            renderer.set_page_layout(1,4);
            for (int i = 0; i < 4; i++)
            {    
                
                set_params(init_diop + (i * 2.0), pupil_size, model, age);
                SchematicEye();
                EyeTracer(object_distance, off_axis);
                //EyePlots();
                
                std::cout << "lens (diopters): " << init_diop + (i * 2.0) << "  Age (years): " 
                            << age << "  pupil size (mm):" << pupil_size << std::endl;
                            
                std::cout << "unaccommodated defocus: " << FindOpticalPower(2)
                            << "  best focus (diopter): " << FindOpticalPower(1) << std::endl;
                
                renderer.set_page(i);
                
                sys->enable_single<Sys::Source>(*source_rays);
                tracer->get_trace_result().set_generated_save_state(*source_rays);
                sys->draw_2d_fit(renderer);
                sys->draw_2d(renderer);                
                tracer->trace();
                tracer->get_trace_result().draw_2d(renderer);
                
                
            }
        }    
        
    }
}

void Eye::SpotPlot(int option = 1, float object_distance=10000, 
                    float off_axis=5)
{
    switch (option)
    { 
        case 1:
        {
            const int number (1);
            //  change position of light slightly for a series of plots.
            Io::RendererSvg     renderer("img/spot.svg",  300 * 1, 300 * number, Io::rgb_black);   
            
            sys->get_tracer_params().set_default_distribution(
                                                Trace::Distribution(Trace::HexaPolarDist, 200)); 
            renderer.set_margin_ratio(0.1, 0.1, 0.1, 0.1);
            renderer.set_page_layout(1, number);
            
            for (int i = 0; i < number; i++)
            {
                Analysis::Spot spot(*sys);
                
                renderer.set_page(i);
                spot.draw_diagram(renderer);
                
                //source_point.rotate(0, 0.10, 0);
            }
            break;
        }
        case 2:
        {
            Io::RendererSvg renderer1("img/spot.svg", 300,1200, Io::rgb_black);
            Io::RendererSvg renderer2("img/spot_intensity.svg", 640, 480*4);
            Io::RendererSvg renderer3("img/spotBestFocus.svg", 300,1200, Io::rgb_black);
            Io::RendererSvg renderer4("img/spot_intensityBestFocus.svg", 640, 480*4);

            renderer1.set_margin_ratio(0.1, 0.1, 0.1, 0.1);
            renderer2.set_margin_ratio(0.1, 0.1, 0.1, 0.1);
            renderer3.set_margin_ratio(0.1, 0.1, 0.1, 0.1);
            renderer4.set_margin_ratio(0.1, 0.1, 0.1, 0.1);

            renderer1.set_page_layout(1,4);
            renderer2.set_page_layout(1,4);
            renderer3.set_page_layout(1,4);
            renderer4.set_page_layout(1,4);

            int init_diop = 0;
            for (int i = 0; i < 4; i++)
            {    
                
                set_params(init_diop + (i * 2.0), pupil_size, model, age);
                SchematicEye();
                EyeTracer(object_distance, off_axis);
                EyePlots(1);
                
                renderer1.set_page(i);
                renderer2.set_page(i);
                sys->get_tracer_params().set_default_distribution(
                                            Trace::Distribution(Trace::HexaPolarDist, 100)); 
                
                
                Analysis::Spot spot1(*sys);
                
                spot1.draw_diagram(renderer1);    

                ref<Data::Plot> plot1 = spot1.get_encircled_intensity_plot(100);

                plot1->draw(renderer2);                

        		// now repeat measurements at point of best focus.
                
                renderer3.set_page(i);
                renderer4.set_page(i);
                
                Analysis::Spot spot2(*sys);
                Analysis::Focus focus(*sys);
                image->set_plane(focus.get_best_focus());

                spot2.draw_diagram(renderer3);    

                ref<Data::Plot> plot2 = spot2.get_encircled_intensity_plot(100);

                plot2->draw(renderer4);                


            }
        }    
    }
}

void Eye::SimplePlot(float object_distance, float off_axis,
            std::string model)
{
    set_params(model);
    SchematicEye();
    EyeTracer(object_distance, off_axis);
    EyePlots(1, object_distance, off_axis); 
    Diopters(true);

}

void Eye::AccommodationAnalysis(float object_distance, float off_axis, 
            std::string model)
{
    std::cout << "starting values... " << std::endl;
    std::cout << "  "  << std::endl;
    
    set_params(model);
    EyePlots(2, object_distance, off_axis);
    Intensity(1, object_distance, off_axis);
}

void Eye::LSAanalysis(float object_distance, float off_axis,
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
                
                Eye::Eye     eye;
                eye.set_params( LensAccomm, PupilSize, model, AGE );
                eye.SchematicEye();
                eye.EyeTracer(object_distance, off_axis);
                
                AccommOptPower = eye.FindOpticalPower(1);
                RelaxedOptPower = eye.FindOpticalPower(2);
                Defocus_diopters = eye.Diopters(print=true);
                
                outputfile << LensAccomm << "," << PupilSize << "," << AGE << "," << 
                            AccommOptPower << "," << RelaxedOptPower << "," << 
                            Defocus_diopters << std::endl;
            }
        }
    }
    
    outputfile.close();
}
