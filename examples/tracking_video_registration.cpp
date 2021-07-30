/*
* @Author: Jose Tascon
* @Date:   2019-11-18 13:30:52
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2021-07-28 22:37:45
*/

// std libs
#include <iostream>
#include <memory>
#include <sstream>
#include <iomanip>
#include <filesystem>
#include <thread>

// boost libs
#include <boost/program_options.hpp>

// local libs
#include "../src/image.h"
#include "../src/dfield.h"
#include "../src/demons.h"
#include "../src/gradient_descent.h"
#include "../src/registration.h"
#include "../src/viewer_track.h"

#include "../src/utils/timer.h"

using namespace imart;
namespace po = boost::program_options;
namespace fs = std::filesystem;

std::vector<std::string> list_directory(std::string path);

int main(int argc, char *argv[])
{
    using type = float;
    // using type = double;
    using image_type = image_cpu<type>;
    
    // if( argc < 3 )
    // {
    //     std::cerr << "Usage: " << std::endl;
    //     std::cerr << argv[0] << " input_folder output_folder" << std::endl;
    //     return EXIT_FAILURE;
    // }
    // std::string path_input = argv[1];
    // std::string path_output = argv[2];

    std::string path_input;
    std::string path_output = "";
    bool verbose, plot, write;
    double fps, tf;
    double tolerance, step;
    double sigmaf, sigmae;
    std::vector<int> level_scales({8,4,2});
    std::vector<int> level_iterations({120,100,80});

    // Program options
    po::options_description desc("Tracking 2D CineMR images with deformanble image registration. Options");
    desc.add_options()
    ("help,h", "Help message")
    ("input,i", po::value<std::string>(&path_input), "Input folder")
    ("output,o", po::value<std::string>(&path_output), "Output folder")
    ("verbose,v", po::bool_switch(&verbose), "Enable verbose")
    ("plot,p", po::bool_switch(&plot), "Enable plot")
    ("frame_per_second,f", po::value<double>(&fps)->default_value(4.0), "Images frames per second")
    ("opt-step,l", po::value<double>(&step)->default_value(1.0), "Optimizer step")
    ("opt-level-scales,s", po::value<std::vector<int> >(&level_scales)->multitoken(), "Optimizer multilevel scales")
    ("opt-level-iterations,k", po::value<std::vector<int> >(&level_iterations)->multitoken(), "Optimizer multilevel iterations")
    ("opt-tolerance,t", po::value<double>(&tolerance)->default_value(1e-6), "Optimizer tolerance")
    ("sigma-fluid,u", po::value<double>(&sigmaf)->default_value(0.0), "Deformation field sigma fluid")
    ("sigma-elastic,e", po::value<double>(&sigmae)->default_value(2.5), "Deformation field sigma elastic");

    // Parse command line options
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")){ std::cerr << desc << std::endl; return 1; };

    if (level_scales.size() != level_iterations.size()) {std::cerr << "Multilevel scales size is different to multilevel iterations" << std::endl; return 1; };

    if (vm.count("input")) std::cout << "Folder input: " << path_input << std::endl << std::endl;
    else {std::cerr << "Please provide input folder with images" << std::endl; return 1; };

    if (path_output == "") write = false;
    else write = true;

    tf = 1000/fps; // sampling time in milliseconds

    // ============================================
    //          Testing video registration
    // ============================================

    std::string path_image;
    std::vector<std::string> path_organs;

    // Read folders in directory
    auto list_folders = list_directory(path_input);

    for (auto const& folder : list_folders)
    {
        std::string relative = fs::path(folder).filename();
        if (relative == "boundary")
            continue;
        else if (relative.find("image") != std::string::npos)
            path_image = folder;
        else if ( not fs::is_directory(folder) )
            ;
        else
            path_organs.push_back(folder);
    }
    
    // Images
    auto img_fixed = image_cpu<type>::new_pointer();
    auto img_input = image_cpu<type>::new_pointer();
    std::vector<image_cpu<type>::pointer> img_organs(path_organs.size());
    std::vector<image_cpu<type>::pointer> img_organs_warped(path_organs.size());
    for (int k = 0; k < path_organs.size(); k++)
    {
        img_organs[k] = image_cpu<type>::new_pointer();
        img_organs_warped[k] = image_cpu<type>::new_pointer();
    }
    
    // Read fixed image
    auto list_files_images = list_directory(path_image);
    std::string file_fixed = list_files_images.front();
    img_fixed->read(file_fixed);
    std::cout << "Read fixed: " << file_fixed << std::endl;

    for(int k = 0; k < path_organs.size(); k++)
    {
        auto list_files_organ = list_directory(path_organs[k]);
        img_organs[k]->read(list_files_organ.front());
        std::cout << "Read " << std::string(fs::path(path_organs[k]).filename()) 
            << ": " << list_files_organ.front() << std::endl;
    }
    std::cout << std::endl;
    
    // Transform
    auto trfm = dfield<type,vector_cpu<type>>::new_pointer(img_fixed);
    trfm->set_sigma_fluid(sigmaf);
    trfm->set_sigma_elastic(sigmae);

    // Metric
    auto metric_demons = demons<type,vector_cpu<type>>::new_pointer(img_fixed, img_input, trfm);
    auto opt = gradient_descent<type,vector_cpu<type>>::new_pointer();
    opt->set_step(step);
    opt->set_tolerance(tolerance);
    opt->set_unchanged_times(10);
    opt->set_verbose(false);

    // Registration
    auto registro = registration<type,vector_cpu<type>>::new_pointer(img_fixed, img_input, trfm);
    registro->set_levels(3);
    registro->set_levels_scales(level_scales);
    registro->set_levels_iterations(level_iterations);
    registro->set_metric(metric_demons);
    registro->set_optimizer(opt);
    registro->set_normalize(true);
    registro->set_verbose(false);
    auto xi = grid_cpu<type>::new_pointer(img_fixed); // Post registration

    // Plot
    auto view = viewer_track<image_type>::new_pointer();
    auto img_view_input = img_fixed->copy();
    std::vector<image_cpu<type>::pointer> img_view_organs(path_organs.size());
    for(int k = 0; k < path_organs.size(); k++)
        img_view_organs[k] = img_organs[k]->copy();
    
    // Setup viewer
    if (plot)
    {
        view->add_image(img_view_input);
        for(int k = 0; k < path_organs.size(); k++) 
            view->add_image(img_view_organs[k]);
        view->setup();
    }

    timer t("ms");
    t.start();
    double now_time = 0.0;
    double sum_time = 0.0;

    // Output folder
    std::vector<std::string> folders_output(path_organs.size());
    if (write)
    {
        for(int k = 0; k < path_organs.size(); k++) 
        {
            std::string folder_output = path_output + "/" + std::string(fs::path(path_organs[k]).filename());
            std::cout << "Create output folder: " << folder_output << std::endl;
            fs::create_directories(folder_output);
            folders_output[k] = folder_output;
        }
        std::cout << std::endl;
    }

    // Main Loop
    for(size_t i = 1; i < list_files_images.size(); i++ )
    {
        // Read new input
        std::string file_input = list_files_images[i];
        img_input->read(file_input);
        if (verbose) std::cout << "Read input: " << file_input << std::endl;
        std::string num = fs::path(file_input).stem();
        num = num.substr(num.size()-4,num.size());
        std::string ext = fs::path(file_input).extension();
        // std::cout << "Number: " << num << std::endl;

        if (verbose)
        {
            t.lap(); 
            now_time = t.get_elapsed();
            sum_time += now_time;
            printf("Read time: \t%5.2f [ms]\n", now_time);
        }

        // Registration
        trfm = dfield<type,vector_cpu<type>>::new_pointer(img_fixed);
        trfm->set_sigma_fluid(sigmaf);
        trfm->set_sigma_elastic(sigmae);
        registro->set_moving(img_input);
        registro->set_transform(trfm);
        registro->apply();

        if (verbose)
        {
            t.lap(); 
            now_time = t.get_elapsed();
            sum_time += now_time;
            printf("Register time: \t%5.2f [ms]\n", now_time);
        }

        // Warp output images
        auto transformation = registro->get_transform()->inverse();
        for(size_t k = 0; k < img_organs.size(); k++)
        {
            auto interpol_organ = ilinear_cpu<type>::new_pointer(img_organs[k]);
            img_organs_warped[k] = interpol_organ->apply(transformation->apply(xi));
            // std::cout << k << " " << std::endl;
            std::string base = fs::path(path_organs[k]).filename();
            std::string file_output = folders_output[k] + "/" + base + "_" + num + ext;
            if (write)
            {
                std::cout << "Write output: " << file_output << std::endl;
                img_organs_warped[k]->write(file_output);
            }
        }

        if (verbose)
        {
            t.lap(); 
            now_time = t.get_elapsed();
            sum_time += now_time;
            printf("Warp time: \t%5.2f [ms]\n", now_time);
        }

        // Plot
        if (plot)
        {
            img_view_input->equal(*img_input);
            view->update(0);
            
            for(int k = 0; k < path_organs.size(); k++) 
            {
                img_view_organs[k]->equal(*img_organs_warped[k]);
                view->update(k+1);
            }
            
            view->render();

            if (verbose)
            {
                t.lap();
                now_time = t.get_elapsed();
                sum_time += now_time;
                printf("Viewer time: \t%5.2f [ms]\n", now_time);
            }
        }

        if (verbose)
        {
            t.lap(); 
            now_time = t.get_elapsed();
            sum_time += now_time;
            printf("Total time: \t%5.2f [ms]\n", sum_time);
        }

        if (sum_time < tf)
            std::this_thread::sleep_for(std::chrono::milliseconds( int(tf-sum_time) ));

        // update for next iteration
        t.lap();
        sum_time = 0.0;
    };

    return 0;
};

std::vector<std::string> list_directory(std::string path)
{
    std::list<std::string> list_files;
    for (const auto & entry : fs::directory_iterator(path))
        list_files.push_back(entry.path());
    list_files.sort();
    return std::vector<std::string>(list_files.begin(), list_files.end());
};