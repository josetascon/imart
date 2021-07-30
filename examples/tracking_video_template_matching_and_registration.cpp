/*
* @Author: Jose Tascon
* @Date:   2019-11-18 13:30:52
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2021-07-29 17:16:40
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
#include "../src/template_matching.h"
#include "../src/viewer_track.h"

#include "../src/grid.h"
#include "../src/ilinear.h"
#include "../src/dfield.h"
#include "../src/demons.h"
#include "../src/gradient_descent.h"
#include "../src/registration.h"

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

    std::string path_input;
    std::string path_output = "";
    std::string file_mask;
    bool verbose, plot, write;
    double fps, tf;
    int slide, extra_pixels;
    double tolerance, step;
    double sigmaf, sigmae;
    std::vector<int> level_scales({4,2,1});
    std::vector<int> level_iterations({120,100,80});

    // Program options
    po::options_description desc("Tracking 2D CineMR images with template matching. Options");
    desc.add_options()
    ("help,h", "Help message")
    ("input,i", po::value<std::string>(&path_input), "Input folder")
    ("mask,m", po::value<std::string>(&file_mask), "Mask image used for tracking")
    ("output,o", po::value<std::string>(&path_output), "Output folder")
    ("verbose,v", po::bool_switch(&verbose), "Enable verbose")
    ("plot,p", po::bool_switch(&plot), "Enable plot")
    ("frame_per_second,f", po::value<double>(&fps)->default_value(4.0), "Images frames per second")
    ("slide,w", po::value<int>(&slide)->default_value(10), "Sliding window search in pixels")
    ("extra-pixels,x", po::value<int>(&extra_pixels)->default_value(4), "Extra pixels of bounding box")
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

    if (vm.count("input")) std::cout << "Folder input: " << path_input << std::endl << std::endl;
    else {std::cerr << "Please provide input folder with images" << std::endl; return 1; };

    if (vm.count("mask")) std::cout << "File mask: " << path_input << std::endl << std::endl;
    else {std::cerr << "Please provide image mask or template" << std::endl; return 1; };

    if (path_output == "") write = false;
    else write = true;

    tf = 1000/fps; // sampling time in milliseconds

    // ============================================
    //       Testing video template matching
    // ============================================
    
    // Images
    auto img_fixed = image_cpu<type>::new_pointer();
    auto img_input = image_cpu<type>::new_pointer();
    auto img_mask = image_cpu<type>::new_pointer();
    auto img_mask_warped = image_cpu<type>::new_pointer();
    
    // Read fixed image
    auto list_files_images = list_directory(path_input);
    std::string file_fixed = list_files_images.front();
    img_fixed->read(file_fixed);
    std::cout << "Read fixed: " << file_fixed << std::endl;

    img_mask->read(file_mask);
    std::cout << "Read tracking mask: " << file_mask << std::endl << std::endl;
    
    // Template matching
    auto bbox_fixed = bounding_box(img_mask, extra_pixels);
    auto img_mask_region = img_mask->region(bbox_fixed[0], bbox_fixed[1]);
    auto img_fixed_region = img_fixed->region(bbox_fixed[0],bbox_fixed[1]);
    auto tm = template_matching<type, vector_cpu<type>>::new_pointer(img_fixed, bbox_fixed);
    tm->set_slide(std::vector<int>{slide,slide});

    // Registration
    auto trfm = dfield<type,vector_cpu<type>>::new_pointer(img_fixed);
    trfm->set_sigma_fluid(sigmaf);
    trfm->set_sigma_elastic(sigmae);

    auto metric_demons = demons<type,vector_cpu<type>>::new_pointer(img_fixed_region, img_fixed_region, trfm);
    auto opt = gradient_descent<type,vector_cpu<type>>::new_pointer();
    opt->set_step(step);
    opt->set_tolerance(tolerance);
    opt->set_unchanged_times(10);
    opt->set_verbose(false);

    auto registro = registration<type,vector_cpu<type>>::new_pointer(img_fixed_region, img_fixed_region, trfm);
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
    auto img_view_mask = img_mask->copy();
    
    // Setup viewer
    if (plot)
    {
        view->add_image(img_view_input);
        view->add_image(img_view_mask);
        view->setup();
    }

    timer t("ms");
    t.start();
    double now_time = 0.0;
    double sum_time = 0.0;

    // Output folder
    std::string folder_output = path_output + "/";
    if (write)
    {
        std::cout << "Create output folder: " << folder_output << std::endl;
        std::cout << std::endl;
        fs::create_directories(folder_output);
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

        // Template matching
        auto bbox_moving = tm->apply(img_input);
        tm->set_fixed(img_input->clone());
        tm->set_box_fixed(bbox_moving);

        if (verbose)
        {
            t.lap(); 
            now_time = t.get_elapsed();
            sum_time += now_time;
            printf("Template time: \t%5.2f [ms]\n", now_time);
        }

        // Registration
        auto img_moving_region = img_input->region(bbox_moving[0],bbox_moving[1]);
        registro->set_moving(img_moving_region);
        registro->apply();

        if (verbose)
        {
            t.lap(); 
            now_time = t.get_elapsed();
            sum_time += now_time;
            printf("Register time: \t%5.2f [ms]\n", now_time);
        }

        // Warp with registration
        std::vector<int> pre = bbox_moving[0];
        std::vector<int> post(bbox_moving[1].size(),0);
        for (int k = 0; k < img_fixed->get_dimension(); k++)
        {
            post[k] = img_fixed->get_size()[k] - bbox_moving[0][k] - bbox_moving[1][k];
        }

        auto transformation = registro->get_transform()->inverse();
        auto interpolation = ilinear_cpu<type>::new_pointer(img_mask_region);
        auto x1 = grid_cpu<type>::new_pointer(img_moving_region);
        auto img_mask_region_warped = interpolation->apply(transformation->apply(x1));

        auto img_mask_warped = pad(img_mask_region_warped, pre, post);
        img_mask_warped->set_spacing(img_mask->get_spacing());
        img_mask_warped->set_origin(img_mask->get_origin());
        img_mask_warped->set_direction(img_mask->get_direction());

        // std::string base = fs::path(path_organs[k]).filename();
        std::string base = "test";
        std::string file_output = folder_output + "/" + base + "_" + num + ext;
        if (write)
        {
            std::cout << "Write output: " << file_output << std::endl;
            img_mask_warped->write(file_output);
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
            img_view_mask->equal(*img_mask_warped);
            view->update(0);
            view->update(1);
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