/*
* @Author: Jose Tascon
* @Date:   2019-11-18 13:30:52
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2021-03-11 00:19:23
*/

// std libs
#include <iostream>
#include <memory>
#include <sstream>
#include <iomanip>
#include <filesystem>
#include <thread>

// local libs
#include "../src/image.h"
#include "../src/dfield.h"
#include "../src/demons.h"
#include "../src/gradient_descent.h"
#include "../src/registration.h"

#include "../src/utils/timer.h"
#include "../src/viewer_track.h"

using namespace imart;

std::string to_zero_lead(const int value, const unsigned precision)
{
     std::ostringstream oss;
     oss << std::setw(precision) << std::setfill('0') << value;
     return oss.str();
}

int main(int argc, char *argv[])
{
    using type = float;
    // using type = double;
    using image_type = image_cpu<type>;
    
    if( argc < 2 )
    {
        std::cerr << "Usage: " << std::endl;
        std::cerr << argv[0] << " input_folder" << std::endl;
        return EXIT_FAILURE;
    }
    std::string input_path = argv[1];
    // std::string output_path = argv[2];

    // ============================================
    //          Testing video registration
    // ============================================

    // string variables
    size_t init_image = 1;
    size_t num_images = 400;
    
    std::string ext = ".nii";
    std::string fix = "0000";
    std::string num = "0001";
    // std::string file_input = "";

    // Images
    auto img_fixed = image_cpu<type>::new_pointer();
    auto img_input = image_cpu<type>::new_pointer();
    auto img_tumor = image_cpu<type>::new_pointer();
    auto img_liver = image_cpu<type>::new_pointer();
    auto img_lung = image_cpu<type>::new_pointer();
    auto img_tumor_warped = image_cpu<type>::new_pointer();
    auto img_lung_warped = image_cpu<type>::new_pointer();
    auto img_liver_warped = image_cpu<type>::new_pointer();
    auto img_view_input  = image_cpu<type>::new_pointer();
    auto img_view_tumor  = image_cpu<type>::new_pointer();
    auto img_view_liver  = image_cpu<type>::new_pointer();
    auto img_view_lung  = image_cpu<type>::new_pointer();

    // Read fixed image
    std::string file_fixed = input_path + "/images/patient01_" + fix + ext;
    img_fixed->read(file_fixed);
    std::cout << "Read fixed: " << file_fixed << std::endl;
    img_view_input = img_fixed->copy();

    std::string file_input = input_path + "/images/patient01_" + num + ext;
    img_input->read(file_input);
    std::cout << "Read input: " << file_input << std::endl;
    
    std::string file_tumor = input_path + "/tumor/" + num + ext;
    img_tumor->read(file_tumor);
    std::cout << "Read mask: " << file_tumor << std::endl;
    img_view_tumor = img_tumor->copy();

    std::string file_liver = input_path + "/liver/" + num + ext;
    img_liver->read(file_liver);
    std::cout << "Read mask: " << file_liver << std::endl;
    img_view_liver = img_liver->copy();

    std::string file_lung = input_path + "/lung/" + num + ext;
    img_lung->read(file_lung);
    std::cout << "Read mask: " << file_lung << std::endl;
    img_view_lung = img_lung->copy();

    // Transform
    double sigmaf = 0.0;
    double sigmae = 4.5;
    auto trfm = dfield<type,vector_cpu<type>>::new_pointer(img_fixed);
    trfm->set_sigma_fluid(sigmaf);
    trfm->set_sigma_elastic(sigmae);

    // Metric
    auto demons1 = demons<type,vector_cpu<type>>::new_pointer(img_fixed, img_input, trfm);
    auto opt = gradient_descent<type,vector_cpu<type>>::new_pointer();
    opt->set_step(0.8);
    opt->set_tolerance(1e-6);
    opt->set_unchanged_times(10);
    opt->set_verbose(false);

    // Registration
    auto registro = registration<type,vector_cpu<type>>::new_pointer(img_fixed, img_input, trfm);
    registro->set_levels(3);
    registro->set_levels_scales(std::vector<int>{8,4,2});
    registro->set_levels_iterations(std::vector<int>{120,100,80});
    registro->set_metric(demons1);
    registro->set_optimizer(opt);
    registro->set_normalize(true);
    registro->set_verbose(false);
    // registro->apply();

    // Post registration
    auto interpol_tumor = ilinear_cpu<type>::new_pointer(img_tumor);
    auto interpol_liver = ilinear_cpu<type>::new_pointer(img_liver);
    auto interpol_lung = ilinear_cpu<type>::new_pointer(img_lung);
    auto xi = grid_cpu<type>::new_pointer(img_fixed);

    // Setup viewer
    viewer_track<image_type> view;
    view.add_image(img_view_input);
    view.add_image(img_view_tumor);
    view.add_image(img_view_liver);
    view.add_image(img_view_lung);
    view.setup();

    timer t("ms");
    t.start();
    double now_time = 0.0;
    double sum_time = 0.0;

    for(size_t i = init_image; i < num_images; i++ )
    {
        file_input = input_path + "/images/patient01_" + num + ext;
        img_input->read(file_input);
        std::cout << "Read input: " << file_input << std::endl;
        // img_view_input->equal(*img_input);
        // view.update(0);

        t.lap(); 
        now_time = t.get_elapsed();
        sum_time += now_time;
        printf("Read time: \t%5.2f [ms]\n", now_time);

        trfm = dfield<type,vector_cpu<type>>::new_pointer(img_fixed);
        trfm->set_sigma_fluid(sigmaf);
        trfm->set_sigma_elastic(sigmae);

        registro->set_moving(img_input);
        registro->set_transform(trfm);
        registro->apply();

        t.lap(); 
        now_time = t.get_elapsed();
        sum_time += now_time;
        printf("Register time: \t%5.2f [ms]\n", now_time);

        auto transformation = registro->get_transform()->inverse();
        // auto transformation = registro->get_transform();
        auto img_tumor_warped = interpol_tumor->apply(transformation->apply(xi));
        auto img_liver_warped = interpol_liver->apply(transformation->apply(xi));
        auto img_lung_warped = interpol_lung->apply(transformation->apply(xi));

        t.lap(); 
        now_time = t.get_elapsed();
        sum_time += now_time;
        printf("Mask time: \t%5.2f [ms]\n", now_time);

        img_view_input->equal(*img_input);
        img_view_tumor->equal(*img_tumor_warped);
        img_view_liver->equal(*img_liver_warped);
        img_view_lung->equal(*img_lung_warped);

        t.lap(); 
        now_time = t.get_elapsed();
        sum_time += now_time;
        printf("Equal time: \t%5.2f [ms]\n", now_time);

        view.update(0);
        view.update(1);
        view.update(2);
        view.update(3);

        t.lap(); 
        now_time = t.get_elapsed();
        sum_time += now_time;
        printf("Modify time: \t%5.2f [ms]\n", now_time);
        
        view.render();

        t.lap();
        now_time = t.get_elapsed();
        sum_time += now_time;
        printf("Render time: \t%5.2f [ms]\n", now_time);
        printf("Total time: \t%5.2f [ms]\n", sum_time);

        if (sum_time < 250)
            std::this_thread::sleep_for(std::chrono::milliseconds(250-int(sum_time)));

        // update for next iteration
        t.lap();
        sum_time = 0.0;
        num = to_zero_lead(i+1,4);
        // std::cout << "num " << num << std::endl;
    };
        


    return 0;
};