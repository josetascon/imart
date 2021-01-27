/*
* @Author: Jose Tascon
* @Date:   2019-11-18 13:30:52
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2021-01-27 06:47:08
*/

// std libs
#include <iostream>
#include <memory>
#include <sstream>
#include <iomanip>
#include <filesystem>

// local libs
#include "../src/image.h"
#include "../src/dfield.h"
#include "../src/demons.h"
#include "../src/gradient_descent.h"
#include "../src/registration.h"

using namespace imart;

std::string to_zero_lead(const int value, const unsigned precision)
{
     std::ostringstream oss;
     oss << std::setw(precision) << std::setfill('0') << value;
     return oss.str();
}

int main(int argc, char *argv[])
{
    // using type = float;
    using type = double;
    
    if( argc < 3 )
    {
        std::cerr << "Usage: " << std::endl;
        std::cerr << argv[0] << " input_folder output_folder" << std::endl;
        return EXIT_FAILURE;
    }
    std::string input_path = argv[1];
    std::string output_path = argv[2];

    // ============================================
    //          Testing video registration
    // ============================================

    // string variables
    size_t num_images = 80;
    // std::string input_path = "/home/jose/Public/workspace/medical_imaging/liver/scripts/video_out/patient06/video1/";
    std::string ext = ".nii";
    std::string fix = "0000";
    std::string mov = "0001";

    // Images
    auto img_fixed = image_cpu<type>::new_pointer();
    auto img_moving = image_cpu<type>::new_pointer();
    auto img_tumor = image_cpu<type>::new_pointer();
    auto img_liver = image_cpu<type>::new_pointer();

    // Read Reference
    std::string file_fixed = input_path + "/image/image_" + "0000" + ext;
    img_fixed->read(file_fixed);
    std::cout << "Read reference: " << file_fixed << std::endl;

    std::string file_tumor = input_path + "/struct00/structure_" + "0000" + ext;
    img_tumor->read(file_tumor);
    std::cout << "Read tumor: " << file_tumor << std::endl;

    std::string file_liver = input_path + "/struct01/structure_" + "0000" + ext;
    img_liver->read(file_liver);
    std::cout << "Read liver: " << file_liver << std::endl;    

    // Transform
    auto trfm = dfield<type,vector_cpu<type>>::new_pointer();
    double sigmaf = 1.0;
    double sigmae = 3.0;
    trfm->set_sigma_fluid(sigmaf);
    trfm->set_sigma_elastic(sigmae);

    // Registration
    auto demonsreg = demons<type,vector_cpu<type>>::new_pointer(img_fixed, img_moving, trfm);
    auto opt = gradient_descent<type,vector_cpu<type>>::new_pointer();
    opt->set_tolerance(1e-7);

    auto registro = registration<type,vector_cpu<type>>::new_pointer(img_fixed, img_moving, trfm);
    registro->set_metric(demonsreg);
    registro->set_optimizer(opt);
    registro->set_levels(4);
    registro->set_levels_scales(std::vector<int>{8,4,2,1});
    registro->set_levels_iterations(std::vector<int>{600,600,600,600});

    // Output folders
    std::filesystem::create_directories(output_path + "/image/");
    std::filesystem::create_directories(output_path + "/tumor/");
    std::filesystem::create_directories(output_path + "/liver/");

    for(size_t i = 1; i < num_images; i++ )
    {

        std::string file_moving = input_path + "/image/image_" + mov + ext;
        img_moving->read(file_moving);
        std::cout << "Read moving: " << file_moving << std::endl;

        demonsreg->get_fixed()->print("Fixed");
        demonsreg->get_moving()->print("Moving");

        trfm = dfield<type,vector_cpu<type>>::new_pointer(img_fixed);
        trfm->set_sigma_fluid(sigmaf);
        trfm->set_sigma_elastic(sigmae);

        registro->set_moving(img_moving);
        registro->set_transform(trfm);

        registro->apply();
        auto moving_warped = demonsreg->warped_moving();
        moving_warped->write(output_path + "/image/" + fix + "to" + mov + ext);

        auto transformation = registro->get_transform()->inverse();
        // transformation->print("transform inverse");
        auto interpolation = ilinear_cpu<type>::new_pointer(img_tumor);
        auto x1 = grid_cpu<type>::new_pointer(img_moving);
        auto tumor_warped = interpolation->apply(transformation->apply(x1));

        auto interpolation2 = ilinear_cpu<type>::new_pointer(img_liver);
        auto liver_warped = interpolation2->apply(transformation->apply(x1));

        // write the images
        tumor_warped->write(output_path + "/tumor/" + fix + "to" + mov + ext);
        liver_warped->write(output_path + "/liver/" + fix + "to" + mov + ext);

        // update for next iteration
        // trfm = std::static_pointer_cast<dfield<type,vector_cpu<type>>>(transformation);
        mov = to_zero_lead(i,4);
    };
        


    return 0;
};