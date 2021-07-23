/*
* @Author: Jose Tascon
* @Date:   2019-11-18 13:30:52
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2021-03-08 17:33:39
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
// #include "../src/dfield.h"
// #include "../src/demons.h"
// #include "../src/gradient_descent.h"
// #include "../src/registration.h"

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
    // using type = float;
    using type = double;
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
    size_t num_images = 100;
    // std::string input_path = "/home/jose/Public/workspace/medical_imaging/liver/scripts/video_out/patient06/video1/";
    std::string ext = ".nii";
    std::string fix = "0000";
    std::string num = "0000";

    // Images
    auto img_fixed = image_cpu<type>::new_pointer();
    auto img_input = image_cpu<type>::new_pointer();
    auto img_view  = image_cpu<type>::new_pointer();
    auto img_tumor = image_cpu<type>::new_pointer();
    auto img_liver = image_cpu<type>::new_pointer();
    auto img_lung = image_cpu<type>::new_pointer();

    // Read Reference
    std::string file_input = input_path + "/images/patient01_" + num + ext;
    img_input->read(file_input);
    std::cout << "Read input: " << file_input << std::endl;
    img_view = img_input->copy();

    std::string file_tumor = input_path + "/tumor/" + num + ext;
    img_tumor->read(file_tumor);
    std::cout << "Read input: " << file_tumor << std::endl;

    std::string file_liver = input_path + "/liver/" + num + ext;
    img_liver->read(file_liver);
    std::cout << "Read input: " << file_liver << std::endl;

    std::string file_lung = input_path + "/lung/" + num + ext;
    img_lung->read(file_lung);
    std::cout << "Read input: " << file_lung << std::endl;

    // Setup viewer
    viewer_track<image_type> view;
    view.add_image(img_view);
    view.add_image(img_tumor);
    view.add_image(img_liver);
    view.add_image(img_lung);
    view.setup();

    for(size_t i = 1; i < num_images; i++ )
    {

        file_input = input_path + "/images/patient01_" + num + ext;
        img_input->read(file_input);
        std::cout << "Read input: " << file_input << std::endl;
        img_view->equal(*img_input);
        view.update(0);

        file_tumor = input_path + "/tumor/" + num + ext;
        img_input->read(file_tumor);
        std::cout << "Read input: " << file_tumor << std::endl;
        img_tumor->equal(*img_input);
        view.update(1);

        file_liver = input_path + "/liver/" + num + ext;
        img_input->read(file_liver);
        std::cout << "Read input: " << file_liver << std::endl;
        img_liver->equal(*img_input);
        view.update(2);

        file_lung = input_path + "/lung/" + num + ext;
        img_input->read(file_lung);
        std::cout << "Read input: " << file_lung << std::endl;
        img_lung->equal(*img_input);
        view.update(3);
        
        view.render();

        // update for next iteration
        num = to_zero_lead(i,4);

        std::this_thread::sleep_for(std::chrono::milliseconds(300));
    };
        


    return 0;
};