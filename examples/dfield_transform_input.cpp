/*
* @Author: Jose Tascon
* @Date:   2019-11-18 17:17:46
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-12-07 18:55:56
*/

// std libs
#include <iostream>
#include <memory>
#include <vector>

// boost libs
#include <boost/program_options.hpp>

// local libs
#include "../src/image.h"
#include "../src/grid.h"
#include "../src/dfield.h"
#include "../src/ilinear.h"

using namespace imart;
namespace po = boost::program_options;

int main(int argc, char *argv[])
{
    // ============================================
    //              Testing resolution
    // ============================================
    // using type = float;
    using type = double;

    // Variables
    bool verbose;
    unsigned int dim;
    std::string file_transform;
    std::string file_fixed;
    std::string file_moving;
    std::string file_output;

    // Program description
    po::options_description desc("Transform an image with a deformation field. Options");
    desc.add_options()
    ("help,h", "Help message")
    ("dimension,d", po::value<unsigned int>(&dim)->default_value(2), "Dimension")
    ("transform,t", po::value<std::string>(&file_transform), "Transform file name")
    ("reference,r", po::value<std::string>(&file_fixed), "Reference image")
    ("input,i", po::value<std::string>(&file_moving), "Input image")
    ("output,o", po::value<std::string>(&file_output), "Output file name")
    ("verbose,v", po::bool_switch(&verbose), "Enable verbose");

    // Parse command line options
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")){ std::cout << desc << std::endl; return 1; };

    if (vm.count("reference")) std::cout << "Reference image: " << file_fixed << std::endl;
    else {std::cout << "Please provide reference image" << std::endl; return 1; };

    if (vm.count("input")) std::cout << "Input image: " << file_moving << std::endl;
    else {std::cout << "Please provide moving image" << std::endl; return 1; };

    if (vm.count("transform")) std::cout << "Transform: " << file_transform << std::endl;
    else {std::cout << "Please provide moving image" << std::endl; return 1; };

    // ============================================
    //                  Code
    // ============================================

    // Read Images
    auto img_fixed = image_cpu<type>::new_pointer(dim);
    img_fixed->read(file_fixed);

    auto img_moving = image_cpu<type>::new_pointer(dim);
    img_moving->read(file_moving);

    // Transform
    auto transformation = dfield<type,vector_cpu<type>>::new_pointer(dim);
    transformation->read(file_transform);
    auto interpolation = ilinear_cpu<type>::new_pointer(img_moving);
    auto x0 = grid_cpu<type>::new_pointer(img_fixed);
    auto moving_warped = interpolation->apply(transformation->apply(x0));

    moving_warped->write(file_output);

    return 0;
};