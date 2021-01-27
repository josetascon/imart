/*
* @Author: jose
* @Date:   2019-11-07 10:12:34
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2021-01-27 07:09:19
*/

// std libs
#include <iostream>
#include <memory>

// boost libs
#include <boost/program_options.hpp>

// local libs
#include "../src/image.h"
#include "../src/dfield.h"
#include "../src/lddmm.h"
#include "../src/gradient_descent.h"
#include "../src/registration.h"
#include "../src/viewer.h"

using namespace imart;
namespace po = boost::program_options;

int main(int argc, char *argv[])
{
    // ============================================
    //              Parse Options
    // ============================================
    // using type = float;
    using type = double;

    // Variables
    std::string file_input;
    std::string file_output;
    unsigned int dim;

    // Program description
    po::options_description desc("Read an image a write output. Options");
    desc.add_options()
    ("help,h", "Help message")
    ("dimension,d", po::value<unsigned int>(&dim)->default_value(2), "Dimension")
    ("input,i", po::value<std::string>(&file_input), "File name of input image")
    ("output,o", po::value<std::string>(&file_output)->default_value("./image_output.nrrd"), "File name of output image");

    // Parse command line options
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")){ std::cout << desc << std::endl; return 1; };

    if (vm.count("input")) std::cout << "Input image: " << file_input << std::endl;
    else {std::cout << "Please provide input image" << std::endl << desc << std::endl; return 1; };

    if (vm.count("output")) std::cout << "Output image: " << file_output << std::endl;
    else {std::cout << "Please provide output image" << std::endl << desc << std::endl; return 1; };

    // ============================================
    //                  Testing
    // ============================================
    auto img = image_cpu<type>::new_pointer(dim);
    img->read(file_input);

    img->write(file_output);
    
    return 0;
};