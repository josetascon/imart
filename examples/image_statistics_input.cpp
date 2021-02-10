/*
* @Author: Jose Tascon
* @Date:   2019-11-18 17:17:46
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2021-01-28 08:40:39
*/

// std libs
#include <iostream>
#include <memory>
#include <vector>

// boost libs
#include <boost/program_options.hpp>

// local libs
#include "../src/image.h"
#include "../src/dfield.h"

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
    std::string file_image;

    // Program description
    po::options_description desc("Jacobian stats of deformation field. Options");
    desc.add_options()
    ("help,h", "Help message")
    ("dimension,d", po::value<unsigned int>(&dim)->default_value(2), "Dimension")
    ("image,i", po::value<std::string>(&file_image), "Transform file name")
    ("verbose,v", po::bool_switch(&verbose), "Enable verbose");

    // Parse command line options
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")){ std::cout << desc << std::endl; return 1; };

    if (vm.count("image")) std::cout << "Image: " << file_image << std::endl;
    else {std::cout << "Please provide image file" << std::endl << desc << std::endl; return 1; };

    // ============================================
    //                  Code
    // ============================================

    // Transform
    auto img = image<type,vector_cpu<type>>::new_pointer(dim);
    img->read(file_image);

    img->print();
    printf("\nImage Stats\n");
    printf("Min:     \t%f\n", img->min());
    printf("Max:     \t%f\n", img->max());
    printf("Mean:    \t%f\n", img->mean());
    printf("Std dev: \t%f\n", img->std_dev());

    return 0;
};