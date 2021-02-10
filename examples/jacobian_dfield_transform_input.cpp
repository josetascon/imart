/*
* @Author: Jose Tascon
* @Date:   2019-11-18 17:17:46
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2021-01-28 10:32:11
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
    std::string file_transform;
    std::string file_output;

    // Program description
    po::options_description desc("Jacobian stats of deformation field. Options");
    desc.add_options()
    ("help,h", "Help message")
    ("dimension,d", po::value<unsigned int>(&dim)->default_value(2), "Dimension")
    ("transform,t", po::value<std::string>(&file_transform), "Transform file name")
    ("output,o", po::value<std::string>(&file_output)->default_value("jac.nrrd"), "Output file name")
    ("verbose,v", po::bool_switch(&verbose), "Enable verbose");

    // Parse command line options
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")){ std::cout << desc << std::endl; return 1; };

    if (vm.count("transform")) std::cout << "Transform: " << file_transform << std::endl;
    else {std::cout << "Please provide dfield transform" << std::endl << desc << std::endl; return 1; };

    // ============================================
    //                  Code
    // ============================================

    // Transform
    auto transformation = dfield<type,vector_cpu<type>>::new_pointer(dim);
    transformation->read(file_transform);

    typename image<type,vector_cpu<type>>::pointer img;
    img = jacobian(transformation->get_parameters_vector(), true);

    transformation->print();
    printf("\nJacobian Determinant Stats\n");
    printf("Min:     \t%f\n", img->min());
    printf("Max:     \t%f\n", img->max());
    printf("Mean:    \t%f\n", img->mean());
    printf("Std dev: \t%f\n", img->std_dev());

    img->write(file_output);

    return 0;
};