/*
* @Author: Jose Tascon
* @Date:   2019-11-18 13:30:52
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2021-01-27 06:44:29
*/

// std libs
#include <iostream>
#include <memory>

// boost libs
#include <boost/program_options.hpp>

// local libs
#include "../src/image.h"
#include "../src/dfield.h"
#include "../src/demons.h"
#include "../src/demons_diffeomorphic.h"
#include "../src/gradient_descent.h"
#include "../src/registration.h"
#include "../src/viewer.h"

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
    std::string file_fixed;
    std::string file_moving;
    std::string file_output = "./moving_demons_";
    bool verbose;
    bool plot;
    unsigned int dim; int iter;
    double tolerance, step;
    double sigma_f, sigma_e;

    // Program description
    po::options_description desc("Register two images with demons algorithm. Options");
    desc.add_options()
    ("help,h", "Help message")
    ("dimension,d", po::value<unsigned int>(&dim)->default_value(2), "Dimension")
    ("fixed,f", po::value<std::string>(&file_fixed), "Fixed image")
    ("moving,m", po::value<std::string>(&file_moving), "Moving image")
    ("output,o", po::value<std::string>(&file_output), "Output prefix name")
    ("verbose,v", po::bool_switch(&verbose), "Enable verbose")
    ("plot,p", po::bool_switch(&plot), "Enable plot")
    ("opt-step,l", po::value<double>(&step)->default_value(1.0), "Optimizer step")
    ("opt-iterations,i", po::value<int>(&iter)->default_value(150), "Optimizer iterations")
    ("opt-tolerance,t", po::value<double>(&tolerance)->default_value(1e-7), "Optimizer tolerance")
    ("sigma_fluid,u", po::value<double>(&sigma_f)->default_value(0.0), "Deformation field sigma fluid")
    ("sigma_elastic,e", po::value<double>(&sigma_e)->default_value(3.0), "Deformation field sigma elastic");

    // Parse command line options
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help")){ std::cout << desc << std::endl; return 1; };

    if (vm.count("fixed")) std::cout << "Fixed image: " << file_fixed << std::endl;
    else {std::cout << "Please provide fixed image" << std::endl; return 1; };

    if (vm.count("moving")) std::cout << "Moving image: " << file_moving << std::endl;
    else {std::cout << "Please provide moving image" << std::endl; return 1; };

    // ============================================
    //          Testing metric
    // ============================================
    // Create small imame
    // std::string file_fixed = argv[1];
    auto img_fixed = image_cuda<type>::new_pointer(dim);
    img_fixed->read(file_fixed);

    // std::string file_moving = argv[2];
    auto img_moving = image_cuda<type>::new_pointer(dim);
    img_moving->read(file_moving);

    img_fixed->print("Fixed Image");
    // std::cout << "min: " << img_fixed->min() << std::endl;
    // std::cout << "max: " << img_fixed->max() << std::endl;

    img_moving->print("Moving Image");
    // std::cout << "min: " << img_moving->min() << std::endl;
    // std::cout << "max: " << img_moving->max() << std::endl;

    // *img_fixed = normalize<type>(*img_fixed);
    // *img_moving = normalize<type>(*img_moving);

    // img_fixed->print();
    // img_fixed->print_data();

    // img_moving->print();
    // img_moving->print_data();

    // auto mg = gradient(img_moving);

    auto trfm = dfield<type,vector_cuda<type>>::new_pointer(img_fixed);
    trfm->set_sigma_fluid(sigma_f);
    trfm->set_sigma_elastic(sigma_e);
    // trfm->print();
    // trfm->print_data();

    auto demons1 = demons<type,vector_cuda<type>>::new_pointer(img_fixed, img_moving, trfm);
    auto opt = gradient_descent<type,vector_cuda<type>>::new_pointer();
    opt->set_step(step);
    opt->set_tolerance(tolerance);
    opt->set_unchanged_times(15);

    auto registro = registration<type,vector_cuda<type>>::new_pointer(img_fixed, img_moving, trfm);
    registro->set_levels(5);
    registro->set_levels_scales(std::vector<int>{10,6,4,2,1});
    registro->set_levels_iterations(std::vector<int>{iter,iter,iter,iter,iter});
    registro->set_metric(demons1);
    registro->set_optimizer(opt);
    registro->set_normalize(true);
    // registro->set_padding(true);
    registro->apply();


    auto transformation = registro->get_transform();
    auto interpolation = ilinear_cuda<type>::new_pointer(img_moving);
    auto x0 = grid_cuda<type>::new_pointer(img_fixed);
    auto moving_warped = interpolation->apply(transformation->apply(x0));

    // auto moving_warped = demons1->warped_moving();
    // moving_warped->print();
    // moving_warped->print_data();

    // auto moving_cast = image_cpu<intype>::new_pointer();
    // cast((*moving_warped)*(type(255)), *moving_cast);

    // auto output = image_cpu<intype>::new_pointer();
    // *output = normalize<intype>(*moving_cast, min, max);

    // auto moving = (*moving_warped)*(type(255));
    // moving.print_data();
    
    // cast(*moving, *output);
    // cast(*moving_warped, *output);
    // output->print_data();
    // moving_warped->write("./moving_warped.nii");
    // output->write("./demons_warped.png");

    moving_warped->write(file_output + "warped.nrrd");

    transformation->write(file_output + "trfm_warp.nrrd");

    // if (plot)
    // {
    //     auto view = viewer<image_cpu<type>>::new_pointer();
    //     view->size(1000,400);
    //     view->subplot(1,3);
    //     view->add_image(img_moving);
    //     view->add_image(img_fixed);
    //     view->add_image(moving_warped);
    //     view->setup();
    //     // view->visualize();
    //     view->show();
    // };

    return 0;
};