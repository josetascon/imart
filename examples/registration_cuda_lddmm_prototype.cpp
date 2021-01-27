/*
* @Author: Jose Tascon
* @Date:   2019-11-18 13:30:52
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2021-01-27 06:50:36
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
#include "../src/utils/timer.h"

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
    std::string file_output = "./moving_lddmm_";
    bool verbose;
    bool plot;
    unsigned int dim; int iter, tsteps;
    double tolerance, step;
    double sigma_f, sigma_e, sigma_l, alpha, gamma;

    // Program description
    po::options_description desc("Register two images with lddmm algorithm. Options");
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
    ("sigma_elastic,e", po::value<double>(&sigma_e)->default_value(3.0), "Deformation field sigma elastic")
    ("sigma_lddmm,s", po::value<double>(&sigma_l)->default_value(0.1), "LDDMM sigma")
    ("alpha_lddmm,a", po::value<double>(&alpha)->default_value(10.0), "LDDMM alpha")
    ("gamma_lddmm,g", po::value<double>(&gamma)->default_value(1.0), "LDDMM gamma")
    ("time_steps,n", po::value<int>(&tsteps)->default_value(32), "LDDMM time steps");

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
    img_fixed = normalize(img_fixed);

    // std::string file_moving = argv[2];
    auto img_moving = image_cuda<type>::new_pointer(dim);
    img_moving->read(file_moving);
    img_moving = normalize(img_moving);

    img_fixed->print("Fixed Image");
    std::cout << "Min: " << img_fixed->min() << std::endl;
    std::cout << "Max: " << img_fixed->max() << std::endl;

    img_moving->print("Moving Image");
    std::cout << "Min: " << img_moving->min() << std::endl;
    std::cout << "Max: " << img_moving->max() << std::endl;

    // *img_fixed = normalize<type>(*img_fixed);
    // *img_moving = normalize<type>(*img_moving);

    // img_fixed->print();
    // img_fixed->print_data();

    // img_moving->print();
    // img_moving->print_data();

    // auto mg = gradient(img_moving);

    auto trfm = dfield<type,vector_cuda<type>>::new_pointer(img_fixed);
    auto trfm_out = transform<type,vector_cuda<type>>::new_pointer(dim);
    trfm->set_sigma_fluid(sigma_f);
    trfm->set_sigma_elastic(sigma_e);
    // trfm->print();
    // trfm->print_data();

    auto lddmm1 = lddmm<type,vector_cuda<type>>::new_pointer(img_fixed, img_moving, trfm);
    lddmm1->set_learning_rate(step);
    lddmm1->set_time_steps(tsteps);
    lddmm1->set_sigma(sigma_l);
    lddmm1->set_alpha_gamma(alpha, gamma);

    type diff = 0.0;
    type current_cost = 0.0;
    type previous_cost = 1e100;
    timer t("ms");
    t.start();

    std::string termination = "iterations";
    
    for (int k = 0; k < iter; k++)
    {
        trfm_out = lddmm1->derivative();

        current_cost = lddmm1->cost();

        diff = abs(previous_cost - current_cost);
        t.lap();
        printf( "iteration: %4d  cost: %7.3e  diff: %7.3e  time: %7.3f\n", k, current_cost, diff, t.get_elapsed() );

        if (lddmm1->fault())
        {
            termination = lddmm1->fault_info();
            // break;q
            // std::cout << "determinant error" << std::endl;
        };

        
        // if (plot)
        // {
        //     auto view = viewer<image_cpu<type>>::new_pointer();
        //     view->size(1600,800);
        //     view->subplot(2,5);

        //     view->add_image(normalize(lddmm1->j1->at(0),0.0,255.0));
        //     view->add_image(normalize(lddmm1->j1->at(7),0.0,255.0));
        //     view->add_image(normalize(lddmm1->j1->at(15),0.0,255.0));
        //     view->add_image(normalize(lddmm1->j1->at(23),0.0,255.0));
        //     view->add_image(normalize(lddmm1->j1->at(31),0.0,255.0));

        //     view->add_image(normalize(lddmm1->j0->at(0),0.0,255.0));
        //     view->add_image(normalize(lddmm1->j0->at(7),0.0,255.0));
        //     view->add_image(normalize(lddmm1->j0->at(15),0.0,255.0));
        //     view->add_image(normalize(lddmm1->j0->at(23),0.0,255.0));
        //     view->add_image(normalize(lddmm1->j0->at(31),0.0,255.0));

        //     // view->add_image(lddmm1->j1->at(0));
        //     // view->add_image(lddmm1->j1->at(7));
        //     // view->add_image(lddmm1->j1->at(15));
        //     // view->add_image(lddmm1->j1->at(23));
        //     // view->add_image(lddmm1->j1->at(31));

        //     // view->add_image(lddmm1->j0->at(0));
        //     // view->add_image(lddmm1->j0->at(7));
        //     // view->add_image(lddmm1->j0->at(15));
        //     // view->add_image(lddmm1->j0->at(23));
        //     // view->add_image(lddmm1->j0->at(31));
            
        //     view->setup();
        //     // view->visualize();
        //     view->show();
        // };

        previous_cost = current_cost;
    };

    std::cout << "Termination:\t" << termination << std::endl;

    // Save velocity
    std::string vfile0 = "lddmm_v0_";
    std::string vfile1 = "lddmm_v1_";
    for (int d = 0; d < dim; d++)
    {
        lddmm1->v->at(0)[d]->write(vfile0 + std::to_string(d) + ".nrrd");
        lddmm1->v->at(tsteps-1)[d]->write(vfile1 + std::to_string(d) + ".nrrd");
    }

    // Save displacements
    std::string pfile0 = "lddmm_phi0_";
    std::string pfile1 = "lddmm_phi1_";
    // for (int t = 0; t < tsteps; t++)
    // {
    for (int d = 0; d < dim; d++)
    {
        lddmm1->phi0->at(0)[d]->write(pfile0 + "0_" + std::to_string(d) + ".nrrd");
        lddmm1->phi0->at(tsteps-1)[d]->write(pfile0 + "1_" + std::to_string(d) + ".nrrd");

        lddmm1->phi1->at(1)[d]->write(pfile1 + "0_" + std::to_string(d) + ".nrrd");
        lddmm1->phi1->at(tsteps-1)[d]->write(pfile1 + "1_" + std::to_string(d) + ".nrrd");
    }
    // }

    auto transformation = trfm_out;
    auto interpolation = ilinear_cuda<type>::new_pointer(img_moving);
    auto x0 = grid_cuda<type>::new_pointer(img_fixed);
    auto moving_warped = interpolation->apply(transformation->apply(x0));

    // auto moving_warped = lddmm1->warped_moving();
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
    // output->write("./lddmm_warped.png");

    moving_warped->write(file_output + "warped.nrrd");

    transformation->write(file_output + "trfm_warp.nrrd");

    // if (plot)
    // {
    //     auto view = viewer<image_cpu<type>>::new_pointer();
    //     view->size(1000,400);
    //     view->subplot(1,3);

    //     view->add_image( normalize(img_moving,0.0,255.0)    );
    //     view->add_image( normalize(img_fixed,0.0,255.0)     );
    //     view->add_image( normalize(moving_warped,0.0,255.0) );
    //     // view->add_image(img_moving);
    //     // view->add_image(img_fixed);
    //     // view->add_image(moving_warped);
    //     view->setup();
    //     // view->visualize();
    //     view->show();
    // };

    return 0;
};