/*
* @Author: Jose Tascon
* @Date:   2019-11-18 13:30:52
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2021-02-10 13:09:08
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
    std::string file_output = "./moving_lddmm_proto_";
    bool verbose;
    bool plot, plot_per_iter;
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
    ("plot_per_iter,q", po::bool_switch(&plot_per_iter), "Enable plot per iteration")
    ("opt-step,l", po::value<double>(&step)->default_value(1.0), "Optimizer step")
    ("opt-iterations,i", po::value<int>(&iter)->default_value(150), "Optimizer iterations")
    ("opt-tolerance,t", po::value<double>(&tolerance)->default_value(1e-7), "Optimizer tolerance")
    ("sigma_fluid,u", po::value<double>(&sigma_f)->default_value(0.0), "Deformation field sigma fluid")
    ("sigma_elastic,e", po::value<double>(&sigma_e)->default_value(1.0), "Deformation field sigma elastic")
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
    auto img_fixed = image_cpu<type>::new_pointer(dim);
    img_fixed->read(file_fixed);
    img_fixed = normalize(img_fixed);

    // std::string file_moving = argv[2];
    auto img_moving = image_cpu<type>::new_pointer(dim);
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

    auto trfm = dfield<type,vector_cpu<type>>::new_pointer(img_fixed);
    trfm->set_sigma_fluid(sigma_f);
    trfm->set_sigma_elastic(sigma_e);
    // trfm->print();
    // trfm->print_data();

    auto lddmm1 = lddmm<type,vector_cpu<type>>::new_pointer(img_fixed, img_moving, trfm);
    lddmm1->set_learning_rate(step);
    lddmm1->set_time_steps(tsteps);
    lddmm1->set_sigma(sigma_l);
    lddmm1->set_alpha_gamma(alpha, gamma);

    type diff = 0.0;
    type current_cost = 0.0;
    type previous_cost = 1e100;
    timer t("ms");
    t.start();

    typename transform<type,vector_cpu<type>>::pointer tr_base = lddmm1->get_transform();
    typename transform<type,vector_cpu<type>>::pointer tr_derive = tr_base->mimic();

    std::string termination = "iterations";

    auto view = viewer<image_cpu<type>>::new_pointer();
    if (plot_per_iter)
    {
        view->size(1600,800);
        view->subplot(2,5);
        for (int i=0; i<5; i++)
        {
            int p = int(tsteps/4)*i - 1;
            if (i == 0) p = 0;
            view->add_image(normalize(lddmm1->j1->at(p),0.0,255.0));
        }
        for (int i=0; i<5; i++)
        {
            int p = int(tsteps/4)*i - 1;
            if (i == 0) p = 0;
            view->add_image(normalize(lddmm1->j0->at(p),0.0,255.0));
        }
        view->setup();
        view->visualize();
    };
    
    for (int k = 0; k < iter; k++)
    {
        tr_derive = lddmm1->derivative();

        current_cost = lddmm1->cost();

        *tr_base = *tr_base - (*tr_derive)*step;

        diff = abs(previous_cost - current_cost);
        t.lap();
        printf( "iteration: %4d  cost: %7.3e  diff: %7.3e  time: %7.3f\n", k, current_cost, diff, t.get_elapsed() );

        if (lddmm1->fault())
        {
            termination = lddmm1->fault_info();
            // break;q
            // std::cout << "determinant error" << std::endl;
        };

        // TODO: VIEWER ONLY WORK WITH DOUBLE. FIX
        int c = 0;
        if (plot_per_iter)
        {
            // auto view = viewer<image_cpu<type>>::new_pointer();
            // view->size(1600,800);
            // view->subplot(2,5);
            for (int i=0; i<5; i++)
            {
                int p = int(tsteps/4)*i - 1;
                if (i == 0) p = 0;
                // view->add_image(normalize(lddmm1->j1->at(p),0.0,255.0));
                view->update_image(normalize(lddmm1->j1->at(p),0.0,255.0), c);
                c++;
            }
            for (int i=0; i<5; i++)
            {
                int p = int(tsteps/4)*i - 1;
                if (i == 0) p = 0;
                view->update_image(normalize(lddmm1->j0->at(p),0.0,255.0), c);
                c++;
            }
            // view->setup();
            view->visualize();
        };

        
        // if (plot_per_iter)
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
    std::vector<std::string> xx({"x","y","z"});

    std::string vfile0 = "lddmm_v_";
    std::string dvfile0 = "lddmm_dv_";
    for (int d = 0; d < dim; d++)
    {
        lddmm1->v->at(0)[d]->write(vfile0 + "tinit_" + xx[d] + ".nrrd");
        lddmm1->v->at(tsteps-1)[d]->write(vfile0 + "tend_" + xx[d] + ".nrrd");

        lddmm1->dv->at(0)[d]->write(dvfile0 + "tinit_" + xx[d] + ".nrrd");
        lddmm1->dv->at(tsteps-1)[d]->write(dvfile0 + "tend_" + xx[d] + ".nrrd");
    }

    // Save displacements
    std::string pfile0 = "lddmm_phi_zero_";
    std::string pfile1 = "lddmm_phi_one_";
    // for (int t = 0; t < tsteps; t++)
    // {
    for (int d = 0; d < dim; d++)
    {
        lddmm1->phi0->at(0)[d]->write(pfile0 + "tinit_" + xx[d] + ".nrrd");
        lddmm1->phi0->at(tsteps-1)[d]->write(pfile0 + "tend_" + xx[d] + ".nrrd");

        lddmm1->phi1->at(0)[d]->write(pfile1 + "tinit_" + xx[d] + ".nrrd");
        lddmm1->phi1->at(tsteps-1)[d]->write(pfile1 + "tend_" + xx[d] + ".nrrd");
    }
    // }

     // Save jacobian
    std::string jfile0 = "lddmm_jac_";

    lddmm1->det_phi1->at(0)->write(jfile0 + "tinit.nrrd");
    lddmm1->det_phi1->at(tsteps-1)->write(jfile0 + "tend.nrrd");

    auto transformation = tr_base->clone();
    auto interpolation = ilinear_cpu<type>::new_pointer(img_moving);
    auto x0 = grid_cpu<type>::new_pointer(img_fixed);
    auto moving_warped = interpolation->apply(transformation->apply(x0));
    
    transformation->print();
    x0->print();
    moving_warped->print();


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

    if (plot)
    {
        auto view = viewer<image_cpu<type>>::new_pointer();
        view->size(1000,400);
        view->subplot(1,3);
        type max_gray = 255.0;

        view->add_image( normalize(img_moving,    0.0, max_gray) );
        view->add_image( normalize(img_fixed,     0.0, max_gray) );
        view->add_image( normalize(moving_warped, 0.0, max_gray) );
        // view->add_image(img_moving);
        // view->add_image(img_fixed);
        // view->add_image(moving_warped);
        view->setup();
        // view->visualize();
        view->show();
    };

    return 0;
};