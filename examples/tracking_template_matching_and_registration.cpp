/*
* @Author: Jose Tascon
* @Date:   2019-11-18 13:30:52
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2021-07-22 13:09:15
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
#include "../src/image_utils.h"
#include "../src/template_matching.h"

#include "../src/grid.h"
#include "../src/ilinear.h"
#include "../src/dfield.h"
#include "../src/demons.h"
#include "../src/gradient_descent.h"
#include "../src/registration.h"
#include "../src/viewer.h"

using namespace imart;

int main(int argc, char *argv[])
{
    // using type = float;
    using type = double;
    using image_type = image_cpu<type>;
    
    if( argc < 4 )
    {
        std::cerr << "Usage: " << std::endl;
        std::cerr << argv[0] << " <fixed> <template> <moving>" << std::endl;
        return EXIT_FAILURE;
    }
    // Files
    std::string file_fixed = argv[1];
    std::string file_template = argv[2];
    std::string file_moving = argv[3];

    // ============================================
    //      Testing tracking with two images
    // ============================================

    // Images
    auto img_fixed = image_cpu<type>::new_pointer();
    auto img_moving = image_cpu<type>::new_pointer();
    auto img_template = image_cpu<type>::new_pointer();
    
    // Read
    img_fixed->read(file_fixed);
    img_moving->read(file_moving);
    img_template->read(file_template);

    // Setup for template matching
    int extra_pixels_bbox = 4;
    auto bbox_fixed = bounding_box(img_template, extra_pixels_bbox);
    auto img_fixed_region = img_fixed->region(bbox_fixed[0],bbox_fixed[1]);

    // Setup for registration
    auto trfm = dfield<type,vector_cpu<type>>::new_pointer(img_fixed_region);
    double sigmaf = 0.0;
    double sigmae = 2.5;
    trfm->set_sigma_fluid(sigmaf);
    trfm->set_sigma_elastic(sigmae);
    
    auto demonsreg = demons<type,vector_cpu<type>>::new_pointer(img_fixed, img_moving, trfm);
    auto opt = gradient_descent<type,vector_cpu<type>>::new_pointer();
    opt->set_tolerance(1e-6);
    
    auto registro = registration<type,vector_cpu<type>>::new_pointer(img_fixed_region, img_fixed_region, trfm);
    registro->set_metric(demonsreg);
    registro->set_optimizer(opt);
    registro->set_levels(3);
    registro->set_levels_scales(std::vector<int>{4,2,1});
    registro->set_levels_iterations(std::vector<int>{120,100,80});
    registro->set_transform(trfm);

    timer t("ms");

    // Template Matching
    t.start();
    auto tm = template_matching<type, vector_cpu<type>>::new_pointer(img_fixed, bbox_fixed);
    tm->set_slide(std::vector<int>{30,30});
    auto bbox_moving = tm->apply(img_moving);

    t.finish();
    printf("Total time: \t%5.2f [ms]\n", t.get_elapsed());
    
    printf("Template bounding box: (%d, %d, %d, %d)\n",
        bbox_fixed[0][0], bbox_fixed[0][1], bbox_fixed[1][0], bbox_fixed[1][1]);
    printf("Matched bounding box: (%d, %d, %d, %d)\n",
        bbox_moving[0][0], bbox_moving[0][1], bbox_moving[1][0], bbox_moving[1][1]);


    // Registration
    t.start();
    auto img_moving_region = img_moving->region(bbox_moving[0],bbox_moving[1]);
    registro->set_moving(img_moving_region);
    registro->apply();
    t.finish();
    printf("Total time: \t%5.2f [ms]\n", t.get_elapsed());

    auto img_template_region = img_template->region(bbox_fixed[0], bbox_fixed[1]);

    // Get template transform
    auto transformation = registro->get_transform()->inverse();
    auto interpolation = ilinear_cpu<type>::new_pointer(img_template_region);
    auto x1 = grid_cpu<type>::new_pointer(img_moving_region);
    auto img_template_region_warped = interpolation->apply(transformation->apply(x1));


    std::vector<int> pre = bbox_moving[0];
    std::vector<int> post(bbox_moving[1].size(),0);

    for (int k = 0; k < img_fixed->get_dimension(); k++)
    {
        post[k] = img_fixed->get_size()[k] - bbox_moving[0][k] - bbox_moving[1][k];
    }

    auto img_template_warped = pad(img_template_region_warped, pre, post);
    img_template_warped->set_spacing(img_template->get_spacing());
    img_template_warped->set_origin(img_template->get_origin());
    img_template_warped->set_direction(img_template->get_direction());

    img_template->print();
    img_template_warped->print();

    img_template->write("template.nii");
    img_template_warped->write("output.nii");


    return 0;
};