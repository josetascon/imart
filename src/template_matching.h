/*
* @Author: jose
* @Date:   2019-12-09 00:00:00
* @Last Modified by:   jose
* @Last Modified time: 2020-03-18 00:00:00
*/

#ifndef __template_matching_H__
#define __template_matching_H__

#include "tracker.h"

namespace imart
{

// Class tracker
template <typename type, typename container=vector_cpu<type>>
class template_matching: public inherit<template_matching<type,container>, tracker<type,container>>
{
public:
    //Type definitions
    using self    = template_matching;
    using pointer = std::shared_ptr<self>;
    using vector  = std::vector<self::pointer>;

    // Inherited variables
    using tracker<type,container>::fixed; //Parent class properties to be used here (avoid use "this->")
    using tracker<type,container>::moving;
    using tracker<type,container>::current;
    using tracker<type,container>::box_fixed;
    using tracker<type,container>::box_moving;
    using tracker<type,container>::box_current;

    using inherit<template_matching<type,container>, tracker<type,container>>::inherit;

    std::vector<int> slide;

    // ===========================================
    // Functions
    // ===========================================
    std::vector<int> sliding_window(typename image<type,container>::pointer img_box,
                                    std::vector<int> corner);

public:
    // ===========================================
    // Create Functions
    // ===========================================
    // Constructors
    template_matching() : inherit<template_matching<type,container>, tracker<type,container>>()
          { this->class_name = "template_matching"; slide = std::vector<int>(2,10); };

    template_matching(int d) : inherit<template_matching<type,container>, tracker<type,container>>(d)
               { this->class_name = "template_matching"; slide = std::vector<int>(d,10); };

    template_matching(typename image<type,container>::pointer fixed_image, 
        std::vector<std::vector<int>> bbox_fixed)
        : inherit<template_matching<type,container>, tracker<type,container>>(fixed_image, bbox_fixed)
        { this->class_name = "template_matching"; slide = std::vector<int>(fixed_image->get_dimension(),10); };

    // ===========================================
    // Set/Get Functions
    // ===========================================
    std::vector<int> get_slide() const { return slide; };
    void set_slide(std::vector<int> ss) { slide = ss; };

    // ===========================================
    // Functions
    // ===========================================
    // !apply
    std::vector<std::vector<int>> apply(typename image<type,container>::pointer moving_image);
    
};


// ===========================================
//      Functions of Class tracker
// ===========================================

// ===========================================
// Functions
// ===========================================
template <typename type, typename container>
std::vector<std::vector<int>> template_matching<type,container>::apply(typename image<type,container>::pointer moving_image)
{
    moving = moving_image;
    typename image<type,container>::pointer img_box;
    std::vector<int> corner, rsize;
    
    bool run_fixed = true;
    if (this->current_mode)
    {
        run_fixed = false;

        // type ds1 = sum_squared_difference(fixed, moving);
        // type ds2 = sum_squared_difference(fixed, current);
        // printf("f-m %f, f-c %f\n",ds1, ds2);

        bool diff = true;
        type ccfm = cross_correlation(fixed, moving);
        // printf("cc %f\n",ccfm);
        if (ccfm > 0.9989){ diff = false; };

        if (diff)
        {
            // printf("Current\n");
            corner = box_current[0];
            rsize = box_current[1];
            img_box = current->region(corner, rsize);
        }
        else
            run_fixed = true;
    }
    if (run_fixed)
    {
        // printf("Fixed\n");
        corner = box_fixed[0];
        rsize = box_fixed[1];
        img_box = fixed->region(corner, rsize);
    }
    
    auto new_corner = sliding_window(img_box, corner);
    box_moving = std::vector<std::vector<int>>{new_corner, rsize};

    // update current image and box
    current = moving_image->clone();
    box_current = box_moving;

    return box_moving;
};

template <typename type, typename container>
std::vector<int> template_matching<type,container>::sliding_window(typename image<type,container>::pointer img_box,
                                std::vector<int> corner)
{
    int si = slide[0];
    int sj = slide[1];

    // Verify if sliding window is within the image coordinates
    if (corner[0] - si < 0) si = corner[0];
    if (corner[1] - sj < 0) sj = corner[1];
    // if (corner[0] + si > img_box->get_width() ) si = img_box->get_width()  - corner[0];
    // if (corner[1] + sj > img_box->get_height()) sj = img_box->get_height() - corner[1];

    int w = 2*si + 1;
    int h = 2*sj + 1;

    // fix for cpu
    auto corr = image<type,vector_cpu<type>>::new_pointer(w,h);
    type * p = corr->get_data()->data();

    #pragma omp parallel for collapse(2)
    for(int j = 0; j < h; j++)
    {
        for(int i = 0; i < w; i++)
        {
            int x = i + corner[0] - si;
            int y = j + corner[1] - sj;
            
            auto moving_region = moving->region(std::vector<int>{x,y}, img_box->get_size());
            p[i + j*w] = cross_correlation(moving_region, img_box);
        }
    }

    int max_idx = corr->argmax();
    int y = (int)max_idx / w;
    int x = max_idx % w;

    return std::vector<int>{x - si + corner[0], y - sj + corner[1]};
};

}; //end namespace

#endif