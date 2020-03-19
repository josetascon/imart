/*
* @Author: jose
* @Date:   2019-12-09 00:00:00
* @Last Modified by:   jose
* @Last Modified time: 2020-03-18 00:00:00
*/

#ifndef __SSD_H__
#define __SSD_H__

#include "metric.h"

template <typename pixel_type>
class ssd: public metric<pixel_type>
{
public:
    //Type definitions
    using self    = ssd;
    using pointer = std::shared_ptr<self>;
    using vector  = std::vector<self::pointer>;

    //Parent class properties to be used here (avoid use "this->")
    using metric<pixel_type>::cost_value;
    using metric<pixel_type>::fixed;
    using metric<pixel_type>::moving;
    using metric<pixel_type>::transform;
    using metric<pixel_type>::x0;
    using metric<pixel_type>::x1;
    using metric<pixel_type>::interpolator0;
    using metric<pixel_type>::interpolator1;


protected:
    // ===========================================
    // Functions
    // ===========================================
    void init( typename image_base<pixel_type>::pointer fixed_image, 
               typename image_base<pixel_type>::pointer moving_image,
               typename transform_base<pixel_type>::pointer transform_reg );

public:
    // ===========================================
    // Create Functions
    // ===========================================
    // Constructors
    ssd(typename image_base<pixel_type>::pointer fixed_image, 
        typename image_base<pixel_type>::pointer moving_image, 
        typename transform_base<pixel_type>::pointer transform_reg) 
        : metric<pixel_type> (fixed_image, moving_image, transform_reg) 
        { this->class_name = "ssd"; };

    template<typename... ARGS>
    static pointer new_pointer(const ARGS&... args);

    // ===========================================
    // Functions
    // ===========================================
    // !compute the cost
    pixel_type cost();
    // !calculate derivative
    transform_base<pixel_type> derivative();
};


// ===========================================
//      Functions of Class metric
// ===========================================

template <typename pixel_type>
void ssd<pixel_type>::init( typename image_base<pixel_type>::pointer fixed_image, 
                            typename image_base<pixel_type>::pointer moving_image,
                            typename transform_base<pixel_type>::pointer transform_reg)
{
    assert(fixed_image->get_dimension() == moving_image->get_dimension());
    metric<pixel_type>::init(fixed_image, moving_image, transform_reg);
    this->class_name = "ssd";
};

template <typename pixel_type>
template <typename ... ARGS>
typename ssd<pixel_type>::pointer ssd<pixel_type>::new_pointer(const ARGS&... args)
{
    return std::make_shared< ssd<pixel_type> >(args...); // not working for inherited classes
};


template <typename pixel_type>
pixel_type ssd<pixel_type>::cost()
{
    // this->x0.print_data();
    // this->x1.print_data();

    // interpolate<pixel_type> image0_p(fixed, x0);
    // interpolate<pixel_type> image1_p(moving, x1);

    // transform->print_data();
    
    // grid<pixel_type> xx = transform->transform(*x1);
    // xx.print_data();
    // image_base<pixel_type> tmp = image0_p*(transform->transform(x1));
    // tmp.print_data();
    // image_base<pixel_type> ssd_ = moving - image0_p*(transform->transform(x1));
    // // ssd_.print_data();
    // cost_value = ssd_.sum();

    image_base<pixel_type> ssd_ = *moving - interpolator0->linear(transform->transform(*x1));
    cost_value = ssd_.sum();

    return cost_value;
};

template <typename pixel_type>
transform_base<pixel_type> ssd<pixel_type>::derivative()
{
    return *(transform);
};

#endif