/*
* @Author: jose
* @Date:   2019-12-09 00:00:00
* @Last Modified by:   jose
* @Last Modified time: 2020-03-18 00:00:00
*/

#ifndef __SSD_H__
#define __SSD_H__

#include "metric.h"

namespace imart
{

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
    typename transform_base<pixel_type>::pointer derivative();
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
    // std::cout << "ssd init" << std::endl;
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

    image_base<pixel_type> ssd_ = *fixed - interpolator1->linear(transform->transform(*x0));
    // image_base<pixel_type> ssd_ = *moving - interpolator0->linear(transform->transform(*x1));
    cost_value = ssd_.sum();

    return cost_value;
};

template <typename pixel_type>
typename transform_base<pixel_type>::pointer ssd<pixel_type>::derivative()
{
    if (transform->get_name() == "affine")
    {
        int d = fixed->get_dimension();
        int N = fixed->get_total_elements();

        image_base<pixel_type> x = (x0->ptr())[0];
        image_base<pixel_type> y = (x0->ptr())[1];
        // x.print_data();

        auto dm_di = image_base<pixel_type>::new_pointer(d);
        auto di_dp = image_base<pixel_type>::new_pointer(d);

        grid<pixel_type> x_new = transform->transform(*x0);
        // x_new.print_data();

        *dm_di = interpolator1->linear(x_new) - *fixed;
        // dm_di->print_data();

        // Computing gradient
        // std::cout << "\tGradient" << std::endl;
        typename image_base<pixel_type>::vector grad(d);
        grad = gradient(*moving);
        // grad[0]->print_data();
        // std::cout << "image max: " << moving->max() << std::endl;
        // std::cout << "grad x max: " << grad[0]->max() << std::endl;


        auto gradx = image<pixel_type>::new_pointer(*grad[0]); // image_base casting to image
        auto grady = image<pixel_type>::new_pointer(*grad[1]);

        // interpolate<pixel_type> grad_interp_x(grad[0], x1);
        // interpolate<pixel_type> grad_interp_y(grad[1], x1);
        // auto grad_interp_x = interpolate<pixel_type>::new_pointer(grad[0], x1);
        // auto grad_interp_y = interpolate<pixel_type>::new_pointer(grad[1], x1);
        auto grad_interp_x = interpolate<pixel_type>::new_pointer(gradx, x1);
        auto grad_interp_y = interpolate<pixel_type>::new_pointer(grady, x1);

        auto gx = grad_interp_x->linear(x_new);
        auto gy = grad_interp_y->linear(x_new);
        // gx.print_data();

        // Affine parameters
        // pixel_type p0 = 0, p1 = 0, p2 = 0;
        // pixel_type p3 = 0, p4 = 0, p5 = 0;
        // std::cout << "\tValues" << std::endl;
        auto param = image_base<pixel_type>::new_pointer(6,1);
        pixel_type * p = param->ptr();
        // param->print_data();

        // std::cout<< dm_di->dot(gx*x);
        // std::cout << ( dm_di->dot(gx*x) )*(1.0/N);

        p[0] = ( dm_di->dot(gx*x) )*(1.0/N);
        p[1] = (1.0/N)*( dm_di->dot(gx*y) );
        p[2] = (1.0/N)*( dm_di->dot(gy*x) );
        p[3] = (1.0/N)*( dm_di->dot(gy*y) );
        p[4] = (1.0/N)*( dm_di->dot(gx) );
        p[5] = (1.0/N)*( dm_di->dot(gy) );

        // std::cout<< N;
        // std::cout<< p[0];
        
        // std::cout << "\tTransform" << std::endl;
        auto trfm = affine<pixel_type>::new_pointer();
        trfm->imitate(*transform);
        trfm->set_parameters( param );
        return trfm;
    };

    return transform;
};

}; //end namespace

#endif