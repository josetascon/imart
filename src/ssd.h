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

// Class metric
template <typename type, typename container=vector_cpu<type>>
class ssd: public inherit<ssd<type,container>, metric<type,container>>
{
public:
    //Type definitions
    using self    = ssd;
    using pointer = std::shared_ptr<self>;
    using vector  = std::vector<self::pointer>;

    // Inherited variables
    using metric<type,container>::cost_value; //Parent class properties to be used here (avoid use "this->")
    using metric<type,container>::fixed;
    using metric<type,container>::moving;
    using metric<type,container>::transformation;
    using metric<type,container>::x0;
    using metric<type,container>::x1;
    using metric<type,container>::cost;
    using metric<type,container>::fixed_prime;
    using metric<type,container>::moving_prime;

    using inherit<ssd<type,container>, metric<type,container>>::inherit;

public:
    // ===========================================
    // Create Functions
    // ===========================================
    // Constructors
    ssd() : inherit<ssd<type,container>, metric<type,container>>()
          { this->class_name = "ssd"; };

    ssd(int d) : inherit<ssd<type,container>, metric<type,container>>(d)
               { this->class_name = "ssd"; };

    ssd(typename image<type,container>::pointer fixed_image, 
        typename image<type,container>::pointer moving_image)
        : inherit<ssd<type,container>, metric<type,container>>(fixed_image, moving_image)
        { this->class_name = "ssd"; };

    ssd(typename image<type,container>::pointer fixed_image, 
        typename image<type,container>::pointer moving_image,
        typename transform<type,container>::pointer transformd)
        : inherit<ssd<type,container>, metric<type,container>>(fixed_image, moving_image, transformd)
        { this->class_name = "ssd"; };

    // ===========================================
    // Functions
    // ===========================================
    // !compute the cost
    type cost();
    // !calculate derivative
    typename transform<type,container>::pointer derivative();
};


// ===========================================
//      Functions of Class metric
// ===========================================

// ===========================================
// Functions
// ===========================================
template <typename type, typename container>
type ssd<type,container>::cost()
{
    this->warped_moving(); // update moving_prime
    // auto moving_prime = this->warped_moving();
    // moving_prime->print_data("moving prime");
    type N = (type)fixed->get_total_elements();
    image<type,container> ssd_ = *fixed - *moving_prime;
    cost_value = (0.5/N)*( (ssd_^2).sum() );
    return cost_value;
};

template <typename type, typename container>
typename transform<type,container>::pointer ssd<type,container>::derivative()
{
    int d = fixed->get_dimension();
    type N = (type)fixed->get_total_elements();

    auto trfm = transformation->mimic();
    auto param = transformation->get_parameters()->mimic();
    
    // moving prime already computed
    // auto moving_prime = this->warped_moving();        // consider store this to avoid compute again
    auto dm_di = image<type,container>::new_pointer();
    *dm_di = (*moving_prime - *fixed);                // **** mising product 1/N

    if (transformation->get_name() == "affine")
    {
        // Computing gradient
        auto grad = gradient(moving_prime);

        // dimensional data
        auto x = *(x0->ptr()[0]);
        auto y = *(x0->ptr()[1]);
        auto gx = *grad[0];
        auto gy = *grad[1];
        
        if (d == 2)
        {
            // compute derivatives
            std::vector<type> p(6);
            p[0] = (1.0/N)*( dm_di->dot(gx*x) );
            p[1] = (1.0/N)*( dm_di->dot(gx*y) );
            p[2] = (1.0/N)*( dm_di->dot(gy*x) );
            p[3] = (1.0/N)*( dm_di->dot(gy*y) );
            p[4] = (1.0/N)*( dm_di->dot(gx) );
            p[5] = (1.0/N)*( dm_di->dot(gy) );
            // write to parameters
            param->get_data()->read_ram(p.data(),p.size());
        }
        else if (d == 3)
        {
            auto z = *(x0->ptr()[2]);
            auto gz = *grad[2];
            // compute derivatives
            std::vector<type> p(12);
            p[0]  = (1.0/N)*( dm_di->dot(gx*x) );
            p[1]  = (1.0/N)*( dm_di->dot(gx*y) );
            p[2]  = (1.0/N)*( dm_di->dot(gx*z) );
            p[3]  = (1.0/N)*( dm_di->dot(gy*x) );
            p[4]  = (1.0/N)*( dm_di->dot(gy*y) );
            p[5]  = (1.0/N)*( dm_di->dot(gy*z) );
            p[6]  = (1.0/N)*( dm_di->dot(gz*x) );
            p[7]  = (1.0/N)*( dm_di->dot(gz*y) );
            p[8]  = (1.0/N)*( dm_di->dot(gz*z) );
            p[9]  = (1.0/N)*( dm_di->dot(gx) );
            p[10] = (1.0/N)*( dm_di->dot(gy) );
            p[11] = (1.0/N)*( dm_di->dot(gz) );
            // write to parameters
            param->get_data()->read_ram(p.data(),p.size());
        };
    };

    trfm->set_parameters( param );
    return trfm;

};

}; //end namespace

#endif