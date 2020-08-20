/*
* @Author: jose
* @Date:   2019-12-09 00:00:00
* @Last Modified by:   jose
* @Last Modified time: 2020-03-18 00:00:00
*/

#ifndef __DEMONS_H__
#define __DEMONS_H__

#include "metric.h"

namespace imart
{

// Class metric
template <typename type, typename container=vector_cpu<type>>
class demons: public inherit<demons<type,container>, metric<type,container>>
{
public:
    //Type definitions
    using self    = demons;
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

    using inherit<demons<type,container>, metric<type,container>>::inherit;

public:
    // ===========================================
    // Create Functions
    // ===========================================
    // Constructors
    demons() : inherit<demons<type,container>, metric<type,container>>()
          { this->class_name = "demons"; };

    demons(int d) : inherit<demons<type,container>, metric<type,container>>(d)
               { this->class_name = "demons"; };

    demons(typename image<type,container>::pointer fixed_image, 
        typename image<type,container>::pointer moving_image)
        : inherit<demons<type,container>, metric<type,container>>(fixed_image, moving_image)
        { this->class_name = "demons"; };

    demons(typename image<type,container>::pointer fixed_image, 
        typename image<type,container>::pointer moving_image,
        typename transform<type,container>::pointer transformd)
        : inherit<demons<type,container>, metric<type,container>>(fixed_image, moving_image, transformd)
        { this->class_name = "demons"; };

    // ===========================================
    // Functions
    // ===========================================
    // !compute the cost
    type cost();
    // !calculate derivative
    typename transform<type,container>::pointer derivative();
};

template<typename type>
using demons_cpu = demons<type,vector_cpu<type>>;

template<typename type>
using demons_gpu = demons<type,vector_ocl<type>>;


// ===========================================
//      Functions of Class metric
// ===========================================

// ===========================================
// Functions
// ===========================================
template <typename type, typename container>
type demons<type,container>::cost()
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
typename transform<type,container>::pointer demons<type,container>::derivative()
{
    int d = fixed->get_dimension();
    // transform
    typename transform<type,container>::pointer trfm = transformation->mimic();
    
    // mimic parameters
    auto param = std::make_shared< std::vector< typename image<type,container>::pointer >>(d);
    for(int i = 0; i < d; i++)
    {
        (*param)[i] = image<type,container>::new_pointer();
        (*param)[i]->mimic_(*(transformation->get_parameters(i)));
    }
    
    // moving prime already computed
    // auto moving_prime = this->warped_moving();        // consider store this to avoid compute again
    type sigma = 2.0;
    type low = 1e-6;

    auto dif = (*moving_prime - *fixed);
    auto dif_sq = dif*dif;

    if (transformation->get_name() == "dfield")
    {
        // std::cout << "der 1" << std::endl;
        // Computing gradient
        auto grad = gradient(moving_prime);

        // dimensional data
        auto gx = *grad[0];
        auto gy = *grad[1];
        
        if (d == 2)
        {
            // compute mag
            auto mag_sq = (gx*gx + gy*gy);

            // write parameters
            *((*param)[0]) = (dif*gx)/(mag_sq + dif_sq + low);
            *((*param)[1]) = (dif*gy)/(mag_sq + dif_sq + low);

            // type sigma = 2.0;
            // (*param)[0] = gaussian_filter((*param)[0],sigma,3);
            // (*param)[1] = gaussian_filter((*param)[1],sigma,3);
        }
        else if (d == 3)
        {
            // z
            auto gz = *grad[2];
            // compute mag
            auto mag_sq = (gx*gx + gy*gy + gz*gz);
            
            // write parameters
            *((*param)[0]) = dif*gx/(mag_sq + dif_sq + low);
            *((*param)[1]) = dif*gy/(mag_sq + dif_sq + low);
            *((*param)[2]) = dif*gz/(mag_sq + dif_sq + low);
        };
    };

    trfm->set_parameters_vector( param );
    return trfm;
};

}; //end namespace

#endif