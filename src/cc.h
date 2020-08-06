/*
* @Author: jose
* @Date:   2020-06-21 00:00:00
* @Last Modified by:   jose
* @Last Modified time: 2020-06-21 00:00:00
*/

#ifndef __CC_H__
#define __CC_H__

#include "metric.h"

namespace imart
{

// Class metric
template <typename type, typename container=vector_cpu<type>>
class cc: public inherit<cc<type,container>, metric<type,container>>
{
public:
    //Type definitions
    using self    = cc;
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

    using inherit<cc<type,container>, metric<type,container>>::inherit;

public:
    // ===========================================
    // Create Functions
    // ===========================================
    // Constructors
    cc() : inherit<cc<type,container>, metric<type,container>>()
         { this->class_name = "cc"; };

    cc(int d) : inherit<cc<type,container>, metric<type,container>>(d)
              { this->class_name = "cc"; };

    cc(typename image<type,container>::pointer fixed_image, 
       typename image<type,container>::pointer moving_image)
       : inherit<cc<type,container>, metric<type,container>>(fixed_image, moving_image)
       { this->class_name = "cc"; };

    cc(typename image<type,container>::pointer fixed_image, 
       typename image<type,container>::pointer moving_image,
       typename transform<type,container>::pointer transformd)
       : inherit<cc<type,container>, metric<type,container>>(fixed_image, moving_image, transformd)
       { this->class_name = "cc"; };

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
type cc<type,container>::cost()
{
    auto moving_prime = this->warped_moving();
    
    // using elastix definition
    type N = (type)fixed->get_total_elements();
    type f_mean = fixed->sum()/N;
    type m_mean = moving_prime->sum()/N;
    // std::cout << "f mean: " << f_mean << std::endl;
    // std::cout << "m mean: " << m_mean << std::endl;

    auto f_dif = *fixed - f_mean;
    auto m_dif = *moving_prime - m_mean;

    type num = (f_dif*m_dif).sum();
    type den = sqrt(((f_dif^2.0).sum())*((m_dif^2.0).sum()));
    // std::cout << "numerator: " << num << std::endl;
    // std::cout << "denominator: " << den << std::endl;

    cost_value = num/den;
    return cost_value;
};

template <typename type, typename container>
typename transform<type,container>::pointer cc<type,container>::derivative()
{
    return transformation->mimic();
};

}; //end namespace

#endif