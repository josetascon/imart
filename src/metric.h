/*
* @Author: jose
* @Date:   2019-12-09 00:00:00
* @Last Modified by:   jose
* @Last Modified time: 2020-03-18 00:00:00
*/

#ifndef __METRIC_H__
#define __METRIC_H__

#include "image_base.h"
#include "transform_base.h"

template <typename pixel_type>
class metric: public object<pixel_type>
{
public:
    //Type definitions
    using self    = metric;
    using pointer = std::shared_ptr<self>;
    using vector  = std::vector<self::pointer>;
    
protected:
    // ===========================================
    // Internal Variables
    // ===========================================
    pixel_type cost_value;
    typename image_base<pixel_type>::pointer fixed;
    typename image_base<pixel_type>::pointer moving;
    typename transform_base<pixel_type>::pointer transform;
    typename grid<pixel_type>::pointer x0;
    typename grid<pixel_type>::pointer x1;
    typename interpolate<pixel_type>::pointer interpolator0;
    typename interpolate<pixel_type>::pointer interpolator1;

    // ===========================================
    // Functions
    // ===========================================
    virtual void init( typename image_base<pixel_type>::pointer fixed_image, 
                       typename image_base<pixel_type>::pointer moving_image,
                       typename transform_base<pixel_type>::pointer transform_reg );

public:
    // ===========================================
    // Create Functions
    // ===========================================
    // Constructors
    metric( typename image_base<pixel_type>::pointer fixed_image, 
            typename image_base<pixel_type>::pointer moving_image, 
            typename transform_base<pixel_type>::pointer transform_reg );

    template<typename... ARGS>
    static pointer new_pointer(const ARGS&... args);

    // ===========================================
    // Get Functions
    // ===========================================
    pixel_type get_cost() const;
    typename image_base<pixel_type>::pointer get_fixed() const;
    typename image_base<pixel_type>::pointer get_moving() const;
    typename transform_base<pixel_type>::pointer get_transform() const;

    // ===========================================
    // Print Functions
    // ===========================================
    std::string info(std::string msg);
    std::string info_data(std::string msg);

    // ===========================================
    // Functions
    // ===========================================
    // !compute the cost
    virtual pixel_type cost();
    // !calculate derivative
    virtual transform_base<pixel_type> derivative();
};


// ===========================================
//      Functions of Class metric
// ===========================================

// ===========================================
// Create Functions
// ===========================================
// Constructor
template <typename pixel_type>
metric<pixel_type>::metric( typename image_base<pixel_type>::pointer fixed_image, 
                            typename image_base<pixel_type>::pointer moving_image, 
                            typename transform_base<pixel_type>::pointer transform_reg)
{
    // std::cout << "metric\n";
    init(fixed_image, moving_image, transform_reg);
};

template <typename pixel_type>
void metric<pixel_type>::init( typename image_base<pixel_type>::pointer fixed_image, 
                               typename image_base<pixel_type>::pointer moving_image,
                               typename transform_base<pixel_type>::pointer transform_reg)
{
    assert(fixed_image->get_dimension() == moving_image->get_dimension());
    this->class_name = "metric";
    cost_value = 1e40;
    fixed = fixed_image;
    moving = moving_image;
    transform = transform_reg;
    x0 = grid<pixel_type>::new_pointer(*fixed);
    x1 = grid<pixel_type>::new_pointer(*moving);
    interpolator0 = interpolate<pixel_type>::new_pointer(*fixed, *x0);
    interpolator1 = interpolate<pixel_type>::new_pointer(*moving, *x1);

};

template <typename pixel_type>
template <typename ... ARGS>
typename metric<pixel_type>::pointer metric<pixel_type>::new_pointer(const ARGS&... args)
{
    return std::make_shared< metric<pixel_type> >(args...); // not working for inherited classes
};

// ===========================================
// Get Functions
// ===========================================
template <typename pixel_type>
pixel_type metric<pixel_type>::get_cost() const
{
    return cost_value;
};

template <typename pixel_type>
typename image_base<pixel_type>::pointer metric<pixel_type>::get_fixed() const
{
    return fixed;
};

template <typename pixel_type>
typename image_base<pixel_type>::pointer metric<pixel_type>::get_moving() const
{
    return moving;
};

template <typename pixel_type>
typename transform_base<pixel_type>::pointer metric<pixel_type>::get_transform() const
{
    return transform;
};

// ===========================================
// Print Functions
// ===========================================
template <typename pixel_type>
std::string metric<pixel_type>::info(std::string msg)
{
    std::stringstream ss;
    std::string title = "Metric Information";
    if (msg != "") { title = msg; };
    // Summary of the metric information
    ss << object<pixel_type>::info(title);
    ss << std::endl;

    return ss.str();
};

template <typename pixel_type>
std::string metric<pixel_type>::info_data(std::string msg)
{
    std::stringstream ss;
    if (msg != "") { ss << msg << std::endl; }
    else { ss << "Metric parameters:\n"; };
    ss << "Fixed Image: \t\t";
    ss << fixed->ptr() << std::endl;
    ss << "Moving Image: \t\t";
    ss << moving->ptr() << std::endl;
    ss << "Transform: \t\t";
    ss << transform->get_parameters() << std::endl;

    return ss.str();
};

template <typename pixel_type>
pixel_type metric<pixel_type>::cost()
{
    return cost_value;
};

template <typename pixel_type>
transform_base<pixel_type> metric<pixel_type>::derivative()
{
    return *transform;
};

#endif