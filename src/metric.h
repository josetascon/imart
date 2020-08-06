/*
* @Author: jose
* @Date:   2019-12-09 00:00:00
* @Last Modified by:   jose
* @Last Modified time: 2020-03-18 00:00:00
*/

#ifndef __METRIC_H__
#define __METRIC_H__

#include "process_object.h"
#include "image.h"
#include "grid.h"
#include "transform.h"
#include "interpolator.h"
#include "ilinear.h"

namespace imart
{

// Class metric
template <typename type, typename container=vector_cpu<type>>
class metric: public inherit<metric<type,container>, process_object>
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
    type cost_value;
    typename image<type,container>::pointer fixed;
    typename image<type,container>::pointer moving;
    typename transform<type,container>::pointer transformation;
    typename grid<type,container>::pointer x0;
    typename grid<type,container>::pointer x1;
    typename interpolator<type,container>::pointer interpolate_fixed;
    typename interpolator<type,container>::pointer interpolate_moving;
    typename image<type,container>::pointer fixed_prime;
    typename image<type,container>::pointer moving_prime;

    // ===========================================
    // Functions
    // ===========================================
    void init(int d);
    virtual void init( typename image<type,container>::pointer fixed_image, 
                       typename image<type,container>::pointer moving_image,
                       typename transform<type,container>::pointer transformd);

public:
    // ===========================================
    // Create Functions
    // ===========================================
    // Constructors
    metric();
    metric(int d);
    metric(typename image<type,container>::pointer fixed_image, 
           typename image<type,container>::pointer moving_image);
    metric(typename image<type,container>::pointer fixed_image,
           typename image<type,container>::pointer moving_image,
           typename transform<type,container>::pointer transformd);

    // ===========================================
    // Get Functions
    // ===========================================
    type get_cost() const;
    typename image<type,container>::pointer get_fixed() const;
    typename image<type,container>::pointer get_moving() const;
    typename transform<type,container>::pointer get_transform() const;
    typename grid<type,container>::pointer get_grid_fixed() const;
    typename grid<type,container>::pointer get_grid_moving() const;

    // ===========================================
    // Set Functions
    // ===========================================
    void set_fixed(typename image<type,container>::pointer fixed_image);
    void set_moving(typename image<type,container>::pointer moving_image);
    void set_transform(typename transform<type,container>::pointer transformd);

    // ===========================================
    // Print Functions
    // ===========================================
    std::string info(std::string msg);

    // ===========================================
    // Functions
    // ===========================================
    // !compute the cost
    virtual type cost();
    // !update transform and compute cost
    type cost(typename transform<type,container>::pointer transformd);
    // !update fixed and moving image and compute cost
    type cost(typename image<type,container>::pointer fixed_image, 
              typename image<type,container>::pointer moving_image,
              typename transform<type,container>::pointer transformd);

    // !calculate derivative with regards to transform
    virtual typename transform<type,container>::pointer derivative();

    typename grid<type,container>::pointer grid_moving();
    typename grid<type,container>::pointer grid_fixed();
    typename image<type,container>::pointer warped_moving();
    typename image<type,container>::pointer warped_fixed();
    
};


// ===========================================
//      Functions of Class metric
// ===========================================

// ===========================================
// Create Functions
// ===========================================
// Constructor empty
template <typename type, typename container>
metric<type,container>::metric()
{
    this->class_name = "metric";
    init(2);
};

template <typename type, typename container>
metric<type,container>::metric(int d)
{
    this->class_name = "metric";
    init(d);
};

template <typename type, typename container>
metric<type,container>::metric( typename image<type,container>::pointer fixed_image, 
                                typename image<type,container>::pointer moving_image)
{
    this->class_name = "metric";
    auto transformd = transform<type,container>::new_pointer(fixed_image->get_dimension());
    init(fixed_image, moving_image, transformd);
};

template <typename type, typename container>
metric<type,container>::metric( typename image<type,container>::pointer fixed_image, 
                                typename image<type,container>::pointer moving_image,
                                typename transform<type,container>::pointer transformd)
{
    this->class_name = "metric";
    init(fixed_image, moving_image, transformd);
};

template <typename type, typename container>
void metric<type,container>::init(int d)
{
    auto fixed_image = image<type,container>::new_pointer(d);
    auto moving_image = image<type,container>::new_pointer(d);
    auto transformd = transform<type,container>::new_pointer(d);
    init(fixed_image, moving_image, transformd);
}

template <typename type, typename container>
void metric<type,container>::init( typename image<type,container>::pointer fixed_image, 
                                   typename image<type,container>::pointer moving_image,
                                   typename transform<type,container>::pointer transformd)
{
    // std::cout << "metric init" << std::endl;
    assert(fixed_image->get_dimension() == moving_image->get_dimension());
    cost_value = 1e40;
    fixed = fixed_image;
    moving = moving_image;
    transformation = transformd;
    grid_fixed(); // x0 = grid<type,container>::new_pointer(fixed);
    // x1 = grid<type,container>::new_pointer(moving);
    interpolate_fixed = ilinear<type,container>::new_pointer(fixed);
    interpolate_moving = ilinear<type,container>::new_pointer(moving);
    
    this->set_total_inputs(3);      //process_object::init
    this->set_total_outputs(0);     //process_object::init
    this->setup_input(fixed, moving, transformation);
    this->setup_output();           // output is scalar
    // std::cout << "metric init end" << std::endl;
    return;
};

// ===========================================
// Get Functions
// ===========================================
template <typename type, typename container>
type metric<type,container>::get_cost() const
{
    return cost_value;
};

template <typename type, typename container>
typename image<type,container>::pointer metric<type,container>::get_fixed() const
{
    return fixed;
};

template <typename type, typename container>
typename image<type,container>::pointer metric<type,container>::get_moving() const
{
    return moving;
};

template <typename type, typename container>
typename transform<type,container>::pointer metric<type,container>::get_transform() const
{
    return transformation;
};



template <typename type, typename container>
typename grid<type,container>::pointer metric<type,container>::get_grid_fixed() const
{
    return x0; 
};

template <typename type, typename container>
typename grid<type,container>::pointer metric<type,container>::get_grid_moving() const
{
    return x1;
};

// ===========================================
// Set Functions
// ===========================================
template <typename type, typename container>
void metric<type,container>::set_fixed(typename image<type,container>::pointer fixed_image)
{
    fixed = fixed_image;
    interpolate_fixed = ilinear<type,container>::new_pointer(fixed);
};

template <typename type, typename container>
void metric<type,container>::set_moving(typename image<type,container>::pointer moving_image)
{
    moving = moving_image;
    interpolate_moving = ilinear<type,container>::new_pointer(moving);
};

template <typename type, typename container>
void metric<type,container>::set_transform(typename transform<type,container>::pointer transformd)
{
    transformation = transformd;
};

// ===========================================
// Print Functions
// ===========================================
template <typename type, typename container>
std::string metric<type,container>::info(std::string msg)
{
    std::stringstream ss;
    std::string title = "Metric Information";
    if (msg != "") { title = msg; };

    ss << object::info(title);
    ss << process_object::info("");
    return ss.str();
};

// ===========================================
// Functions
// ===========================================
template <typename type, typename container>
type metric<type,container>::cost()
{
    return cost_value;
};

template <typename type, typename container>
type metric<type,container>::cost(typename transform<type,container>::pointer transformd)
{
    transformation = transformd;
    return cost();
};

template <typename type, typename container>
type metric<type,container>::cost( typename image<type,container>::pointer fixed_image, 
                                   typename image<type,container>::pointer moving_image,
                                   typename transform<type,container>::pointer transformd)
{
    init(fixed_image, moving_image, transformd);
    return cost();
};

template <typename type, typename container>
typename transform<type,container>::pointer metric<type,container>::derivative()
{
    return transformation->mimic(); // derivate
};

template <typename type, typename container>
typename grid<type,container>::pointer metric<type,container>::grid_moving()
{
    x1 = grid<type,container>::new_pointer(moving);
    return x1;
};

template <typename type, typename container>
typename grid<type,container>::pointer metric<type,container>::grid_fixed()
{
    x0 = grid<type,container>::new_pointer(fixed); // computed during init
    return x0;
};

template <typename type, typename container>
typename image<type,container>::pointer metric<type,container>::warped_moving()
{
    moving_prime = interpolate_moving->apply(transformation->apply(x0));
    return moving_prime;
};

template <typename type, typename container>
typename image<type,container>::pointer metric<type,container>::warped_fixed()
{
    fixed_prime = interpolate_fixed->apply(transformation->apply(x1));
    return fixed_prime;
};

}; //end namespace

#endif