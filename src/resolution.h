/*
* @Author: jose
* @Date:   2020-06-21 00:00:00
* @Last Modified by:   jose
* @Last Modified time: 2020-06-21 00:00:00
*/

#ifndef __RESOLUTION_H__
#define __RESOLUTION_H__

#include "process_object.h"
#include "image.h"
#include "transform.h"
#include "interpolator.h"
#include "ilinear.h"

namespace imart
{

// Class resolution
template <typename type, typename container=vector_cpu<type>>
class resolution: public inherit<resolution<type,container>, process_object>
{

public:
    //Type definitions
    using self    = resolution;
    using pointer = std::shared_ptr<self>;
    using vector  = std::vector<self::pointer>;
    
    using inherit<resolution<type,container>, process_object>::inherit;
    
protected:
    // ===========================================
    // Internal Variables
    // ===========================================
    typename image<type,container>::pointer in;
    typename image<type,container>::pointer out;
    typename interpolator<type,container>::pointer interpolation;
    double scale;

public:
    // ===========================================
    // Create Functions
    // ===========================================
    // Constructors
    resolution(typename image<type,container>::pointer input);

    // ===========================================
    // Get Functions
    // ===========================================
    double get_scale() const;
    typename image<type,container>::pointer get_image() const;
    typename interpolator<type,container>::pointer get_interpolator() const;

    // ===========================================
    // Set Functions
    // ===========================================
    void set_scale(double scalar);
    void set_image(typename image<type,container>::pointer input);
    void set_interpolator(typename interpolator<type,container>::pointer interpolationd);

    // ===========================================
    // Print Functions
    // ===========================================
    std::string info(std::string msg);

    // ===========================================
    // Functions
    // ===========================================
    // !apply the scale
    typename image<type,container>::pointer apply();
    typename image<type,container>::pointer apply(double scalar);
    typename transform<type,container>::pointer apply(double scalar, typename transform<type,container>::pointer trfm);
    typename transform<type,container>::pointer apply(std::vector<int> sz, std::vector<double> space, typename transform<type,container>::pointer trfm);
};


// ===========================================
//      Functions of Class resolution
// ===========================================

// ===========================================
// Create Functions
// ===========================================
// Constructor empty
template <typename type, typename container>
resolution<type,container>::resolution(typename image<type,container>::pointer input)
{
    this->class_name = "resolution";
    in = input;
    out = image<type,container>::new_pointer(input->get_dimension());
    interpolation = ilinear<type,container>::new_pointer(in);
    scale = 1.0;
    
    // setup process
    this->set_total_inputs(2);      //process_object::init
    this->set_total_outputs(0);     //process_object::init
    this->setup_input(in, out);
    this->setup_output();           // output is scalar
};

// ===========================================
// Get Functions
// ===========================================
template <typename type, typename container>
double resolution<type,container>::get_scale() const
{
    return scale;
};

template <typename type, typename container>
typename image<type,container>::pointer resolution<type,container>::get_image() const
{
    return in;
};

template <typename type, typename container>
typename interpolator<type,container>::pointer resolution<type,container>::get_interpolator() const
{
    return interpolation;
};

// ===========================================
// Set Functions
// ===========================================
template <typename type, typename container>
void resolution<type,container>::set_scale(double scalar)
{
    scale = scalar;
};

template <typename type, typename container>
void resolution<type,container>::set_image(typename image<type,container>::pointer input)
{
    in = input;
    interpolation->set_reference(in);
    // interpolation = ilinear<type,container>::new_pointer(in);
};

template <typename type, typename container>
void resolution<type,container>::set_interpolator(typename interpolator<type,container>::pointer interpolationd)
{
    interpolation = interpolationd->mimic();
    interpolation->set_reference(in);
};

// ===========================================
// Print Functions
// ===========================================
template <typename type, typename container>
std::string resolution<type,container>::info(std::string msg)
{
    std::stringstream ss;
    std::string title = "Resolution Information";
    if (msg != "") { title = msg; };

    ss << object::info(title);
    ss << process_object::info("");
    ss << "Scale: \t\t" << scale << std::endl;
    return ss.str();
};

// ===========================================
// Functions
// ===========================================
template <typename type, typename container>
typename image<type,container>::pointer resolution<type,container>::apply()
{
    std::vector<int> sz = in->get_size();
    std::vector<double> space = in->get_spacing();

    for(int i = 0; i < sz.size(); i++)
    {
        // sz[i] = 1 + (sz[i]-1)/scale;
        sz[i] = std::round(sz[i]/scale);
        space[i] = space[i]*scale;
    };

    out = image<type,container>::new_pointer(sz);
    out->set_spacing(space);
    out->set_origin(in->get_origin());
    out->set_direction(in->get_direction());
    // out->print();

    auto x = grid<type,container>::new_pointer(out);
    // x->ptr()[0]->print_data();
    // std::cout << "interpolation" << std::endl;
    return interpolation->apply(x);
};

template <typename type, typename container>
typename image<type,container>::pointer resolution<type,container>::apply(double scalar)
{
    set_scale(scalar);
    return apply();
};

template <typename type, typename container>
typename transform<type,container>::pointer resolution<type,container>::apply(double scalar, typename transform<type,container>::pointer trfm)
{
    set_scale(scalar);

    std::vector<int> sz = trfm->get_size();
    std::vector<double> space = trfm->get_spacing();

    for(int i = 0; i < sz.size(); i++)
    {
        // sz[i] = 1 + (sz[i]-1)/scale;
        sz[i] = std::round(sz[i]/scale);
        space[i] = space[i]*scale;
    };
    return apply(sz,space,trfm);
};

template <typename type, typename container>
typename transform<type,container>::pointer resolution<type,container>::apply(std::vector<int> sz, std::vector<double> space, typename transform<type,container>::pointer trfm)
{
    auto out_trfm = trfm->mimic();
    out_trfm->change_size(sz);
    out_trfm->set_spacing(space);
    out_trfm->set_origin(trfm->get_origin());
    out_trfm->set_direction(trfm->get_direction());

    auto ref = image<type,container>::new_pointer(sz);
    ref->set_spacing(space);
    ref->set_origin(in->get_origin());
    ref->set_direction(in->get_direction());

    auto intertrfm = interpolation->mimic();

    for(int i = 0; i < sz.size(); i ++)
    {
        auto x = grid<type,container>::new_pointer(ref);
        intertrfm->set_reference(trfm->get_parameters(i));
        out_trfm->set_parameters(intertrfm->apply(x),i);
    };
    return out_trfm;
};

}; //end namespace

#endif