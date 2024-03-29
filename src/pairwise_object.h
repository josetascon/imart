/*
* @Author: jose
* @Date:   2020-08-20 00:00:00
* @Last Modified by:   jose
* @Last Modified time: 2020-08-20 00:00:00
*/

#ifndef __PAIRWISE_OBJECT_H__
#define __PAIRWISE_OBJECT_H__

#include "process_object.h"
#include "image.h"
#include "transform.h"
#include "interpolator.h"
// #include "viewer.h"

namespace imart
{

// Class pairwise_object
template <typename type, typename container=vector_cpu<type>>
class pairwise_object: public inherit<pairwise_object<type,container>, process_object>
{

public:
    //Type definitions
    using self    = pairwise_object;
    using pointer = std::shared_ptr<self>;
    using vector  = std::vector<self::pointer>;
    
protected:
    // ===========================================
    // Internal Variables
    // ===========================================
    bool plot;
    typename image<type,container>::pointer fixed;
    typename image<type,container>::pointer moving;
    typename transform<type,container>::pointer transformation;
    typename interpolator<type,container>::pointer interpolation;
    // typename viewer<image<type,container>>::pointer view;

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
    pairwise_object();
    pairwise_object(int d);
    pairwise_object(typename image<type,container>::pointer fixed_image, 
           typename image<type,container>::pointer moving_image);
    pairwise_object(typename image<type,container>::pointer fixed_image,
           typename image<type,container>::pointer moving_image,
           typename transform<type,container>::pointer transformd);

    // ===========================================
    // Get Functions
    // ===========================================
    bool get_plot() const;
    typename image<type,container>::pointer get_fixed() const;
    typename image<type,container>::pointer get_moving() const;
    typename transform<type,container>::pointer get_transform() const;
    typename interpolator<type,container>::pointer get_interpolator() const;
    // typename viewer<image<type,container>>::pointer get_viewer() const;

    // ===========================================
    // Set Functions
    // ===========================================
    void set_plot(bool p);
    virtual void set_fixed(typename image<type,container>::pointer fixed_image);
    virtual void set_moving(typename image<type,container>::pointer moving_image);
    virtual void set_transform(typename transform<type,container>::pointer transformd);
    virtual void set_interpolator(typename interpolator<type,container>::pointer interpolationd);
    // virtual void set_viewer(typename viewer<image<type,container>>::pointer viewd);

    // ===========================================
    // Print Functions
    // ===========================================
    std::string info(std::string msg);
};


// ===========================================
//      Functions of Class pairwise_object
// ===========================================

// ===========================================
// Create Functions
// ===========================================
// Constructor empty
template <typename type, typename container>
pairwise_object<type,container>::pairwise_object()
{
    this->class_name = "pairwise_object";
    init(2);
};

template <typename type, typename container>
pairwise_object<type,container>::pairwise_object(int d)
{
    this->class_name = "pairwise_object";
    init(d);
};

template <typename type, typename container>
pairwise_object<type,container>::pairwise_object( typename image<type,container>::pointer fixed_image, 
                                typename image<type,container>::pointer moving_image)
{
    this->class_name = "pairwise_object";
    auto transformd = transform<type,container>::new_pointer(moving_image->get_dimension());
    init(fixed_image, moving_image, transformd);
};

template <typename type, typename container>
pairwise_object<type,container>::pairwise_object( typename image<type,container>::pointer fixed_image, 
                                typename image<type,container>::pointer moving_image,
                                typename transform<type,container>::pointer transformd)
{
    this->class_name = "pairwise_object";
    init(fixed_image, moving_image, transformd);
};

template <typename type, typename container>
void pairwise_object<type,container>::init(int d)
{
    // std::cout << "init pairwise_object" << std::endl;
    auto fixed_image = image<type,container>::new_pointer(d);
    auto moving_image = image<type,container>::new_pointer(d);
    auto transformd = transform<type,container>::new_pointer(d);
    init(fixed_image, moving_image, transformd);
    // std::cout << "end init pairwise_object" << std::endl;
}

template <typename type, typename container>
void pairwise_object<type,container>::init( typename image<type,container>::pointer fixed_image, 
                                   typename image<type,container>::pointer moving_image,
                                   typename transform<type,container>::pointer transformd)
{
    // std::cout << "pairwise_object init" << std::endl;
    assert(fixed_image->get_dimension() == moving_image->get_dimension());
    assert(fixed_image->get_dimension() == transformd->get_dimension());
    
    fixed = fixed_image;
    moving = moving_image;
    transformation = transformd;
    interpolation = interpolator<type,container>::new_pointer(moving->get_dimension());
    // view = std::make_shared<viewer<image<type,container>>>();
    
    this->set_total_inputs(3);      //process_object::init
    this->set_total_outputs(0);     //process_object::init
    this->setup_input(fixed, moving, transformation);
    this->setup_output();           // no defined output, this object is a base class
    // std::cout << "pairwise_object init end" << std::endl;
    return;
};

// ===========================================
// Get Functions
// ===========================================
template <typename type, typename container>
bool pairwise_object<type,container>::get_plot() const
{
    return plot;
};

template <typename type, typename container>
typename image<type,container>::pointer pairwise_object<type,container>::get_fixed() const
{
    return fixed;
};

template <typename type, typename container>
typename image<type,container>::pointer pairwise_object<type,container>::get_moving() const
{
    return moving;
};

template <typename type, typename container>
typename transform<type,container>::pointer pairwise_object<type,container>::get_transform() const
{
    return transformation;
};

template <typename type, typename container>
typename interpolator<type,container>::pointer pairwise_object<type,container>::get_interpolator() const
{
    return interpolation; // single interpolator
};

// template <typename type, typename container>
// typename viewer<image<type,container>>::pointer pairwise_object<type,container>::get_viewer() const
// {
//     return view;
// };

// ===========================================
// Set Functions
// ===========================================
template <typename type, typename container>
void pairwise_object<type,container>::set_plot(bool p)
{
    plot = p;
};

template <typename type, typename container>
void pairwise_object<type,container>::set_fixed(typename image<type,container>::pointer fixed_image)
{
    fixed = fixed_image;
    this->set_input(0, fixed); // update process input
};

template <typename type, typename container>
void pairwise_object<type,container>::set_moving(typename image<type,container>::pointer moving_image)
{
    moving = moving_image;
    this->set_input(1, moving); // update process input
};

template <typename type, typename container>
void pairwise_object<type,container>::set_transform(typename transform<type,container>::pointer transformd)
{
    transformation = transformd;
    this->set_input(2, transformation); // update process input
};

template <typename type, typename container>
void pairwise_object<type,container>::set_interpolator(typename interpolator<type,container>::pointer interpolationd)
{
    interpolation = interpolationd;
};

// template <typename type, typename container>
// void pairwise_object<type,container>::set_viewer(typename viewer<image<type,container>>::pointer viewd)
// {
//     view = viewd;
// };

// ===========================================
// Print Functions
// ===========================================
template <typename type, typename container>
std::string pairwise_object<type,container>::info(std::string msg)
{
    std::stringstream ss;
    std::string title = "Pairwise Object Information";
    if (msg != "") { title = msg; };

    ss << object::info(title);
    ss << process_object::info("");
    return ss.str();
};

}; //end namespace

#endif