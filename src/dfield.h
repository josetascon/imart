/*
* @Author: jose
* @Date:   2019-11-15 13:55:00
* @Last Modified by:   jose
* @Last Modified time: 2019-11-05 13:55:00
*/

#ifndef __DFIELD_H__
#define __DFIELD_H__

// local libs
#include "transform.h"

namespace imart
{

// Class dfield
template <typename type, typename container=vector_cpu<type>>
class dfield: public inherit<dfield<type,container>, transform<type,container>>
{
public:
    //Type definitions
    using self    = dfield;
    using pointer = std::shared_ptr<self>;
    using vector  = std::vector<self::pointer>;

    // Inherited variables
    using transform<type,container>::parameters;
    using transform<type,container>::inverse_parameters;
    using transform<type,container>::allocate;
    using transform<type,container>::get_parameters;
    using transform<type,container>::get_inverse_parameters;
    using transform<type,container>::get_parameters_vector;
    using transform<type,container>::get_inverse_parameters_vector;

    // using transform<type,container>::operator=;
    // using transform<type,container>::operator+;

    using inherit<dfield<type,container>, transform<type,container>>::inherit;

protected:
    // ===========================================
    // Functions
    // ===========================================
    void init(int d);
    void init(std::vector<int> sz);
    void init(const space_object & input);
    void inverse_();

    std::vector<type> transform_2d(std::vector<type> & point);
    std::vector<type> transform_3d(std::vector<type> & point);
    
    typename grid<type,container>::pointer transform_2d(const typename grid<type,container>::pointer input);
    typename grid<type,container>::pointer transform_3d(const typename grid<type,container>::pointer input);

public:
    // ===========================================
    // Constructor Functions
    // ===========================================
    dfield();
    dfield(int d);
    dfield(int d, typename transform<type,container>::ptr_vector_image params);
    dfield(std::vector<int> sz);
    dfield(typename image<type,container>::pointer input);
    dfield(typename grid<type,container>::pointer input);
    dfield(const dfield<type,container> & input);

    // ===========================================
    // Functions
    // ===========================================
    void change_size(std::vector<int> sz); // special function for multiresolution
    void identity();
    std::vector<type> apply(std::vector<type> & point);
    typename grid<type,container>::pointer apply(const typename grid<type,container>::pointer input);
};

template<typename type>
using dfield_cpu = dfield<type,vector_cpu<type>>;

template<typename type>
using dfield_gpu = dfield<type,vector_ocl<type>>;

template<typename type>
using dfield_cuda = dfield<type,vector_cuda<type>>;


// ===========================================
//          Functions of Class dfield
// ===========================================

// ===========================================
// Constructor Functions
// ===========================================
template <typename type, typename container>
dfield<type,container>::dfield()
{
    this->class_name = "dfield";
    init(2);
    identity();
};

template <typename type, typename container>
dfield<type,container>::dfield(int d)
{
    // std::cout << "dfield" << std::endl;
    this->class_name = "dfield";
    init(d);
    identity();
    // std::cout << "created" << std::endl;
};

template <typename type, typename container>
dfield<type,container>::dfield(std::vector<int> sz)
{
    this->class_name = "dfield";
    init(sz.size());
    init(sz);
    identity();
};

template <typename type, typename container>
dfield<type,container>::dfield(int d, typename transform<type,container>::ptr_vector_image params)
{
    // assert(params->get_total_elements()==(d*d+d));
    this->class_name = "dfield";
    init(d);
    parameters = params;
    inverse_();
};

template <typename type, typename container>
dfield<type,container>::dfield(typename image<type,container>::pointer input)
{
    this->class_name = "dfield";
    init(*input); 
    identity();
};

template <typename type, typename container>
dfield<type,container>::dfield(typename grid<type,container>::pointer input)
{
    this->class_name = "dfield";
    init(*input); 
    identity();
};

template <typename type, typename container>
dfield<type,container>::dfield(const dfield<type,container> & input)
{
    this->class_name = "dfield";
    this->clone_(input);
};

template <typename type, typename container>
void dfield<type,container>::init(int d)
{
    space_object::init(d);
    allocate(d);
};

template <typename type, typename container>
void dfield<type,container>::init(std::vector<int> sz)
{
    // initialize size
    this->set_size(sz);
    // allocated parameters, change the size
    for(int i = 0; i < this->get_dimension(); i++)
    {
        (*parameters)[i] = image<type,container>::new_pointer(sz);        // parameters are 2d
        (*inverse_parameters)[i] = image<type,container>::new_pointer(sz);// parameters are 2d
    };
};

template <typename type, typename container>
void dfield<type,container>::init(const space_object & input)
{
    this->copy_space(input);
    allocate(input.get_dimension());
    init(input.get_size());
};

// ===========================================
// Functions
// ===========================================
template <typename type, typename container>
void dfield<type,container>::change_size(std::vector<int> sz)
{
    space_object::set_size(sz);
};

template <typename type, typename container>
void dfield<type,container>::identity()
{
    // std::cout << "init identity" << std::endl;
    for (int i = 0; i < this->dim; i++)
    {
        // std::cout << "zeros" << std::endl;
        // get_parameters(i)->print();
        get_parameters(i)->zeros();
        // std::cout << "inv" << std::endl;
        get_inverse_parameters(i)->zeros();
    };
    // std::cout << "end identity" << std::endl;
};

template <typename type, typename container>
void dfield<type,container>::inverse_()
{
    *inverse_parameters = (*parameters)*((type)-1.0);
};

//Transform point
template <typename type, typename container>
std::vector<type> dfield<type,container>::apply(std::vector<type> & point)
{
    if (this->dim == 2){ return transform_2d(point); }
    else if (this->dim == 3){ return transform_3d(point); };
    return point;
};

//Transform grid
template <typename type, typename container>
typename grid<type,container>::pointer dfield<type,container>::apply(typename grid<type,container>::pointer input)
{
    if (this->dim == 2){ return transform_2d(input); }
    else if (this->dim == 3){ return transform_3d(input); };
    return input;
};

// ===========================================
// Functions 2d
// ===========================================
// Transform point
template <typename type, typename container>
std::vector<type> dfield<type,container>::transform_2d(std::vector<type> & point)
{
    // TODO: consider if the point uses the inverse or the direct transform. I think direct.
    assert(this->get_dimension() == point.size());
    std::vector<type> out(point.size());
    
    // interpolation required!

    // std::vector<type> va = parameters->get_data()->std_vector();
    // type * a = va.data();
    // out[0] = a[0]*point[0] + a[1]*point[1] + a[4];
    // out[1] = a[2]*point[0] + a[3]*point[1] + a[5];
    return out;
};

//Transform grid
template <typename type, typename container>
typename grid<type,container>::pointer dfield<type,container>::transform_2d(typename grid<type,container>::pointer input)
{
    assert(this->get_dimension() == input->get_dimension());
    auto output = input->mimic();

    typename image<type,container>::pointer * xin = input->ptr();
    typename image<type,container>::pointer * xout = output->ptr();

    container::dfield_2d(xin[0]->get_data(), xin[1]->get_data(),
                         get_parameters(0)->get_data(), get_parameters(1)->get_data(),
                         xout[0]->get_data(), xout[1]->get_data());
    return output;
};

// ===========================================
// Functions 3d
// ===========================================
// Transform point
template <typename type, typename container>
std::vector<type> dfield<type,container>::transform_3d(std::vector<type> & point)
{
    // TODO: consider if the point uses the inverse or the direct transform. I think direct.
    assert(this->get_dimension() == point.size());
    std::vector<type> out(point.size());

    // interpolation required!

    // std::vector<type> va = parameters->get_data()->std_vector();
    // type * a = va.data();
    // out[0] = a[0]*point[0] + a[1]*point[1] + a[2]*point[2] + a[9];
    // out[1] = a[3]*point[0] + a[4]*point[1] + a[5]*point[2] + a[10];
    // out[2] = a[6]*point[0] + a[7]*point[1] + a[8]*point[2] + a[11];
    return out;
};

//Transform grid
template <typename type, typename container>
typename grid<type,container>::pointer dfield<type,container>::transform_3d(typename grid<type,container>::pointer input)
{
    assert(this->get_dimension() == input->get_dimension());
    auto output = input->mimic();
    
    typename image<type,container>::pointer * xin = input->ptr();
    typename image<type,container>::pointer * xout = output->ptr();

    container::dfield_3d(xin[0]->get_data(), xin[1]->get_data(), xin[2]->get_data(),
                         get_parameters(0)->get_data(), get_parameters(1)->get_data(), get_parameters(2)->get_data(),
                         xout[0]->get_data(), xout[1]->get_data(), xout[2]->get_data());
    return output;
};

}; //end namespace

#endif