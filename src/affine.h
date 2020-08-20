/*
* @Author: jose
* @Date:   2019-11-15 13:55:00
* @Last Modified by:   jose
* @Last Modified time: 2019-11-05 13:55:00
*/

#ifndef __AFFINE_H__
#define __AFFINE_H__

// local libs
#include "transform.h"

namespace imart
{

// Class affine
template <typename type, typename container=vector_cpu<type>>
class affine: public inherit<affine<type,container>, transform<type,container>>
{
public:
    //Type definitions
    using self    = affine;
    using pointer = std::shared_ptr<self>;
    using vector  = std::vector<self::pointer>;

    // Inherited variables
    using transform<type,container>::parameters;
    using transform<type,container>::inverse_parameters;
    using transform<type,container>::get_parameters;
    using transform<type,container>::get_inverse_parameters;

    using transform<type,container>::operator=;
    // using transform<type,container>::operator+;

    using inherit<affine<type,container>, transform<type,container>>::inherit;

protected:
    // ===========================================
    // Functions
    // ===========================================
    void init(int d);
    void inverse_();
    void inverse_2d();
    void inverse_3d();

    std::vector<type> transform_2d(std::vector<type> & point);
    std::vector<type> transform_3d(std::vector<type> & point);
    
    typename grid<type,container>::pointer transform_2d(const typename grid<type,container>::pointer input);
    typename grid<type,container>::pointer transform_3d(const typename grid<type,container>::pointer input);

public:
    // ===========================================
    // Constructor Functions
    // ===========================================
    affine();
    affine(int d);
    affine(int d, typename image<type,container>::pointer params);

    // ===========================================
    // Overloading Functions
    // ===========================================
    // affine<type,container> & operator = (const affine<type,container> & input);
    // affine<type,container> operator + (const affine<type,container> & input);
    // affine<type,container> operator - (const affine<type,container> & input);
    // affine<type,container> operator * (const image<type,container> & input);
    // affine<type,container> operator * (type scalar);

    // ===========================================
    // Functions
    // ===========================================
    void identity();
    std::vector<type> apply(std::vector<type> & point);
    typename grid<type,container>::pointer apply(const typename grid<type,container>::pointer input);
    // void transform(image<type,container> & image); // not going to implement this
};

template<typename type>
using affine_cpu = affine<type,vector_cpu<type>>;

template<typename type>
using affine_gpu = affine<type,vector_ocl<type>>;


// ===========================================
//          Functions of Class affine
// ===========================================

// ===========================================
// Constructor Functions
// ===========================================
template <typename type, typename container>
affine<type,container>::affine()
{
    this->class_name = "affine";
    init(2);
    identity();
};

template <typename type, typename container>
affine<type,container>::affine(int d)
{
    this->class_name = "affine";
    init(d);
    identity();
};

template <typename type, typename container>
affine<type,container>::affine(int d, typename image<type,container>::pointer params)
{
    assert(params->get_total_elements()==(d*d+d));
    this->class_name = "affine";
    init(d);
    (*parameters)[0] = params;
    inverse_();
};

template <typename type, typename container>
void affine<type,container>::init(int d)
{
    space_object::init(d);
    int n = d*d + d;
    
    parameters = std::make_shared< std::vector< typename image<type,container>::pointer >>(1);
    inverse_parameters = std::make_shared< std::vector< typename image<type,container>::pointer >>(1);
    
    (*parameters)[0] = image<type,container>::new_pointer(n,1);          // parameters are 2d
    (*inverse_parameters)[0] = image<type,container>::new_pointer(n,1);  // parameters are 2d
};

// ===========================================
// Functions
// ===========================================
template <typename type, typename container>
void affine<type,container>::identity()
{
    if (this->dim == 2)
    { 
        std::vector<type> p(6);
        p[0] = 1; p[1] = 0; 
        p[2] = 0; p[3] = 1;
        p[4] = 0; p[5] = 0;
        get_parameters()->get_data()->read_ram(p.data(),p.size());
        get_inverse_parameters()->get_data()->read_ram(p.data(),p.size());
    }
    else if (this->dim == 3)
    {
        std::vector<type> p(12);
        p[0] = 1; p[1]  = 0; p[2]  = 0; 
        p[3] = 0; p[4]  = 1; p[5]  = 0;
        p[6] = 0; p[7]  = 0; p[8]  = 1;
        p[9] = 0; p[10] = 0; p[11] = 0;
        get_parameters()->get_data()->read_ram(p.data(),p.size());
        get_inverse_parameters()->get_data()->read_ram(p.data(),p.size());
    };
};

template <typename type, typename container>
void affine<type,container>::inverse_()
{
    if (this->dim == 2){ inverse_2d(); }
    else if (this->dim == 3){ inverse_3d(); };
};


template <typename type, typename container>
void affine<type,container>::inverse_2d()
{
    std::vector<type> vp(6);
    type * p = vp.data();

    std::vector<type> va = get_parameters()->get_data()->std_vector();
    type * a = va.data();

    p[0] = a[3]/(a[0]*a[3] - a[1]*a[2]);
    p[1] = -a[1]/(a[0]*a[3] - a[1]*a[2]);
    p[2] = -a[2]/(a[0]*a[3] - a[1]*a[2]);
    p[3] = a[0]/(a[0]*a[3] - a[1]*a[2]);
    p[4] = (a[1]*(a[0]*a[5] - a[2]*a[4]) - a[4]*(a[0]*a[3] - a[1]*a[2]))/(a[0]*(a[0]*a[3] - a[1]*a[2]));
    p[5] = (-a[0]*a[5] + a[2]*a[4])/(a[0]*a[3] - a[1]*a[2]);

    auto inv = image<type,container>::new_pointer(6,1);
    auto vec = container::new_pointer(6);
    vec->read_ram(vp.data(),6);
    inv->set_data(vec);
    (*inverse_parameters)[0] = inv;
};

template <typename type, typename container>
void affine<type,container>::inverse_3d()
{
    std::vector<type> vp(12);
    type * p = vp.data();

    std::vector<type> va = get_parameters()->get_data()->std_vector();
    type * a = va.data();
    
    p[0]  = (a[0]*a[4]*((a[0]*a[4] - a[1]*a[3])*(a[0]*a[8] - a[2]*a[6]) - (a[0]*a[5] - a[2]*a[3])*(a[0]*a[7] - a[1]*a[6])) - (-a[1]*(a[0]*a[5] - a[2]*a[3]) + a[2]*(a[0]*a[4] - a[1]*a[3]))*(a[3]*(a[0]*a[7] - a[1]*a[6]) - a[6]*(a[0]*a[4] - a[1]*a[3])))/(a[0]*(a[0]*a[4] - a[1]*a[3])*((a[0]*a[4] - a[1]*a[3])*(a[0]*a[8] - a[2]*a[6]) - (a[0]*a[5] - a[2]*a[3])*(a[0]*a[7] - a[1]*a[6])));
    p[1]  = (-a[0]*a[1]*((a[0]*a[4] - a[1]*a[3])*(a[0]*a[8] - a[2]*a[6]) - (a[0]*a[5] - a[2]*a[3])*(a[0]*a[7] - a[1]*a[6])) + a[0]*(a[0]*a[7] - a[1]*a[6])*(-a[1]*(a[0]*a[5] - a[2]*a[3]) + a[2]*(a[0]*a[4] - a[1]*a[3])))/(a[0]*(a[0]*a[4] - a[1]*a[3])*((a[0]*a[4] - a[1]*a[3])*(a[0]*a[8] - a[2]*a[6]) - (a[0]*a[5] - a[2]*a[3])*(a[0]*a[7] - a[1]*a[6])));
    p[2]  = -(-a[1]*(a[0]*a[5] - a[2]*a[3]) + a[2]*(a[0]*a[4] - a[1]*a[3]))/((a[0]*a[4] - a[1]*a[3])*(a[0]*a[8] - a[2]*a[6]) - (a[0]*a[5] - a[2]*a[3])*(a[0]*a[7] - a[1]*a[6]));
    p[3]  = (-a[3]*((a[0]*a[4] - a[1]*a[3])*(a[0]*a[8] - a[2]*a[6]) - (a[0]*a[5] - a[2]*a[3])*(a[0]*a[7] - a[1]*a[6])) - (a[0]*a[5] - a[2]*a[3])*(a[3]*(a[0]*a[7] - a[1]*a[6]) - a[6]*(a[0]*a[4] - a[1]*a[3])))/((a[0]*a[4] - a[1]*a[3])*((a[0]*a[4] - a[1]*a[3])*(a[0]*a[8] - a[2]*a[6]) - (a[0]*a[5] - a[2]*a[3])*(a[0]*a[7] - a[1]*a[6])));
    p[4]  = (a[0]*(a[0]*a[5] - a[2]*a[3])*(a[0]*a[7] - a[1]*a[6]) + a[0]*((a[0]*a[4] - a[1]*a[3])*(a[0]*a[8] - a[2]*a[6]) - (a[0]*a[5] - a[2]*a[3])*(a[0]*a[7] - a[1]*a[6])))/((a[0]*a[4] - a[1]*a[3])*((a[0]*a[4] - a[1]*a[3])*(a[0]*a[8] - a[2]*a[6]) - (a[0]*a[5] - a[2]*a[3])*(a[0]*a[7] - a[1]*a[6])));
    p[5]  = -a[0]*(a[0]*a[5] - a[2]*a[3])/((a[0]*a[4] - a[1]*a[3])*(a[0]*a[8] - a[2]*a[6]) - (a[0]*a[5] - a[2]*a[3])*(a[0]*a[7] - a[1]*a[6]));
    p[6]  = (a[3]*(a[0]*a[7] - a[1]*a[6]) - a[6]*(a[0]*a[4] - a[1]*a[3]))/((a[0]*a[4] - a[1]*a[3])*(a[0]*a[8] - a[2]*a[6]) - (a[0]*a[5] - a[2]*a[3])*(a[0]*a[7] - a[1]*a[6]));
    p[7]  = -a[0]*(a[0]*a[7] - a[1]*a[6])/((a[0]*a[4] - a[1]*a[3])*(a[0]*a[8] - a[2]*a[6]) - (a[0]*a[5] - a[2]*a[3])*(a[0]*a[7] - a[1]*a[6]));
    p[8]  = a[0]*(a[0]*a[4] - a[1]*a[3])/((a[0]*a[4] - a[1]*a[3])*(a[0]*a[8] - a[2]*a[6]) - (a[0]*a[5] - a[2]*a[3])*(a[0]*a[7] - a[1]*a[6]));
    p[9]  = (-(-a[1]*(a[0]*a[10] - a[3]*a[9]) + a[9]*(a[0]*a[4] - a[1]*a[3]))*((a[0]*a[4] - a[1]*a[3])*(a[0]*a[8] - a[2]*a[6]) - (a[0]*a[5] - a[2]*a[3])*(a[0]*a[7] - a[1]*a[6])) + (-a[1]*(a[0]*a[5] - a[2]*a[3]) + a[2]*(a[0]*a[4] - a[1]*a[3]))*(-(a[0]*a[10] - a[3]*a[9])*(a[0]*a[7] - a[1]*a[6]) + (a[0]*a[11] - a[6]*a[9])*(a[0]*a[4] - a[1]*a[3])))/(a[0]*(a[0]*a[4] - a[1]*a[3])*((a[0]*a[4] - a[1]*a[3])*(a[0]*a[8] - a[2]*a[6]) - (a[0]*a[5] - a[2]*a[3])*(a[0]*a[7] - a[1]*a[6])));
    p[10] = (-(a[0]*a[10] - a[3]*a[9])*((a[0]*a[4] - a[1]*a[3])*(a[0]*a[8] - a[2]*a[6]) - (a[0]*a[5] - a[2]*a[3])*(a[0]*a[7] - a[1]*a[6])) + (a[0]*a[5] - a[2]*a[3])*(-(a[0]*a[10] - a[3]*a[9])*(a[0]*a[7] - a[1]*a[6]) + (a[0]*a[11] - a[6]*a[9])*(a[0]*a[4] - a[1]*a[3])))/((a[0]*a[4] - a[1]*a[3])*((a[0]*a[4] - a[1]*a[3])*(a[0]*a[8] - a[2]*a[6]) - (a[0]*a[5] - a[2]*a[3])*(a[0]*a[7] - a[1]*a[6])));
    p[11] = ((a[0]*a[10] - a[3]*a[9])*(a[0]*a[7] - a[1]*a[6]) - (a[0]*a[11] - a[6]*a[9])*(a[0]*a[4] - a[1]*a[3]))/((a[0]*a[4] - a[1]*a[3])*(a[0]*a[8] - a[2]*a[6]) - (a[0]*a[5] - a[2]*a[3])*(a[0]*a[7] - a[1]*a[6]));

    auto inv = image<type,container>::new_pointer(12,1);
    auto vec = container::new_pointer(12);
    vec->read_ram(vp.data(),12);
    inv->set_data(vec);
    (*inverse_parameters)[0] = inv;
};

//Transform point
template <typename type, typename container>
std::vector<type> affine<type,container>::apply(std::vector<type> & point)
{
    if (this->dim == 2){ return transform_2d(point); }
    else if (this->dim == 3){ return transform_3d(point); };
    return point;
};

//Transform grid
template <typename type, typename container>
typename grid<type,container>::pointer affine<type,container>::apply(typename grid<type,container>::pointer input)
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
std::vector<type> affine<type,container>::transform_2d(std::vector<type> & point)
{
    // TODO: consider if the point uses the inverse or the direct transform. I think direct.
    assert(this->get_dimension() == point.size());
    std::vector<type> out(point.size());
    std::vector<type> va = get_parameters()->get_data()->std_vector();
    type * a = va.data();

    out[0] = a[0]*point[0] + a[1]*point[1] + a[4];
    out[1] = a[2]*point[0] + a[3]*point[1] + a[5];

    return out;
};

//Transform grid
template <typename type, typename container>
typename grid<type,container>::pointer affine<type,container>::transform_2d(typename grid<type,container>::pointer input)
{
    assert(this->get_dimension() == input->get_dimension());
    auto output = input->mimic();
    // std::vector<type> va = parameters->get_data()->std_vector();
    // type * a = va.data();
    typename image<type,container>::pointer * xin = input->ptr();
    typename image<type,container>::pointer * xout = output->ptr();

    container::affine_2d(xin[0]->get_data(), xin[1]->get_data(),
                         xout[0]->get_data(), xout[1]->get_data(),
                         get_parameters()->get_data());
    // *(xout[0]) = (*xin[0])*a[0] + (*xin[1])*a[1] + a[4];
    // *(xout[1]) = (*xin[0])*a[2] + (*xin[1])*a[3] + a[5];
    return output;
};

// ===========================================
// Functions 3d
// ===========================================
// Transform point
template <typename type, typename container>
std::vector<type> affine<type,container>::transform_3d(std::vector<type> & point)
{
    // TODO: consider if the point uses the inverse or the direct transform. I think direct.
    assert(this->get_dimension() == point.size());
    std::vector<type> out(point.size());
    std::vector<type> va = get_parameters()->get_data()->std_vector();
    type * a = va.data();
    out[0] = a[0]*point[0] + a[1]*point[1] + a[2]*point[2] + a[9];
    out[1] = a[3]*point[0] + a[4]*point[1] + a[5]*point[2] + a[10];
    out[2] = a[6]*point[0] + a[7]*point[1] + a[8]*point[2] + a[11];
    return out;
};

//Transform grid
template <typename type, typename container>
typename grid<type,container>::pointer affine<type,container>::transform_3d(typename grid<type,container>::pointer input)
{
    assert(this->get_dimension() == input->get_dimension());
    auto output = input->mimic();
    // std::vector<type> va = parameters->get_data()->std_vector();
    // type * a = va.data();
    typename image<type,container>::pointer * xin = input->ptr();
    typename image<type,container>::pointer * xout = output->ptr();

    container::affine_3d(xin[0]->get_data(), xin[1]->get_data(), xin[2]->get_data(),
                         xout[0]->get_data(), xout[1]->get_data(), xout[2]->get_data(),
                         get_parameters()->get_data());
    // (*xout[0]) = (*xin[0])*a[0] + (*xin[1])*a[1] + (*xin[2])*a[2] + a[9];
    // (*xout[1]) = (*xin[0])*a[3] + (*xin[1])*a[4] + (*xin[2])*a[5] + a[10];
    // (*xout[2]) = (*xin[0])*a[6] + (*xin[1])*a[7] + (*xin[2])*a[8] + a[11];
    return output;
};

// ===========================================
// Overloading Functions
// ===========================================
// Equal
// template <typename type, typename container>
// affine<type,container> & affine<type,container>::operator = (const affine<type,container> & input)
// {
//     // delete &data;
//     this->copy_(input);
//     return *this;
// };

// // Transform to Transform
// template <typename type, typename container>
// affine<type,container> affine<type,container>::operator + (const affine<type,container> & input)
// {
//     auto pp = image<type,container>::new_pointer();
//     *pp = *parameters + *(input.get_parameters());

//     affine<type,container> output;
//     output.mimic_(*this);
//     output.set_parameters( pp );
//     return output;
//     // pointer output = this->mimic();
//     // output->set_parameters( pp );
//     // return *output;
// };

// template <typename type, typename container>
// affine<type,container> affine<type,container>::operator - (const affine<type,container> & input)
// {
//     auto pp = image<type,container>::new_pointer();
//     *pp = *parameters - *(input.get_parameters());

//     affine<type,container> output;
//     output.mimic_(*this);
//     output.set_parameters( pp );
//     return output;
// };

// template <typename type, typename container>
// affine<type,container> affine<type,container>::operator * (const image<type,container> & input)
// {
//     auto pp = image<type,container>::new_pointer();
//     *pp = (*parameters)*input;

//     affine<type,container> output;
//     output.mimic_(*this);
//     output.set_parameters( pp );
//     return output;
// };

// // Scalar
// template <typename type, typename container>
// affine<type,container> affine<type,container>::operator * (type scalar)
// {
//     auto pp = image<type,container>::new_pointer();
//     *pp = scalar*(*parameters);

//     affine<type,container> output;
//     output.mimic_(*this);
//     output.set_parameters( pp );
//     return output;
// };

}; //end namespace

#endif