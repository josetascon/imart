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

    using transform<type,container>::init;

protected:
    // ===========================================
    // Functions
    // ===========================================
    // void init(int d);
    void inverse_();
    void inverse_2d();
    void inverse_3d();

    std::vector<type> transform_2d(std::vector<type> & point);
    std::vector<type> transform_3d(std::vector<type> & point);
    template <typename gcontainer>
    typename grid<type,gcontainer>::pointer transform_2d(const typename grid<type,gcontainer>::pointer input);
    template <typename gcontainer>
    typename grid<type,gcontainer>::pointer transform_3d(const typename grid<type,gcontainer>::pointer input);

public:
    // ===========================================
    // Constructor Functions
    // ===========================================
    affine();
    affine(int d);
    affine(int d, typename image<type,container>::pointer params);

    // ===========================================
    // Functions
    // ===========================================
    std::vector<type> apply(std::vector<type> & point);
    
    template <typename gcontainer>
    typename grid<type,gcontainer>::pointer apply(const typename grid<type,gcontainer>::pointer input);
    // void transform(image<type,container> & image); // not going to implement this
};


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
    this->init(2);
};

template <typename type, typename container>
affine<type,container>::affine(int d)
{
    this->class_name = "affine";
    this->init(d);
};

template <typename type, typename container>
affine<type,container>::affine(int d, typename image<type,container>::pointer params)
{
    assert(params->get_total_elements()==(d*d+d));
    this->class_name = "affine";
    this->init(d);
    this->parameters = params;
    inverse_();
};

// template <typename type, typename container>
// void affine<type,container>::init(int d)
// {
//     transform<type,container>::init(d);
// };

// ===========================================
// Functions
// ===========================================
template <typename type, typename container>
void affine<type,container>::inverse_()
{
    if (this->dim == 2){ inverse_2d(); };
    if (this->dim == 3){ inverse_3d(); };
};


template <typename type, typename container>
void affine<type,container>::inverse_2d()
{
    auto inv = image<type,container>::new_pointer(6,1);
    // inv->mimic_(*(this->parameters));

    type * a = this->parameters->ptr();
    type * p = inv->ptr();

    p[0] = a[3]/(a[0]*a[3] - a[1]*a[2]);
    p[1] = -a[1]/(a[0]*a[3] - a[1]*a[2]);
    p[2] = -a[2]/(a[0]*a[3] - a[1]*a[2]);
    p[3] = a[0]/(a[0]*a[3] - a[1]*a[2]);
    p[4] = (a[1]*(a[0]*a[5] - a[2]*a[4]) - a[4]*(a[0]*a[3] - a[1]*a[2]))/(a[0]*(a[0]*a[3] - a[1]*a[2]));
    p[5] = (-a[0]*a[5] + a[2]*a[4])/(a[0]*a[3] - a[1]*a[2]);

    this->inverse_parameters = inv;
};

template <typename type, typename container>
void affine<type,container>::inverse_3d()
{
    auto inv = image<type,container>::new_pointer(12,1);
    // inv->mimic_(*this->parameters);
    type * a = this->parameters->ptr();
    type * p = inv->ptr();

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

    this->inverse_parameters = inv;
};

//Transform point
template <typename type, typename container>
std::vector<type> affine<type,container>::apply(std::vector<type> & point)
{
    if (this->dim == 2){ return transform_2d(point); };
    if (this->dim == 3){ return transform_3d(point); };
    return point;
};

//Transform grid
template <typename type, typename container>
template <typename gcontainer>
typename grid<type,gcontainer>::pointer affine<type,container>::apply(typename grid<type,gcontainer>::pointer input)
{
    if (this->dim == 2){ return transform_2d<gcontainer>(input); };
    if (this->dim == 3){ return transform_3d<gcontainer>(input); };
    return input;
};

// ===========================================
// Functions 2d
// ===========================================
// Transform point
template <typename type, typename container>
std::vector<type> affine<type,container>::transform_2d(std::vector<type> & point)
{
    // TODO: assert*****
    // TODO: consider if the point uses the inverse or the direct transform. I think direct.
    std::vector<type> out(point.size());
    type * a = this->parameters->ptr();

    out[0] = a[0]*point[0] + a[1]*point[1] + a[4];
    out[1] = a[2]*point[0] + a[3]*point[1] + a[5];

    return out;
};

//Transform grid
template <typename type, typename container>
template <typename gcontainer>
typename grid<type,gcontainer>::pointer affine<type,container>::transform_2d(typename grid<type,gcontainer>::pointer input)
{
    // TODO: assert*****
    auto output = grid<type,gcontainer>::new_pointer(input->get_dimension());
    output->mimic_(*input);

    type * a = this->parameters->ptr();
    typename image<type,gcontainer>::pointer * xin = input->ptr();
    typename image<type,gcontainer>::pointer * xout = output->ptr();

    *(xout[0]) = (*xin[0])*a[0] + (*xin[1])*a[1] + a[4];
    *(xout[1]) = (*xin[0])*a[2] + (*xin[1])*a[3] + a[5];
    return output;
};

// ===========================================
// Functions 3d
// ===========================================
// Transform point
template <typename type, typename container>
std::vector<type> affine<type,container>::transform_3d(std::vector<type> & point)
{
    // TODO: assert*****
    // TODO: consider if the point uses the inverse or the direct transform. I think direct.
    std::vector<type> out(point.size());
    type * a = this->parameters->ptr();
    out[0] = a[0]*point[0] + a[1]*point[1] + a[2]*point[2] + a[9];
    out[1] = a[3]*point[0] + a[4]*point[1] + a[5]*point[2] + a[10];
    out[2] = a[6]*point[0] + a[7]*point[1] + a[8]*point[2] + a[11];
    return out;
};

//Transform grid
template <typename type, typename container>
template <typename gcontainer>
typename grid<type,gcontainer>::pointer affine<type,container>::transform_3d(typename grid<type,gcontainer>::pointer input)
{
    // TODO: assert*****
    // typename grid<type,gcontainer>::pointer output;
    auto output = grid<type,gcontainer>::new_pointer(input->get_dimension());
    output->mimic_(*input);
    type * a = this->parameters->ptr();
    typename image<type,gcontainer>::pointer * xin = input->ptr();
    typename image<type,gcontainer>::pointer * xout = output->ptr();
    (*xout[0]) = (*xin[0])*a[0] + (*xin[1])*a[1] + (*xin[2])*a[2] + a[9];
    (*xout[1]) = (*xin[0])*a[3] + (*xin[1])*a[4] + (*xin[2])*a[5] + a[10];
    (*xout[2]) = (*xin[0])*a[6] + (*xin[1])*a[7] + (*xin[2])*a[8] + a[11];
    return output;
};

}; //end namespace

#endif