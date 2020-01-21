/*
* @Author: jose
* @Date:   2019-11-26 00:00:00
* @Last Modified by:   jose
* @Last Modified time: 2019-11-26 00:00:00
*/

#ifndef __AFFINE_3D_H__
#define __AFFINE_3D_H__

// images 
#include "transform_base.h"
#include "image_3d.h"

// Class image_base_2d
template <typename pixel_type>
class affine_3d: public transform_base<pixel_type>
{
protected:
    void init(int d);

public:
    affine_3d();
    affine_3d(image_2d<pixel_type> & params);

    // ===========================================
    // Functions
    // ===========================================
    void inverse();

    std::vector<pixel_type> transform(std::vector<pixel_type> & point);
    grid<pixel_type> transform(grid<pixel_type> & input);
    // void transform(image_3d<pixel_type> & image); // not going to implement this

};


// ===========================================
//          Functions of Class affine_2d
// ===========================================


// ===========================================
// Create Functions
// ===========================================
template <typename pixel_type>
affine_3d<pixel_type>::affine_3d()
{
    init(3);
    this->class_name = "affine 3d";
    this->parameters = image_base<pixel_type>();
};

template <typename pixel_type>
affine_3d<pixel_type>::affine_3d(image_2d<pixel_type> & params)
{
    init(3);
    this->class_name = "affine 3d";
    this->parameters = params;
    this->inverse();
};

template <typename pixel_type>
void affine_3d<pixel_type>::init(int d)
{
    transform_base<pixel_type>::init(d);
};


// template <typename pixel_type>
// std::ostream & operator << (std::ostream & os, transform_base<pixel_type> & input)
// {
//     os << input.info("");
//     return os;
// };


template <typename pixel_type>
void affine_3d<pixel_type>::inverse()
{
    
    image_base<pixel_type> inv(this->parameters);
    pixel_type * a = this->parameters.ptr();
    pixel_type * p = inv.ptr();

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


// ===========================================
// Functions
// ===========================================
// Transform point
template <typename pixel_type>
std::vector<pixel_type> affine_3d<pixel_type>::transform(std::vector<pixel_type> & point)
{
    // TODO: assert*****
    // TODO: consider if the point uses 
    // the inverse or the direct transform. I think direct.
    std::vector<pixel_type> out(this->dim);
    pixel_type * a = this->parameters.ptr();
    
    out[0] = a[0]*point[0] + a[1]*point[1] + a[2]*point[2] + a[9];
    out[1] = a[3]*point[0] + a[4]*point[1] + a[5]*point[2] + a[10];
    out[2] = a[6]*point[0] + a[7]*point[1] + a[8]*point[2] + a[11];
    return out;
};

//Transform grid
template <typename pixel_type>
grid<pixel_type> affine_3d<pixel_type>::transform(grid<pixel_type> & input)
{
    // TODO: assert*****
    grid<pixel_type> output(input);
    pixel_type * a = this->parameters.ptr();

    image_base<pixel_type> * xin = input.ptr();
    image_base<pixel_type> * xout = output.ptr();
    
    xout[0] = a[0]*xin[0] + a[1]*xin[1] + a[2]*xin[2] + a[9];
    xout[1] = a[3]*xin[0] + a[4]*xin[1] + a[5]*xin[2] + a[10];
    xout[2] = a[6]*xin[0] + a[7]*xin[1] + a[8]*xin[2] + a[11];

    return output;
};


#endif