/*
* @Author: jose
* @Date:   2019-11-15 13:55:00
* @Last Modified by:   jose
* @Last Modified time: 2019-11-05 13:55:00
*/

#ifndef __affine_H__
#define __affine_H__

// local libs
#include "image.h"
#include "transform_base.h"

// Class affine
template <typename pixel_type>
class affine: public transform_base<pixel_type>
{
protected:
    // image<pixel_type> parameters;
    void init(int d);

    void inverse_();
    void inverse_2d();
    void inverse_3d();

    std::vector<pixel_type> transform_2d(std::vector<pixel_type> & point);
    std::vector<pixel_type> transform_3d(std::vector<pixel_type> & point);
    grid<pixel_type> transform_2d(grid<pixel_type> & input);
    grid<pixel_type> transform_3d(grid<pixel_type> & input);

public:
    affine();
    affine(int d);
    affine(int d, image<pixel_type> & params);

    // using transform_base<pixel_type>::print;

    // template <typename pixel_type_>
    // friend std::ostream & operator << (std::ostream & os, affine<pixel_type_> & input);

    // ===========================================
    // Functions
    // ===========================================
    // transform_base<pixel_type> inverse();

    std::vector<pixel_type> transform(std::vector<pixel_type> & point);
    grid<pixel_type> transform(grid<pixel_type> & input);
    // void transform(image<pixel_type> & image); // not going to implement this

};


// ===========================================
//          Functions of Class affine
// ===========================================


// ===========================================
// Create Functions
// ===========================================
template <typename pixel_type>
affine<pixel_type>::affine()
{
    this->class_name = "affine";
    this->parameters = image<pixel_type>();
    init(2);  
};

template <typename pixel_type>
affine<pixel_type>::affine(int d)
{
    this->class_name = "affine";
    this->parameters = image<pixel_type>();
    init(d);  
};

template <typename pixel_type>
affine<pixel_type>::affine(int d, image<pixel_type> & params)
{
    this->class_name = "affine";
    this->parameters = params;
    init(d);
};

template <typename pixel_type>
void affine<pixel_type>::init(int d)
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
void affine<pixel_type>::inverse_()
{
    if (this->dim == 2){ inverse_2d(); };
    if (this->dim == 3){ inverse_3d(); };
};


template <typename pixel_type>
void affine<pixel_type>::inverse_2d()
{
    // image<pixel_type> inv(this->parameters.get_width(),this->parameters.get_height());
    // int w = this->parameters.get_width();
    // int h = this->parameters.get_height();
    // this->inverse_parameters = image<pixel_type>(w,h);
    image<pixel_type> inv(this->parameters);
    pixel_type * a = this->parameters.ptr();
    pixel_type * p = inv.ptr();
    // pixel_type * p = this->inverse_parameters.ptr();
    p[0] = a[3]/(a[0]*a[3] - a[1]*a[2]);
    p[1] = -a[1]/(a[0]*a[3] - a[1]*a[2]);
    p[2] = -a[2]/(a[0]*a[3] - a[1]*a[2]);
    p[3] = a[0]/(a[0]*a[3] - a[1]*a[2]);
    p[4] = (a[1]*(a[0]*a[5] - a[2]*a[4]) - a[4]*(a[0]*a[3] - a[1]*a[2]))/(a[0]*(a[0]*a[3] - a[1]*a[2]));
    p[5] = (-a[0]*a[5] + a[2]*a[4])/(a[0]*a[3] - a[1]*a[2]);

    // inv.print_data();
    this->inverse_parameters = inv;
    // this->inverse_parameters.print_data();
};

template <typename pixel_type>
void affine<pixel_type>::inverse_3d()
{
    
    image<pixel_type> inv(this->parameters);
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

//Transform point
template <typename pixel_type>
std::vector<pixel_type> affine<pixel_type>::transform(std::vector<pixel_type> & point)
{
    std::vector<pixel_type> out(this->dim);
    if (this->dim == 2){ out = transform_2d(point); };
    if (this->dim == 3){ out = transform_3d(point); };
    return out;
};

//Transform grid
template <typename pixel_type>
grid<pixel_type> affine<pixel_type>::transform(grid<pixel_type> & input)
{
    
    grid<pixel_type> output(input);
    if (this->dim == 2){ output = transform_2d(input); };
    if (this->dim == 3){ output = transform_3d(input); };
    return output;
};

// ===========================================
// Functions 2d
// ===========================================
// Transform point
template <typename pixel_type>
std::vector<pixel_type> affine<pixel_type>::transform_2d(std::vector<pixel_type> & point)
{
    // TODO: assert*****
    // TODO: consider if the point uses 
    // the inverse or the direct transform. I think direct.
    std::vector<pixel_type> out(this->dim);
    // pixel_type * a = this->inverse_parameters.ptr();
    pixel_type * a = this->parameters.ptr();
    out[0] = a[0]*point[0] + a[1]*point[1] + a[4];
    out[1] = a[2]*point[0] + a[3]*point[1] + a[5];
    return out;
};

//Transform grid
template <typename pixel_type>
grid<pixel_type> affine<pixel_type>::transform_2d(grid<pixel_type> & input)
{
    // TODO: assert*****
    grid<pixel_type> output(input);
    // pixel_type * a = this->inverse_parameters.ptr();
    pixel_type * a = this->parameters.ptr();

    image<pixel_type> * xin = input.ptr();
    image<pixel_type> * xout = output.ptr();
    
    xout[0] = a[0]*xin[0] + a[1]*xin[1] + a[4];
    xout[1] = a[2]*xin[0] + a[3]*xin[1] + a[5];
    return output;
};

// ===========================================
// Functions 3d
// ===========================================
// Transform point
template <typename pixel_type>
std::vector<pixel_type> affine<pixel_type>::transform_3d(std::vector<pixel_type> & point)
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
grid<pixel_type> affine<pixel_type>::transform_3d(grid<pixel_type> & input)
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