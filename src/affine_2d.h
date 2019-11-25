/*
* @Author: jose
* @Date:   2019-11-15 13:55:00
* @Last Modified by:   jose
* @Last Modified time: 2019-11-05 13:55:00
*/

#ifndef __AFFINE_2D_H__
#define __AFFINE_2D_H__

// images 
#include "transform_base.h"
#include "image_2d.h"

// parallel
// openmp
// opencl


// Definitions
// typedef float pixel_type;

// Class image_base_2d
template <typename pixel_type>
class affine_2d: public transform_base<pixel_type>
{
protected:
    // image_2d<pixel_type> parameters;
    void init(int d);

public:
    affine_2d();
    affine_2d(image_2d<pixel_type> & params);

    // using transform_base<pixel_type>::print;

    // template <typename pixel_type_>
    // friend std::ostream & operator << (std::ostream & os, affine_2d<pixel_type_> & input);

    // ===========================================
    // Functions
    // ===========================================
    void inverse();

    std::vector<pixel_type> transform(std::vector<pixel_type> & point);
    grid<pixel_type> transform(grid<pixel_type> & input);
    // void transform(image_2d<pixel_type> & image); // not going to implement this

};


// ===========================================
//          Functions of Class affine_2d
// ===========================================


// ===========================================
// Create Functions
// ===========================================
template <typename pixel_type>
affine_2d<pixel_type>::affine_2d()
{
    init(2);
    this->class_name = "affine 2d";
    this->parameters = image_base<pixel_type>();
};

template <typename pixel_type>
affine_2d<pixel_type>::affine_2d(image_2d<pixel_type> & params)
{
    init(2);
    this->class_name = "affine 2d";
    this->parameters = params;
    this->inverse();
};

template <typename pixel_type>
void affine_2d<pixel_type>::init(int d)
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
void affine_2d<pixel_type>::inverse()
{
    // image_base<pixel_type> inv(this->parameters.get_width(),this->parameters.get_height());
    // int w = this->parameters.get_width();
    // int h = this->parameters.get_height();
    // this->inverse_parameters = image_base<pixel_type>(w,h);
    image_base<pixel_type> inv(this->parameters);
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


// ===========================================
// Functions
// ===========================================
// Transform point
template <typename pixel_type>
std::vector<pixel_type> affine_2d<pixel_type>::transform(std::vector<pixel_type> & point)
{
    // TODO: assert*****
    // TODO: consider if the point uses 
    // the inverse or the direct transform. I think direct.
    std::vector<pixel_type> out(this->dim);
    pixel_type * a = this->inverse_parameters.ptr();
    out[0] = a[0]*point[0] + a[1]*point[1] + a[4];
    out[1] = a[2]*point[0] + a[3]*point[1] + a[5];
    return out;
};

//Transform grid
template <typename pixel_type>
grid<pixel_type> affine_2d<pixel_type>::transform(grid<pixel_type> & input)
{
    // TODO: assert*****
    grid<pixel_type> output(input);
    pixel_type * a = this->inverse_parameters.ptr();

    image_base<pixel_type> * xin = input.ptr();
    image_base<pixel_type> * xout = output.ptr();
    
    xout[0] = a[0]*xin[0] + a[1]*xin[1] + a[4];
    xout[1] = a[2]*xin[0] + a[3]*xin[1] + a[5];
    return output;
};


// This is NOT going to be implemented
// Transform image. 
// template <typename pixel_type>
// void affine_2d<pixel_type>::transform(image_2d<pixel_type> & image)
// {
//     using vector_image = std::vector<image_2d<pixel_type>>;
//     std::shared_ptr<vector_image> xy = image.get_grid();

//     // CONTINUE HERE
//     // image_2d<pixel_type> x = parameters(0)*(*xy)[0] + parameters(1)*(*xy)[0] + parameters(4);
//     // image_2d<pixel_type> y = parameters(2)*(*xy)[1] + parameters(3)*(*xy)[1] + parameters(5);
// };

#endif