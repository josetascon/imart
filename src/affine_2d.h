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

public:
    affine_2d();
    affine_2d(image_2d<pixel_type> & params);

    // using transform_base<pixel_type>::print;

    // template <typename pixel_type_>
    // friend std::ostream & operator << (std::ostream & os, affine_2d<pixel_type_> & input);

    // ===========================================
    // Functions
    // ===========================================
    std::vector<pixel_type> transform(std::vector<pixel_type> & point);
    void transform(image_2d<pixel_type> & image);

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
    this->dim = 2;
    this->type = "affine";
    this->parameters = image_base<pixel_type>();
};

template <typename pixel_type>
affine_2d<pixel_type>::affine_2d(image_2d<pixel_type> & params)
{
    this->dim = 2;
    this->type = "affine";
    this->parameters = params;
};


// template <typename pixel_type>
// std::ostream & operator << (std::ostream & os, transform_base<pixel_type> & input)
// {
//     os << input.info("");
//     return os;
// };


// ===========================================
// Functions
// ===========================================
// Transform point
template <typename pixel_type>
std::vector<pixel_type> affine_2d<pixel_type>::transform(std::vector<pixel_type> & point)
{
    std::vector<pixel_type> out(this->dim);
    image_base<pixel_type> param = this->parameters;
    out[0] = param(0)*point[0] + param(1)*point[1] + param(4);
    out[1] = param(2)*point[0] + param(3)*point[1] + param(5);
    return out;
};

// Transform image
template <typename pixel_type>
void affine_2d<pixel_type>::transform(image_2d<pixel_type> & image)
{
    using vector_image = std::vector<image_2d<pixel_type>>;
    std::shared_ptr<vector_image> xy = image.get_grid();

    // CONTINUE HERE
    // image_2d<pixel_type> x = parameters(0)*(*xy)[0] + parameters(1)*(*xy)[0] + parameters(4);
    // image_2d<pixel_type> y = parameters(2)*(*xy)[1] + parameters(3)*(*xy)[1] + parameters(5);
};

#endif