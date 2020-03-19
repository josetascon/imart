/*
* @Author: jose
* @Date:   2019-11-15 13:55:00
* @Last Modified by:   jose
* @Last Modified time: 2019-11-05 13:55:00
*/

#ifndef __AFFINE_H__
#define __AFFINE_H__

// local libs
#include "image.h"
#include "transform_base.h"

// Class affine
template <typename pixel_type>
class affine: public transform_base<pixel_type>
{
public:
    //Type definitions
    using self    = affine;
    using pointer = std::shared_ptr<self>;
    using vector  = std::vector<self::pointer>;

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
    affine(int d, typename image<pixel_type>::pointer params);

    template<typename... ARGS>
    static pointer new_pointer(const ARGS&... args);

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
    init(2);
};

template <typename pixel_type>
affine<pixel_type>::affine(int d)
{
    this->class_name = "affine";
    init(d);
};

template <typename pixel_type>
affine<pixel_type>::affine(int d, typename image<pixel_type>::pointer params)
{
    assert(params->get_total_elements()==(d*d+d));
    this->class_name = "affine";
    init(d);
    this->parameters = params;
    inverse_();
};

template <typename pixel_type>
void affine<pixel_type>::init(int d)
{
    transform_base<pixel_type>::init(d);
};

template <typename pixel_type>
template <typename ... ARGS>
typename affine<pixel_type>::pointer affine<pixel_type>::new_pointer(const ARGS&... args)
{
    return std::make_shared< affine<pixel_type> >(args...); // not working for inherited classes
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
    typename image<pixel_type>::pointer inv(new image<pixel_type>(6,1));
    // inv->imitate(*(this->parameters));

    pixel_type * a = this->parameters->ptr();
    pixel_type * p = inv->ptr();

    p[0] = a[3]/(a[0]*a[3] - a[1]*a[2]);
    p[1] = -a[1]/(a[0]*a[3] - a[1]*a[2]);
    p[2] = -a[2]/(a[0]*a[3] - a[1]*a[2]);
    p[3] = a[0]/(a[0]*a[3] - a[1]*a[2]);
    p[4] = (a[1]*(a[0]*a[5] - a[2]*a[4]) - a[4]*(a[0]*a[3] - a[1]*a[2]))/(a[0]*(a[0]*a[3] - a[1]*a[2]));
    p[5] = (-a[0]*a[5] + a[2]*a[4])/(a[0]*a[3] - a[1]*a[2]);

    this->inverse_parameters = inv;
};

template <typename pixel_type>
void affine<pixel_type>::inverse_3d()
{
    typename image<pixel_type>::pointer inv;
    inv->imitate(*this->parameters);
    pixel_type * a = this->parameters->ptr();
    pixel_type * p = inv->ptr();

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
    std::vector<pixel_type> out(point.size());
    if (this->dim == 2){ out = transform_2d(point); };
    if (this->dim == 3){ out = transform_3d(point); };
    return out;
};

//Transform grid
template <typename pixel_type>
grid<pixel_type> affine<pixel_type>::transform(grid<pixel_type> & input)
{
    grid<pixel_type> output(input.get_dimension());
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
    std::vector<pixel_type> out(point.size());
    pixel_type * a = this->parameters->ptr();

    out[0] = a[0]*point[0] + a[1]*point[1] + a[4];
    out[1] = a[2]*point[0] + a[3]*point[1] + a[5];

    return out;
};

//Transform grid
template <typename pixel_type>
grid<pixel_type> affine<pixel_type>::transform_2d(grid<pixel_type> & input)
{
    // TODO: assert*****
    grid<pixel_type> output;
    output.imitate(input);
    pixel_type * a = this->parameters->ptr();

    image<pixel_type> * xin = input.ptr();
    image<pixel_type> * xout = output.ptr();
    
    xout[0] = xin[0]*a[0] + a[1]*xin[1] + a[4];
    xout[1] = xin[0]*a[2] + a[3]*xin[1] + a[5];   

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
    std::vector<pixel_type> out(point.size());
    pixel_type * a = this->parameters->ptr();
    
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
    grid<pixel_type> output;
    output.imitate(input);
    pixel_type * a = this->parameters->ptr();

    image_base<pixel_type> * xin = input.ptr();
    image_base<pixel_type> * xout = output.ptr();
    
    xout[0] = a[0]*xin[0] + a[1]*xin[1] + a[2]*xin[2] + a[9];
    xout[1] = a[3]*xin[0] + a[4]*xin[1] + a[5]*xin[2] + a[10];
    xout[2] = a[6]*xin[0] + a[7]*xin[1] + a[8]*xin[2] + a[11];

    return output;
};

#endif