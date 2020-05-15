/*
* @Author: jose
* @Date:   2019-11-19 13:55:00
* @Last Modified by:   jose
* @Last Modified time: 2019-11-19 13:55:00
*/

#ifndef __GRID_H__
#define __GRID_H__

// std libs
#include <iostream>     // std::cout
#include <memory>       // smart pointers
#include <vector>       // std::vector
#include <cassert>      // assert

// images 
// #include "image_base.h"
#include "image.h"
// #include "image_2d.h"
// #include "image_3d.h"

namespace imart
{

template <typename pixel_type>
class grid: public object<pixel_type>
{
public:
    //Type definitions
    using self    = grid;
    using pointer = std::shared_ptr<self>;
    using vector  = std::vector<self::pointer>;
    
protected:
    //Type definitions
    using vector_image     = std::vector<image<pixel_type>>;
    using ptr_vector_image = std::shared_ptr<vector_image>;

    // ===========================================
    // Internal Variables
    // ===========================================
    ptr_vector_image xyz;

    // ===========================================
    // Functions
    // ===========================================
    void init(int dim);
    void copy_properties(const grid<pixel_type> & input);
    void copy_properties(const image_base<pixel_type> & input);
    // void update(const image_2d<pixel_type> & input);
    // void update(const image_3d<pixel_type> & input);
    void meshgrid_2d();
    void meshgrid_3d();

    void allocate(int dimension);

public:
    // ===========================================
    // Create Functions
    // ===========================================
    grid();
    grid(int dim);
    grid(const grid<pixel_type> & input);
    grid(const image_base<pixel_type> & input);
    // grid(const typename image_base<pixel_type>::pointer input);

    // grid(const image_2d<pixel_type> & input);
    // grid(const image_3d<pixel_type> & input);

    ~grid();

    // static pointer new_pointer();
    template<typename... ARGS>
    static pointer new_pointer(const ARGS&... args);
    // {
    //     return std::make_shared<grid>(args...);
    // };

    void copy(const grid<pixel_type> & input);
    void duplicate(const grid<pixel_type> & input);
    void imitate(const grid<pixel_type> & input);
    
    // ===========================================
    // Get Functions
    // ===========================================
    image<pixel_type> * ptr() const;
    ptr_vector_image get_grid() const;

    // ===========================================
    // Print Functions
    // ===========================================
    // void print(std::string msg = "");
    std::string info(std::string msg = "");

    // void print_data(std::string msg = "");
    std::string info_data(std::string msg = "");

    // ===========================================
    // Overloading Operators
    // ===========================================
    grid<pixel_type> & operator = (const grid<pixel_type> & input);

    // ===========================================
    // Functions
    // ===========================================
    void meshgrid();
    void meshgrid(image_base<pixel_type> & input);
    // void meshgrid(image_2d<pixel_type> & input);
    // void meshgrid(image_3d<pixel_type> & input);
};


// ===========================================
//      Functions of Class grid
// ===========================================

// ===========================================
// Create Functions
// ===========================================
template <typename pixel_type>
grid<pixel_type>::grid()
{
    init(2);
};

template <typename pixel_type>
grid<pixel_type>::grid(int d)
{
    init(d);
}

template <typename pixel_type>
grid<pixel_type>::grid(const grid<pixel_type> & input)
{
    init(input.get_dimension());
    this->class_name = "grid";
    copy(input);
}

template <typename pixel_type>
grid<pixel_type>::grid(const image_base<pixel_type> & input)
{
    copy_properties(input);
    meshgrid();
};

template <typename pixel_type>
grid<pixel_type>::~grid()
{
    xyz.reset();
};

template <typename pixel_type>
void grid<pixel_type>::init(int d)
{
    this->class_name = "grid";
    // xyz.reset();
    // xyz = std::make_shared<vector_image>(d);
    allocate(d);
    object<pixel_type>::init(d);
};

template <typename pixel_type>
template <typename ... ARGS>
typename grid<pixel_type>::pointer grid<pixel_type>::new_pointer(const ARGS&... args)
{
    return std::make_shared< grid<pixel_type> >(args...); // not working for inherited classes
};


template <typename pixel_type>
void grid<pixel_type>::allocate(int dimension)
{
    xyz.reset();
    xyz = std::make_shared< std::vector< image<pixel_type >>>(dimension);
    for (int i = 0; i < dimension; i++)
    {
        (*xyz)[i] = image_base<pixel_type>();
    };

};

template <typename pixel_type>
void grid<pixel_type>::copy(const grid<pixel_type> & input)
{
    int d = input.get_dimension();
    image<pixel_type> * p1 = input.ptr();
    image<pixel_type> * p2 = ptr();
    
    copy_properties(input);
    for (int i = 0; i < d; i++){ p2[i].copy(p1[i]); }; // Copy xyz
};

template <typename pixel_type>
void grid<pixel_type>::duplicate(const grid<pixel_type> & input)
{
    copy_properties(input);
    xyz.reset();
    xyz = input.get_grid();
};

template <typename pixel_type>
void grid<pixel_type>::imitate(const grid<pixel_type> & input)
{
    copy_properties(input);
    xyz.reset();
    xyz = std::make_shared<vector_image>(input.get_dimension());
};

template <typename pixel_type>
void grid<pixel_type>::copy_properties(const grid<pixel_type> & input)
{
    int d = input.get_dimension();

    xyz.reset();
    xyz = std::make_shared<vector_image>(d);
    object<pixel_type>::copy_properties(input);
};

template <typename pixel_type>
void grid<pixel_type>::copy_properties(const image_base<pixel_type> & input)
{
    int d = input.get_dimension();

    xyz.reset();
    xyz = std::make_shared<vector_image>(d);
    object<pixel_type>::copy_properties(input);
};

// template <typename pixel_type>
// void grid<pixel_type>::update(const image_2d<pixel_type> & input)
// {
//     int d = input.get_dimension();

//     xyz.reset();
//     xyz = std::make_shared<vector_image>(d);
//     object<pixel_type>::copy_properties(input);
    
//     // this->dim = input.get_dimension();
//     // this->size = input.get_size();
//     // this->spacing = input.get_spacing();
//     // this->origin = input.get_origin();
//     // this->direction = input.get_direction();
// }

// template <typename pixel_type>
// void grid<pixel_type>::update(const image_3d<pixel_type> & input)
// {
//     int d = input.get_dimension();

//     xyz.reset();
//     xyz = std::make_shared<vector_image>(d);
//     object<pixel_type>::copy_properties(input);
// }

// ===========================================
// Get Functions
// ===========================================
template <typename pixel_type>
typename grid<pixel_type>::ptr_vector_image grid<pixel_type>::get_grid() const
{
    return xyz;
};

template <typename pixel_type>
image<pixel_type> * grid<pixel_type>::ptr() const
{
    return xyz.get()->data();
};

// ===========================================
// Overloading Functions
// ===========================================
// Equal
template <typename pixel_type>
grid<pixel_type> & grid<pixel_type>::operator = (const grid<pixel_type> & input)
{
    duplicate(input);
    return *this;
};

// ===========================================
// Print Functions
// ===========================================
template <typename pixel_type>
std::string grid<pixel_type>::info(std::string msg)
{
    std::stringstream ss;
    std::string title = "Grid Information";
    if (msg != "") { title = msg; };

    ss << object<pixel_type>::info(title);
    
    ss << "Size: \t\t\t[ ";
    for(int i = 0; i < this->size.size(); i++) { ss << this->size[i] << " "; };
    ss << "]" << std::endl;
    
    ss << "Length (mm): \t\t[ ";
    for(int i = 0; i < this->spacing.size(); i++) { ss << this->spacing[i] << " "; };
    ss << "]" << std::endl;

    ss << "Origin (mm): \t\t[ ";
    for(int i = 0; i < this->origin.size(); i++) { ss << this->origin[i] << " "; };
    ss << "]" << std::endl;

    return ss.str();
};

template <typename pixel_type>
std::string grid<pixel_type>::info_data(std::string msg)
{
    std::stringstream ss;
    if (msg != "") { ss << msg << std::endl; };

    std::vector<std::string> abc{"x:","y:","z:","u:","v:","w:"};

    assert(this->dim < 7);
    
    for(int i = 0; i < this->dim; i++)
    {
        ss << (*xyz)[i].info_data(abc[i]);
    }
    
    return ss.str();
};



// ===========================================
// Functions
// ===========================================
template <typename pixel_type>
void grid<pixel_type>::meshgrid()
{   
    if (this->dim == 2) { meshgrid_2d(); };
    if (this->dim == 3) { meshgrid_3d(); };
};

template <typename pixel_type>
void grid<pixel_type>::meshgrid(image_base<pixel_type> & input)
{   
    copy_properties(input);
    meshgrid();
};

template <typename pixel_type>
void grid<pixel_type>::meshgrid_2d()
{
    // TODO: DIRECTION IS DISABLED ****
    // object<pixel_type>::print("internal"); // debug
    int w = this->get_size()[0];
    int h = this->get_size()[1];
    int elements = w*h;

    std::vector<double> s = this->get_spacing();
    std::vector<double> o = this->get_origin();
    std::vector<double> d = this->get_direction();

    image<pixel_type> x(w,h);
    image<pixel_type> y(w,h);
    
    pixel_type * px = x.ptr();
    pixel_type * py = y.ptr();

    // #pragma omp parallel for
    for(int j = 0; j < h; j++)
    {
        for(int i = 0; i < w; i++)
        {
            px[i+j*w] = s[0]*i + o[0];
            py[i+j*w] = s[1]*j + o[1];
            // with direction
            // px[i+j*w] = d[0]*s[0]*i + d[1]*s[1]*j + o[0];
            // py[i+j*w] = d[2]*s[0]*i + d[3]*s[1]*j + o[1];
        }
    }

    (*xyz)[0] = x;
    (*xyz)[1] = y;
};

template <typename pixel_type>
void grid<pixel_type>::meshgrid_3d()
{
    // TODO: DIRECTION IS DISABLED ****
    int w = this->get_size()[0];
    int h = this->get_size()[1];
    int l = this->get_size()[2];
    int elements = w*h*l;

    std::vector<double> s = this->get_spacing();
    std::vector<double> o = this->get_origin();
    std::vector<double> d = this->get_direction();

    image<pixel_type> x(w,h,l);
    image<pixel_type> y(w,h,l);
    image<pixel_type> z(w,h,l);
    
    pixel_type * px = x.ptr();
    pixel_type * py = y.ptr();
    pixel_type * pz = z.ptr();

    // #pragma omp parallel for
    for(int k = 0; k < l; k++)
    {
        for(int j = 0; j < h; j++)
        {
            for(int i = 0; i < w; i++)
            {
                px[i + j*w + k*w*h] = s[0]*i + o[0];
                py[i + j*w + k*w*h] = s[1]*j + o[1];
                pz[i + j*w + k*w*h] = s[2]*k + o[2];
            };
        };
    };

    (*xyz)[0] = x;
    (*xyz)[1] = y;
    (*xyz)[2] = z;
};

// template <typename pixel_type>
// void grid<pixel_type>::meshgrid(image_2d<pixel_type> & input)
// {   
//     update(input);
//     meshgrid_2d();
// };

// template <typename pixel_type>
// void grid<pixel_type>::meshgrid(image_3d<pixel_type> & input)
// {   
//     update(input);
//     meshgrid_3d();
// };

}; //end namespace

#endif