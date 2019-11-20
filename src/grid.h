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
#include "image_base.h"
#include "image_2d.h"
// #include "image_3d.h"

template <typename pixel_type>
class grid 
{
    //Type definitions
    template<typename pixel_t> using image_vector = std::vector<image_base<pixel_t>>;
    template<typename pixel_t> using ptr_image_vector = std::shared_ptr<image_vector<pixel_t>>;
    
protected:
    // ===========================================
    // Internal Variables
    // ===========================================
    int dim;
    std::vector<int> size;
    std::vector<pixel_type> spacing;
    std::vector<pixel_type> origin;
    std::vector<pixel_type> direction;

    ptr_image_vector<pixel_type> xyz;

    // ===========================================
    // Functions
    // ===========================================
    void init(int dim);
    void update(grid<pixel_type> & input);
    void update(image_base<pixel_type> & input);

public:
    // ===========================================
    // Create Functions
    // ===========================================
    grid();
    grid(int dim);
    ~grid();
    
    // ===========================================
    // Get Functions
    // ===========================================
    int get_dimension();
    std::vector<int> get_size();
    std::vector<pixel_type> get_spacing();
    std::vector<pixel_type> get_origin();
    std::vector<pixel_type> get_direction();

    // ===========================================
    // Print Functions
    // ===========================================
    void print(std::string msg = "");
    std::string info(std::string msg = "");

    void print_data(std::string msg = "");
    std::string info_data(std::string msg = "");

    // ===========================================
    // Overloading Operators
    // ===========================================
    grid<pixel_type> & operator = (grid<pixel_type> input);

    // ===========================================
    // Functions
    // ===========================================
    void meshgrid(image_2d<pixel_type> & input);
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
grid<pixel_type>::~grid()
{
    ;
};

template <typename pixel_type>
void grid<pixel_type>::init(int d)
{
    dim = d;
    xyz = std::make_shared<image_vector<pixel_type>>(dim);
    size = std::vector<int>(dim, 0);
    spacing = std::vector<pixel_type>(dim, 1.0);
    origin = std::vector<pixel_type>(dim, 0.0);
    direction = std::vector<pixel_type>(dim);
}

template <typename pixel_type>
void grid<pixel_type>::update(grid<pixel_type> & input)
{
    size = input.get_size();
    spacing = input.get_spacing();
    origin = input.get_origin();
    direction = input.get_direction();
}

template <typename pixel_type>
void grid<pixel_type>::update(image_base<pixel_type> & input)
{
    size = input.get_size();
    spacing = input.get_spacing();
    origin = input.get_origin();
    direction = input.get_direction();
}

// ===========================================
// Overloading Functions
// ===========================================
// Equal
template <typename pixel_type>
grid<pixel_type> & grid<pixel_type>::operator = (grid<pixel_type> input)
{
    // delete &data;
    this->xyz.reset();
    this->xyz = input.get_grid();
    this->update(input);
    return *this;
};


// ===========================================
// Print Functions
// ===========================================
template <typename pixel_type>
void grid<pixel_type>::print(std::string msg)
{
    std::cout << this->info(msg);
};

template <typename pixel_type>
std::string grid<pixel_type>::info(std::string msg)
{
    std::stringstream ss;
    std::string title = "Grid Information";
    if (msg != "") { title = msg; };

    // Summary of the image information
    ss << "\n===== " << title << " =====\n";
    
    ss << "Pixel type: \t\t" << typeid(spacing[0]).name() << std::endl;
    ss << "Dimensions: \t\t" << dim << std::endl;
    
    ss << "Size: \t\t\t[ ";
    for(int i = 0; i < size.size(); i++) { ss << size[i] << " "; };
    ss << "]" << std::endl;
    
    ss << "Length (mm): \t\t[ ";
    for(int i = 0; i < spacing.size(); i++) { ss << spacing[i] << " "; };
    ss << "]" << std::endl;

    ss << "Origin (mm): \t\t[ ";
    for(int i = 0; i < origin.size(); i++) { ss << origin[i] << " "; };
    ss << "]" << std::endl;

    return ss.str();
};

template <typename pixel_type>
void grid<pixel_type>::print_data(std::string msg)
{
    std::cout << this->info_data(msg) << std::endl;
};

template <typename pixel_type>
std::string grid<pixel_type>::info_data(std::string msg)
{
    std::stringstream ss;
    if (msg != "") { ss << msg << std::endl; };

    std::vector<std::string> abc{"x:","y:","z:","u:","v:","w:"};

    assert(dim < 7);
    
    for(int i = 0; i < dim; i++)
    {
        ss << (*(this->xyz))[i].info_data(abc[i]);
    }
    
    return ss.str();
};



// ===========================================
// Functions
// ===========================================
template <typename pixel_type>
void grid<pixel_type>::meshgrid(image_2d<pixel_type> & input)
{
    assert(this->dim == input.get_dimension());
    // using vector_image = std::vector<image_base<pixel_type>>;
    int w = input.get_width();
    int h = input.get_height();
    int elements = w*h;

    std::vector<pixel_type> s = input.get_spacing();
    std::vector<pixel_type> o = input.get_origin();
    std::vector<pixel_type> d = input.get_direction();

    image_2d<pixel_type> x(w,h);
    image_2d<pixel_type> y(w,h);
    // std::unique_ptr<vec> y = std::make_unique<vec>(elements);
    std::shared_ptr<pixel_type[]> px = x.get_data();
    std::shared_ptr<pixel_type[]> py = y.get_data();

    pixel_type u;
    pixel_type v;

    #pragma omp parallel for
    for(int j = 0; j < h; j++)
    {
        for(int i = 0; i < w; i++)
        {
            u = (pixel_type)i;
            v = (pixel_type)j;
            px[i+j*w] = d[0]*s[0]*u + d[1]*s[1]*v + o[0];
            py[i+j*w] = d[2]*s[0]*u + d[3]*s[1]*v + o[1];
        }
    }

    (*(this->xyz))[0] = x;
    (*(this->xyz))[1] = y;
    update(input);
    
    // this->grid = std::make_shared<vector_image>(2);
    // (*(this->grid))[0] = x;
    // (*(this->grid))[1] = y;
};


#endif