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

// local libs
#include "image.h"

namespace imart
{

template <typename type, typename container=vector_cpu<type>>
class grid: public inherit<grid<type,container>,space_object>
{
public:
    //Type definitions
    using self    = grid;
    using pointer = std::shared_ptr<self>;
    using vector  = std::vector<self::pointer>;
    
protected:
    //Type definitions
    using vector_image     = std::vector<typename image<type,container>::pointer>;
    using ptr_vector_image = std::shared_ptr<vector_image>;

    // ===========================================
    // Internal Variables
    // ===========================================
    ptr_vector_image xyz;

    // ===========================================
    // Functions
    // ===========================================
    void init(int dim);
    void allocate(int dimension);
    void copy_properties(const grid<type,container> & input);
    void copy_properties(const image<type,container> & input);
    
    void meshgrid_2d();
    void meshgrid_3d();

public:
    // ===========================================
    // Constructor Functions
    // ===========================================
    grid();
    grid(int dim);
    grid(const grid<type,container> & input);
    grid(const image<type,container> & input);
    // ~grid();

    // ===========================================
    // Create Functions
    // ===========================================
    void clone_(const grid<type,container> & input);
    void copy_(const grid<type,container> & input);
    void mimic_(const grid<type,container> & input);
    
    // ===========================================
    // Get Functions
    // ===========================================
    typename image<type,container>::pointer * ptr() const;
    ptr_vector_image get_grid() const;

    // ===========================================
    // Print Functions
    // ===========================================
    std::string info(std::string msg = "");
    std::string info_data(std::string msg = "");

    // ===========================================
    // Overloading Operators
    // ===========================================
    grid<type,container> & operator = (const grid<type,container> & input);

    // ===========================================
    // Functions
    // ===========================================
    void meshgrid();
    void meshgrid(image<type,container> & input);
};

template<typename type>
using grid_cpu = grid<type,vector_cpu<type>>;

template<typename type>
using grid_gpu = grid<type,vector_ocl<type>>;


// ===========================================
//      Functions of Class grid
// ===========================================

// ===========================================
// Create Functions
// ===========================================
template <typename type, typename container>
grid<type,container>::grid()
{
    this->class_name = "grid";
    init(2);
};

template <typename type, typename container>
grid<type,container>::grid(int d)
{
    this->class_name = "grid";
    init(d);
}

template <typename type, typename container>
grid<type,container>::grid(const grid<type,container> & input)
{
    this->class_name = "grid";
    allocate(input.get_dimension());
    clone_(input);
}

template <typename type, typename container>
grid<type,container>::grid(const image<type,container> & input)
{
    this->class_name = "grid";
    allocate(input.get_dimension());
    copy_properties(input);
    meshgrid();
};

// template <typename type, typename container>
// grid<type,container>::~grid()
// {
//     ;// xyz.reset();
// };

template <typename type, typename container>
void grid<type,container>::init(int d)
{
    allocate(d);
    space_object::init(d);
};

template <typename type, typename container>
void grid<type,container>::allocate(int dimension)
{
    xyz.reset();
    xyz = std::make_shared< std::vector< typename image<type,container>::pointer >>(dimension);
    for (int i = 0; i < dimension; i++)
    {
        (*xyz)[i] = image<type,container>::new_pointer(this->get_size(),1);
    };
};

template <typename type, typename container>
void grid<type,container>::clone_(const grid<type,container> & input)
{
    copy_properties(input);
    typename image<type,container>::pointer * p1 = input.ptr();
    typename image<type,container>::pointer * p2 = ptr();
    int d = input.get_dimension();
    for (int i = 0; i < d; i++){ p2[i] = p1[i]->clone(); }; // Copy xyz
};

template <typename type, typename container>
void grid<type,container>::copy_(const grid<type,container> & input)
{
    space_object::mimic_(input);
    xyz.reset();
    xyz = input.get_grid();
};

template <typename type, typename container>
void grid<type,container>::mimic_(const grid<type,container> & input)
{
    space_object::mimic_(input);
    allocate(input.get_dimension());
};

template <typename type, typename container>
void grid<type,container>::copy_properties(const grid<type,container> & input)
{
    space_object::mimic_(input);
};

template <typename type, typename container>
void grid<type,container>::copy_properties(const image<type,container> & input)
{
    space_object::mimic_(input);
};

// ===========================================
// Get Functions
// ===========================================
template <typename type, typename container>
typename grid<type,container>::ptr_vector_image grid<type,container>::get_grid() const
{
    return xyz;
};

template <typename type, typename container>
typename image<type,container>::pointer * grid<type,container>::ptr() const
{
    return xyz->data();
};

// ===========================================
// Overloading Functions
// ===========================================
// Equal
template <typename type, typename container>
grid<type,container> & grid<type,container>::operator = (const grid<type,container> & input)
{
    copy_(input);
    return *this;
};

// ===========================================
// Print Functions
// ===========================================
template <typename type, typename container>
std::string grid<type,container>::info(std::string msg)
{
    std::stringstream ss;
    std::string title = "Grid Information";
    if (msg != "") { title = msg; };

    ss << space_object::info(title);
    return ss.str();
};

template <typename type, typename container>
std::string grid<type,container>::info_data(std::string msg)
{
    std::stringstream ss;
    if (msg != "") { ss << msg << std::endl; };
    std::vector<std::string> abc{"x:","y:","z:","u:","v:","w:"};
    assert(this->dim < 7);
    for(int i = 0; i < this->dim; i++)
    {
        ss << (*xyz)[i]->info_data(abc[i]);
    }
    
    return ss.str();
};

// ===========================================
// Functions
// ===========================================
template <typename type, typename container>
void grid<type,container>::meshgrid()
{   
    if (this->dim == 2) { meshgrid_2d(); };
    if (this->dim == 3) { meshgrid_3d(); };
};

template <typename type, typename container>
void grid<type,container>::meshgrid(image<type,container> & input)
{
    allocate(input.get_dimension());
    copy_properties(input);
    meshgrid();
};

template <typename type, typename container>
void grid<type,container>::meshgrid_2d()
{
    int w = this->get_size()[0];
    int h = this->get_size()[1];

    std::vector<double> p = this->get_sod_parameters();

    auto imgx = image<type,container>::new_pointer(w,h);
    auto imgy = image<type,container>::new_pointer(w,h);

    std::vector<typename container::pointer> gxy;
    gxy = container::grid_2d(w,h,p);    // grid computation

    imgx->set_data(gxy[0]);
    imgy->set_data(gxy[1]);

    (*xyz)[0] = imgx;
    (*xyz)[1] = imgy;
};

template <typename type, typename container>
void grid<type,container>::meshgrid_3d()
{
    int w = this->get_size()[0];
    int h = this->get_size()[1];
    int l = this->get_size()[2];

    std::vector<double> p = this->get_sod_parameters();

    auto imgx = image<type,container>::new_pointer(w,h,l);
    auto imgy = image<type,container>::new_pointer(w,h,l);
    auto imgz = image<type,container>::new_pointer(w,h,l);
    
    std::vector<typename container::pointer> gxyz;
    gxyz = container::grid_3d(w,h,l,p); // grid computation

    imgx->set_data(gxyz[0]);
    imgy->set_data(gxyz[1]);
    imgz->set_data(gxyz[2]);

    (*xyz)[0] = imgx;
    (*xyz)[1] = imgy;
    (*xyz)[2] = imgz;
};

}; //end namespace

#endif