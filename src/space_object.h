/*
* @Author: jose
* @Date:   2019-11-20 00:00:00
* @Last Modified by:   jose
* @Last Modified time: 2019-11-20 00:00:00
*/

#ifndef __SPACE_OBJECT_H__
#define __SPACE_OBJECT_H__

// std libs
#include <iostream>     // std::cout
#include <cassert>      // assert

// local libs
#include "inherit.h"
#include "object.h"

namespace imart
{

// Class object
class space_object: public inherit<space_object, object>
{
public:
    //Type definitions
    using self    = space_object;
    using pointer = std::shared_ptr<self>;
    using vector  = std::vector<self::pointer>;

protected:
    // ===========================================
    // Internal Variables
    // ===========================================
    int dim;                        // dimension
    std::vector<int>    size;       // object size
    std::vector<double> spacing;    // spacing between elements
    std::vector<double> origin;     // origin of coordinates
    std::vector<double> direction;  // direction of elements

    // ===========================================
    // Functions
    // ===========================================
    virtual void init(int d);                                       // init default properties
    // virtual void copy_properties(const space_object & input);    // copy only properties
    virtual std::string info(std::string msg);

public:
    // ===========================================
    // Constructor Functions
    // ===========================================
    space_object();                                 // constructor empty
    space_object(int d);                            // constructor with dimension
    space_object(const space_object & input);       // constructor with same type    
    ~space_object();                                // destructor empty

    // ===========================================
    // Create Functions
    // ===========================================
    virtual void clone_(const space_object & input);// copy everything
    virtual void copy_ (const space_object & input);// share data
    virtual void mimic_(const space_object & input);// copy meta data

    // ===========================================
    // Get Functions
    // ===========================================
    int get_dimension() const;
    std::vector<int> get_size() const;
    std::vector<double> get_spacing() const;
    std::vector<double> get_origin() const;
    std::vector<double> get_direction() const;

    // ===========================================
    // Set Functions
    // ===========================================
    void set_spacing(std::vector<double> s);
    void set_origin(std::vector<double> o);
    void set_direction(std::vector<double> d);
};


// ===========================================
//          Functions of Class object
// ===========================================

// ===========================================
// Create Functions
// ===========================================
//! Constructor empty
space_object::space_object()
{
    this->class_name = "space_object";
    init(1);
};

//! Constructor with number of dimensions
space_object::space_object(int d)
{
    this->class_name = "space_object";
    init(d);
};

//! Constructor to clone
space_object::space_object(const space_object & input)
{
    clone_(input);                // call the virtual
    // copy_properties(input);    // in object class everything is copied with this method)
};

// Destructor
space_object::~space_object()
{
    ;
};

// ===========================================
// Create Functions
// ===========================================
// Initialization function
void space_object::init(int d)
{
    // Attributes initialization
    dim = d;
    size = std::vector<int>(dim, 0);
    spacing = std::vector<double>(dim, 1.0);
    origin = std::vector<double>(dim, 0.0);
    direction = std::vector<double>(dim*dim);

    // Initialize direction, identity matrix
    int den = dim + 1;
    for(int i=0; i < dim*dim; i++){ if((i%den)==0) { direction[i] = 1.0; }; };
};

// template <typename type>
// void space_object::copy_properties(const space_object & input)
// {

// };

void space_object::clone_(const space_object & input)
{
    mimic_(input);
};

void space_object::copy_(const space_object & input)
{
    mimic_(input);
};

void space_object::mimic_(const space_object & input)
{
    dim = input.get_dimension();
    size = input.get_size();
    spacing = input.get_spacing();
    origin = input.get_origin();
    direction = input.get_direction();
};

// ===========================================
// Get Functions
// ===========================================
int space_object::get_dimension() const
{
    return dim;
};

std::vector<int> space_object::get_size() const
{
    return size;
};

std::vector<double> space_object::get_spacing() const
{
    return spacing;
};

std::vector<double> space_object::get_origin() const
{
    return origin;
};

std::vector<double> space_object::get_direction() const
{
    return direction;
};

// ===========================================
// Set Functions
// ===========================================
void space_object::set_spacing(std::vector<double> s)
{
    assert(dim == s.size());
    spacing = s;
};

void space_object::set_origin(std::vector<double> o)
{
    assert(dim == o.size());
    origin = o;
};

void space_object::set_direction(std::vector<double> d)
{
    assert(dim*dim == d.size());
    direction = d;
};

// ===========================================
// Print Functions
// ===========================================
std::string space_object::info(std::string msg)
{
    std::stringstream ss;
    std::string title = "Space Object Information";
    if (msg != "") { title = msg; };

    // Summary of the object information
    ss << object::info(title);
    ss << "Dimensions: \t\t" << get_dimension() << std::endl;
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

}; //end namespace

#endif