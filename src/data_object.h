/*
* @Author: jose
* @Date:   2019-11-20 00:00:00
* @Last Modified by:   jose
* @Last Modified time: 2019-11-20 00:00:00
*/

#ifndef __DATA_OBJECT_H__
#define __DATA_OBJECT_H__

// std libs
#include <iostream>     // std::cout

// local libs
#include "inherit.h"
#include "object.h"

namespace imart
{

// Class object
template <typename type, typename container>
class data_object: public inherit_multiple<data_object<type,container>, virtualize<object>>
{
public:
    //Type definitions
    using self    = data_object;
    using pointer = std::shared_ptr<self>;
    using vector  = std::vector<self::pointer>;

    //Abbreviations
    using container_pointer = std::shared_ptr<container>;

protected:
    // ===========================================
    // Internal Variables
    // ===========================================
    int num_elements;                               // number of elements
    container_pointer data;                         // data storage

    // ===========================================
    // Functions
    // ===========================================
    virtual void init(int n);                                       // init default properties
    virtual void allocate(int total_elements);
    // virtual void copy_properties(const data_object & input);    // copy only properties
    virtual std::string info(std::string msg);
    virtual std::string info_data(std::string msg);

public:
    // ===========================================
    // Constructor Functions
    // ===========================================
    data_object();                                 // constructor empty
    data_object(int n);                            // constructor with dimension
    data_object(const data_object & input);        // constructor with same type    
    ~data_object();                                // destructor empty

    // ===========================================
    // Create Functions
    // ===========================================
    virtual void clone_(const data_object & input);// copy everything
    virtual void copy_ (const data_object & input);// share data
    virtual void mimic_(const data_object & input);// copy meta data

    // ===========================================
    // Get Functions
    // ===========================================
    std::string get_type() const;
    int get_total_elements() const;
    container_pointer get_data() const;
    type * ptr() const;

    // ===========================================
    // Set Functions
    // ===========================================
    void set_data(container_pointer d);

    // ===========================================
    // Functions
    // ===========================================
    // Free memory. Use when temporary data objects are created and want to force memory cleaning before destruction
    void clear();
};


// ===========================================
//          Functions of Class object
// ===========================================

// ===========================================
// Create Functions
// ===========================================
//! Constructor empty
template <typename type, typename container>
data_object<type,container>::data_object()
{
    this->class_name = "data_object";
    init(0);
};

//! Constructor with number of dimensions
template <typename type, typename container>
data_object<type,container>::data_object(int total_elements)
{
    this->class_name = "data_object";
    init(total_elements);
};

//! Constructor to clone
template <typename type, typename container>
data_object<type,container>::data_object(const data_object & input)
{
    this->class_name = "data_object";
    clone_(input);                // call the virtual
};

// Destructor
template <typename type, typename container>
data_object<type,container>::~data_object()
{
    ;
};

// Initialization function
template <typename type, typename container>
void data_object<type,container>::init(int total_elements)
{
    // Attributes initialization
    num_elements = total_elements;
    allocate(total_elements);
};

template <typename type, typename container>
void data_object<type,container>::allocate(int total_elements)
{
    data.reset();
    data = std::make_shared<container>(total_elements);
};

// Free memory
template <typename type, typename container>
void data_object<type,container>::clear()
{
    data.reset();
};

template <typename type, typename container>
void data_object<type,container>::clone_(const data_object & input)
{
    num_elements = input.get_total_elements();
    container_pointer ptr_data = input.get_data();
    data.reset();
    data = std::make_shared<container>(*ptr_data); // using the default clone constructor of container
};

template <typename type, typename container>
void data_object<type,container>::copy_(const data_object & input)
{
    num_elements = input.get_total_elements();
    data.reset();
    data = input.get_data();
};

template <typename type, typename container>
void data_object<type,container>::mimic_(const data_object & input)
{
    init(input.get_total_elements());
};

// ===========================================
// Get Functions
// ===========================================
template <typename type, typename container>
std::string data_object<type,container>::get_type() const
{   
    // type v;
    // return typeid(v).name();
    return string_type<type>();
};

template <typename type, typename container>
int data_object<type,container>::get_total_elements() const
{
    return num_elements;
};

template <typename type, typename container>
std::shared_ptr<container> data_object<type,container>::get_data() const
{
    return data;
};

template <typename type, typename container>
type * data_object<type,container>::ptr() const
{
    return data->data();
};

template <typename type, typename container>
void data_object<type,container>::set_data(std::shared_ptr<container> d)
{
    data.reset();
    data = d;
};

// ===========================================
// Print Functions
// ===========================================
template <typename type, typename container>
std::string data_object<type,container>::info(std::string msg)
{
    std::stringstream ss;
    std::string title = "Data Object Information";
    if (msg != "") { title = msg; };

    // Summary of the object information
    // ss << object::info(title);
    ss << "Total elements: \t" << num_elements << std::endl;    //Get the total number of pixels
    ss << "Data type: \t\t" << get_type() << std::endl;
    if (data != nullptr) ss << "Container type: \t" << data->get_name() << std::endl;
    else ss << "Container type: \t" << "-" << std::endl;
    ss << "Container pointer: \t" << data << std::endl;         //Print the pointer of container
    return ss.str();
};

template <typename type, typename container>
std::string data_object<type,container>::info_data(std::string msg)
{
    std::stringstream ss;
    ss << data->info_data(msg);
    return ss.str();
};

}; //end namespace

#endif