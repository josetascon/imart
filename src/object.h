/*
* @Author: jose
* @Date:   2019-11-20 00:00:00
* @Last Modified by:   jose
* @Last Modified time: 2019-11-20 00:00:00
*/

#ifndef __OBJECT_H__
#define __OBJECT_H__

// std libs
#include <vector>       // std::vector
#include <sstream>      // std::stringstream
#include <typeinfo>     // operator typeids

// local libs
#include "inherit.h"

namespace imart
{

//! Class object
template <typename type>
class object : public inherit<object<type>, base>
{
public:
    // Type definitions
    using self    = object;
    using pointer = std::shared_ptr<self>;
    using vector  = std::vector<self::pointer>;

protected:
    // ===========================================
    // Internal Variables
    // ===========================================
    std::string class_name;                     // class string name
    type value;                                 // single value of type

    // ===========================================
    // Print Functions
    // ===========================================
    virtual std::string info(std::string msg);
    virtual std::string info_data(std::string msg);

public:
    // ===========================================
    // Constructor Functions
    // ===========================================
    object();                                   // constructor empty
    object(const object & input);               // constructor clone
    ~object();                                  // destructor empty

    // ===========================================
    // Create Functions
    // ===========================================
    virtual void clone_(const object & input);  // copy everything
    virtual void copy_ (const object & input);  // copy properties and share data
    virtual void mimic_(const object & input);  // copy only properties

    // ===========================================
    // Get Functions
    // ===========================================
    std::string get_name() const;
    std::string get_type() const;

    // ===========================================
    // Print Functions
    // ===========================================
    void print(std::string msg = "");
    void print_data(std::string msg = "");

    template<typename pixel_t>
    friend std::ostream & operator << (std::ostream & os, object<pixel_t> & input);

    // ===========================================
    // Overloading Functions
    // ===========================================
    virtual object<type> & operator = (const object<type> & input);
};


// ===========================================
//          Functions of Class object
// ===========================================

// ===========================================
// Create Functions
// ===========================================
//! Constructor empty
template <typename type>
object<type>::object()
{
    class_name = "object";
};

//! Constructor to clone
template <typename type>
object<type>::object(const object<type> & input)
{
    clone_(input);                // call the virtual function clone_
};

//! Destructor
template <typename type>
object<type>::~object()
{
    ;
};

template <typename type>
void object<type>::clone_(const object & input)
{
    ;
};

template <typename type>
void object<type>::copy_(const object & input)
{
    ;
};

template <typename type>
void object<type>::mimic_(const object & input)
{
    ;
};

// ===========================================
// Get Functions
// ===========================================
template <typename type>
std::string object<type>::get_name() const
{
    return class_name;
};

template <typename type>
std::string object<type>::get_type() const
{
    return typeid(value).name();
};

// ===========================================
// Print Functions
// ===========================================
template <typename type>
void object<type>::print(std::string msg)
{
    std::string ss = info(msg);
    std::cout << ss;
};

template <typename type>
void object<type>::print_data(std::string msg)
{
    std::string ss = info_data(msg);
    std::cout << ss;
};

template <typename type>
std::ostream & operator << (std::ostream & os, object<type> & input)
{
    os << input.info("");
    return os;
};

template <typename type>
std::string object<type>::info(std::string msg)
{
    std::stringstream ss;
    std::string title = "Object Information";
    if (msg != "") { title = msg; };
    
    // Summary of the object information
    ss << "\n===== " << title << " =====\n";
    ss << "Pointer: \t\t" << this << std::endl;
    ss << "Class name: \t\t" << get_name() << std::endl;
    ss << "Data type: \t\t" << get_type() << std::endl;

    return ss.str();
};

template <typename type>
std::string object<type>::info_data(std::string msg)
{
    // Totally override when implemented in inherited classes
    std::stringstream ss;
    ss << "This method of " << get_name();
    ss << " is not implemented" << std::endl;
    return ss.str();
};

// ===========================================
// Overloading Functions
// ===========================================
template <typename type>
object<type> & object<type>::operator = (const object<type> & input)
{
    copy_(input);
    return *this;
};

}; //end namespace

#endif