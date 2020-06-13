/*
* @Author: jose
* @Date:   2019-11-20 00:00:00
* @Last Modified by:   jose
* @Last Modified time: 2019-11-20 00:00:00
*/

#ifndef __OBJECT_H__
#define __OBJECT_H__

// std libs
#include <iostream>     // std::cout
#include <vector>       // std::vector
#include <sstream>      // std::stringstream
#include <typeinfo>     // operator typeids

// local libs
#include "inherit.h"

namespace imart
{

//! Class object
class object : public inherit<object, base>
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
    // type value;                              // single value of type

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
    // std::string get_type() const;

    // ===========================================
    // Print Functions
    // ===========================================
    void print(std::string msg = "");
    void print_data(std::string msg = "");

    friend std::ostream & operator << (std::ostream & os, object & input);

    // ===========================================
    // Overloading Functions
    // ===========================================
    virtual object & operator = (const object & input);
};


// ===========================================
//          Functions of Class object
// ===========================================

// ===========================================
// Constructor Functions
// ===========================================
//! Constructor empty
object::object()
{
    class_name = "object";
};

//! Constructor to clone
object::object(const object & input)
{
    clone_(input);                // call the virtual function clone_
};

//! Destructor
object::~object()
{
    ;
};

// ===========================================
// Create Functions
// ===========================================
void object::clone_(const object & input)
{
    ;
};

void object::copy_(const object & input)
{
    ;
};

void object::mimic_(const object & input)
{
    ;
};

// ===========================================
// Get Functions
// ===========================================
std::string object::get_name() const
{
    return class_name;
};

// template <typename type>
// std::string object::get_type() const
// {
//     return typeid(value).name();
// };

// ===========================================
// Print Functions
// ===========================================
void object::print(std::string msg)
{
    std::string ss = info(msg);
    std::cout << ss;
};

void object::print_data(std::string msg)
{
    std::string ss = info_data(msg);
    std::cout << ss;
};

std::ostream & operator << (std::ostream & os, object & input)
{
    os << input.info("");
    return os;
};

std::string object::info(std::string msg)
{
    std::stringstream ss;
    std::string title = "Object Information";
    if (msg != "") { title = msg; };
    
    // Summary of the object information
    ss << "\n===== " << title << " =====\n";
    ss << "Pointer: \t\t" << this << std::endl;
    ss << "Class name: \t\t" << get_name() << std::endl;
    // ss << "Data type: \t\t" << get_type() << std::endl;
    return ss.str();
};

std::string object::info_data(std::string msg)
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
object & object::operator = (const object & input)
{
    copy_(input);
    return *this;
};

}; //end namespace

#endif