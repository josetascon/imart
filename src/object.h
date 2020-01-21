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
#include <sstream>      // stringstream
#include <vector>       // std::vector
#include <typeinfo>     // operator typeids

// Class object
template <typename pixel_type>
class object
{
protected:
    // ===========================================
    // Internal Variables
    // ===========================================
    std::string class_name;             // class string name
    pixel_type value;                   // single value of pixel_type

    int dim;                            // dimension
    std::vector<int> size;              // object size
    std::vector<pixel_type> spacing;    // spacing between elements
    std::vector<pixel_type> origin;     // origin of coordinates
    std::vector<pixel_type> direction;  // direction of elements

    // ===========================================
    // Functions
    // ===========================================
    virtual void init(int d);                   // init default properties
    virtual void update(const object & input);  // copy only properties

public:
    // ===========================================
    // Create Functions
    // ===========================================
    object();                           // constructor empty
    object(int d);                      // constructor with dimension
    object(const object & input);       // constructor with same type
    
    ~object();                          // destructor empty

    // ===========================================
    // Get Functions
    // ===========================================
    std::string get_name() const;
    std::string get_type() const;
    int get_dimension() const;
    std::vector<int> get_size() const;
    std::vector<pixel_type> get_spacing() const;
    std::vector<pixel_type> get_origin() const;
    std::vector<pixel_type> get_direction() const;

    // ===========================================
    // Set Functions
    // ===========================================
    void set_spacing(std::vector<pixel_type> s) const;
    void set_origin(std::vector<pixel_type> o) const;
    void set_direction(std::vector<pixel_type> d) const;
    
    // ===========================================
    // Print Functions
    // ===========================================
    void print(std::string msg = "");
    void print_data(std::string msg = "");

    template<typename pixel_t>
    friend std::ostream & operator << (std::ostream & os, object<pixel_t> & input);

    virtual std::string info(std::string msg);
    virtual std::string info_data(std::string msg);
};







// ===========================================
//          Functions of Class object
// ===========================================

// ===========================================
// Create Functions
// ===========================================
template <typename pixel_type>
object<pixel_type>::object()
{
    class_name = "object";
    init(1);
};

template <typename pixel_type>
object<pixel_type>::object(int d)
{
    class_name = "object";
    init(d);
};

template <typename pixel_type>
object<pixel_type>::object(const object<pixel_type> & input)
{
    class_name = "object";
    update(input);
};

// Destructor
template <typename pixel_type>
object<pixel_type>::~object()
{
    ;
};

template <typename pixel_type>
void object<pixel_type>::init(int d)
{
    dim = d;
    size = std::vector<int>{0, 0};
    spacing = std::vector<pixel_type>(dim, 1.0);
    origin = std::vector<pixel_type>(dim, 0.0);
    direction = std::vector<pixel_type>(dim*dim);

    // initialize direction, identity matrix
    int den = dim + 1;
    for(int i=0; i < dim*dim; i++){ if((i%den)==0) { direction[i] = 1.0; }; };
};

template <typename pixel_type>
void object<pixel_type>::update(const object<pixel_type> & input)
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
template <typename pixel_type>
std::string object<pixel_type>::get_name() const
{
    return class_name;
};

template <typename pixel_type>
std::string object<pixel_type>::get_type() const
{
    return typeid(value).name();
};

template <typename pixel_type>
int object<pixel_type>::get_dimension() const
{
    return dim;
};

template <typename pixel_type>
std::vector<int> object<pixel_type>::get_size() const
{
    return size;
};

template <typename pixel_type>
std::vector<pixel_type> object<pixel_type>::get_spacing() const
{
    return spacing;
};

template <typename pixel_type>
std::vector<pixel_type> object<pixel_type>::get_origin() const
{
    return origin;
};

template <typename pixel_type>
std::vector<pixel_type> object<pixel_type>::get_direction() const
{
    return direction;
};

// ===========================================
// Set Functions
// ===========================================
template <typename pixel_type>
void object<pixel_type>::set_spacing(std::vector<pixel_type> s) const
{
    assert(dim == s.size());
    spacing = s;
};

template <typename pixel_type>
void object<pixel_type>::set_origin(std::vector<pixel_type> o) const
{
    assert(dim == o.size());
    origin = o;
};

template <typename pixel_type>
void object<pixel_type>::set_direction(std::vector<pixel_type> d) const
{
    assert(dim == d.size());
    direction = d;
};

// ===========================================
// Print Functions
// ===========================================
template <typename pixel_type>
void object<pixel_type>::print(std::string msg)
{
    std::string ss = info(msg);
    std::cout << ss;
};

template <typename pixel_type>
void object<pixel_type>::print_data(std::string msg)
{
    std::string ss = info_data(msg);
    std::cout << ss;
};

template <typename pixel_type>
std::ostream & operator << (std::ostream & os, object<pixel_type> & input)
{
    os << input.info("");
    return os;
};

template <typename pixel_type>
std::string object<pixel_type>::info(std::string msg)
{
    std::stringstream ss;
    std::string title = "Object Information";
    if (msg != "") { title = msg; };
    // Summary of the object information
    ss << "\n===== " << title << " =====\n";
    ss << "Class name: \t\t" << get_name() << std::endl;
    ss << "Data type: \t\t" << get_type() << std::endl;
    ss << "Dimensions: \t\t" << get_dimension() << std::endl;
    return ss.str();
};

template <typename pixel_type>
std::string object<pixel_type>::info_data(std::string msg)
{
    // Totally override when implemented in inherited classes
    std::stringstream ss;
    // if (msg != "") { ss << msg << std::endl; };
    ss << "This method of " << get_name();
    ss << " is not implemented" << std::endl;
    return ss.str();
};


#endif