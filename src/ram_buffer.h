/*
* @Author: jose
* @Date:   2019-11-20 00:00:00
* @Last Modified by:   jose
* @Last Modified time: 2019-11-20 00:00:00
*/

#ifndef __RAM_BUFFER_H__
#define __RAM_BUFFER_H__

// std libs
#include <iostream>     // std::cout
#include <vector>       // std::vector
#include <sstream>      // std::stringstream

// local libs
#include "object.h"

namespace imart
{

//! Class ram_buffer
template <typename type>
class ram_buffer : public inherit<ram_buffer<type>, object>
{
public:
    // Type definitions
    using self    = ram_buffer;
    using pointer = std::shared_ptr<self>;
    using vector  = std::vector<self::pointer>;

protected:
    // ===========================================
    // Internal Variables
    // ===========================================
    size_t _size_;
    type * mem;

    // ===========================================
    // Functions
    // ===========================================
    void allocate_memory();
    void free_memory();

public:
    // ===========================================
    // Constructor Functions
    // ===========================================
    ram_buffer();                                   // constructor empty
    ram_buffer(size_t s);                           // constructor empty
    ram_buffer(size_t s, type value);               // constructor empty
    ram_buffer(std::initializer_list<type> list);
    ram_buffer(const ram_buffer & input);           // constructor clone
    ~ram_buffer();                                  // destructor empty

    // ===========================================
    // Create Functions
    // ===========================================
    virtual void clone_(const ram_buffer & input);  // copy everything
    virtual void copy_ (const ram_buffer & input);  // copy properties and share data
    virtual void mimic_(const ram_buffer & input);  // copy only properties

    // ===========================================
    // Get Functions
    // ===========================================
    int size() const;
    type * data() const;

    // ===========================================
    // Print Functions
    // ===========================================
    virtual std::string info(std::string msg);

    // ===========================================
    // Functions
    // ===========================================
    void resize(size_t sz);
    void clear();

    // ===========================================
    // Overloading Functions
    // ===========================================
    virtual ram_buffer & operator = (const ram_buffer & input);
    type & operator [] (size_t idx);
};


// ===========================================
//          Functions of Class ram_buffer
// ===========================================

// ===========================================
// Constructor Functions
// ===========================================
//! Constructor empty
template <typename type>
ram_buffer<type>::ram_buffer()
{
    this->class_name = "ram_buffer";
    _size_ = 0;
    allocate_memory();
};

template <typename type>
ram_buffer<type>::ram_buffer(size_t s)
{
    this->class_name = "ram_buffer";
    _size_ = size_t(s);
    allocate_memory();
};

template <typename type>
ram_buffer<type>::ram_buffer(std::initializer_list<type> list)
{
    // printf("constructor init list\n");
    this->class_name = "ram_buffer";
    _size_ = list.size();
    allocate_memory();
    if (_size_ > 0)
    {
        const type * d = list.begin();
        for(size_t i=0; i<_size_; i++)
        {
            *(mem+i) = d[i];
        };
    };
};

//! Constructor to clone
template <typename type>
ram_buffer<type>::ram_buffer(const ram_buffer & input)
{
    this->class_name = "ram_buffer";
    clone_(input);                // call the virtual function clone_
};

//! Destructor
template <typename type>
ram_buffer<type>::~ram_buffer()
{
    free_memory();
};

template <typename type>
void ram_buffer<type>::allocate_memory()
{
    // if (_size_ > 0) mem = (type *)malloc(_size_*sizeof(type));
    if (_size_ > 0) mem = new type[_size_];
    else mem = nullptr;
};

template <typename type>
void ram_buffer<type>::free_memory()
{
    // if (_size_ > 0) free(mem);
    if (_size_ > 0) delete[] mem;
};

// ===========================================
// Create Functions
// ===========================================
template <typename type>
void ram_buffer<type>::clone_(const ram_buffer & input)
{
    mimic_(input);
    //full copy here
};

template <typename type>
void ram_buffer<type>::copy_(const ram_buffer & input)
{
    mimic_(input);
};

template <typename type>
void ram_buffer<type>::mimic_(const ram_buffer & input)
{
    free_memory();
    _size_ = input.size();
    allocate_memory();
};

// ===========================================
// Get Functions
// ===========================================
template <typename type>
int ram_buffer<type>::size() const
{
    return (int)_size_;
};

template <typename type>
type * ram_buffer<type>::data() const
{
    return mem;
};

// ===========================================
// Print Functions
// ===========================================
template <typename type>
std::string ram_buffer<type>::info(std::string msg)
{
    std::stringstream ss;
    std::string title = "RAM Buffer Information";
    if (msg != "") { title = msg; };

    // Summary of the object information
    ss << object::info(title);
    ss << "Size: \t\t\t" << this->size() << std::endl;
    return ss.str();
};

// ===========================================
// Functions
// ===========================================
template <typename type>
void ram_buffer<type>::resize(size_t sz)
{
    if(_size_ == 0)
    {
        _size_ = sz;
        allocate_memory();
    }
    else
    {
        free_memory();
        _size_ = sz;
        allocate_memory();
    };
};

template <typename type>
void ram_buffer<type>::clear()
{
    free_memory();
    _size_ = 0;
    mem = nullptr;
};

// ===========================================
// Overloading Functions
// ===========================================
template <typename type>
ram_buffer<type> & ram_buffer<type>::operator = (const ram_buffer & input)
{
    copy_(input);
    return *this;
};

template <typename type>
type & ram_buffer<type>::operator [] (size_t idx)
{   
    // if (_size_ > 0) return mem[idx];
    return mem[idx];
};

}; //end namespace

#endif