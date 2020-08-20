/*
* @Author: jose
* @Date:   2020-08-20 00:00:00
* @Last Modified by:   jose
* @Last Modified time: 2020-08-20 00:00:00
*/

#ifndef __REGISTRATION_H__
#define __REGISTRATION_H__

#include "pairwise_object.h"
#include "metric.h"
#include "optimizer.h"

namespace imart
{

template <typename type, typename container=vector_cpu<type>>
class registration: public inherit<registration<type,container>, pairwise_object<type,container>>
{
public:
    //Type definitions
    using self    = registration;
    using pointer = std::shared_ptr<self>;
    using vector  = std::vector<self::pointer>;

protected:
    // ===========================================
    // Internal Variables
    // ===========================================
    int levels;
    std::vector<int> levels_scales;
    std::vector<int> levels_sigmas;
    std::vector<int> levels_iterations;

    using pairwise_object<type,container>::fixed;
    using pairwise_object<type,container>::moving;
    using pairwise_object<type,container>::transformation;
    using pairwise_object<type,container>::interpolation;
    typename metric<type,container>::pointer method;
    typename optimizer<type,container>::pointer optimization;

    // ===========================================
    // Functions
    // ===========================================
    void init();

public:
    // ===========================================
    // Create Functions
    // ===========================================
    // Constructors
    registration();
    
    // ===========================================
    // Get Functions
    // ===========================================
    int get_levels() const;
    std::vector<int> get_levels_scales() const;
    std::vector<int> get_levels_sigmas() const;
    std::vector<int> get_levels_iterations() const;
    
    // ===========================================
    // Set Functions
    // ===========================================
    void set_levels(int lvls);
    void set_levels_scales(std::vector<int> scales);
    void set_levels_sigmas(std::vector<int> sigmas);
    void set_levels_iterations(std::vector<int> iterations);

    // ===========================================
    // Print Functions
    // ===========================================
    std::string info(std::string msg);

    // ===========================================
    // Functions
    // ===========================================
    // optimization
    virtual void apply(typename metric<type,container>::pointer method);
};


// ===========================================
//      Functions of Class metric
// ===========================================

// ===========================================
// Create Functions
// ===========================================
// Constructor
template <typename type, typename container>
registration<type,container>::registration()
{
    this->class_name = "registration";
    init();
};

template <typename type, typename container>
void registration<type,container>::init()
{
    // Initilize control variables
    
};

// ===========================================
// Get Functions
// ===========================================
template <typename type, typename container>
int registration<type,container>::get_levels() const
{
    return levels;
};

template <typename type, typename container>
std::vector<int> registration<type,container>::get_levels_scales() const
{
    return levels_scales;
};

template <typename type, typename container>
std::vector<int> registration<type,container>::get_levels_sigmas() const
{
    return levels_sigmas;
};

template <typename type, typename container>
std::vector<int> registration<type,container>::get_levels_iterations() const
{
    return levels_iterations;
};

// ===========================================
// Set Functions
// ===========================================
template <typename type, typename container>
void registration<type,container>::set_levels(int lvls)
{
    levels = lvls;
};

template <typename type, typename container>
void registration<type,container>::set_levels_scales(std::vector<int> scales)
{
    assert(levels == scales.size());
    levels_scales = scales;
};

template <typename type, typename container>
void registration<type,container>::set_levels_sigmas(std::vector<int> sigmas)
{
    assert(levels == sigmas.size());
    levels_sigmas = sigmas;
};

template <typename type, typename container>
void registration<type,container>::set_levels_iterations(std::vector<int> iterations)
{
    assert(levels == iterations.size());
    levels_iterations = iterations;
};

// ===========================================
// Print Functions
// ===========================================
template <typename type, typename container>
std::string registration<type,container>::info(std::string msg)
{
    std::stringstream ss;
    std::string title = "Registration Information";
    if (msg != "") { title = msg; };

    ss << object::info(title);
    ss << process_object::info("");
    return ss.str();
};

// ===========================================
// Functions
// ===========================================
template <typename type, typename container>
void registration<type,container>::apply(typename metric<type,container>::pointer method)
{
    ;
};

}; //end namespace

#endif