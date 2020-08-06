/*
* @Author: jose
* @Date:   2020-01-27 00:00:00
* @Last Modified by:   jose
* @Last Modified time: 2020-01-27 00:00:00
*/

#ifndef __OPTIMIZER_H__
#define __OPTIMIZER_H__

#include "image.h"
#include "transform.h"
#include "metric.h"

namespace imart
{

template <typename type, typename container=vector_cpu<type>>
class optimizer: public inherit<optimizer<type,container>, process_object>
{
public:
    //Type definitions
    using self    = optimizer;
    using pointer = std::shared_ptr<self>;
    using vector  = std::vector<self::pointer>;

protected:
    // ===========================================
    // Internal Variables
    // ===========================================
    int iterations;
    int max_iterations;
    
    int unchanged_times;
    int max_unchanged_times;
    
    double tolerance;
    double current_cost;
    double previous_cost;
    double lowest_cost;

    std::string termination;

    // ===========================================
    // Functions
    // ===========================================
    void init();

public:
    // ===========================================
    // Create Functions
    // ===========================================
    // Constructors
    optimizer();
    
    // ===========================================
    // Get Functions
    // ===========================================
    int get_iterations() const;
    int get_unchanged_times() const;
    double get_tolerance() const;
    std::string get_termination() const;

    // ===========================================
    // Set Functions
    // ===========================================
    void set_iterations(int i);
    void set_unchanged_times(int u);
    void set_tolerance(double t);

    // ===========================================
    // Print Functions
    // ===========================================
    std::string info(std::string msg);

    // ===========================================
    // Functions
    // ===========================================
    // optimization
    virtual void optimize(typename metric<type,container>::pointer method);
};


// ===========================================
//      Functions of Class metric
// ===========================================

// ===========================================
// Create Functions
// ===========================================
// Constructor
template <typename type, typename container>
optimizer<type,container>::optimizer()
{
    this->class_name = "optimizer";
    init();
};

template <typename type, typename container>
void optimizer<type,container>::init()
{
    // Initilize control variables
    iterations = 0;
    max_iterations = 300;
    
    unchanged_times = 0;
    max_unchanged_times = 15;
    
    tolerance = 1e-6;
    previous_cost = 1e41;
    current_cost = 1e40;
    lowest_cost = 1e41;
    termination = "none";
};

// ===========================================
// Get Functions
// ===========================================
template <typename type, typename container>
int optimizer<type,container>::get_iterations() const
{
    return max_iterations;
};

template <typename type, typename container>
int optimizer<type,container>::get_unchanged_times() const
{
    return max_unchanged_times;
};

template <typename type, typename container>
double optimizer<type,container>::get_tolerance() const
{
    return tolerance;
};

template <typename type, typename container>
std::string optimizer<type,container>::get_termination() const
{
    return termination;
};

// ===========================================
// Print Functions
// ===========================================
template <typename type, typename container>
std::string optimizer<type,container>::info(std::string msg)
{
    std::stringstream ss;
    std::string title = "Optimizer Information";
    if (msg != "") { title = msg; };

    ss << object::info(title);
    ss << process_object::info("");
    return ss.str();
};

// ===========================================
// Functions
// ===========================================
template <typename type, typename container>
void optimizer<type,container>::optimize(typename metric<type,container>::pointer method)
{
    ;
};

}; //end namespace

#endif