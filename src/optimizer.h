/*
* @Author: jose
* @Date:   2020-01-27 00:00:00
* @Last Modified by:   jose
* @Last Modified time: 2020-01-27 00:00:00
*/

#ifndef __OPTIMIZER_H__
#define __OPTIMIZER_H__

#include "image.h"
#include "image_utils.h"
#include "transform.h"
#include "metric.h"
#include "utils/timer.h"

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

    typename metric<type,container>::pointer _method_;

    // ===========================================
    // Functions
    // ===========================================
    void init();
    void defaults();

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
    // std::cout << "init optimizer" << std::endl;
    defaults();
    tolerance = 1e-6;
    max_unchanged_times = 15;

    this->set_total_inputs(1);      //process_object::init
    this->set_total_outputs(0);     //process_object::init

    _method_ = metric<type,container>::new_pointer();
    this->setup_input(_method_);
    this->setup_output();           // output is transformation
    // std::cout << "end init optimizer" << std::endl;
};

template <typename type, typename container>
void optimizer<type,container>::defaults()
{
    // Initilize control variables
    iterations = 0;
    max_iterations = 300;
    unchanged_times = 0;
    
    // tolerance = 1e-5;
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
// Set Functions
// ===========================================
template <typename type, typename container>
void optimizer<type,container>::set_iterations(int i)
{   
    defaults(); // reset actual iterations
    max_iterations = i;
};

template <typename type, typename container>
void optimizer<type,container>::set_unchanged_times(int u)
{
    max_unchanged_times = u;
};

template <typename type, typename container>
void optimizer<type,container>::set_tolerance(double t)
{
    tolerance = t;
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
    _method_ = method;
};

}; //end namespace

#endif