/*
* @Author: jose
* @Date:   2020-01-27 00:00:00
* @Last Modified by:   jose
* @Last Modified time: 2020-01-27 00:00:00
*/

#ifndef __OPTIMIZER_H__
#define __OPTIMIZER_H__

#include "image_base.h"
#include "transform_base.h"

namespace imart
{

template <typename pixel_type>
class optimizer: public object<pixel_type>
{
public:
    //Type definitions
    ;
protected:
    // ===========================================
    // Internal Variables
    // ===========================================
    int iteration;
    int num_iterations;
    
    int unchanged_times;
    int unchanged_max_times;
    
    double tolerance;
    double previous_cost;
    double current_cost;

    pixel_type step;

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

    // ===========================================
    // Print Functions
    // ===========================================
    virtual std::string info(std::string msg);
    virtual std::string info_data(std::string msg);

    // ===========================================
    // Functions
    // ===========================================
    // compute the cost
    void optimize(metric<pixel_type> & method);

};






// ===========================================
//      Functions of Class metric
// ===========================================


// ===========================================
// Create Functions
// ===========================================
// Constructor
template <typename pixel_type>
optimizer<pixel_type>::optimizer()
{
    init();
};

template <typename pixel_type>
void optimizer<pixel_type>::init()
{
    // Initilize control variables
    iteration = 0;
    num_iterations = 100;
    step = 1.0;
    
    unchanged_times = 0;
    unchanged_max_times = 15;
    
    tolerance = 1e-4;
    previous_cost = 1e41;
    current_cost = 1e40;
    termination = "none";
};


template <typename pixel_type>
std::string optimizer<pixel_type>::info(std::string msg)
{
    std::stringstream ss;
    std::string title = "Optimizer Information";
    if (msg != "") { title = msg; };
    // Summary of the optimizer information
    ss << object<pixel_type>::info(title);

    return ss.str();
};

template <typename pixel_type>
std::string optimizer<pixel_type>::info_data(std::string msg)
{
    std::stringstream ss;
    if (msg != "") { ss << msg << std::endl; }
    else { ss << "Optimizer parameters:\n"; };

    return ss.str();
};

template <typename pixel_type>
void optimizer<pixel_type>::optimize(metric<pixel_type> & method)
{
    std::cout << "Init optimization" << std::endl;

    typename transform_base<pixel_type>::pointer tr_base = method.get_transform();
    auto tr_derive = transform_base<pixel_type>::new_pointer();
    // typename transform_base<pixel_type>::pointer tr_derive;

    double diff = 0.0;

    auto scales = image_base<pixel_type>::new_pointer();
    scales->imitate(*tr_base->get_parameters());
    // *** CONTINUE here to include scales

    while( iteration < num_iterations )
    {
        // std::cout << "Cost:" << std::endl;
        current_cost = method.cost();
        // std::cout << "Derivative:" << std::endl;
        tr_derive = method.derivative();
        // tr_derive->print_data("gradient");
        // std::cout << "Update:" << std::endl;
        *tr_base = *tr_base - (*tr_derive)*step;
        // tr_base->print_data("transform update");

        diff = abs(previous_cost - current_cost);
        previous_cost = current_cost;

        if( (diff) < tolerance )
        {
            termination = "tolerance";
            std::cout << "Optimization terminated: tolerance" << std::endl;
            break;
        }

        std::cout << "iteration: " << iteration << "\tcost: " << current_cost << "\tdiff: " << diff << std::endl;

        iteration++;
    };

    if (iteration >= num_iterations) 
    {
        termination = "iterations";
        std::cout << "Optimization terminated: max iterations" << std::endl;
    };

};

}; //end namespace

#endif