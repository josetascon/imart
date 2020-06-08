/*
* @Author: jose
* @Date:   2020-01-27 00:00:00
* @Last Modified by:   jose
* @Last Modified time: 2020-01-27 00:00:00
*/

#ifndef __OPTIMIZER_H__
#define __OPTIMIZER_H__

#include "image.h"

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
    void optimize();

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
    ;
};


template <typename pixel_type>
std::string optimizer::<pixel_type>::info(std::string msg)
{
    std::stringstream ss;
    std::string title = "Optimizer Information";
    if (msg != "") { title = msg; };
    // Summary of the optimizer information
    ss << object<pixel_type>::info(title);
};

}; //end namespace

#endif