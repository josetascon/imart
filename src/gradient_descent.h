/*
* @Author: jose
* @Date:   2020-01-27 00:00:00
* @Last Modified by:   jose
* @Last Modified time: 2020-01-27 00:00:00
*/

#ifndef __GRADIENT_DESCENT_H__
#define __GRADIENT_DESCENT_H__

#include "image.h"
#include "transform.h"
#include "metric.h"
#include "optimizer.h"
#include "utils/timer.h"

namespace imart
{

template <typename type, typename container=vector_cpu<type>>
class gradient_descent: public inherit<gradient_descent<type,container>, optimizer<type,container>>
{
public:
    //Type definitions
    using self    = gradient_descent;
    using pointer = std::shared_ptr<self>;
    using vector  = std::vector<self::pointer>;

    // Inherited variables
    using optimizer<type,container>::iterations;
    using optimizer<type,container>::max_iterations;
    using optimizer<type,container>::unchanged_times;
    using optimizer<type,container>::max_unchanged_times;
    using optimizer<type,container>::tolerance;
    using optimizer<type,container>::current_cost;
    using optimizer<type,container>::previous_cost;
    using optimizer<type,container>::lowest_cost;
    using optimizer<type,container>::termination;
    using optimizer<type,container>::_method_;

    using optimizer<type,container>::init;

protected:
    type step;

public:
    // ===========================================
    // Create Functions
    // ===========================================
    // Constructors
    gradient_descent();

    // ===========================================
    // Get Functions
    // ===========================================
    type get_step();

    // ===========================================
    // Get Functions
    // ===========================================
    void set_step(type s);

    // ===========================================
    // Functions
    // ===========================================
    // optimization
    void optimize(typename metric<type,container>::pointer method);
};


// ===========================================
//      Functions of Class metric
// ===========================================

// ===========================================
// Create Functions
// ===========================================
// Constructor
template <typename type, typename container>
gradient_descent<type,container>::gradient_descent()
{
    this->class_name = "gradient_descent";
    init();
    step = 1.0;
};

// ===========================================
// Get Functions
// ===========================================
template <typename type, typename container>
type gradient_descent<type,container>::get_step()
{
    return step;
};

// ===========================================
// Set Functions
// ===========================================
template <typename type, typename container>
void gradient_descent<type,container>::set_step(type s)
{
    step = s;
};

// ===========================================
// Functions
// ===========================================
template <typename type, typename container>
void gradient_descent<type,container>::optimize(typename metric<type,container>::pointer method)
{
    _method_ = method;
    std::cout << "Init optimization" << std::endl;
    // Affine
    // typename affine<type,container>::pointer tr_base = std::static_pointer_cast<affine<type,container>>(method->get_transform());
    // typename affine<type,container>::pointer tr_derive = tr_base->mimic();
    // typename affine<type,container>::pointer tr_best = tr_base->mimic();

    typename transform<type,container>::pointer tr_base = method->get_transform();
    typename transform<type,container>::pointer tr_derive = tr_base->mimic();
    typename transform<type,container>::pointer tr_best = tr_base->mimic();


    // tr_base->print("transform base");
    // tr_base->print_data("transform base");
    // tr_derive->print("transform derive");
    // tr_derive->print_data("transform derive");
    // *tr_base = *tr_base - *tr_derive;

    int d = method->get_fixed()->get_dimension();
    termination = "max iterations";
    double diff = 0.0;

    // Affine transform
    // typename transform<type,container>::pointer scales = tr_base->mimic();
    // std::vector<type> p(d*d + d);
    // if (d == 2)
    // {
    //     p[0] = 0.12; p[1] = 0.12; p[2] = 0.12; p[3] = 0.12;
    //     p[4] = 2000; p[5] = 2000;
    // }
    // else if (d == 3)
    // {
    //     p[0] = 0.12; p[1] = 0.12; p[2] = 0.12; 
    //     p[3] = 0.12; p[4] = 0.12; p[5] = 0.12; 
    //     p[6] = 0.12; p[7] = 0.12; p[8] = 0.12;
    //     p[9] = 2000; p[10] = 2000; p[11] = 2000;
    // }
    // scales->get_parameters()->get_data()->read_ram(p.data(),p.size());

    timer t("ms");
    t.start();

    while( iterations < max_iterations )
    {
        // std::cout << "Cost:" << std::endl;
        current_cost = method->cost();

        // auto output = image<unsigned char>::new_pointer();
        // auto tmp = method->warped_moving();
        // cast((*tmp)*(type(255)),*output);
        // output->write("./out" + std::to_string(iterations) + ".png");

        
        // std::cout << "Derivative:" << std::endl;
        tr_derive = method->derivative();
        // tr_derive->print_data("gradient");
        
        // std::cout << "Update:" << std::endl;
        // *tr_base = *tr_base - (*tr_derive)*(*scales)*step;
        *tr_base = *tr_base - (*tr_derive)*step;
        // tr_base->print_data("transform update");

        type sigma = 2.0;
        tr_base->set_parameters( gaussian_filter(tr_base->get_parameters(0),sigma,3), 0 );
        tr_base->set_parameters( gaussian_filter(tr_base->get_parameters(1),sigma,3), 1 );

        diff = abs(previous_cost - current_cost);
        t.lap();
        std::cout << std::showpoint << std::setprecision(4);
        std::cout << "iteration: " << std::setw(4) << iterations << "  cost: " << std::setw(4) <<  current_cost;
        std::cout << "  diff: " << std::setw(4) << diff << "  time:" << std::setw(4) << t.get_elapsed() << std::endl;
        iterations++;

        if (diff > 0.0 and diff < tolerance)
        {
            lowest_cost = current_cost;
            tr_best = tr_base->copy();
            termination = "convergence";
            break;
        };

        if (current_cost < lowest_cost)
        {
            lowest_cost = current_cost;
            tr_best = tr_base->copy();
            unchanged_times = 0;
        }
        else unchanged_times += 1;

        if (unchanged_times >= max_unchanged_times)
        {
            termination = "max unchaged";
            break;
        };

        previous_cost = current_cost;
    };

    t.finish();
    std::cout << "Optimization Summary" << std::endl;
    std::cout << "Iterations:\t" << iterations << std::endl;
    std::cout << "Lowest cost:\t" << lowest_cost << std::endl;
    std::cout << "Termination:\t" << termination << std::endl;
    // tr_best->print_data();
    method->set_transform(tr_best);

};

}; //end namespace

#endif