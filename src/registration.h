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
#include "ssd.h"
#include "optimizer.h"
#include "resolution.h"
#include "resolution.h"
#include "utils/timer.h"

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
    bool normalize;
    bool padding;
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
    typename resolution<type,container>::pointer multiresolution;

    // ===========================================
    // Functions
    // ===========================================
    void init();
    void init(typename image<type,container>::pointer fixed_image, 
              typename image<type,container>::pointer moving_image,
              typename transform<type,container>::pointer transformd);

public:
    // ===========================================
    // Create Functions
    // ===========================================
    // Constructors
    registration();
    registration( typename image<type,container>::pointer fixed_image, 
                  typename image<type,container>::pointer moving_image,
                  typename transform<type,container>::pointer transformd );
    
    // ===========================================
    // Get Functions
    // ===========================================
    bool get_normalize() const;
    bool get_padding() const;

    int get_levels() const;
    std::vector<int> get_levels_scales() const;
    std::vector<int> get_levels_sigmas() const;
    std::vector<int> get_levels_iterations() const;

    typename metric<type,container>::pointer get_metric() const;
    typename optimizer<type,container>::pointer get_optimizer() const;
    
    // ===========================================
    // Set Functions
    // ===========================================
    void set_normalize(bool norm);
    void set_padding(bool p);

    void set_levels(int lvls);
    void set_levels_scales(std::vector<int> scales);
    void set_levels_sigmas(std::vector<int> sigmas);
    void set_levels_iterations(std::vector<int> iterations);

    void set_fixed(typename image<type,container>::pointer fixed_image);
    void set_moving(typename image<type,container>::pointer moving_image);
    void set_transform(typename transform<type,container>::pointer transformd);
    void set_interpolator(typename interpolator<type,container>::pointer interpolationd);
    void set_metric(typename metric<type,container>::pointer methode);
    void set_optimizer(typename optimizer<type,container>::pointer optim);

    // ===========================================
    // Print Functions
    // ===========================================
    std::string info(std::string msg);

    // ===========================================
    // Functions
    // ===========================================
    // optimization
    virtual void apply();
    void operator() ();
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
registration<type,container>::registration( typename image<type,container>::pointer fixed_image, 
                                            typename image<type,container>::pointer moving_image,
                                            typename transform<type,container>::pointer transformd)
{
    this->class_name = "registration";
    init(fixed_image, moving_image, transformd);
};

template <typename type, typename container>
void registration<type,container>::init()
{
    fixed = image<type,container>::new_pointer();
    moving = image<type,container>::new_pointer();
    transformation = transform<type,container>::new_pointer();
    init(fixed, moving, transformation);
};

template <typename type, typename container>
void registration<type,container>::init(typename image<type,container>::pointer fixed_image, 
                                        typename image<type,container>::pointer moving_image,
                                        typename transform<type,container>::pointer transformd)
{
    // std::cout << "registration init" << std::endl;
    assert(fixed_image->get_dimension() == moving_image->get_dimension());
    assert(fixed_image->get_dimension() == transformd->get_dimension());

    // Initilize control variables
    normalize = true;
    padding = false;

    levels = 3;
    levels_scales = std::vector<int>({4,2,1});
    levels_sigmas = std::vector<int>({2,2,1});
    levels_iterations = std::vector<int>({150,100,50});
    
    fixed = fixed_image;
    moving = moving_image;
    transformation = transformd;
    interpolation = ilinear<type,container>::new_pointer(moving->get_dimension());
    method = ssd<type,container>::new_pointer(fixed,moving,transformation);
    optimization = gradient_descent<type,container>::new_pointer();
    multiresolution = resolution<type,container>::new_pointer(fixed);

    method->set_interpolator(interpolation);
    
    this->set_total_inputs(5);      //process_object::init
    this->set_total_outputs(0);     //process_object::init
    this->setup_input(fixed, moving, transformation, method, optimization);
    this->setup_output();
};

// ===========================================
// Get Functions
// ===========================================
template <typename type, typename container>
bool registration<type,container>::get_normalize() const
{
    return normalize;
};

template <typename type, typename container>
bool registration<type,container>::get_padding() const
{
    return padding;
};

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
void registration<type,container>::set_normalize(bool norm)
{
    normalize = norm;
};

template <typename type, typename container>
void registration<type,container>::set_padding(bool p)
{
    padding = p;
};

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

template <typename type, typename container>
void registration<type,container>::set_fixed(typename image<type,container>::pointer fixed_image)
{
    fixed = fixed_image;
    method->set_fixed(fixed);
};

template <typename type, typename container>
void registration<type,container>::set_moving(typename image<type,container>::pointer moving_image)
{
    moving = moving_image;
    method->set_moving(moving);
};

template <typename type, typename container>
void registration<type,container>::set_transform(typename transform<type,container>::pointer transformd)
{
    transformation = transformd;
    method->set_transform(transformation);
};

template <typename type, typename container>
void registration<type,container>::set_interpolator(typename interpolator<type,container>::pointer interpolationd)
{
    interpolation = interpolationd;
    method->set_interpolator(interpolationd);
};

template <typename type, typename container>
void registration<type,container>::set_metric(typename metric<type,container>::pointer methode)
{
    method = methode;
    method->set_fixed(fixed);
    method->set_moving(moving);
    method->set_transform(transformation);
    method->set_interpolator(interpolation);
};

template <typename type, typename container>
void registration<type,container>::set_optimizer(typename optimizer<type,container>::pointer optim)
{
    optimization = optim;
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
void registration<type,container>::apply()
{
    std::cout << "============================" << std::endl;
    std::cout << "Multiresolution Registration" << std::endl;
    std::cout << "============================" << std::endl;

    timer t("ms");
    t.start();

    auto tmp_fixed = image<type,container>::new_pointer(fixed->get_dimension());
    auto tmp_moving = image<type,container>::new_pointer(moving->get_dimension());
    auto tmp_transform = transformation->copy();

    // init transformation 
    if(transformation->get_name() == "dfield")
    {
        tmp_transform = multiresolution->apply(levels_scales[0], method->get_transform());
        // tmp_transform->print("tmp transform");
        method->set_transform(tmp_transform);
    };

    // multiresolution
    for(int i = 0; i < levels; i++)
    {
        // method->print("Demons");

        // std::cout << "Fixed scale" << std::endl;
        multiresolution->set_image(fixed);
        tmp_fixed = multiresolution->apply(levels_scales[i]);

        // std::cout << "Moving scale" << std::endl;
        multiresolution->set_image(moving);
        tmp_moving = multiresolution->apply(levels_scales[i]);

        method->set_fixed(tmp_fixed);
        method->set_moving(tmp_moving);

        // method->print("Demons");
        // tmp_fixed->print("Fixed");
        // tmp_moving->print("Moving");
        // method->get_transform()->print("Transform");
        // method->get_transform()->get_parameters(0)->print("x");
        // method->get_transform()->get_parameters(1)->print("y");

        // Info
        std::cout << "Level:\t\t" << i << std::endl;
        std::cout << "Scale:\t\t" << levels_scales[i] << std::endl;
        std::cout << "Iters:\t\t" << levels_iterations[i] << std::endl;
        std::cout << "Image:\t\t[ ";
        for(int k = 0; k < tmp_fixed->get_size().size(); k++) { std::cout << tmp_fixed->get_size()[k] << " "; };
        std::cout << "]" << std::endl;

        // optimization
        optimization->set_iterations(levels_iterations[i]);
        optimization->optimize(method);

        // increase resolution of transform;
        if (i < levels - 1 && transformation->get_name() == "dfield")
        {
            tmp_transform = multiresolution->apply(double(levels_scales[i+1])/double(levels_scales[i]), method->get_transform());
            method->set_transform(tmp_transform);
        };
        // levels_sigmas[i]
    };

    // set_transform;
    set_transform(tmp_transform);

    t.lap();
    std::cout << "Total registration time: " << std::setw(4) << t.get_elapsed() << std::endl;
};

template <typename type, typename container>
void registration<type,container>::operator() ()
{
    apply(method);
};

}; //end namespace

#endif