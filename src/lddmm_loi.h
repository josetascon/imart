/*
* @Author: jose
* @Date:   2019-12-09 00:00:00
* @Last Modified by:   jose
* @Last Modified time: 2020-03-18 00:00:00
*/

#ifndef __LDDMM_LOI_H__
#define __LDDMM_LOI_H__

#include "metric.h"

namespace imart
{

// Class demons local orderless information
template <typename type, typename container=vector_cpu<type>>
class lddmm_loi: public inherit<lddmm_loi<type,container>, lddmm<type,container>>
{
public:
    //Type definitions
    using self    = lddmm_loi;
    using pointer = std::shared_ptr<self>;
    using vector  = std::vector<self::pointer>;

    // Inherited variables
    using metric<type,container>::cost_value; //Parent class properties to be used here (avoid use "this->")
    using pairwise_object<type,container>::fixed;
    using pairwise_object<type,container>::moving;
    using pairwise_object<type,container>::transformation;
    using metric<type,container>::x0;
    using metric<type,container>::x1;
    using metric<type,container>::cost;
    using metric<type,container>::fixed_prime;
    using metric<type,container>::moving_prime;

    using inherit<lddmm_loi<type,container>, lddmm<type,container>>::inherit;

    double alpha;
    double sigma_gradient;

    typename lddmm<type,container>::pointer lddmm_gradients;
    typename image<type,container>::pointer fixed_dir_gradient;
    typename image<type,container>::pointer moving_dir_gradient;



protected:
    void init(){ lddmm_gradients};
    void init_vars(){ alpha = 0.5; sigma_gradient = 2.0;};

public:
    // ===========================================
    // Create Functions
    // ===========================================
    // Constructors
    lddmm_loi() : inherit<lddmm_loi<type,container>, lddmm<type,container>>()
          { this->class_name = "lddmm_loi"; };//init_sigma(); };

    lddmm_loi(int d) : inherit<lddmm_loi<type,container>, lddmm<type,container>>(d)
               { this->class_name = "lddmm_loi"; };//init_sigma(); };

    lddmm_loi(typename image<type,container>::pointer fixed_image, 
        typename image<type,container>::pointer moving_image)
        : inherit<lddmm_loi<type,container>, lddmm<type,container>>(fixed_image, moving_image)
        { this->class_name = "lddmm_loi"; init_vars(); };//init_sigma(); };

    lddmm_loi(typename image<type,container>::pointer fixed_image, 
        typename image<type,container>::pointer moving_image,
        typename transform<type,container>::pointer transformd)
        : inherit<lddmm_loi<type,container>, lddmm<type,container>>(fixed_image, moving_image, transformd)
        { this->class_name = "lddmm_loi"; init_vars(); };//init_sigma(); };

    // ===========================================
    // Get Functions
    // ===========================================
    double get_alpha() const;
    double get_sigma_gradient() const;

    // ===========================================
    // Set Functions
    // ===========================================
    void set_alpha(double a);
    void set_sigma_gradient(double s);

    // ===========================================
    // Functions
    // ===========================================
    // !compute the cost
    // type cost();
    // !calculate derivative
    typename transform<type,container>::pointer derivative();
    // !regularize
    // void regularize();
    // // !max scale
    // void max_scale();

    typename image<type,container>::pointer moving_directional_gradient();
    typename image<type,container>::pointer fixed_directional_gradient();
};

template<typename type>
using lddmm_loi_cpu = lddmm_loi<type,vector_cpu<type>>;

// template<typename type>
// using lddmm_loi_gpu = lddmm_loi<type,vector_ocl<type>>;

#ifdef IMART_WITH_OPENCL
template<typename type>
using lddmm_loi_ocl = lddmm_loi<type,vector_ocl<type>>;
#endif

#ifdef IMART_WITH_CUDA
template<typename type>
using lddmm_loi_cuda = lddmm_loi<type,vector_cuda<type>>;
#endif

// ===========================================
//      Functions of Class lddmm_loi
// ===========================================

// ===========================================
// Get Functions
// ===========================================
template <typename type, typename container>
double lddmm_loi<type,container>::get_alpha() const
{
    return alpha;
};

template <typename type, typename container>
double lddmm_loi<type,container>::get_sigma_gradient() const
{
    return sigma_gradient;
};

// ===========================================
// Set Functions
// ===========================================
template <typename type, typename container>
void lddmm_loi<type,container>::set_alpha(double a)
{
    alpha = a;
};

template <typename type, typename container>
void lddmm_loi<type,container>::set_sigma_gradient(double s)
{
    sigma_gradient = s;
};

// ===========================================
// Functions
// ===========================================
template <typename type, typename container>
typename transform<type,container>::pointer lddmm_loi<type,container>::derivative()
{
    int d = fixed->get_dimension();
    // transform
    typename transform<type,container>::pointer trfm = transformation->mimic();
    
    // mimic parameters
    auto param = std::make_shared< std::vector< typename image<type,container>::pointer >>(d);
    for(int i = 0; i < d; i++)
    {
        (*param)[i] = image<type,container>::new_pointer();
        (*param)[i]->mimic_(*(transformation->get_parameters(i)));
    };
    // moving prime already computed
    // auto moving_prime = this->warped_moving();        // consider store this to avoid compute again

    // parameters derivative
    auto param_der = std::make_shared< std::vector< typename image<type,container>::pointer >>(d);
    for(int i = 0; i < d; i++)
    {
        (*param_der)[i] = image<type,container>::new_pointer();
        (*param_der)[i]->mimic_(*(transformation->get_parameters(i)));
    };
    
    if (transformation->get_name() == "dfield")
    {
        // this->optical_flow(fixed, moving_prime, param);

        // moving_directional_gradient();
        // fixed_directional_gradient();

        // if (alpha > 0.0)
        // {
        //     this->optical_flow(fixed_dir_gradient, moving_dir_gradient, param_der);
        // };
    };

    if (alpha > 0.0)
    {
        *param      = (*param)*(1.0-alpha);
        *param_der  = (*param_der)*alpha;
        *param      = (*param) + (*param_der);
    };

    trfm->set_parameters_vector( param );
    return trfm;
};

template <typename type, typename container>
typename image<type,container>::pointer lddmm_loi<type,container>::moving_directional_gradient()
{
    moving_dir_gradient = gradient_direction_gaussian( moving_prime, sigma_gradient );
    return moving_dir_gradient;
};

template <typename type, typename container>
typename image<type,container>::pointer lddmm_loi<type,container>::fixed_directional_gradient()
{
    fixed_dir_gradient = gradient_direction_gaussian( fixed, sigma_gradient );
    return fixed_dir_gradient;
};

}; //end namespace

#endif