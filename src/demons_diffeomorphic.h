/*
* @Author: jose
* @Date:   2019-12-09 00:00:00
* @Last Modified by:   jose
* @Last Modified time: 2020-03-18 00:00:00
*/

#ifndef __DEMONS_DIFFEOMORPHIC_H__
#define __DEMONS_DIFFEOMORPHIC_H__

#include "metric.h"

namespace imart
{

// Class metric
template <typename type, typename container=vector_cpu<type>>
class demons_diffeomorphic: public inherit<demons_diffeomorphic<type,container>, metric<type,container>>
{
public:
    //Type definitions
    using self    = demons_diffeomorphic;
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

    using inherit<demons_diffeomorphic<type,container>, metric<type,container>>::inherit;

    double weight;
    double step;
    void init_w();

public:
    // ===========================================
    // Create Functions
    // ===========================================
    // Constructors
    demons_diffeomorphic() : inherit<demons_diffeomorphic<type,container>, metric<type,container>>()
          { this->class_name = "demons_diffeomorphic"; init_w(); };

    demons_diffeomorphic(int d) : inherit<demons_diffeomorphic<type,container>, metric<type,container>>(d)
               { this->class_name = "demons_diffeomorphic"; init_w(); };

    demons_diffeomorphic(typename image<type,container>::pointer fixed_image, 
        typename image<type,container>::pointer moving_image)
        : inherit<demons_diffeomorphic<type,container>, metric<type,container>>(fixed_image, moving_image)
        { this->class_name = "demons_diffeomorphic"; init_w(); };

    demons_diffeomorphic(typename image<type,container>::pointer fixed_image, 
        typename image<type,container>::pointer moving_image,
        typename transform<type,container>::pointer transformd)
        : inherit<demons_diffeomorphic<type,container>, metric<type,container>>(fixed_image, moving_image, transformd)
        { this->class_name = "demons_diffeomorphic"; init_w(); };

    void set_step(double s);

    // ===========================================
    // Functions
    // ===========================================
    // !compute the cost
    type cost();
    // !calculate derivative
    typename transform<type,container>::pointer derivative();
    // !regularize
    // void regularize();
    // // !max scale
    void max_scale(typename transform<type,container>::pointer trfm);

};

template<typename type>
using demons_diffeomorphic_cpu = demons_diffeomorphic<type,vector_cpu<type>>;

template<typename type>
using demons_diffeomorphic_gpu = demons_diffeomorphic<type,vector_ocl<type>>;


// ===========================================
//      Functions of Class metric
// ===========================================

// ===========================================
// Functions
// ===========================================
template <typename type, typename container>
void demons_diffeomorphic<type,container>::init_w()
{
    step = 1.0;
    weight = 1.0;
};

template <typename type, typename container>
void demons_diffeomorphic<type,container>::set_step(double s)
{
    step = s;
};

// ===========================================
// Functions
// ===========================================
template <typename type, typename container>
type demons_diffeomorphic<type,container>::cost()
{
    // std::cout << "init cost" << std::endl;
    this->warped_moving(); // update moving_prime
    // auto moving_prime = this->warped_moving();
    // moving_prime->print_data("moving prime");
    // std::cout << "dif cost" << std::endl;
    // fixed->print("dem fix");
    // moving_prime->print("dem mov");
    type N = (type)fixed->get_total_elements();

    // e_reg
    image<type,container> ssd_ = *fixed - *moving_prime;
    type e_sim = (1/N)*( (ssd_*ssd_).sum() );
    // e_sim
    auto jac = jacobian(transformation->get_parameters_vector()); // vector field
    type e_reg = ( ((*jac)*(*jac)).sum() )/N;
    // cost
    cost_value = e_sim + ((weight*weight)/(step*step)) * e_reg;
    
    // std::cout << "end cost" << std::endl;
    return cost_value;
};

template <typename type, typename container>
typename transform<type,container>::pointer demons_diffeomorphic<type,container>::derivative()
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
    
    type low = 1e-8;    // value to avoid 0 division

    auto dif = (*moving_prime - *fixed);
    auto dif_sq = dif*dif;

    if (transformation->get_name() == "dfield")
    {
        // std::cout << "der 1" << std::endl;
        // Computing gradient
        auto grad = gradient(moving_prime);

        // dimensional data
        auto gx = *grad[0];
        auto gy = *grad[1];
        
        if (d == 2)
        {
            // compute mag
            auto mag_sq = (gx*gx + gy*gy);

            // write parameters
            // auto den = ( mag_sq + dif_sq*((weight*weight)/(step*step)) );
            auto scale = dif/( mag_sq + dif_sq*((weight*weight)/(step*step)) + low );
            scale.replace(mag_sq <= low, 0.0);
            scale.replace(dif_sq <= low, 0.0);
            // dif.replace(den <= 1e-10, 0);
            // den.replace(den <= 1e-10, 1);
            // printf("den {min: %e, max: %e}\n",den.min(), den.max());
            // auto scale = dif/den;

            auto ux = scale*gx*fixed->get_spacing()[0];
            auto uy = scale*gy*fixed->get_spacing()[1];

            ux.replace(*fixed <= low, 0.0);
            ux.replace(*moving_prime <= low, 0.0);
            uy.replace(*fixed <= low, 0.0);
            uy.replace(*moving_prime <= low, 0.0);
            // printf("ux {min: %e, max: %e}\n",ux.min(), ux.max());
            // printf("uy {min: %e, max: %e}\n",uy.min(), uy.max());


            *((*param)[0]) = ux;
            *((*param)[1]) = uy;
            // *((*param)[0]) = scale*gx;
            // *((*param)[1]) = scale*gy;
        }
        else if (d == 3)
        {
            // z
            auto gz = *grad[2];
            // compute mag
            auto mag_sq = (gx*gx + gy*gy + gz*gz);
            
            // write parameters
            auto scale = dif/( mag_sq + dif_sq*((weight*weight)/(step*step)) );
            scale.replace(mag_sq <= 0,0);
            scale.replace(dif_sq <= 0,0);
            *((*param)[0]) = scale*gx;
            *((*param)[1]) = scale*gy;
            *((*param)[2]) = scale*gz;
        };
    };

    // regularize update (gaussian filter with sigma fluid)

    // type sigma = 2.0;
    // (*param)[0] = gaussian_filter((*param)[0],sigma,3);
    // (*param)[1] = gaussian_filter((*param)[1],sigma,3);


    trfm->set_parameters_vector( param );
    return trfm;
};


template <typename type, typename container>
void demons_diffeomorphic<type,container>::max_scale(typename transform<type,container>::pointer trfm)
{
    type pixel = 0.4;
    type vmax = 0.0;
    auto param = trfm->get_parameters_vector();

    for (size_t i = 0; i < param->size(); i ++)
    {
        type pmax = ((*param)[i])->max();
        if (pmax > vmax) vmax = pmax;
    };
    printf("vmax = %f\n", vmax);
    for (size_t i = 0; i < param->size(); i ++)
    {
        *((*param)[i]) = (*((*param)[i]))*(pixel/vmax);
    };
};

}; //end namespace

#endif