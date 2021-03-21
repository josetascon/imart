/*
* @Author: jose
* @Date:   2019-12-09 00:00:00
* @Last Modified by:   jose
* @Last Modified time: 2020-03-18 00:00:00
*/

#ifndef __DEMONS_H__
#define __DEMONS_H__

#include "metric.h"

namespace imart
{

// Class metric
template <typename type, typename container=vector_cpu<type>>
class demons: public inherit<demons<type,container>, metric<type,container>>
{
public:
    //Type definitions
    using self    = demons;
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

    using inherit<demons<type,container>, metric<type,container>>::inherit;

    // double sigma;
    // int kernel_width;
    // void init_sigma();

protected:
    void optical_flow(std::shared_ptr<image<type,container>> img_fixed, 
                      std::shared_ptr<image<type,container>> img_moving, 
                      std::shared_ptr< std::vector< std::shared_ptr < image<type,container> >>> param );

public:
    // ===========================================
    // Create Functions
    // ===========================================
    // Constructors
    demons() : inherit<demons<type,container>, metric<type,container>>()
          { this->class_name = "demons"; };//init_sigma(); };

    demons(int d) : inherit<demons<type,container>, metric<type,container>>(d)
               { this->class_name = "demons"; };//init_sigma(); };

    demons(typename image<type,container>::pointer fixed_image, 
        typename image<type,container>::pointer moving_image)
        : inherit<demons<type,container>, metric<type,container>>(fixed_image, moving_image)
        { this->class_name = "demons"; };//init_sigma(); };

    demons(typename image<type,container>::pointer fixed_image, 
        typename image<type,container>::pointer moving_image,
        typename transform<type,container>::pointer transformd)
        : inherit<demons<type,container>, metric<type,container>>(fixed_image, moving_image, transformd)
        { this->class_name = "demons"; };//init_sigma(); };

    // ===========================================
    // Get Functions
    // ===========================================
    // double get_sigma() const;
    // int get_kernel_width() const;

    // ===========================================
    // Set Functions
    // ===========================================
    // void set_sigma(double s);
    // void set_kernel_width(int k); // automatically assign based on sigma, this is manual modification

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
    // void max_scale();
};

template<typename type>
using demons_cpu = demons<type,vector_cpu<type>>;

// template<typename type>
// using demons_gpu = demons<type,vector_opencl<type>>;

#ifdef IMART_WITH_OPENCL
template<typename type>
using demons_opencl = demons<type,vector_opencl<type>>;
#endif

#ifdef IMART_WITH_CUDA
template<typename type>
using demons_cuda = demons<type,vector_cuda<type>>;
#endif

// ===========================================
//      Functions of Class demons
// ===========================================

// ===========================================
// Functions
// ===========================================
template <typename type, typename container>
type demons<type,container>::cost()
{
    // std::cout << "init cost" << std::endl;
    this->warped_moving(); // update moving_prime
    // auto moving_prime = this->warped_moving();
    // moving_prime->print_data("moving prime");
    // std::cout << "dif cost" << std::endl;
    // fixed->print("dem fix");
    // moving_prime->print("dem mov");
    type N = (type)fixed->get_total_elements();
    image<type,container> ssd_ = *fixed - *moving_prime;
    cost_value = (0.5/N)*( (ssd_*ssd_).sum() );
    // printf("cost: %7.3e\n", cost_value);
    // std::cout << "end cost" << std::endl;
    return cost_value;
};

template <typename type, typename container>
typename transform<type,container>::pointer demons<type,container>::derivative()
{
    int d = fixed->get_dimension();
    // transform
    typename transform<type,container>::pointer trfm = transformation->mimic();
    
    // mimic parameters
    // std::cout << "Params:" << std::endl;
    auto param = std::make_shared< std::vector< typename image<type,container>::pointer >>(d);
    for(int i = 0; i < d; i++)
    {
        (*param)[i] = image<type,container>::new_pointer(d);
        (*param)[i]->mimic_(*(transformation->get_parameters(i)));
    };
    // moving prime already computed (at initialization or when cost function)
    // auto moving_prime = this->warped_moving();        // consider store this to avoid compute again
    
    // std::cout << "Optical Flow:" << std::endl;
    if (transformation->get_name() == "dfield")
    {
        optical_flow(fixed, moving_prime, param);
    };

    trfm->set_parameters_vector( param );
    return trfm;
};

template <typename type, typename container>
void demons<type,container>::optical_flow(std::shared_ptr<image<type,container>> img_fixed, 
                                          std::shared_ptr<image<type,container>> img_moving, 
                                          std::shared_ptr< std::vector< std::shared_ptr < image<type,container> >>> param )
{
    // std::cout << "Init flow:" << std::endl;
    int d = img_fixed->get_dimension();
    type low = 1e-7;    // value considered as 0

    auto dif = (*img_moving - *img_fixed);
    auto dif_sq = dif*dif;
    
    // Computing gradient
    // std::cout << "Gradient:" << std::endl;
    auto grad = gradient(img_moving); //moving prime

    // dimensional data
    auto gx = *grad[0];
    auto gy = *grad[1];
    
    if (d == 2)
    {
        // compute mag
        auto mag_sq = (gx*gx + gy*gy);

        // scale parameter based on optical flow
        auto scale = dif/(mag_sq + dif_sq);
        scale.replace(mag_sq <= low,0.0);   // remove 0 division
        scale.replace(dif_sq <= low,0.0);   // remove 0 division

        // update depends of space and orientation
        // type spacex = img_fixed->get_spacing()[0]*img_fixed->get_direction()[0];
        // type spacey = img_fixed->get_spacing()[1]*img_fixed->get_direction()[3];

        // demons derivative of dfield
        auto ux = scale*gx;//*spacex;
        auto uy = scale*gy;//*spacey;

        // no update where image is empty
        ux.replace(*img_fixed <= low, 0.0);
        ux.replace(*img_moving <= low, 0.0);
        uy.replace(*img_fixed <= low, 0.0);
        uy.replace(*img_moving <= low, 0.0);

        // return parameters
        *((*param)[0]) = ux;
        *((*param)[1]) = uy;
    }
    else if (d == 3)
    {
        // z
        auto gz = *grad[2];
        // compute mag
        auto mag_sq = (gx*gx + gy*gy + gz*gz);

        // scale parameter based on optical flow
        auto scale = dif/(mag_sq + dif_sq);
        scale.replace(mag_sq <= low,0.0);   // remove 0 division
        scale.replace(dif_sq <= low,0.0);   // remove 0 division

        // update depends of space and orientation
        // type spacex = img_fixed->get_spacing()[0]*img_fixed->get_direction()[0];
        // type spacey = img_fixed->get_spacing()[1]*img_fixed->get_direction()[4];
        // type spacez = img_fixed->get_spacing()[2]*img_fixed->get_direction()[8];

        // demons derivative of dfield
        // std::cout << "Scale:" << std::endl;
        auto ux = scale*gx;//*spacex;
        auto uy = scale*gy;//*spacey;
        auto uz = scale*gz;//*spacez;

        // no update where image is empty
        ux.replace(*img_fixed <= low, 0.0);
        ux.replace(*img_moving <= low, 0.0);
        uy.replace(*img_fixed <= low, 0.0);
        uy.replace(*img_moving <= low, 0.0);
        uz.replace(*img_fixed <= low, 0.0);
        uz.replace(*img_moving <= low, 0.0);

        // return parameters
        *((*param)[0]) = ux;
        *((*param)[1]) = uy;
        *((*param)[2]) = uz;

        // std::cout << "End flow:" << std::endl;
    };
    return;
};

// template <typename type, typename container>
// typename transform<type,container>::pointer demons<type,container>::derivative()
// {
//     int d = fixed->get_dimension();
//     // transform
//     typename transform<type,container>::pointer trfm = transformation->mimic();
    
//     // mimic parameters
//     auto param = std::make_shared< std::vector< typename image<type,container>::pointer >>(d);
//     for(int i = 0; i < d; i++)
//     {
//         (*param)[i] = image<type,container>::new_pointer();
//         (*param)[i]->mimic_(*(transformation->get_parameters(i)));
//     };
//     // moving prime already computed
//     // auto moving_prime = this->warped_moving();        // consider store this to avoid compute again
    
//     type low = 1e-7;    // value considered as 0

//     auto dif = (*moving_prime - *fixed);
//     auto dif_sq = dif*dif;

//     if (transformation->get_name() == "dfield")
//     {
//         // Computing gradient
//         auto grad = gradient(moving_prime);

//         // dimensional data
//         auto gx = *grad[0];
//         auto gy = *grad[1];
        
//         if (d == 2)
//         {
//             // compute mag
//             auto mag_sq = (gx*gx + gy*gy);

//             // scale parameter based on optical flow
//             auto scale = dif/(mag_sq + dif_sq);
//             scale.replace(mag_sq <= low,0.0);   // remove 0 division
//             scale.replace(dif_sq <= low,0.0);   // remove 0 division

//             // update depends of space and orientation
//             type spacex = fixed->get_spacing()[0]*fixed->get_direction()[0];
//             type spacey = fixed->get_spacing()[1]*fixed->get_direction()[3];

//             // demons derivative of dfield
//             auto ux = scale*gx*spacex;
//             auto uy = scale*gy*spacey;

//             // no update where image is empty
//             ux.replace(*fixed <= low, 0.0);
//             ux.replace(*moving_prime <= low, 0.0);
//             uy.replace(*fixed <= low, 0.0);
//             uy.replace(*moving_prime <= low, 0.0);

//             // return parameters
//             *((*param)[0]) = ux;
//             *((*param)[1]) = uy;
//         }
//         else if (d == 3)
//         {
//             // z
//             auto gz = *grad[2];
//             // compute mag
//             auto mag_sq = (gx*gx + gy*gy + gz*gz);

//             // scale parameter based on optical flow
//             auto scale = dif/(mag_sq + dif_sq);
//             scale.replace(mag_sq <= low,0.0);   // remove 0 division
//             scale.replace(dif_sq <= low,0.0);   // remove 0 division

//             // update depends of space and orientation
//             type spacex = fixed->get_spacing()[0]*fixed->get_direction()[0];
//             type spacey = fixed->get_spacing()[1]*fixed->get_direction()[4];
//             type spacez = fixed->get_spacing()[2]*fixed->get_direction()[8];

//             // demons derivative of dfield
//             auto ux = scale*gx*spacex;
//             auto uy = scale*gy*spacey;
//             auto uz = scale*gz*spacez;

//             // no update where image is empty
//             ux.replace(*fixed <= low, 0.0);
//             ux.replace(*moving_prime <= low, 0.0);
//             uy.replace(*fixed <= low, 0.0);
//             uy.replace(*moving_prime <= low, 0.0);
//             uz.replace(*fixed <= low, 0.0);
//             uz.replace(*moving_prime <= low, 0.0);

//             // return parameters
//             *((*param)[0]) = ux;
//             *((*param)[1]) = uy;
//             *((*param)[2]) = uz;
//         };
//     };

//     trfm->set_parameters_vector( param );
//     return trfm;
// };

}; //end namespace

#endif