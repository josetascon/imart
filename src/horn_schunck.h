/*
* @Author: jose
* @Date:   2019-12-09 00:00:00
* @Last Modified by:   jose
* @Last Modified time: 2020-03-18 00:00:00
*/

#ifndef __HORN_SCHUNCK_H__
#define __HORN_SCHUNCK_H__

#include "metric.h"

namespace imart
{

// Class metric
template <typename type, typename container=vector_cpu<type>>
class horn_schunck: public inherit<horn_schunck<type,container>, metric<type,container>>
{
public:
    //Type definitions
    using self    = horn_schunck;
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

    using inherit<horn_schunck<type,container>, metric<type,container>>::inherit;

    // double sigma;
    // int kernel_width;
    // void init_sigma();

    typename image<type,container>::pointer kernel_x;
    typename image<type,container>::pointer kernel_y;
    typename image<type,container>::pointer kernel_z;
    typename image<type,container>::pointer kernel_t;

    typename image<type,container>::pointer fx;
    typename image<type,container>::pointer fy;
    typename image<type,container>::pointer fz;
    typename image<type,container>::pointer ft;

protected:
    void optical_flow(std::shared_ptr<image<type,container>> img_fixed, 
                      std::shared_ptr<image<type,container>> img_moving, 
                      std::shared_ptr< std::vector< std::shared_ptr < image<type,container> >>> param );

public:
    // ===========================================
    // Create Functions
    // ===========================================
    // Constructors
    horn_schunck() : inherit<horn_schunck<type,container>, metric<type,container>>()
          { this->class_name = "horn_schunck"; };//init_sigma(); };

    horn_schunck(int d) : inherit<horn_schunck<type,container>, metric<type,container>>(d)
               { this->class_name = "horn_schunck"; };//init_sigma(); };

    horn_schunck(typename image<type,container>::pointer fixed_image, 
        typename image<type,container>::pointer moving_image)
        : inherit<horn_schunck<type,container>, metric<type,container>>(fixed_image, moving_image)
        { this->class_name = "horn_schunck"; };//init_sigma(); };

    horn_schunck(typename image<type,container>::pointer fixed_image, 
        typename image<type,container>::pointer moving_image,
        typename transform<type,container>::pointer transformd)
        : inherit<horn_schunck<type,container>, metric<type,container>>(fixed_image, moving_image, transformd)
        { this->class_name = "horn_schunck"; };//init_sigma(); };

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
    void resolution_update();
};

template<typename type>
using horn_schunck_cpu = horn_schunck<type,vector_cpu<type>>;

// template<typename type>
// using horn_schunck_gpu = horn_schunck<type,vector_opencl<type>>;

#ifdef IMART_WITH_OPENCL
template<typename type>
using horn_schunck_opencl = horn_schunck<type,vector_opencl<type>>;
#endif

#ifdef IMART_WITH_CUDA
template<typename type>
using horn_schunck_cuda = horn_schunck<type,vector_cuda<type>>;
#endif

// ===========================================
//      Functions of Class horn_schunck
// ===========================================

// ===========================================
// Init Functions
// ===========================================
template <typename type, typename container>
type horn_schunck<type,container>::init()
{
    // HSKERN = np.array([[1/12, 1/6, 1/12],
    //                [1/6,    0, 1/6],
    //                [1/12, 1/6, 1/12]], float)

    // kernel_x = np.array([[-1, 1],
    //                     [-1, 1]]) * .25  # kernel for computing d/dx

    // kernel_y = np.array([[-1, -1],
    //                     [1, 1]]) * .25  # kernel for computing d/dy

    // kernel_t = np.ones((2, 2))*.25

    int d = fixed->get_dimension();

    if (d == 2)
    {

        auto a = container<type>::new_pointer(std::initializer_list{-1,1,-1,1});
        auto b = container<type>::new_pointer(std::initializer_list{-1,-1,1,1});

        kernel_x = image<type,container>::new_pointer(a,2,2);
        kernel_y = image<type,container>::new_pointer(b,2,2);
        kernel_t = image<type,container>::new_pointer(2,2);
        kernel_t->ones();

        *kernel_x = (*kernel_x)*0.25;
        *kernel_y = (*kernel_y)*0.25;
        *kernel_t = (*kernel_t)*0.25;
    }
};


// ===========================================
// Functions
// ===========================================
template <typename type, typename container>
type horn_schunck<type,container>::cost()
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
typename transform<type,container>::pointer horn_schunck<type,container>::derivative()
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
void horn_schunck<type,container>::optical_flow(std::shared_ptr<image<type,container>> img_fixed, 
                                          std::shared_ptr<image<type,container>> img_moving, 
                                          std::shared_ptr< std::vector< std::shared_ptr < image<type,container> >>> param )
{
    // std::cout << "Init flow:" << std::endl;
    int d = img_fixed->get_dimension();
    
    if (d == 2)
    {
        auto ux_last = *transformation->get_parameters(0);
        auto uy_last = *transformation->get_parameters(1);

        auto num = (*fx)*(*ux_last) + (*fy)*(*uy_last) + *ft;
        auto den = (*fx)^2 + (*fy)^2 + alpha^2;
        auto der =  num / den;

        // horn_schunck derivative of dfield
        auto ux = der*gx;
        auto uy = der*gy;

        // return parameters
        *((*param)[0]) = ux;
        *((*param)[1]) = uy;
        // std::cout << "End flow:" << std::endl;
    }
    else if (d == 3)
    {
        // horn_schunck derivative of dfield
        // auto ux = der*gx;
        // auto uy = der*gy;
        // auto uz = der*gz;

        // return parameters
        // *((*param)[0]) = ux;
        // *((*param)[1]) = uy;
        // *((*param)[2]) = uz;
        // std::cout << "End flow:" << std::endl;
    };
    return;
};

template <typename type, typename container>
void horn_schunck<type,container>::resolution_update()
{
    // compute derivatives
    auto im1x = convolve(fixed, kernel_x);
    auto im1y = convolve(fixed, kernel_y);
    auto im1t = convolve(fixed, kernel_t);

    auto im2x = convolve(moving, kernel_x);
    auto im2y = convolve(moving, kernel_y);
    auto im2t = convolve(moving, kernel_t);

    *fx = *im1x + *im2x;
    *fy = *im1y + *im2y;
    *ft = *im1t + *im2t;
};

}; //end namespace

#endif