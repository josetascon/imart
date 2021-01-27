/*
* @Author: jose
* @Date:   2019-12-09 00:00:00
* @Last Modified by:   jose
* @Last Modified time: 2020-03-18 00:00:00
*/

#ifndef __LDDMM_H__
#define __LDDMM_H__

#include "metric.h"
#include "regularizer.h"

namespace imart
{

// NOTE: Doing this I changed image_utils operations input1 to reference (&). Before input1 copy, maybe error??

// Class metric
template <typename type, typename container=vector_cpu<type>>
class lddmm: public inherit<lddmm<type,container>, metric<type,container>>
{
public:
    //Type definitions
    using self    = lddmm;
    using pointer = std::shared_ptr<self>;
    using vector  = std::vector<self::pointer>;

    // Alias
    using ptr_image = typename image<type,container>::pointer;
    using image_group = typename image<type,container>::vector;
    using field_group = std::vector< image_group >;
    using ptr_image_group = std::shared_ptr< image_group >;
    using ptr_field_group = std::shared_ptr< field_group >;
    using ptr_trfm = typename transform<type,container>::pointer;
    using ptr_grid = typename grid<type,container>::pointer;

    ptr_field_group v;
    ptr_field_group dv;

    ptr_image_group j0;
    ptr_image_group j1;

    ptr_field_group phi0;
    ptr_field_group phi1;
    
protected:
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

    using inherit<lddmm<type,container>, metric<type,container>>::inherit;

    // Variables
    int tsteps;
    int integration_steps;
    double learning_rate;
    double sigma;
    
    regularizer<type,container> regularize;

    // Init Functions
    void init();
    void init_fields();
    // NOTE: Consider init functions without zero initialization of images
    ptr_field_group init_field_group(ptr_image img);
    ptr_image_group init_image_group(ptr_image img, int n, bool zero = true);
    ptr_image sample(ptr_image img, ptr_grid xr, ptr_image_group xo);
    type norm(std::vector<std::shared_ptr<image<type,container>>> & a);

    // Functions
    ptr_image_group update_lddmm();
    ptr_field_group integrate_velocity_backward();
    ptr_field_group integrate_velocity_forward();
    ptr_image backward_alpha(ptr_image v_td, ptr_grid x_vtd, ptr_image_group x);
    ptr_image forward_alpha(ptr_image v_td, ptr_grid x_vtd, ptr_image_group x);
    ptr_image_group push_forward(ptr_field_group phi0);
    ptr_image_group pull_backward(ptr_field_group phi1);
    ptr_field_group image_gradients(ptr_image_group j0);
    ptr_image_group jacobian_determinant(ptr_field_group phi);
    bool is_injectivity_violated(ptr_image_group det_phi);


public:
    // ===========================================
    // Create Functions
    // ===========================================
    // Constructors
    lddmm() : inherit<lddmm<type,container>, metric<type,container>>()
          { this->class_name = "lddmm"; init(); };

    lddmm(int d) : inherit<lddmm<type,container>, metric<type,container>>(d)
               { this->class_name = "lddmm"; init(); };

    lddmm(ptr_image fixed_image, ptr_image moving_image)
        : inherit<lddmm<type,container>, metric<type,container>>(fixed_image, moving_image)
        { this->class_name = "lddmm"; init(); };

    lddmm(ptr_image fixed_image, ptr_image moving_image, ptr_trfm transformd)
        : inherit<lddmm<type,container>, metric<type,container>>(fixed_image, moving_image, transformd)
        { this->class_name = "lddmm"; init(); };

    // ===========================================
    // Get Functions
    // ===========================================
    int get_time_steps() const;
    int get_integration_steps() const;
    double get_learning_rate() const;
    double get_alpha() const;
    double get_gamma() const;
    double get_sigma() const;

    // ===========================================
    // Set Functions
    // ===========================================
    void set_time_steps(int ts);
    void set_integration_steps(int is);
    void set_learning_rate(double l);
    void set_alpha_gamma(double a, double g);
    void set_sigma(double g);

    // ===========================================
    // Functions
    // ===========================================
    // !compute the cost
    type cost();
    // !calculate derivative
    ptr_trfm derivative();
    // !regularize
    // void regularize();
    // !resolution
    void resolution_update();
    // // !max scale
    // void max_scale();
};

template<typename type>
using lddmm_cpu = lddmm<type,vector_cpu<type>>;

template<typename type>
using lddmm_gpu = lddmm<type,vector_ocl<type>>;

#ifdef IMART_WITH_OPENCL
template<typename type>
using lddmm_ocl = lddmm<type,vector_ocl<type>>;
#endif

#ifdef IMART_WITH_CUDA
template<typename type>
using lddmm_cuda = lddmm<type,vector_cuda<type>>;
#endif

// ===========================================
//      Functions of Class lddmm
// ===========================================


// ===========================================
// Get Functions
// ===========================================
template <typename type, typename container>
int lddmm<type,container>::get_time_steps() const
{
    return tsteps;
};

template <typename type, typename container>
int lddmm<type,container>::get_integration_steps() const
{
    return integration_steps;
};

template <typename type, typename container>
double lddmm<type,container>::get_learning_rate() const
{
    return learning_rate;
};

template <typename type, typename container>
double lddmm<type,container>::get_alpha() const
{
    return regularize->get_alpha();
};

template <typename type, typename container>
double lddmm<type,container>::get_gamma() const
{
    return regularize->get_gamma();
};

template <typename type, typename container>
double lddmm<type,container>::get_sigma() const
{
    return sigma;
};

// ===========================================
// Set Functions
// ===========================================
template <typename type, typename container>
void lddmm<type,container>::set_time_steps(int ts)
{
    tsteps = ts;
    init_fields();
};

template <typename type, typename container>
void lddmm<type,container>::set_integration_steps(int is)
{
    integration_steps = is;
};

template <typename type, typename container>
void lddmm<type,container>::set_learning_rate(double l)
{
    learning_rate = l;
};

template <typename type, typename container>
void lddmm<type,container>::set_alpha_gamma(double a, double g)
{
    regularize.set_alpha(a);
    regularize.set_gamma(g);
    regularize.update_a(fixed);
};

template <typename type, typename container>
void lddmm<type,container>::set_sigma(double s)
{
    sigma = s;
};

// ===========================================
// Init Functions
// ===========================================
template <typename type, typename container>
void lddmm<type,container>::init()
{
    integration_steps = 5;
    learning_rate = 0.01;
    tsteps = 16;
    
    set_sigma(0.1);
    set_alpha_gamma(1.0,1.0);

    init_fields();

    // j0 = init_image_group(fixed, tsteps);
    // j0->at(tsteps-1) = fixed->clone();
};

template <typename type, typename container>
void lddmm<type,container>::init_fields()
{
    // std::cout << "init fields" << std::endl;
    v = init_field_group(fixed);
    dv = init_field_group(fixed);

    // initialize to be able to call cost function
    j0 = init_image_group(fixed, tsteps, false);
    j0->at(tsteps-1) = fixed->clone();

    // learning_rate = 0.001;
    // sigma = 0.1;
    
    // regularize.set_alpha(100.0);
    // regularize.set_gamma(1.0);
    // regularize.update_a(fixed);
};

template <typename type, typename container>
typename lddmm<type,container>::ptr_field_group lddmm<type,container>::init_field_group(typename lddmm<type,container>::ptr_image img)
{
    // init pointer and outer vector
    auto output = std::make_shared< field_group >(tsteps);
    // dimension
    int dim = img->get_dimension();
    // init transforms vectors
    for(int t = 0; t < tsteps; t++)
    {
        output->at(t) = image_group(dim);
        
        for(int d = 0; d < dim; d++)
        {
            output->at(t)[d] = img->mimic();
            output->at(t)[d]->zeros();
        };
    };
    return output;
};

template <typename type, typename container>
typename lddmm<type,container>::ptr_image_group lddmm<type,container>::init_image_group(typename lddmm<type,container>::ptr_image img, int n, bool zero)
{
    // init pointer and outer vector
    auto output = std::make_shared< image_group >(n);
    // init image vectors
    for(int k = 0; k < n; k++)
    {
        output->at(k) = img->mimic();
        if(zero) { output->at(k)->zeros(); }
    };
    return output;
};

template <typename type, typename container>
typename lddmm<type,container>::ptr_image lddmm<type,container>::sample(typename lddmm<type,container>::ptr_image img, 
    typename lddmm<type,container>::ptr_grid xr, typename lddmm<type,container>::ptr_image_group xo)
{
    // xr used to create a similar grid
    auto xout = xr->mimic();
    xout->set_grid(xo);

    auto itp = ilinear<type,container>::new_pointer(img);
    return itp->apply(xout);
};

template <typename type, typename container>
type lddmm<type,container>::norm(std::vector<std::shared_ptr<image<type,container>>> & a)
{
    type e_d = 0.0;
    for(int i = 0; i < a.size(); i++)
    {
        e_d = e_d + ((*a[i])*(*a[i])).sum();
    };
    return std::sqrt(e_d);
};

// ===========================================
// Functions
// ===========================================
template <typename type, typename container>
type lddmm<type,container>::cost()
{
    // type N = (type)fixed->get_total_elements();
    // image<type,container> ssd_ = *j0->at(tsteps-1) - *moving;
    // cost_value = (0.5/N)*( (ssd_*ssd_).sum() );
    // return cost_value;
    
    int dim = fixed->get_dimension();

    image<type,container> ssd_ = *j0->at(tsteps-1) - *moving;
    type e_intensity = (1.0/(sigma*sigma))*( (ssd_*ssd_).sum() );

    type e_regularize = 0.0;
    for(int t = 0; t < tsteps; t++)
    {
        auto vv = std::make_shared< std::vector<std::shared_ptr<image<type,container>>> >(dim);
        for(int i = 0; i < dim; i++) vv->at(i) = v->at(t)[i];
        auto a = regularize.l(vv);
        e_regularize = e_regularize + norm(a);
    };
    cost_value = e_intensity + e_regularize;
    return cost_value;
};

template <typename type, typename container>
void lddmm<type,container>::resolution_update()
{
    int dim = fixed->get_dimension();

    regularize.update_a(fixed);

    for(int t = 0; t < tsteps; t++)
    {
        for(int d = 0; d < dim; d++)
        {
            auto intepolate_velocity = this->interpolate_fixed->mimic();
            intepolate_velocity->set_reference(v->at(t)[d]);
            v->at(t)[d] = intepolate_velocity->apply(x0);
            dv->at(t)[d] = v->at(t)[d]->mimic();
            dv->at(t)[d]->zeros();
        };
    };
};

template <typename type, typename container>
typename lddmm<type,container>::ptr_trfm lddmm<type,container>::derivative()
{
    //dimension
    int dim = fixed->get_dimension();
    
    // transform
    ptr_trfm trfm = transformation->mimic();
    
    // mimic parameters
    // std::cout << "Params:" << std::endl;
    auto param = std::make_shared< std::vector< ptr_image >>(dim);
    for(int i = 0; i < dim; i++)
    {
        (*param)[i] = image<type,container>::new_pointer(dim);
        (*param)[i]->mimic_(*(transformation->get_parameters(i)));
    };
    // moving prime already computed
    // auto moving_prime = this->warped_moving();        // consider store this to avoid compute again
    
    // std::cout << "Update:" << std::endl;
    if (transformation->get_name() == "dfield")
    {
        auto pu = update_lddmm();
        for(int i = 0; i < dim; i++) param->at(i)->set_data( pu->at(i)->get_data());
    };

    trfm->set_parameters_vector( param );
    return trfm;
};

template <typename type, typename container>
typename lddmm<type,container>::ptr_image_group lddmm<type,container>::update_lddmm()
{
    //dimension
    int dim = fixed->get_dimension();
    
    // (1): Calculate new estimate of velocity
    // std::cout << "(1)" << std::endl;
    for(int t = 0; t < tsteps; t++)
    {
        // v->at(t) =  v->at(t) - (dv->at(t) * learning_rate);
        auto dd = (dv->at(t) * learning_rate);
        v->at(t) = v->at(t) - dd;
    };

    // (3): calculate backward flows
    // std::cout << "(3)" << std::endl;
    // ptr_field_group phi1 = integrate_velocity_backward();
    phi1 = integrate_velocity_backward();

    // (4): calculate forward flows
    // std::cout << "(4)" << std::endl;
    // ptr_field_group phi0 = integrate_velocity_forward();
    phi0 = integrate_velocity_forward();

    // (5): push-forward i0
    // std::cout << "(5)" << std::endl;
    // ptr_image_group j0 = push_forward(phi0);
    j0 = push_forward(phi0);

    // (6): pull back i1
    // std::cout << "(6)" << std::endl;
    // ptr_image_group j1 = pull_backward(phi1);
    j1 = pull_backward(phi1);

    // (7): Calculate image gradient
    // std::cout << "(7)" << std::endl;
    ptr_field_group dj0 = image_gradients(j0);
    
    // (8): Calculate Jacobian determinant of the transformation
    // std::cout << "(8)" << std::endl;
    ptr_image_group det_phi1 = jacobian_determinant(phi1);
    if (is_injectivity_violated(det_phi1))
    {
        this->fail = true;
        this->fail_msg = "injectivity failure";
    };

    // (9): Calculate the gradient
    // std::cout << "(9)" << std::endl;
    ptr_image_group g = init_image_group(j0->at(0), dim, false);
    for(int t = 0; t < tsteps; t++)
    {
        auto s = *(j0->at(t)) - *(j1->at(t));
        auto p = (*(det_phi1->at(t))) * s * (2 / (sigma*sigma));
        for(int d = 0; d < dim; d++) *(g->at(d)) = p * (*(dj0->at(t)[d]));
        auto r = regularize.k(g);
        auto v_t_2 = v->at(t)*2.0;
        dv->at(t) = v_t_2  - r;
        // de = 2 / sigma**2 * (det_phi1[t])[np.newaxis] * dj0[t] * (j0[t] - j1[t])[np.newaxis]
        // dv[t] = 2*v[t] - bregularizer.K(de)
    };

    // (10) calculate norm of the gradient, stop if small
    // std::cout << "(10)" << std::endl;
    // dv_norm = np.linalg.norm(dv)
    // if dv_norm < 0.001:
    //     print(dv_norm)
    //     termination = 'Gradient norm below threshold'
    //     break
    
    // copy parameters
    // transformation based on gradient descent
    auto param = std::make_shared< std::vector< ptr_image >>(dim);
    for(int d = 0; d < dim; d++) 
    {
        // std::cout << d << std::endl;
        // param->at(d) = phi0->at(tsteps-1)[d];
        param->at(d) = phi1->at(0)[d]->clone(); // moving to fixed transformation
        *param->at(d) = (*param->at(d)) - *((x0->get_grid())->at(d));
        *param->at(d) = (*transformation->get_parameters(d) - (*param->at(d)))*(1/learning_rate);
        // param->at(d)->print();
    }

    return param;
};

template <typename type, typename container>
typename lddmm<type,container>::ptr_field_group lddmm<type,container>::integrate_velocity_backward()
{
    // dimension
    int dim = fixed->get_dimension();

    // x
    ptr_image_group x = x0->get_grid();

    // create temporary x pointer image vector
    auto xa = std::make_shared< image_group >(x->size());
    for(int k = 0; k < x->size(); k++)
    {
        xa->at(k) = x->at(k)->mimic();
        // xa->at(k)->zeros();
    }

    // create the fields
    auto phi1 = init_field_group(fixed);

    // reference grid phi1, selecting [0,0], same for all
    ptr_grid x_phi1 = grid<type,container>::new_pointer(phi1->at(0)[0]);

    // reference grid v, selecting [0,0], same for all
    ptr_grid x_vtd = grid<type,container>::new_pointer(v->at(0)[0]);

    // phi1_1 is the identity mapping
    for(int d = 0; d < dim; d++) phi1->at(tsteps - 1)[d] = (*x)[d]->clone();

    for(int t = tsteps-2; t >= 0; t--)
    {
        // std::cout << t << " " << std::endl;
        for(int d = 0; d < dim; d++)
        {
            auto alpha = backward_alpha(v->at(t)[d], x_vtd,  x);

            // operation: xa = x + alpha
            for(int k = 0; k < x->size(); k++) *(xa->at(k)) = *(x->at(k)) + (*alpha);

            phi1->at(t)[d] = sample(phi1->at(t + 1)[d], x_phi1, xa);
        };

    };

    // for t in range(tsteps-2, -1, -1):
    //     for d in range(dim):
    //         alpha = backward_alpha(v[t,d], x)
    //         phi1[t,d] = sample_2d(x + alpha, phi1[t + 1,d])

    return phi1;
};

template <typename type, typename container>
typename lddmm<type,container>::ptr_image lddmm<type,container>::backward_alpha(typename lddmm<type,container>::ptr_image v_td, 
        typename lddmm<type,container>::ptr_grid x_vtd, typename lddmm<type,container>::ptr_image_group x)
{
    // pointer image with zeros
    ptr_image alpha = v_td->mimic();
    alpha->zeros();

    // create temporary pointer image vector
    auto xa = std::make_shared< image_group >(x->size());
    for(int k = 0; k < x->size(); k++) 
    {
        xa->at(k) = x->at(k)->mimic();
        // xa->at(k)->zeros();
    }

    for(int i = 0; i < integration_steps; i++)
    {   
        // operation: xa = x + 0.5 * alpha
        for(int k = 0; k < x->size(); k++) *(xa->at(k)) = *(x->at(k)) + (*alpha)*0.5;
        // interpolation
        alpha = sample(v_td, x_vtd, xa);
    };
    return alpha;
};

template <typename type, typename container>
typename lddmm<type,container>::ptr_field_group lddmm<type,container>::integrate_velocity_forward()
{
    // dimension
    int dim = fixed->get_dimension();

    // x
    ptr_image_group x = x0->get_grid();

    // create temporary x pointer image vector
    auto xa = std::make_shared< image_group >(x->size());
    for(int k = 0; k < x->size(); k++)
    {
        xa->at(k) = x->at(k)->mimic();
        // xa->at(k)->zeros();
    }

    // create the fields
    auto phi0 = init_field_group(fixed);

    // reference grid phi0, selecting [0,0], same for all
    ptr_grid x_phi0 = grid<type,container>::new_pointer(phi0->at(0)[0]);

    // reference grid v, selecting [0,0], same for all
    ptr_grid x_vtd = grid<type,container>::new_pointer(v->at(0)[0]);

    // phi0_0 is the identity mapping
    for(int d = 0; d < dim; d++) phi0->at(0)[d] = (*x)[d]->clone();

    for(int t = 0; t < tsteps - 1; t++)
    {
        for(int d = 0; d < dim; d++)
        {
            auto alpha = forward_alpha(v->at(t)[d], x_vtd,  x);

            // operation: xa = x + alpha
            for(int k = 0; k < x->size(); k++) *(xa->at(k)) = *(x->at(k)) - (*alpha);

            phi0->at(t+1)[d] = sample(phi0->at(t)[d], x_phi0, xa);
        };

    };

    // for t in range(0, tsteps-1):
    //     for d in range(dim):
    //         alpha = forward_alpha(v[t,d], x)
    //         phi0[t+1,d] = sample_2d(x - alpha, phi0[t,d])
    return phi0;
};

template <typename type, typename container>
typename lddmm<type,container>::ptr_image lddmm<type,container>::forward_alpha(typename lddmm<type,container>::ptr_image v_td, 
        typename lddmm<type,container>::ptr_grid x_vtd, typename lddmm<type,container>::ptr_image_group x)
{
    // pointer image with zeros
    ptr_image alpha = v_td->mimic();
    alpha->zeros();

    // create temporary pointer image vector
    auto xa = std::make_shared< image_group >(x->size());
    for(int k = 0; k < x->size(); k++)
    {
        xa->at(k) = x->at(k)->mimic();
        // xa->at(k)->zeros();
    }

    for(int i = 0; i < integration_steps; i++)
    {   
        // operation: xa = x - 0.5 * alpha
        for(int k = 0; k < x->size(); k++) *(xa->at(k)) = *(x->at(k)) - (*alpha)*0.5;
        // interpolation
        alpha = sample(v_td, x_vtd, xa);
    };
    return alpha;
};

template <typename type, typename container>
typename lddmm<type,container>::ptr_image_group lddmm<type,container>::push_forward(typename lddmm<type,container>::ptr_field_group phi0)
{
    auto j0 = init_image_group(fixed, tsteps, false);
    for(int t = 0; t < tsteps; t++)
    {
        auto p = std::make_shared< image_group >(phi0->at(t).size());
        *p = phi0->at(t);
        j0->at(t) = sample(fixed, x0, p);
    };
    return j0;
};

template <typename type, typename container>
typename lddmm<type,container>::ptr_image_group lddmm<type,container>::pull_backward(typename lddmm<type,container>::ptr_field_group phi1)
{
    auto j1 = init_image_group(fixed, tsteps, false);
    for(int t = tsteps-1; t >= 0 ; t--)
    {
        auto p = std::make_shared< image_group >(phi1->at(t).size());
        *p = phi1->at(t);
        j1->at(t) = sample(moving, x0, p);
    }
    return j1;
};

template <typename type, typename container>
typename lddmm<type,container>::ptr_field_group lddmm<type,container>::image_gradients(typename lddmm<type,container>::ptr_image_group j0)
{
    int dim = fixed->get_dimension();
    auto dj0 = init_field_group(j0->at(0));

    for(int t = 0; t < tsteps; t++)
    {
        dj0->at(t) = gradient(j0->at(t));
    };
    return dj0;
};

template <typename type, typename container>
typename lddmm<type,container>::ptr_image_group lddmm<type,container>::jacobian_determinant(typename lddmm<type,container>::ptr_field_group phi)
{
    auto det_phi = init_image_group(phi->at(0)[0], tsteps, false);

    for(int t = 0; t < tsteps; t++)
    {
        auto p = std::make_shared< image_group >(phi->at(t).size());
        *p = phi->at(t);
        det_phi->at(t) = jacobian(p);
    };
    return det_phi;
};

template <typename type, typename container>
bool lddmm<type,container>::is_injectivity_violated(typename lddmm<type,container>::ptr_image_group det_phi)
{
    for(int t = 0; t < tsteps; t++)
    {
        if (det_phi->at(t)->min() < 0) { std::cout << t << " "; return true; };
    };
    return false;
};


}; //end namespace

#endif