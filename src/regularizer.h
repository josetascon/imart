/*
* @Author: jose
* @Date:   2020-01-27 00:00:00
* @Last Modified by:   jose
* @Last Modified time: 2020-01-27 00:00:00
*/

#ifndef __REGULARIZER_H__
#define __REGULARIZER_H__

#include "image.h"
#include "image_utils.h"
#include "process_object.h"

namespace imart
{

template <typename type, typename container>
class regularizer: public inherit<regularizer<type,container>, process_object>
{
public:
    //Type definitions
    using self    = regularizer;
    using pointer = std::shared_ptr<self>;
    using vector  = std::vector<self::pointer>;

    // template <typename type, typename container>
    // using ptr_image<type,container> = typename image<type,container>::pointer;
    // template <typename type, typename container>
    // using image_group = typename image<type,container>::vector;
    // template <typename type, typename container>
    // using ptr_image_group = std::shared_ptr< image_group >;

    using inherit<regularizer<type,container>, process_object>::inherit;

protected:
    // ===========================================
    // Internal Variables
    // ===========================================
    bool ready;
    double alpha;
    double gamma;
    std::vector<std::shared_ptr<image<type,container>>> A;

    // ===========================================
    // Functions
    // ===========================================
    void init(double a, double g);

public:
    // ===========================================
    // Create Functions
    // ===========================================
    // Constructors
    regularizer();
    regularizer(double a, double g);
    

    // ===========================================
    // Get Functions
    // ===========================================
    double get_alpha() const;
    double get_gamma() const;

    // ===========================================
    // Set Functions
    // ===========================================
    void set_alpha(double a);
    void set_gamma(double g);

    // ===========================================
    // Print Functions
    // ===========================================
    virtual std::string info(std::string msg);

    // ===========================================
    // Functions
    // ===========================================
    // compute
    // template <typename type, typename container>
    std::vector<std::shared_ptr<image<type,container>>> update_a(std::shared_ptr< image<type,container> > sample_img);

    // template <typename type, typename container>
    std::vector<std::shared_ptr<image<type,container>>> k(std::shared_ptr< std::vector<std::shared_ptr<image<type,container>>> > g);
    
    // template <typename type, typename container>
    std::vector<std::shared_ptr<image<type,container>>> l(std::shared_ptr< std::vector<std::shared_ptr<image<type,container>>> > g);

    // template <typename type, typename container>
    std::vector< std::vector<std::shared_ptr<image<type,container>>> > fft_vector(std::shared_ptr< std::vector<std::shared_ptr<image<type,container>>> > input);
    
    // template <typename type, typename container>
    std::vector<std::shared_ptr<image<type,container>>> ifft_vector(std::vector< std::vector<std::shared_ptr<image<type,container>>> > & input);
    
};

// ===========================================
//      Functions of Class metric
// ===========================================

// ===========================================
// Create Functions
// ===========================================
// Constructor
template <typename type, typename container>
regularizer<type,container>::regularizer()
{
    init(1.0, 1.0);
};

template <typename type, typename container>
regularizer<type,container>::regularizer(double a, double g)
{
    init(a, g);
};

template <typename type, typename container>
void regularizer<type,container>::init(double a, double g)
{
    ready = false;
    alpha = a;
    gamma = g;

    this->set_total_inputs(0);      //process_object::init
    this->set_total_outputs(0);     //process_object::init
    this->setup_input();
    this->setup_output();           // no defined output, this object is a base class
};

// ===========================================
// Get Functions
// ===========================================
template <typename type, typename container>
double regularizer<type,container>::get_alpha() const
{
    return alpha;
};

template <typename type, typename container>
double regularizer<type,container>::get_gamma() const
{
    return gamma;
};

// ===========================================
// Set Functions
// ===========================================
template <typename type, typename container>
void regularizer<type,container>::set_alpha(double a)
{
    alpha = a;
};

template <typename type, typename container>
void regularizer<type,container>::set_gamma(double g)
{
    gamma = g;
};

// ===========================================
// Print Functions
// ===========================================
template <typename type, typename container>
std::string regularizer<type,container>::info(std::string msg)
{
    std::stringstream ss;
    std::string title = "Regularizer Information";
    if (msg != "") { title = msg; };

    ss << object::info(title);
    ss << process_object::info("");
    return ss.str();
};

// ===========================================
// Functions
// ===========================================
template <typename type, typename container>
std::vector<std::shared_ptr<image<type,container>>> regularizer<type,container>::update_a(std::shared_ptr< image<type,container> > input)
{
    int dim = input->get_dimension();
    std::shared_ptr<image<type,container>> a = input->mimic();

    int len = 1;
    std::vector<int> sz = a->get_size();
    for (int k = 0; k < sz.size(); k++) len = len*sz[k];    

    auto veca = vector_cpu<type>::new_pointer(len, 0.0);
    type *pa = veca->data();

    if (sz.size() == 2)
    {
        int w = sz[0]; int h = sz[1];

        #pragma omp parallel for collapse(2) if (w*h > IMART_OPENMP_VECTOR_MIN_SIZE)
        for(int j = 0; j < h; j++)
        {
            for(int i = 0; i < w; i++)
            {
                pa[i + j*w] += 2 * alpha * ((1 - cos(2 * M_PI * j / h)) + (1 - cos(2 * M_PI * i / w)));
            };
        };
    }
    else if (sz.size() == 3)
    {
        int w = sz[0]; int h = sz[1]; int l = sz[2];

        #pragma omp parallel for collapse(3) if (w*h*l > IMART_OPENMP_VECTOR_MIN_SIZE)
        for(int k = 0; k < l; k++)
        {
            for(int j = 0; j < h; j++)
            {
                for(int i = 0; i < w; i++)
                {
                    pa[i + j*w + k*w*h] += 2 * alpha * ((1 - cos(2 * M_PI * j / h)) + (1 - cos(2 * M_PI * i / w)) + (1 - cos(2 * M_PI * k / l)));
                };
            };
        };
    };

    auto dat = container::new_pointer(len);
    dat->read_ram(pa, len);
    a->set_data(dat);

    *a = *a + gamma;
    
    this->A = std::vector<std::shared_ptr<image<type,container>>>(dim);

    this->A[0] = a;
    for (int k = 1; k < dim; k++) this->A[k] = a->clone();
    return this->A;
};

template <typename type, typename container>
std::vector<std::shared_ptr<image<type,container>>> regularizer<type,container>::k(std::shared_ptr< std::vector<std::shared_ptr<image<type,container>>> > g)
{
    // A is precomputed, shape [dim, image]

    // transform to fourier domain
    // std::cout << "fft" << std::endl;
    auto G = fft_vector(g); // G shape [dim, real-complex, image]

    // perform operation in fourier space, COMPLEX
    // F = G / (A*A)
    // std::cout << "square" << std::endl;
    std::vector<std::shared_ptr<image<type,container>>> A2(this->A.size());
    for (int k = 0; k < this->A.size(); k++) 
    {
        A2[k] = this->A[k]->mimic();
        *A2[k] = (*this->A[k])*(*this->A[k]);
    }

    // std::cout << "complex division" << std::endl;
    // std::cout << G.size() << std::endl;
    std::vector< std::vector<std::shared_ptr<image<type,container>>> > F(G.size());
    for (int k = 0; k < G.size(); k++)
    {
        auto real = G[k][0]->mimic();
        auto img = G[k][1]->mimic();
        
        *real = (*G[k][0])/(*A2[k]);
        *img  = (*G[k][1])/(*A2[k]);
        F[k] = std::vector<std::shared_ptr<image<type,container>>>{real, img};
    };
    
    // transform back to normal domain
    // std::cout << "ifft" << std::endl;
    auto f = ifft_vector(F);
    return f;
};

template <typename type, typename container>
std::vector<std::shared_ptr<image<type,container>>> regularizer<type,container>::l(std::shared_ptr< std::vector<std::shared_ptr<image<type,container>>> > f)
{
    int dim = f->size();
    std::vector<std::shared_ptr<image<type,container>>> g(dim);

    if (dim == 2)
    {
        auto c = container::new_pointer(std::initializer_list<type>{0.0,1.0,0.0,1.0,-4.0,1.0,0.0,1.0,0.0});
        auto r = image<type,container>::new_pointer(c,3,3);
        
        for (int k = 0; k < dim; k++)
        {
            auto df = convolve(f->at(k), r);
            g[k] = df->mimic();
            *g[k] = (*df)*(-alpha) + (*(f->at(k)))*gamma;
        };
    };

    return g;
    // r = np.array([[0.,  1.,  0.],
    //           [1., -4.,  1.],
    //           [0.,  1.,  0.]])
    // dxx = convolve(f[0, :, :], r)
    // dyy = convolve(f[1, :, :], r)

    // dff = np.array([dxx, dyy])
    // g = - self.alpha * dff + self.gamma * f
};

template <typename type, typename container>
std::vector< std::vector<std::shared_ptr<image<type,container>>> > regularizer<type,container>::fft_vector(std::shared_ptr< std::vector<std::shared_ptr<image<type,container>>> > input)
{
    int n = input->size();
    std::vector< std::vector<std::shared_ptr<image<type,container>>> > output(n);

    for(int k = 0; k < n; k++)
    {
        output[k] = fft(input->at(k));
        // output[k][0]->print_data("real");
        // output[k][1]->print_data("img");
    }
    return output;
};

template <typename type, typename container>
std::vector<std::shared_ptr<image<type,container>>> regularizer<type,container>::ifft_vector(std::vector< std::vector<std::shared_ptr<image<type,container>>> > & input)
{
    int n = input.size();
    std::vector<std::shared_ptr<image<type,container>>> output(n);

    for(int k = 0; k < n; k++)
    {
        // std::cout << k << std::endl;
        output[k] = ifft(input[k]);
    }
    return output;
};

}; //end namespace

#endif