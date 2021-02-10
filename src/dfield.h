/*
* @Author: jose
* @Date:   2019-11-15 13:55:00
* @Last Modified by:   jose
* @Last Modified time: 2019-11-05 13:55:00
*/

#ifndef __DFIELD_H__
#define __DFIELD_H__

// itk headers
#include <itkNiftiImageIO.h>

// local libs
#include "transform.h"

namespace imart
{

// Class dfield
template <typename type, typename container=vector_cpu<type>>
class dfield: public inherit<dfield<type,container>, transform<type,container>>
{
public:
    //Type definitions
    using self    = dfield;
    using pointer = std::shared_ptr<self>;
    using vector  = std::vector<self::pointer>;

    // Inherited variables
    using transform<type,container>::parameters;
    using transform<type,container>::inverse_parameters;
    using transform<type,container>::allocate;
    using transform<type,container>::get_parameters;
    using transform<type,container>::get_inverse_parameters;
    using transform<type,container>::get_parameters_vector;
    using transform<type,container>::get_inverse_parameters_vector;
    using transform<type,container>::sigma;
    using transform<type,container>::kernel;

    // using transform<type,container>::operator=;
    // using transform<type,container>::operator+;

    using inherit<dfield<type,container>, transform<type,container>>::inherit;

protected:
    // ===========================================
    // Functions
    // ===========================================
    void init(int d);
    void init(std::vector<int> sz);
    void init(const space_object & input);
    void inverse_();
    void init_sigma();

    std::vector<type> transform_2d(std::vector<type> & point);
    std::vector<type> transform_3d(std::vector<type> & point);
    
    typename grid<type,container>::pointer transform_2d(const typename grid<type,container>::pointer input);
    typename grid<type,container>::pointer transform_3d(const typename grid<type,container>::pointer input);

    void read_2d(std::string file_name);
    void read_3d(std::string file_name);
    void write_2d(std::string file_name);
    void write_3d(std::string file_name);

public:
    // ===========================================
    // Constructor Functions
    // ===========================================
    dfield();
    dfield(int d);
    dfield(int d, typename transform<type,container>::ptr_vector_image params);
    dfield(std::vector<int> sz);
    dfield(typename image<type,container>::pointer input);
    dfield(typename grid<type,container>::pointer input);
    dfield(const dfield<type,container> & input);

    // ===========================================
    // Get Functions
    // ===========================================
    double get_sigma_elastic() const;
    double get_sigma_fluid() const;
    int get_kernel_elastic() const;
    int get_kernel_fluid() const;

    // ===========================================
    // Set Functions
    // ===========================================
    void set_sigma_elastic(double s);
    void set_sigma_fluid(double s);
    void set_kernel_elastic(int k);   // automatically assign based on sigma, this is manual modification
    void set_kernel_fluid(int k);       // automatically assign based on sigma, this is manual modification

    // ===========================================
    // Functions
    // ===========================================
    void change_size(std::vector<int> sz); // special function for multiresolution
    void identity();
    std::vector<type> apply(std::vector<type> & point);
    typename grid<type,container>::pointer apply(const typename grid<type,container>::pointer input);

    void fluid();
    void elastic();

    void read(std::string file_name);
    void write(std::string file_name);
};

template<typename type>
using dfield_cpu = dfield<type,vector_cpu<type>>;

// template<typename type>
// using dfield_gpu = dfield<type,vector_ocl<type>>;

#ifdef IMART_WITH_OPENCL
template<typename type>
using dfield_ocl = dfield<type,vector_ocl<type>>;
#endif

#ifdef IMART_WITH_CUDA
template<typename type>
using dfield_cuda = dfield<type,vector_cuda<type>>;
#endif


// ===========================================
//          Functions of Class dfield
// ===========================================

// ===========================================
// Constructor Functions
// ===========================================
template <typename type, typename container>
dfield<type,container>::dfield()
{
    this->class_name = "dfield";
    init(2);
    identity();
};

template <typename type, typename container>
dfield<type,container>::dfield(int d)
{
    // std::cout << "dfield" << std::endl;
    this->class_name = "dfield";
    init(d);
    identity();
    // std::cout << "created" << std::endl;
};

template <typename type, typename container>
dfield<type,container>::dfield(std::vector<int> sz)
{
    this->class_name = "dfield";
    init(sz.size());
    init(sz);
    identity();
};

template <typename type, typename container>
dfield<type,container>::dfield(int d, typename transform<type,container>::ptr_vector_image params)
{
    // assert(params->get_total_elements()==(d*d+d));
    this->class_name = "dfield";
    init(d);
    parameters = params;
    inverse_();
};

template <typename type, typename container>
dfield<type,container>::dfield(typename image<type,container>::pointer input)
{
    this->class_name = "dfield";
    init(*input); 
    identity();
};

template <typename type, typename container>
dfield<type,container>::dfield(typename grid<type,container>::pointer input)
{
    this->class_name = "dfield";
    init(*input); 
    identity();
};

template <typename type, typename container>
dfield<type,container>::dfield(const dfield<type,container> & input)
{
    this->class_name = "dfield";
    this->clone_(input);
};

template <typename type, typename container>
void dfield<type,container>::init(int d)
{
    init_sigma();
    space_object::init(d);
    allocate(d);
};

template <typename type, typename container>
void dfield<type,container>::init(std::vector<int> sz)
{
    init_sigma();
    // initialize size
    this->set_size(sz);
    // allocated parameters, change the size
    for(int i = 0; i < this->get_dimension(); i++)
    {
        (*parameters)[i] = image<type,container>::new_pointer(sz);        // parameters are 2d
        (*inverse_parameters)[i] = image<type,container>::new_pointer(sz);// parameters are 2d
    };
};

template <typename type, typename container>
void dfield<type,container>::init(const space_object & input)
{
    this->copy_space(input);
    allocate(input.get_dimension());
    init(input.get_size());
};

template <typename type, typename container>
void dfield<type,container>::init_sigma()
{
    double sigma_elastic = 1.0;
    int kernel_elastic = int(ceil(3*sigma_elastic));
    if (kernel_elastic%2 == 0) kernel_elastic -= 1;

    double sigma_fluid = 0.0;
    int kernel_fluid = int(ceil(3*sigma_fluid));
    if (kernel_fluid%2 == 0) kernel_fluid -= 1;

    sigma = {sigma_elastic, sigma_fluid};
    kernel = {kernel_elastic, kernel_fluid};
};

// ===========================================
// Get Functions
// ===========================================
template <typename type, typename container>
double dfield<type,container>::get_sigma_elastic() const
{
    return sigma[0];
};

template <typename type, typename container>
int dfield<type,container>::get_kernel_elastic() const
{
    return kernel[0];
};

template <typename type, typename container>
double dfield<type,container>::get_sigma_fluid() const
{
    return sigma[1];
};

template <typename type, typename container>
int dfield<type,container>::get_kernel_fluid() const
{
    return kernel[1];
};

// ===========================================
// Set Functions
// ===========================================
template <typename type, typename container>
void dfield<type,container>::set_sigma_elastic(double s)
{
    double sigma_elastic = s;
    int kernel_elastic = int(ceil(3*sigma_elastic));
    if (kernel_elastic%2 == 0) kernel_elastic -= 1;
    sigma[0] = sigma_elastic;
    kernel[0] = kernel_elastic;
};

template <typename type, typename container>
void dfield<type,container>::set_kernel_elastic(int k)
{
    int kernel_elastic = k;
    kernel[0] = kernel_elastic;
};

template <typename type, typename container>
void dfield<type,container>::set_sigma_fluid(double s)
{
    double sigma_fluid = s;
    int kernel_fluid = int(ceil(3*sigma_fluid));
    if (kernel_fluid%2 == 0) kernel_fluid -= 1;
    sigma[1] = sigma_fluid;
    kernel[1] = kernel_fluid;
};

template <typename type, typename container>
void dfield<type,container>::set_kernel_fluid(int k)
{
    int kernel_fluid = k;
    kernel[1] = kernel_fluid;
};

// ===========================================
// Functions
// ===========================================
template <typename type, typename container>
void dfield<type,container>::change_size(std::vector<int> sz)
{
    space_object::set_size(sz);
};

template <typename type, typename container>
void dfield<type,container>::identity()
{
    // std::cout << "init identity" << std::endl;
    for (int i = 0; i < this->dim; i++)
    {
        // std::cout << "zeros" << std::endl;
        // get_parameters(i)->print();
        get_parameters(i)->zeros();
        // std::cout << "inv" << std::endl;
        get_inverse_parameters(i)->zeros();
    };
    // std::cout << "end identity" << std::endl;
};

template <typename type, typename container>
void dfield<type,container>::inverse_()
{
    *inverse_parameters = (*parameters)*((type)-1.0);
};

//Transform point
template <typename type, typename container>
std::vector<type> dfield<type,container>::apply(std::vector<type> & point)
{
    if (this->dim == 2){ return transform_2d(point); }
    else if (this->dim == 3){ return transform_3d(point); };
    return point;
};

//Transform grid
template <typename type, typename container>
typename grid<type,container>::pointer dfield<type,container>::apply(typename grid<type,container>::pointer input)
{
    if (this->dim == 2){ return transform_2d(input); }
    else if (this->dim == 3){ return transform_3d(input); };
    return input;
};

// ===========================================
// Functions 2d
// ===========================================
// Transform point
template <typename type, typename container>
std::vector<type> dfield<type,container>::transform_2d(std::vector<type> & point)
{
    // TODO: consider if the point uses the inverse or the direct transform. I think direct.
    assert(this->get_dimension() == point.size());
    std::vector<type> out(point.size());
    
    // interpolation required!

    // std::vector<type> va = parameters->get_data()->std_vector();
    // type * a = va.data();
    // out[0] = a[0]*point[0] + a[1]*point[1] + a[4];
    // out[1] = a[2]*point[0] + a[3]*point[1] + a[5];
    return out;
};

//Transform grid
template <typename type, typename container>
typename grid<type,container>::pointer dfield<type,container>::transform_2d(typename grid<type,container>::pointer input)
{
    assert(this->get_dimension() == input->get_dimension());
    auto output = input->mimic();

    typename image<type,container>::pointer * xin = input->ptr();
    typename image<type,container>::pointer * xout = output->ptr();

    container::dfield_2d(xin[0]->get_data(), xin[1]->get_data(),
                         get_parameters(0)->get_data(), get_parameters(1)->get_data(),
                         xout[0]->get_data(), xout[1]->get_data());
    return output;
};

// ===========================================
// Functions 3d
// ===========================================
// Transform point
template <typename type, typename container>
std::vector<type> dfield<type,container>::transform_3d(std::vector<type> & point)
{
    // TODO: consider if the point uses the inverse or the direct transform. I think direct.
    assert(this->get_dimension() == point.size());
    std::vector<type> out(point.size());

    // interpolation required!

    // std::vector<type> va = parameters->get_data()->std_vector();
    // type * a = va.data();
    // out[0] = a[0]*point[0] + a[1]*point[1] + a[2]*point[2] + a[9];
    // out[1] = a[3]*point[0] + a[4]*point[1] + a[5]*point[2] + a[10];
    // out[2] = a[6]*point[0] + a[7]*point[1] + a[8]*point[2] + a[11];
    return out;
};

//Transform grid
template <typename type, typename container>
typename grid<type,container>::pointer dfield<type,container>::transform_3d(typename grid<type,container>::pointer input)
{
    assert(this->get_dimension() == input->get_dimension());
    auto output = input->mimic();
    
    typename image<type,container>::pointer * xin = input->ptr();
    typename image<type,container>::pointer * xout = output->ptr();

    container::dfield_3d(xin[0]->get_data(), xin[1]->get_data(), xin[2]->get_data(),
                         get_parameters(0)->get_data(), get_parameters(1)->get_data(), get_parameters(2)->get_data(),
                         xout[0]->get_data(), xout[1]->get_data(), xout[2]->get_data());
    return output;
};

template <typename type, typename container>
void dfield<type,container>::fluid()
{
    // printf("dim = %i\n", this->get_dimension());
    type sigma = get_sigma_fluid();
    int kwidth = get_kernel_fluid();
    if (sigma > 0 )
    {
        // printf("fluid s = %.1f, k = %i\n", sigma, kwidth);
        for(int i = 0; i < this->get_dimension(); i++)
        {
            auto ps = gaussian_filter(this->get_parameters(i),sigma,kwidth);
            this->set_parameters( ps, i );
            // this->set_parameters( gaussian_filter(this->get_parameters(i),sigma,kwidth), i );
        };
    };
};

template <typename type, typename container>
void dfield<type,container>::elastic()
{
    type sigma = get_sigma_elastic();
    int kwidth = get_kernel_elastic();
    if (sigma > 0 )
    {
        // printf("elastic s = %.1f, k = %i\n", sigma, kwidth);
        for(int i = 0; i < this->get_dimension(); i++)
            this->set_parameters( gaussian_filter(this->get_parameters(i),sigma,kwidth), i );
    };
};


template <typename type, typename container>
void dfield<type,container>::read(std::string file_name)
{
    if (this->get_dimension() == 2) { read_2d(file_name); }
    else if (this->get_dimension() == 3) { read_3d(file_name); };
};

template <typename type, typename container>
void dfield<type,container>::read_2d(std::string file_name)
{
    const int d = 2;
    using PixelType = itk::FixedArray<type, d>;
    using itkImageType = itk::Image<PixelType, d>;
    using itkReaderType = itk::ImageFileReader<itkImageType>;
    
    // Initialize itk objects
    typename itkImageType::Pointer image_itk = itkImageType::New();
    typename itkReaderType::Pointer reader = itkReaderType::New();

    // Set the image filename itk
    reader->SetFileName(file_name);
    reader->Update();

    // Read the image from reader
    image_itk = reader->GetOutput();
    // std::cout << image_itk;

    typename itkImageType::RegionType region = image_itk->GetLargestPossibleRegion();
    typename itkImageType::SizeType itksize = region.GetSize();
    PixelType * pitk = image_itk->GetBufferPointer();
    int w = itksize[0];
    int h = itksize[1];

    // Copy all data to std::vector
    std::vector<type> raw(d*w*h);
    type * p = raw.data();
    for(size_t i = 0; i < w*h; i++)
    {
        for (size_t j = 0; j < d; j++)
        {
            p[i+j*w*h] = pitk[i][j];
        };
    };

    // Initilize space
    init(std::vector<int>{w, h});
    int num_elements = w*h;

    // Copy of the image
    (*parameters)[0]->get_data()->read_ram(p, num_elements);        // first channel
    (*parameters)[1]->get_data()->read_ram(p + w*h, num_elements);  // second channel

    typename itkImageType::SpacingType spacing_itk = image_itk->GetSpacing();
    typename itkImageType::PointType origin_itk = image_itk->GetOrigin();
    typename itkImageType::DirectionType direction_itk = image_itk->GetDirection();

    std::vector<double> spacing(d,0), origin(d,0), direction(d*d,0);
    for (int i = 0; i < d; i++)
    {
        spacing[i] = spacing_itk[i];
        origin[i] = origin_itk[i];
    };
    int c = 0;
    for (int i = 0; i < d; i++)
    {
        for (int j = 0; j < d; j++)
        {
            direction[c] = direction_itk[i][j];
            c++;
        };
    };

    this->set_spacing(spacing);
    this->set_origin(origin);
    this->set_direction(direction);
};

template <typename type, typename container>
void dfield<type,container>::read_3d(std::string file_name)
{
    const int d = 3;
    using PixelType = itk::RGBPixel<type>;
    using itkImageType = itk::Image<PixelType, d>;
    using itkReaderType = itk::ImageFileReader<itkImageType>;
    
    // Initialize itk objects
    typename itkImageType::Pointer image_itk = itkImageType::New();
    typename itkReaderType::Pointer reader = itkReaderType::New();

    // Set the image filename itk
    reader->SetFileName(file_name);
    reader->Update();

    // Read the image from reader
    image_itk = reader->GetOutput();
    // std::cout << image_itk;

    typename itkImageType::RegionType region = image_itk->GetLargestPossibleRegion();
    typename itkImageType::SizeType itksize = region.GetSize();
    PixelType * pitk = image_itk->GetBufferPointer();
    int w = itksize[0];
    int h = itksize[1];
    int l = itksize[2];

    // Copy all data to std::vector
    std::vector<type> raw(d*w*h*l);
    type * p = raw.data();
    for(size_t i = 0; i < w*h*l; i++)
    {
        for (size_t j = 0; j < d; j++)
        {
            p[i+j*w*h*l] = pitk[i][j];
            // printf("[%d,%d] = %f",i,j, p[i+j*w*h*l]);
        };
    };

    // for(int k = 0; k < l; k++)
    // {
    //     for(int j = 0; j < h; j++)
    //     {
    //         for (int i = 0; i < w; i++)
    //         {
    //             typename itkImageType::IndexType index = {{i,j,k}};
    //             auto px = image_itk->GetPixel(index);
    //             raw[i+j*w+k*w*h]            = px[0];
    //             raw[i+j*w+k*w*h + w*h*l]    = px[1];
    //             raw[i+j*w+k*w*h + 2*w*h*l]  = px[2];
    //         };
    //     };
    // };

    // Initilize space
    init(std::vector<int>{w, h, l});
    int num_elements = w*h*l;

    // Copy of the image
    (*parameters)[0]->get_data()->read_ram(p, num_elements);            // first channel
    (*parameters)[1]->get_data()->read_ram(p + w*h*l, num_elements);    // second channel
    (*parameters)[2]->get_data()->read_ram(p + 2*w*h*l, num_elements);  // second channel

    typename itkImageType::SpacingType spacing_itk = image_itk->GetSpacing();
    typename itkImageType::PointType origin_itk = image_itk->GetOrigin();
    typename itkImageType::DirectionType direction_itk = image_itk->GetDirection();

    std::vector<double> spacing(d,0), origin(d,0), direction(d*d,0);
    for (int i = 0; i < d; i++)
    {
        spacing[i] = spacing_itk[i];
        origin[i] = origin_itk[i];
    };
    int c = 0;
    for (int i = 0; i < d; i++)
    {
        for (int j = 0; j < d; j++)
        {
            direction[c] = direction_itk[i][j];
            c++;
        };
    };

    this->set_spacing(spacing);
    this->set_origin(origin);
    this->set_direction(direction);
};

template <typename type, typename container>
void dfield<type,container>::write(std::string file_name)
{
    if (this->get_dimension() == 2) { write_2d(file_name); }
    else if (this->get_dimension() == 3) { write_3d(file_name); };
};

template <typename type, typename container>
void dfield<type,container>::write_2d(std::string file_name)
{
    const int d = 2;
    using PixelType = itk::FixedArray<type, d>;
    using itkImageType = itk::Image<PixelType, d>;
    using itkWriterType = itk::ImageFileWriter<itkImageType>;

    typename itkImageType::Pointer image_itk = itkImageType::New();
    typename itkWriterType::Pointer writer = itkWriterType::New();

    typename itkImageType::RegionType region;
    typename itkImageType::IndexType  start;
    typename itkImageType::SizeType size;
    
    for (int i = 0; i < d; i++)
    {
        start[i] = 0;
        size[i] = this->get_size()[i];
    };

    int w = this->get_size()[0];
    int h = this->get_size()[1];
    int num_elements = w*h;

    region.SetSize(size);
    region.SetIndex(start);

    image_itk->SetRegions(region);
    image_itk->Allocate();

    PixelType * pitk = image_itk->GetBufferPointer();

    // Copy dfield parameters to std::vector
    std::vector<type> raw(d*w*h);
    type * p = raw.data();
    for (int i = 0; i < d; i++)
        (*parameters)[i]->get_data()->write_ram(p + i*w*h, num_elements);
    
    // Copy std::vector to itk image
    for(size_t i = 0; i < w*h; i++)
    {
        for (size_t j = 0; j < d; j++)
        {
            pitk[i][j] = p[i+j*w*h];
            // printf("[%d,%d] = %f ",i,j, p[i+j*w*h]);
        };
    };

    // for(size_t j = 0; j < h; j++)
    // {
    //     for (size_t i = 0; i < w; i++)
    //     {
    //         std::array<type,d> a = {raw[i+j*w*h],raw[i+j*w*h + w*h]};
    //         PixelType px(a);
    //         typename itkImageType::IndexType index = {{i,j}};
    //         image_itk->SetPixel(index, px );
    //     };
    // };

    // Writing metadata
    std::vector<double> spacing = this->get_spacing();
    std::vector<double> origin = this->get_origin();
    std::vector<double> direction = this->get_direction();
    
    typename itkImageType::SpacingType spacing_itk;// = image_itk->GetSpacing();
    typename itkImageType::PointType origin_itk;// = image_itk->GetOrigin();
    typename itkImageType::DirectionType direction_itk;// = image_itk->GetDirection();
    
    for (int i = 0; i < d; i++)
    {
        spacing_itk[i] = spacing[i];
        origin_itk[i] = origin[i];
        // std::cout << origin_itk[i];
    };
    int c = 0;
    for (int i = 0; i < d; i++)
    {
        for (int j = 0; j < d; j++)
        {
            direction_itk[i][j] = direction[c];
            c++;
        };
    };

    image_itk->SetSpacing(spacing_itk);
    image_itk->SetOrigin(origin_itk);
    image_itk->SetDirection(direction_itk);

    writer->SetFileName(file_name);
    writer->SetInput(image_itk);

    try
    {
        writer->Update();
    }
    catch (itk::ExceptionObject & error)
    {
        std::cerr << "Error: " << error << std::endl;
    }
};

template <typename type, typename container>
void dfield<type,container>::write_3d(std::string file_name)
{
    const int d = 3;
    using PixelType = itk::FixedArray<type, d>;
    using itkImageType = itk::Image<PixelType, d>;
    using itkWriterType = itk::ImageFileWriter<itkImageType>;

    typename itkImageType::Pointer image_itk = itkImageType::New();
    typename itkWriterType::Pointer writer = itkWriterType::New();

    typename itkImageType::RegionType region;
    typename itkImageType::IndexType  start;
    typename itkImageType::SizeType size;
    
    for (int i = 0; i < d; i++)
    {
        start[i] = 0;
        size[i] = this->get_size()[i];
    };

    int w = this->get_size()[0];
    int h = this->get_size()[1];
    int l = this->get_size()[2];
    int num_elements = w*h*l;

    region.SetSize(size);
    region.SetIndex(start);

    image_itk->SetRegions(region);
    image_itk->Allocate();

    PixelType * pitk = image_itk->GetBufferPointer();

    // Copy dfield parameters to std::vector
    std::vector<type> raw(d*w*h*l);
    type * p = raw.data();
    for (int i = 0; i < d; i++)
        (*parameters)[i]->get_data()->write_ram(p + i*w*h*l, num_elements);

    // Copy std::vector to itk image
    for(size_t i = 0; i < w*h*l; i++)
    {
        for (size_t j = 0; j < d; j++)
        {
            pitk[i][j] = p[i+j*w*h*l];
            // printf("[%d,%d] = %5.1f",i,j, p[i+j*w*h*l]);
        };
    };

    // for(int k = 0; k < l; k++)
    // {
    //     for(int j = 0; j < h; j++)
    //     {
    //         for (int i = 0; i < w; i++)
    //         {
    //             PixelType px;
    //             std::array<type,d> a = {raw[i+j*w+k*w*h],raw[i+j*w+k*w*h + w*h*l],raw[i+j*w+k*w*h + 2*w*h*l]};
                
    //             px[0] = a[0];
    //             px[1] = a[1];
    //             px[2] = a[2];
    //             typename itkImageType::IndexType index = {{i,j,k}};
    //             image_itk->SetPixel(index, px );
    //         };
    //     };
    // };

    // Writing metadata
    std::vector<double> spacing = this->get_spacing();
    std::vector<double> origin = this->get_origin();
    std::vector<double> direction = this->get_direction();
    
    typename itkImageType::SpacingType spacing_itk;// = image_itk->GetSpacing();
    typename itkImageType::PointType origin_itk;// = image_itk->GetOrigin();
    typename itkImageType::DirectionType direction_itk;// = image_itk->GetDirection();
    
    for (int i = 0; i < d; i++)
    {
        spacing_itk[i] = spacing[i];
        origin_itk[i] = origin[i];
        // std::cout << origin_itk[i];
    };
    int c = 0;
    for (int i = 0; i < d; i++)
    {
        for (int j = 0; j < d; j++)
        {
            direction_itk[i][j] = direction[c];
            c++;
        };
    };

    image_itk->SetSpacing(spacing_itk);
    image_itk->SetOrigin(origin_itk);
    image_itk->SetDirection(direction_itk);

    writer->SetFileName(file_name);
    writer->SetInput(image_itk);

    try
    {
        writer->Update();
    }
    catch (itk::ExceptionObject & error)
    {
        std::cerr << "Error: " << error << std::endl;
    }
};

}; //end namespace

#endif