/*
* @Author: jose
* @Date:   2019-11-15 13:55:00
* @Last Modified by:   jose
* @Last Modified time: 2019-11-05 13:55:00
*/

#ifndef __IMAGE_UTILS_H__
#define __IMAGE_UTILS_H__

// std libs
#include <iostream>     // std::cout
#include <memory>       // smart pointers
#include <vector>       // std::vector
#include <cassert>      // assert

// local libs
#include "image.h"
#include "grid.h"

namespace imart
{

// ===========================================
//      Functions Utils
// ===========================================
// template<typename type_, typename container_>
// friend image<type_,container_> normalize(const image<type_,container_> & input, type_ min, type_ max);

// template<typename type_in, typename container_in, typename type_out, typename container_out>
// friend void cast(const image<type_in,container_in> & input, const image<type_out,container_out> & output);

// template<typename type_, typename container_>
// friend image<type_,container_> pad(const image<type_,container_> & input, std::vector<int> pre, std::vector<int> post);

// template<typename type_, typename container_>
// friend image<type_,container_> unpad(const image<type_,container_> & input, std::vector<int> pre, std::vector<int> post);

// // fft
// template<typename type_, typename container_>
// friend typename image<type_,container_>::vector fft(const image<type_,container_> & input);
// template<typename type_, typename container_>
// friend typename image<type_,container_>::vector ifft(const image<type_,container_> & input);
// template<typename type_, typename container_>
// friend typename image<type_,container_>::vector gradient_fft(const image<type_,container_> & input);
// template<typename type_, typename container_>
// friend typename image<type_,container_>::vector gradient(std::shared_ptr<image<type_,container_>> input);
// template<typename type_, typename container_>
// friend typename image<type_,container_>::pointer jacobian(std::shared_ptr<std::vector< std::shared_ptr<image<type_,container_>> >> input);

// friend functions
// template<typename pixel_t>
// friend image<std::complex<pixel_t>> real2complex(const image<pixel_t> & input);
// template<typename pixel_t>
// friend image<pixel_t> complex2real(const image<std::complex<pixel_t>> & input);

template<typename type, typename container>
image<type,container> normalize(const image<type,container> & input, type min = 0.0, type max = 1.0)
{
    typename image<type,container>::pointer output = input.mimic();
    output->set_data((input.get_data()->normalize(min,max)));
    return *output;
};

template<typename type, typename container>
std::shared_ptr<image<type,container>> normalize(std::shared_ptr<image<type,container>> input, type min = 0.0, type max = 1.0)
{
    typename image<type,container>::pointer output = input->mimic();
    output->set_data((input->get_data()->normalize(min,max)));
    return output;
};

template<typename type_in, typename container_in, typename type_out, typename container_out>
void cast(const image<type_in,container_in> & input, image<type_out,container_out> & output)
{
    output = image<type_out,container_out>(input.get_size());
    output.set_data(input.get_data()->template cast<type_out>());
    output.set_sod_parameters(input.get_spacing(), input.get_origin(), input.get_direction());
};

template<typename type_in, typename container_in, typename type_out, typename container_out>
void cast(std::shared_ptr<image<type_in,container_in>> input, std::shared_ptr<image<type_out,container_out>> output)
{
    output = input->mimic();
    output->set_data(input->get_data()->template cast<type_out>());
    output->set_sod_parameters(input->get_spacing(), input->get_origin(), input->get_direction());
};

template<typename type, typename container>
image<type,container> pad(const image<type,container> & input, std::vector<int> pre, std::vector<int> post)
{
    int w = input.get_width();
    int h = input.get_height();
    int l = input.get_length();
    std::vector<int> extra(pre.size(),0);
    for (int i = 0; i < extra.size(); i++){ extra[i] = pre[i]+post[i]; };
    
    typename image<type,container>::pointer output;
    if (input.get_dimension() == 2){ output = image<type,container>::new_pointer(w+extra[0], h+extra[1]); };
    if (input.get_dimension() == 3){ output = image<type,container>::new_pointer(w+extra[0], h+extra[1], l+extra[2]); };
    output->zeros();
    // output->print();

    container::pad(input.get_data(), output->get_data(), input.get_size(), pre, post);
    output->set_sod_parameters(input.get_spacing(), input.get_origin(), input.get_direction());
    return *output;
};

template<typename type, typename container>
std::shared_ptr<image<type,container>> pad(std::shared_ptr<image<type,container>> input, std::vector<int> pre, std::vector<int> post)
{
    int w = input->get_width();
    int h = input->get_height();
    int l = input->get_length();
    std::vector<int> extra(pre.size(),0);
    for (int i = 0; i < extra.size(); i++){ extra[i] = pre[i]+post[i]; };
    
    typename image<type,container>::pointer output;
    if (input->get_dimension() == 2){ output = image<type,container>::new_pointer(w+extra[0], h+extra[1]); };
    if (input->get_dimension() == 3){ output = image<type,container>::new_pointer(w+extra[0], h+extra[1], l+extra[2]); };
    output->zeros();
    // output->print();

    container::pad(input->get_data(), output->get_data(), input->get_size(), pre, post);
    output->set_sod_parameters(input->get_spacing(), input->get_origin(), input->get_direction());
    return output;
};

template<typename type, typename container>
image<type,container> unpad(const image<type,container> & input, std::vector<int> pre, std::vector<int> post)
{
    int w = input.get_width();
    int h = input.get_height();
    int l = input.get_length();
    std::vector<int> extra(pre.size(),0);
    for (int i = 0; i < extra.size(); i++){ extra[i] = pre[i]+post[i]; };
    
    typename image<type,container>::pointer output;
    if (input.get_dimension() == 2){ output = image<type,container>::new_pointer(w-extra[0], h-extra[1]); };
    if (input.get_dimension() == 3){ output = image<type,container>::new_pointer(w-extra[0], h-extra[1], l-extra[2]); };
    
    container::unpad(input.get_data(), output->get_data(), output->get_size(), pre, post);
    output->set_sod_parameters(input.get_spacing(), input.get_origin(), input.get_direction());
    return *output;
};

template<typename type, typename container>
std::shared_ptr<image<type,container>> unpad(std::shared_ptr<image<type,container>> input, std::vector<int> pre, std::vector<int> post)
{
    int w = input->get_width();
    int h = input->get_height();
    int l = input->get_length();
    std::vector<int> extra(pre.size(),0);
    for (int i = 0; i < extra.size(); i++){ extra[i] = pre[i]+post[i]; };
    
    typename image<type,container>::pointer output;
    if (input->get_dimension() == 2){ output = image<type,container>::new_pointer(w-extra[0], h-extra[1]); };
    if (input->get_dimension() == 3){ output = image<type,container>::new_pointer(w-extra[0], h-extra[1], l-extra[2]); };

    container::unpad(input->get_data(), output->get_data(), output->get_size(), pre, post);
    output->set_sod_parameters(input->get_spacing(), input->get_origin(), input->get_direction());
    return output;
};

// template<typename type>
// image<type> complex2real(const image<std::complex<type>> & input)
// {
//     int d = input.get_dimension();
//     int N = input.get_total_elements();
//     std::complex<type> * p1 = input.ptr();

//     image<type> result(input.get_size());
    
//     type * p2 = result.ptr();

//     for (int i = 0; i < N; i++)
//     {
//         p2[i] = std::real(p1[i]);
//     };
//     return result;
// };

template<typename type, typename container>
std::shared_ptr<image<type,container>> gradientx(std::shared_ptr<image<type,container>> input)
{
    auto output = input->mimic();
    container::gradientx(input->get_data(), output->get_data(), input->get_size() );
    *output = (*output) * input->get_spacing()[0] * input->get_direction()[0]; 
    return output;
};

template<typename type, typename container>
std::shared_ptr<image<type,container>> gradienty(std::shared_ptr<image<type,container>> input)
{
    int d = input->get_dimension();
    auto output = input->mimic();
    container::gradienty(input->get_data(), output->get_data(), input->get_size() );
    *output = (*output) * input->get_spacing()[1] * input->get_direction()[1*(d+1)];
    return output;
};

template<typename type, typename container>
std::shared_ptr<image<type,container>> gradientz(std::shared_ptr<image<type,container>> input)
{
    int d = input->get_dimension();
    auto output = input->mimic();
    container::gradientz(input->get_data(), output->get_data(), input->get_size() );
    *output = (*output) * input->get_spacing()[2] * input->get_direction()[2*(d+1)];
    return output;
};

template<typename type, typename container>
typename image<type,container>::vector gradient(std::shared_ptr<image<type,container>> input)
{
    int d = input->get_dimension();
    typename image<type,container>::vector output(d);
    output[0] = gradientx(input);
    output[1] = gradienty(input);
    if (d == 3) output[2] = gradientz(input);
    return output;
};

template<typename type, typename container>
std::shared_ptr<image<type,container>> directional_components(std::shared_ptr<image<type,container>> der1, std::shared_ptr<image<type,container>> der2)
{
    auto output = der1->mimic();
    auto didx = *der1;  // 0 degree
    auto didy = *der2;  // 90 degree
    
    auto did45 = didx*(1/sqrt(2.0)) + didy*(1/sqrt(2.0));
    auto did135 = didx*(-1/sqrt(2.0)) + didy*(1/sqrt(2.0));
    auto did180 = -1.0 * didx;
    auto did225 = didx*(-1/sqrt(2.0)) - didy*(1/sqrt(2.0));
    auto did270 = -1.0 * didy;
    auto did305 = didx*(1/sqrt(2.0)) - didy*(1/sqrt(2.0));
    
    *output = (didx^2.0) + (did45^2.0) + (didy^2.0) + (did135^2.0) + (did180^2.0) + (did225^2.0) + (did270^2.0) + (did305^2.0);
    return output;
};

template<typename type, typename container>
std::shared_ptr<image<type,container>> gradient_direction_gaussian(std::shared_ptr<image<type,container>> input, double sigma = 0.0)
{
    auto output = input->mimic();
    int d = input->get_dimension();
    typename image<type,container>::vector g = gradient(input);

    if (d == 2)
    {
        typename image<type,container>::pointer dsum = directional_components( g[0], g[1] );
        *output = (*dsum)^(0.5);
    }
    else if (d == 3)
    {
        typename image<type,container>::pointer d01 = directional_components( g[0], g[1] );
        typename image<type,container>::pointer d12 = directional_components( g[1], g[2] );
        typename image<type,container>::pointer d20 = directional_components( g[2], g[0] );
        *output = (*d01 + *d12 + *d20)^(0.5);
    }
    else { ; };

    return gaussian_filter(output, sigma, int(3*sigma) );
};

template<typename type, typename container>
typename image<type,container>::pointer jacobian(std::shared_ptr<std::vector< std::shared_ptr<image<type,container>> >> input, bool identity_at_zero = true)
{
    int sz = input->size();
    auto output = image<type,container>::new_pointer();

    auto gx = gradient(input->at(0));
    auto gy = gradient(input->at(1));
    
    auto gx_x = *(gx[0]);
    auto gy_y = *(gy[1]);
    if (identity_at_zero)
    {
        gx_x = gx_x + 1;
        gy_y = gy_y + 1;
    };
    // auto gx_x = *(gx[0]) + 1;
    // auto gy_y = *(gy[1]) + 1;

    // Determinant det_J
    if (sz == 2)
        *output = gx_x*gy_y - ( *(gy[0]) )*( *(gx[1]) );
    else if (sz == 3)
    {
        auto gz = gradient(input->at(2));
        auto gz_z = *(gz[2]);
        if (identity_at_zero) gz_z = gz_z + 1;
        // auto gz_z = *(gz[2]) + 1;

        // Determinant det_J
        *output = gx_x * gy_y * gz_z + \
                  (* gy[0]) * (* gz[1]) * (* gx[2]) + \
                  (* gz[0]) * (* gx[1]) * (* gy[2]) - \
                  (* gz[0]) * gy_y * (* gx[2]) - \
                  (* gy[0]) * (* gx[1]) * gz_z - \
                  gx_x * (* gz[1]) * (* gy[2]);
    };
    return output;
};

template<typename type, typename container>
typename image<type,container>::vector fft(const image<type,container> & input)
{
    // auto in_real = input.copy();
    auto in_img = input.mimic();
    in_img->zeros();
    // input.get_data()->print_data("in fft");
    // in_img->print();

    auto out_real = input.mimic();
    auto out_img = input.mimic();

    typename container::vector vin = {input.get_data(), in_img->get_data()};
    typename container::vector vout = {out_real->get_data(), out_img->get_data()};

    // vin[0]->print_data("vin0"); vin[1]->print_data("vin1");
    container::fft(vin, vout, input.get_size(), true);
    // vout[0]->print_data("vout0"); vout[1]->print_data("vout1");
    
    typename image<type,container>::vector output = {out_real, out_img};
    return output;
};

template<typename type, typename container>
typename image<type,container>::vector fft(std::shared_ptr<image<type,container>> input)
{
    // auto in_real = input.copy();
    // std::cout << "data ";
    auto in_img = input->mimic();
    in_img->zeros();

    auto out_real = input->mimic();
    auto out_img = input->mimic();

    typename container::vector vin = {input->get_data(), in_img->get_data()};
    typename container::vector vout = {out_real->get_data(), out_img->get_data()};

    // vin[0]->print_data("vin0"); vin[1]->print_data("vin1");
    // std::cout << "vin 0:" << (*(vin[0]->get_buffer()))() << std::endl;
    // std::cout << "vin 1:" << (*(vin[1]->get_buffer()))() << std::endl;
    // std::cout << "vin and vout set ";
    container::fft(vin, vout, input->get_size(), true);
    // vout[0]->print_data("vout0"); vout[1]->print_data("vout1");
    // std::cout << "vout 0:" << (*(vout[0]->get_buffer()))() << std::endl;
    // std::cout << "vout 1:" << (*(vout[1]->get_buffer()))() << std::endl;
    
    typename image<type,container>::vector output = {out_real, out_img};
    return output;
};

template<typename type, typename container>
std::shared_ptr< image<type,container> > ifft(const std::vector<std::shared_ptr<image<type,container>>> & input)
{
    const typename image<type,container>::pointer * p = input.data();
    // std::cout << "real" << std::endl;
    auto out_real = p[0]->mimic();
    // std::cout << "img" << std::endl;
    auto out_img = p[1]->mimic();

    typename container::vector vin = {p[0]->get_data(), p[1]->get_data()};
    typename container::vector vout = {out_real->get_data(), out_img->get_data()};

    // vin[0]->print_data("vin0"); vin[1]->print_data("vin1");
    container::fft(vin, vout, p[0]->get_size(), false);
    // vout[0]->print_data("vout0"); vout[1]->print_data("vout1");
    
    // ifft division of total elements
    if(out_real->get_data()->get_name().compare("vector_opencl") != 0)
        *out_real = *out_real/(type)out_real->get_total_elements();
    // *out_img = *out_img/(type)out_img->get_total_elements();
    
    // typename image<type,container>::vector output = {out_real, out_img};
    // return output;
    // return *out_real;
    return out_real;
};

template<typename type, typename container>
typename image<type,container>::vector complex_product(typename image<type,container>::vector a, typename image<type,container>::vector b)
{
    auto real = image<type,container>::new_pointer();
    auto img = image<type,container>::new_pointer();
    *real = (*a[0])*(*b[0]) - (*a[1])*(*b[1]);
    *img = (*a[0])*(*b[1]) + (*a[1])*(*b[0]);
    typename image<type,container>::vector output = {real, img};
    return output;
}


template<typename type, typename container>
typename image<type,container>::vector gradient_fft(const image<type,container> & input)
{
    // padding
    int d = input.get_dimension();
    std::vector<int> none(d);
    std::vector<int> extra(d);
    int a = 0;
    for(int i = 0; i < d; i++)
    {
        none[i] = 0;
        extra[i] = a;
    };
    image<type,container> input_pad = pad(input, none, extra);
    // input_pad.print_data();

    // input fft
    int w = input_pad.get_width();
    int h = input_pad.get_height();
    int l = input_pad.get_length();
    auto fin = fft(input_pad);
    // fin[0]->print_data("fin");
    // fin[1]->print_data("fin");
    input_pad.clear();                      // free memory

    type der[3] = {0.5, 0.0, -0.5};
    typename image<type,container>::vector output(d);
    for(int i = 0; i < d; i++) output[i] = image<type,container>::new_pointer();

    // didx
    auto dx = input_pad.mimic();
    dx->zeros();
    dx->get_data()->read_ram(der, 3, 0);
    // dx->print_data();
    auto didx = ifft(complex_product<type,container>(fin, fft(*dx)));
    dx->clear();
    if(d == 2)          *output[0] = unpad(didx, std::vector<int>{a,0}, std::vector<int>{0,a});
    else if (d == 3)    *output[0] = unpad(didx, std::vector<int>{a,0,0}, std::vector<int>{0,a,a});
    // *output[0] = didx;

    // didy
    auto dy = input_pad.mimic();
    dy->zeros();
    dy->get_data()->read_ram(der, 1, 0);
    dy->get_data()->read_ram(der+1, 1, w);
    dy->get_data()->read_ram(der+2, 1, 2*w);
    auto didy = ifft(complex_product<type,container>(fin, fft(*dy)));
    dy->clear();
    if(d == 2)          *output[1] = unpad(didy, std::vector<int>{0,a}, std::vector<int>{a,0});
    else if (d == 3)    *output[1] = unpad(didy, std::vector<int>{0,a,0}, std::vector<int>{a,0,a});
    // *output[1] = didy;

    if(d == 3)
    {
        auto dz = input_pad.mimic();
        dz->zeros();
        dz->get_data()->read_ram(der, 1, 0);
        dz->get_data()->read_ram(der+1, 1, w*h);
        dz->get_data()->read_ram(der+2, 1, 2*w*h);
        auto didz = ifft(complex_product<type,container>(fin, fft(*dz)));
        dz->clear();
        *output[2] = unpad(didz, std::vector<int>{0,0,2}, std::vector<int>{2,2,0});
        // *output[2] = didz;
    };
    return output;
};

template<typename type, typename container>
typename std::shared_ptr<image<type,container>> convolve(std::shared_ptr<image<type,container>> input, std::shared_ptr<image<type,container>> kernel, bool fft_method = false)
{
    // TODO: select depending on size
    if (fft_method) return fftconvolution(input, kernel);
    return convolution(input, kernel);    
};

template<typename type, typename container>
typename std::shared_ptr<image<type,container>> convolution(std::shared_ptr<image<type,container>> input, std::shared_ptr<image<type,container>> kernel)
{
    auto output = input->mimic();
    output->zeros();
    
    // int kwidth = kernel->get_width(); // kernel image should be squared

    container::convolution(input->get_data(), kernel->get_data(),
                           output->get_data(), input->get_size(), kernel->get_size());
    return output;
};

template<typename type, typename container>
typename std::shared_ptr<image<type,container>> fftconvolution(std::shared_ptr<image<type,container>> input, std::shared_ptr<image<type,container>> kernel)
{
    // padding
    int dim = input->get_dimension();
    std::vector<int> pad_init(dim);
    std::vector<int> pad_extra(dim);
    std::vector<int> pad_kernel(dim);
    std::vector<int> unpad_init(dim);
    std::vector<int> unpad_extra(dim);

    std::vector<int> sz = input->get_size();
    std::vector<int> kw = kernel->get_size();

    // define extra sizes
    for(int i = 0; i < dim; i++)
    {
        pad_init[i] = 0;
        pad_extra[i] = kw[i] >> 1;
        pad_kernel[i] = sz[i] + pad_extra[i] - kw[i];
        unpad_init[i] = pad_extra[i]-1*(1 - kw[i]%2);
        unpad_extra[i] = pad_extra[i] - unpad_init[i];
        // printf("unpad int %d\n", unpad_init[i]);
        // printf("unpad extra %d\n", unpad_extra[i]);
    };
    
    // input padding
    auto input_pad = pad(input, pad_init, pad_extra);
    // input_pad->print_data();

    // kernel padding
    auto kernel_pad = pad(kernel, pad_init, pad_kernel);
    // kernel_pad->print_data();

    auto out_conv = ifft(complex_product<type,container>( fft(input_pad), fft(kernel_pad) ));
    // out_conv->print_data("convolution before unpad");

    // auto output = image<type,container>::new_pointer(dim);
    auto output = unpad(out_conv, unpad_init, unpad_extra);
    // output->print();
    return output;
};

// ===========================================
//      Functions of Gaussian Filter
// ===========================================
template <typename type, typename container>
std::shared_ptr< image<type,container> > gaussian_kernel(int dim = 2, type sigma = 1.0, int width = 3)
{
    assert(width >= 3 and width%2 == 1);
    // Gaussian kernel is small, thats why is computed on cpu
    std::vector<int> sz(dim);
    std::vector<double> origin(dim);
    for(int i = 0; i < dim; i++)
    {
        sz[i] = width;
        origin[i] = -std::floor(width/2.0);
    };
    auto output = image<type,container>::new_pointer(sz);
    output->set_origin(origin);
    auto x = grid<type,container>::new_pointer(output);
    // x->print_data();

    if(dim == 2)
    {
        auto img = ((*x)[0]->operator^(2.0) + (*x)[1]->operator^(2.0))/((type)-2.0*pow(sigma,2.0));
        img = ((type)std::exp(1.0)) ^ ( img );
        *output = img/(img.sum());
    }
    else if(dim == 3)
    {
        auto img = ((*x)[0]->operator^(2.0) + (*x)[1]->operator^(2.0) + (*x)[2]->operator^(2.0))/((type)-2.0*pow(sigma,2.0));
        img = ((type)std::exp(1.0)) ^ ( img );
        *output = img/(img.sum());
    }

    return output;
};

// template<typename type, typename container>
// typename std::shared_ptr<image<type,container>> gaussian_filter(std::shared_ptr<image<type,container>> input, type sigma = 1.0, int kwidth = 3)
// {
//     int d = input->get_dimension();
//     auto output = input->mimic();

//     auto gkernel = gaussian_kernel<type,vector_cpu<type>>(d,sigma,kwidth);

//     if (input->get_data()->get_name() == "vector_cpu")
//     {
//         vector_cpu<type>::convolution(input->get_data(), gkernel->get_data(),
//                            output->get_data(), input->get_size(), kwidth);
//     }
//     else if (input->get_data()->get_name() == "vector_opencl")
//     {
//         // create kernel in gpu (copy from host ram)
//         auto kernel = vector_opencl<type>::new_pointer(gkernel->get_data()->size());
//         kernel->read_ram(gkernel->get_data()->data(),gkernel->get_data()->size());

//         vector_opencl<type>::convolution(input->get_data(), kernel,
//                                output->get_data(), input->get_size(), kwidth);
//     }
//     else ;
//     return output;
// };

template<typename type, typename container>
typename std::shared_ptr<image<type,container>> gaussian_filter(std::shared_ptr<image<type,container>> input, type sigma = 1.0, int kwidth = 3)
{
    if (not (sigma > 0.0)) { return input->copy(); };

    int d = input->get_dimension();
    auto output = input->mimic();
    output->zeros();

    auto gkernel = gaussian_kernel<type,container>(d,sigma,kwidth);
    // gkernel->print_data("kernel");
    // printf("size: %4d\n", input->get_size().size());

    container::convolution(input->get_data(), gkernel->get_data(),
                           output->get_data(), input->get_size(), gkernel->get_size());
    return output;
};

// ===========================================
//      Functions of Vector Image
// ===========================================
template <typename type, typename container>
std::vector< std::shared_ptr< image<type,container> >> operator + (std::vector< std::shared_ptr< image<type,container> >> & input1, std::vector< std::shared_ptr< image<type,container> >> & input2)
{
    assert(input1.size() == input2.size());
    std::vector< std::shared_ptr< image<type,container> >> output(input1.size());

    for(int i = 0; i < input1.size(); i++)
    {
        output[i] = image<type,container>::new_pointer(input1[i]->get_dimension());
        *(output[i]) = *(input1[i]) + *(input2[i]);
    }
    return output;
};

template <typename type, typename container>
std::vector< std::shared_ptr< image<type,container> >> operator - (std::vector< std::shared_ptr< image<type,container> >> & input1, std::vector< std::shared_ptr< image<type,container> >> & input2)
{
    assert(input1.size() == input2.size());
    std::vector< std::shared_ptr< image<type,container> >> output(input1.size());

    for(int i = 0; i < input1.size(); i++)
    {
        output[i] = image<type,container>::new_pointer(input1[i]->get_dimension());
        *(output[i]) = *(input1[i]) - *(input2[i]);
    }
    return output;
};

template <typename type, typename container>
std::vector< std::shared_ptr< image<type,container> >> operator * (std::vector< std::shared_ptr< image<type,container> >> & input1, std::vector< std::shared_ptr< image<type,container> >> & input2)
{
    assert(input1.size() == input2.size());
    std::vector< std::shared_ptr< image<type,container> >> output(input1.size());

    for(int i = 0; i < input1.size(); i++)
    {
        output[i] = image<type,container>::new_pointer(input1[i]->get_dimension());
        *(output[i]) = *(input1[i]) * (*(input2[i]));
    }
    return output;
};

template <typename type, typename container>
std::vector< std::shared_ptr< image<type,container> >> operator / (std::vector< std::shared_ptr< image<type,container> >> & input1, std::vector< std::shared_ptr< image<type,container> >> & input2)
{
    assert(input1.size() == input2.size());
    std::vector< std::shared_ptr< image<type,container> >> output(input1.size());

    for(int i = 0; i < input1.size(); i++)
    {
        output[i] = image<type,container>::new_pointer(input1[i]->get_dimension());
        *(output[i]) = *(input1[i]) / (*(input2[i]));
    }
};

template <typename type, typename container>
std::vector< std::shared_ptr< image<type,container> >> operator * (std::vector< std::shared_ptr< image<type,container> >> & input1, type scalar)
{
    std::vector< std::shared_ptr< image<type,container> >> output(input1.size());

    for(int i = 0; i < input1.size(); i++)
    {
        output[i] = image<type,container>::new_pointer(input1[i]->get_dimension());
        *(output[i]) = (*(input1[i])) * scalar;
    }
    return output;
};

}; //end namespace

#endif