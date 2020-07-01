/*
* @Author: jose
* @Date:   2019-11-19 13:55:00
* @Last Modified by:   jose
* @Last Modified time: 2019-11-19 13:55:00
*/

#ifndef __ILINEAR_H__
#define __ILINEAR_H__

// std libs
#include <iostream>     // std::cout
#include <memory>       // smart pointers
#include <vector>       // std::vector

// images 
#include "interpolator.h"

namespace imart
{

// Class linear interpolator
template <typename type, typename container>
class ilinear: public inherit<ilinear<type,container>, interpolator<type,container>>
{
public:
    //Type definitions
    using self    = ilinear;
    using pointer = std::shared_ptr<self>;
    using vector  = std::vector<self::pointer>;

    using interpolator<type,container>::_default_;
    using interpolator<type,container>::region;
    using interpolator<type,container>::x_reference;
    using interpolator<type,container>::image_reference;

    using inherit<ilinear<type,container>, 
                  interpolator<type,container>>::inherit;

protected:
    
    // ===========================================
    // Functions
    // ===========================================
    void boundaries2();
    void boundaries3();
    int search_(type x, type y);
    int search_(type x, type y, type z);

    type linear_(type r[2], type x); // unitary distance
    type bilinear_(type r[4], type x, type y);
    type trilinear_(type r[8], type x, type y, type z);

    typename image<type,container>::pointer linear2(const typename grid<type,container>::pointer xout);
    typename image<type,container>::pointer linear3(const typename grid<type,container>::pointer xout);

public:
    // ===========================================
    // Create Functions
    // ===========================================
    ilinear(typename image<type,container>::pointer imgref, typename grid<type,container>::pointer xref);

    // ===========================================
    // Functions
    // ===========================================
    typename image<type,container>::pointer apply(const typename grid<type,container>::pointer xout);
};


// ===========================================
//      Functions of Class grid
// ===========================================

// ===========================================
// Constructor
// ===========================================
template <typename type, typename container>
ilinear<type,container>::ilinear(typename image<type,container>::pointer imgref, typename grid<type,container>::pointer xref)
{
    assert(imgref->get_dimension() == xref->get_dimension());
    
    this->class_name = "linear interpolator";
    int d = imgref->get_dimension();
    init(d);
    image_reference = imgref;
    x_reference = xref;
};

// ===========================================
// Functions
// ===========================================
template <typename type, typename container>
typename image<type,container>::pointer ilinear<type,container>::apply(const typename grid<type,container>::pointer xout)
{
    auto a = image<type,container>::new_pointer(xout->get_dimension());
    if(xout->get_dimension() == 2){ a = linear2(xout); };
    if(xout->get_dimension() == 3){ a = linear3(xout); };
    return a;
};


template <typename type, typename container>
typename image<type,container>::pointer ilinear<type,container>::linear2(const typename grid<type,container>::pointer xout)
{
    assert(xout->get_dimension() == image_reference->get_dimension());

    int w = image_reference->get_size()[0];
    int h = image_reference->get_size()[1];

    auto image_out = image<type,container>::new_pointer(w, h);
    image_out->set_sod_parameters(xout->get_spacing(), xout->get_origin(), xout->get_direction());
    image_out->assign(_default_);

    typename image<type,container>::pointer xo = (xout->ptr()[0]);   // raw pointer to
    typename image<type,container>::pointer yo = (xout->ptr()[1]);   // input grid coordinates

    // Future step
    // convert xout to x_reference

    container::linear2(xo->get_data(),yo->get_data(),image_reference->get_data(),
                image_out->get_data(), image_reference->get_size(), xout->get_size() );
    // image_out->set_data(imgo);

    return image_out;
};

template <typename type, typename container>
typename image<type,container>::pointer ilinear<type,container>::linear3(const typename grid<type,container>::pointer xout)
{
    assert(xout->get_dimension() == 3);
    assert(xout->get_dimension() == x_reference->get_dimension());

    // using pixels8 = std::unique_ptr<std::array<type,8>>;

    auto image_out = image<type,container>::new_pointer(xout->get_size()[0],xout->get_size()[1],xout->get_size()[2]);
    // type * pout = image_out->ptr();

    // type * xi = (xout->ptr()[0])->ptr();   // raw pointer to
    // type * yi = (xout->ptr()[1])->ptr();   // input grid coordinates
    // type * zi = (xout->ptr()[2])->ptr();

    // type * xr = (x_reference->ptr()[0])->ptr();  // raw pointer to
    // type * yr = (x_reference->ptr()[1])->ptr();  // ref grid coordinates
    // type * zr = (x_reference->ptr()[2])->ptr();

    // double sx = image_reference->get_spacing()[0];   // scales
    // double sy = image_reference->get_spacing()[1];
    // double sz = image_reference->get_spacing()[2];
    
    // boundaries3(); // call to compute region used in search_
    // int num = xout->ptr()[0]->get_total_elements();
    // for(int k = 0; k < num; k++)
    // {
    //     type xk = xi[k];
    //     type yk = yi[k];
    //     type zk = zi[k];
    //     int idx = search_( xk, yk, zk );
        
    //     type iout = _default_;
        
    //     if (idx != -1)
    //     {
    //         type dx = (xk - xr[idx])/ sx;
    //         type dy = (yk - yr[idx])/ sy;
    //         type dz = (zk - zr[idx])/ sz;
            
    //         pixels8 i8 = image_reference->neighbors8(idx);
    //         iout = trilinear_((*i8).data(), dx, dy, dz);
    //     };

    //     pout[idx] = iout; // TODO: THIS line produce double error
    //                         // check everything in the for loop
    // }
    return image_out;

};

template <typename type, typename container>
void ilinear<type,container>::boundaries2()
{
    int w = x_reference->get_size()[0];
    int h = x_reference->get_size()[1];

    type * xr = (x_reference->ptr()[0])->ptr();  // raw pointer to
    type * yr = (x_reference->ptr()[1])->ptr();  // ref grid coordinates

    region.xmin = xr[0];
    region.xmax = xr[w-1];

    region.ymin = yr[0];
    region.ymax = yr[w*h-1];
}

template <typename type, typename container>
void ilinear<type,container>::boundaries3()
{
    boundaries2();

    int w = x_reference->get_size()[0];
    int h = x_reference->get_size()[1];
    int l = x_reference->get_size()[2];

    type * zr = (x_reference->ptr()[2])->ptr();
    region.zmin = zr[0];
    region.zmax = zr[w*h*l-1];
}

template <typename type, typename container>
int ilinear<type,container>::search_(type x, type y)
{
    // This algorithm assumes a regular fixed grid
    int idx = 0;
    int w = x_reference->get_size()[0];
    // int h = x_reference->get_size()[1];
    // std::cout << "w,h: " << w << " , " << h << "\n";
    // x_reference->print_data();

    // type * xr = (x_reference->ptr()[0]).ptr();  // raw pointer to
    // type * yr = (x_reference->ptr()[1]).ptr();  // ref grid coordinates

    std::vector<double> ss = (x_reference->ptr()[0])->get_spacing();
    double sx = ss[0];
    double sy = ss[1];
    
    // type xmin = xr[0];
    // type xmax = xr[w-1];

    // type ymin = yr[0];
    // type ymax = yr[w*h-1];

    // std::cout << "xmin,xmax: " << xmin << " , " << xmax << "\n";
    // std::cout << "ymin,ymax: " << ymin << " , " << ymax << "\n";
    if(x >= region.xmin and x <= region.xmax and y >= region.ymin and y <= region.ymax)
    {
        int posx = floor((x - region.xmin)/sx);
        int posy = floor((y - region.ymin)/sy);
        idx = posx + posy*w;
    }
    else { idx = -1; };

    return idx;

    // upper left closest intensity value;
    // return closest_image_intensity_index;
};

template <typename type, typename container>
int ilinear<type,container>::search_(type x, type y, type z)
{
    // This algorithm assumes a regular fixed grid
    int idx = 0;
    int w = x_reference->get_size()[0];
    int h = x_reference->get_size()[1];

    std::vector<double> ss = (x_reference->ptr()[0])->get_spacing();
    type sx = ss[0];
    type sy = ss[1];
    type sz = ss[2];

    if(x >= region.xmin and x < region.xmax 
        and y >= region.ymin and y < region.ymax
        and z >= region.zmin and z < region.zmax)
    {
        int posx = floor((x - region.xmin)/sx);
        int posy = floor((y - region.ymin)/sy);
        int posz = floor((z - region.zmin)/sz);

        idx = posx + posy*w + posz*w*h;
    }
    else { idx = -1; };

    return idx;
    // upper left closest intensity value;
    // return closest_image_intensity_index;
};

template <typename type, typename container>
type ilinear<type,container>::linear_(type r[2], type x)
{
    return r[0]*(1-x) + r[1]*(x);
};

template <typename type, typename container>
type ilinear<type,container>::bilinear_(type r[4], type x, type y)
{
    type out[2];
    out[0] = linear_(r, x);
    out[1] = linear_(r+2, x);
    return linear_(out, y);
};

template <typename type, typename container>
type ilinear<type,container>::trilinear_(type r[8], type x, type y, type z)
{
    type out[2];
    out[0] = bilinear_(r, x, y);
    out[1] = bilinear_(r+4, x, y);
    return linear_(out, z);
};

}; //end namespace

#endif