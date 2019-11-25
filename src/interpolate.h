/*
* @Author: jose
* @Date:   2019-11-19 13:55:00
* @Last Modified by:   jose
* @Last Modified time: 2019-11-19 13:55:00
*/

#ifndef __INTERPOLATE_H__
#define __INTERPOLATE_H__

// std libs
#include <iostream>     // std::cout
#include <memory>       // smart pointers
#include <vector>       // std::vector
#include <cassert>      // assert

// images 
#include "image_base.h"
#include "grid.h"


template <typename pixel_type>
struct borders
{
    pixel_type xmin;
    pixel_type xmax;
    pixel_type ymin;
    pixel_type ymax;
    pixel_type zmin;
    pixel_type zmax;
};

// Class interpolate
template <typename pixel_type>
class interpolate: public object<pixel_type>
{
protected:
    // ===========================================
    // Internal Variables
    // ===========================================
    image_2d<pixel_type> image_reference;
    grid<pixel_type> x_reference;
    pixel_type fill;
    borders<pixel_type> region;

    // ===========================================
    // Functions
    // ===========================================
    void boundaries();
    int search_(pixel_type x, pixel_type y);

    pixel_type linear_(pixel_type r[2], pixel_type x); // unitary distance
    pixel_type bilinear_(pixel_type r[4], pixel_type x, pixel_type y);
    pixel_type trilinear_(pixel_type r[8], pixel_type x, pixel_type y, pixel_type z);
    
public:
    // ===========================================
    // Create Functions
    // ===========================================
    // interpolate();
    interpolate(image_2d<pixel_type> & input, grid<pixel_type> & xref);

    // ===========================================
    // Functions
    // ===========================================
    image_2d<pixel_type> linear(grid<pixel_type> & xin);
};






// ===========================================
//      Functions of Class grid
// ===========================================

// ===========================================
// Create Functions
// ===========================================
template <typename pixel_type>
interpolate<pixel_type>::interpolate(image_2d<pixel_type> & input, grid<pixel_type> & xref)
{
    assert(input.get_dimension() == xref.get_dimension());
    image_reference = input;
    x_reference = xref;
    fill = 0;               // fill values that are not interpolated
};


// ===========================================
// Functions
// ===========================================

template <typename pixel_type>
image_2d<pixel_type> interpolate<pixel_type>::linear(grid<pixel_type> & xin)
{
    assert(xin.get_dimension() == 2);
    assert(xin.get_dimension() == x_reference.get_dimension());

    using ptr_pixels4 = std::unique_ptr<std::array<pixel_type,4>>;

    image_2d<pixel_type> image_out(xin.get_size()[0],xin.get_size()[1]);
     pixel_type * pout = image_out.ptr();

    pixel_type * xi = (xin.ptr()[0]).ptr();   // raw pointer to
    pixel_type * yi = (xin.ptr()[1]).ptr();   // input grid coordinates

    pixel_type * xr = (x_reference.ptr()[0]).ptr();  // raw pointer to
    pixel_type * yr = (x_reference.ptr()[1]).ptr();  // ref grid coordinates

    pixel_type sx = image_reference.get_spacing()[0];   // scales
    pixel_type sy = image_reference.get_spacing()[1];
    // image_reference.print_data();    // debug

    // ptr_pixels4 i4;
    // pixel_type iout;
    boundaries(); // call to compute region used in search_
    int num = xin.ptr()[0].get_total_elements();
    for(int k = 0; k < num; k++)
    {
        pixel_type xk = xi[k];
        pixel_type yk = yi[k];
        int idx = search_( xk, yk );
        // print to debug
        // std::cout << "xin (x,y): " << xk << ", " << yk <<"\n";
        // std::cout << "idx: " << idx << "\n";
        
        pixel_type iout = fill;
        
        if (idx != -1)
        {
            pixel_type dx = (xk - xr[idx])/ sx;
            pixel_type dy = (yk - yr[idx])/ sy;
            
            ptr_pixels4 i4 = image_reference.neighbors4(idx);
            iout = bilinear_((*i4).data(), dx, dy);

            // std::cout << "dx,dy: " << dx << ", " << dy <<"\n";
            // for(int i = 0; i<4; i++)
            // {
            //     std::cout << "pixel " << i << ": " << (*i4)[i] << "\n";
            // }
            // std::cout << "iout: " << iout << "\n";
        };
        // else{ iout =  this->fill; };
        // std::cout << "idx: " << idx << "\n";
        // std::cout << "iout: " << iout << "\n";

        pout[idx] = iout; // TODO: THIS line produce double error
                            // check everything in the for loop
    }

    // image_out.print_data(); // debug
    return image_out;

};

template <typename pixel_type>
void interpolate<pixel_type>::boundaries()
{
    int w = x_reference.get_size()[0];
    int h = x_reference.get_size()[1];

    pixel_type * xr = (x_reference.ptr()[0]).ptr();  // raw pointer to
    pixel_type * yr = (x_reference.ptr()[1]).ptr();  // ref grid coordinates

    region.xmin = xr[0];
    region.xmax = xr[w-1];

    region.ymin = yr[0];
    region.ymax = yr[w*h-1];
}

template <typename pixel_type>
int interpolate<pixel_type>::search_(pixel_type x, pixel_type y)
{
    // This algorithm assumes a regular fixed grid
    int idx = 0;
    int w = x_reference.get_size()[0];
    // int h = x_reference.get_size()[1];
    // std::cout << "w,h: " << w << " , " << h << "\n";
    // x_reference.print_data();

    // pixel_type * xr = (x_reference.ptr()[0]).ptr();  // raw pointer to
    // pixel_type * yr = (x_reference.ptr()[1]).ptr();  // ref grid coordinates

    std::vector<pixel_type> ss = (x_reference.ptr()[0]).get_spacing();
    pixel_type sx = ss[0];
    pixel_type sy = ss[1];
    
    // pixel_type xmin = xr[0];
    // pixel_type xmax = xr[w-1];

    // pixel_type ymin = yr[0];
    // pixel_type ymax = yr[w*h-1];

    // std::cout << "xmin,xmax: " << xmin << " , " << xmax << "\n";
    // std::cout << "ymin,ymax: " << ymin << " , " << ymax << "\n";
    if(x >= region.xmin and x < region.xmax and y >= region.ymin and y < region.ymax)
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

template <typename pixel_type>
pixel_type interpolate<pixel_type>::linear_(pixel_type r[2], pixel_type x)
{
    return r[0]*(1-x) + r[1]*(x);
};

template <typename pixel_type>
pixel_type interpolate<pixel_type>::bilinear_(pixel_type r[4], pixel_type x, pixel_type y)
{
    pixel_type out[2];
    out[0] = linear_(r, x);
    out[1] = linear_(r+2, x);
    return linear_(out, y);
};


#endif