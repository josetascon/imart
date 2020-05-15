/*
* @Author: jose
* @Date:   2020-03-18 00:00:00
* @Last Modified by:   jose
* @Last Modified time: 2020-03-18 00:00:00
*/

#ifndef __DFIELD_H__
#define __DFIELD_H__

// local libs
#include "image_base.h"
#include "transform_base.h"

namespace imart
{

// Class affine
template <typename pixel_type>
class dfield: public transform_base<pixel_type>
{
protected:
    void init(int d);
    void inverse_();

public:
    dfield();
    dfield(int d);
    dfield(int d, typename image<pixel_type>::pointer params);

    // ===========================================
    // Functions
    // ===========================================
    std::vector<pixel_type> transform(std::vector<pixel_type> & point);
    grid<pixel_type> transform(grid<pixel_type> & input);
};

}; //end namespace

#endif