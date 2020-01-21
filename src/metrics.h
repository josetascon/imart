/*
* @Author: jose
* @Date:   2019-12-09 00:00:00
* @Last Modified by:   jose
* @Last Modified time: 2019-12-09 00:00:00
*/

#include "image_base.h"
#include "transform_base.h"



template <typename pixel_type>
class metric: public object<pixel_type>
{
public:
    //Type definitions
    ;
protected:
    // ===========================================
    // Internal Variables
    // ===========================================
    pixel_type cost_value;
    image_base<pixel_type> fixed;
    image_base<pixel_type> moving;
    transform_base<pixel_type> transform;

    // ===========================================
    // Functions
    // ===========================================
    void init(image_base<pixel_type> fixed_image, image_base<pixel_type> moving_image, transform_base<pixel_type> transform_reg);

public:
    // ===========================================
    // Create Functions
    // ===========================================
    // Constructors
    metric(image_base<pixel_type> fixed_image, image_base<pixel_type> moving_image, transform_base<pixel_type> transform_reg);

    // ===========================================
    // Get Functions
    // ===========================================
    pixel_type get_cost();

    // ===========================================
    // Print Functions
    // ===========================================
    virtual std::string info(std::string msg);
    virtual std::string info_data(std::string msg);


    // ===========================================
    // Functions
    // ===========================================
    void cost();
    transform_base<pixel_type> derivative();

};






// ===========================================
//      Functions of Class metric
// ===========================================


// ===========================================
// Create Functions
// ===========================================
// Constructor
template <typename pixel_type>
metric<pixel_type>::metric(image_base<pixel_type> fixed_image, image_base<pixel_type> moving_image, transform_base<pixel_type> transform_reg)
{
    init(fixed_image, moving_image, transform_reg)
};



void metric<pixel_type>::init(image_base<pixel_type> fixed_image, image_base<pixel_type> moving_image, transform_base<pixel_type> transform_reg)
{
    fixed = fixed_image;
    moving = moving_image;
    transform = transform_reg;
}



template <typename pixel_type>
std::string metric<pixel_type>::info(std::string msg)
{
    std::stringstream ss;
    std::string title = "Metric Information";
    if (msg != "") { title = msg; };
    // Summary of the metric information
    ss << object<pixel_type>::info(title);
}