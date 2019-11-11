/*
* @Author: jose
* @Date:   2019-11-05 13:55:00
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2019-11-11 15:09:26
*/

#include "image_base.hpp"

// ===========================================
//      Functions of Class image_base_2d
// ===========================================


// ===========================================
// Create Functions
// ===========================================
// Constructor
template <typename pixel_type>
image_base_2d<pixel_type>::image_base_2d()
{
    init(0, 0);
    data.reset();
};

// Constructor
template <typename pixel_type>
image_base_2d<pixel_type>::image_base_2d(int w, int h)
{
    init(w, h);
    data.reset();
    data = std::shared_ptr<pixel_type[]>(new pixel_type[width*height]);
    // data = std::make_shared<float[]>(width*height); // dynamic allocation
};

// Constructor
template <typename pixel_type>
image_base_2d<pixel_type>::image_base_2d(std::shared_ptr<pixel_type[]> buffer, int w, int h)
{
    init(w, h);
    data.reset();
    data = buffer;
};

// Destructor
template <typename pixel_type>
image_base_2d<pixel_type>::~image_base_2d()
{
    data.reset();
};

// Empty init
template <typename pixel_type>
void image_base_2d<pixel_type>::init(int w, int h)
{   
    dim = 2;
    channels = 1;
    width = w;
    height = h;
    num_elements = width*height;

    size = std::vector<int>{width, height};
    spacing = std::vector<pixel_type>(dim, 1.0);
    origin = std::vector<pixel_type>(dim, 0.0);
};

// Copy metadata
template <typename pixel_type>
void image_base_2d<pixel_type>::update(image_base_2d<pixel_type> & input)
{
    width = input.get_width();
    height = input.get_height();
    num_elements = width*height;

    size = input.get_size();
    spacing = input.get_spacing();
    origin = input.get_origin();
};

// ===========================================
// Interface Functions
// ===========================================

// template <size_t type_itk>
// void image_base_2d::read_itk(itk::Image< type_itk, 2 > image_itk)
// {
//     template <size_t type_itk> 
//     using ImageType = itk::Image< type_itk, 2 >;
//     // typedef itk::Image< typename type_itk, 2 > ImageType;

//     ImageType::RegionType region = image_itk->GetLargestPossibleRegion();
//     ImageType::SizeType size = region.GetSize();
//     type * p = image_itk->GetBufferPointer();
//     int w = size[0];
//     int h = size[0];

//     image_base_2d(w,h); // empty current image and create new

//     // Copy of the image
//     for(register int k=0; k<w*h; k++)
//     {
//         *(data+k) = *(p+k);
//     };

// };


// ===========================================
// Get Functions
// ===========================================
template <typename pixel_type>
int image_base_2d<pixel_type>::get_width()
{
    return width;
};

template <typename pixel_type>
int image_base_2d<pixel_type>::get_height()
{
    return height;
};

template <typename pixel_type>
int image_base_2d<pixel_type>::get_total_elements()
{
    return width*height;
};

template <typename pixel_type>
std::vector<int> image_base_2d<pixel_type>::get_size()
{
    return size;
};

template <typename pixel_type>
std::vector<pixel_type> image_base_2d<pixel_type>::get_spacing()
{
    return spacing;
};

template <typename pixel_type>
std::vector<pixel_type> image_base_2d<pixel_type>::get_origin()
{
    return origin;
};

template <typename pixel_type>
std::shared_ptr<pixel_type[]> image_base_2d<pixel_type>::get_data()
{
    return data;
};

template <typename pixel_type>
int image_base_2d<pixel_type>::get_ptr_count()
{
    return data.use_count();
};

// ===========================================
// Print Functions
// ===========================================
template <typename pixel_type>
void image_base_2d<pixel_type>::print(std::string msg)
{
    std::cout << image_base_2d::info(msg);
};

template <typename pixel_type>
void image_base_2d<pixel_type>::print_data()
{
    // std::cout << "Image data:" << std::endl;
    // std::cout << "["
    for(int i=0; i<height; i++)
    {
        for(int j=0; j<width; j++)
        {
            std::cout << data[j+i*width] << " ";
        };
        std::cout << std::endl;
    };
    // std::cout << "]";
    std::cout << std::endl;
};

template <typename pixel_type>
void image_base_2d<pixel_type>::print_ptr_count()
{
    std::cout << "internal image ptr count: ";
    std::cout << data.use_count() << std::endl; // print existing shared pointer
};

template <typename pixel_type>
std::string image_base_2d<pixel_type>::info(std::string msg)
{
    std::stringstream ss;
    std::string title = "Image Information";
    if (msg != "") { title = msg; };

    // Summary of the image information
    ss << "\n===== " << title << " =====\n";
    
    ss << "Pixel type: \t\t" << typeid(data[0]).name() << std::endl;
    ss << "Pixel channels: \t" << channels << std::endl; 
    ss << "Dimensions: \t\t" << dim << std::endl;
    
    ss << "Size: \t\t\t[ ";
    ss << size[0] << " " << size[1] << " ";
    ss << "]" << std::endl;
    
    ss << "Length (mm): \t\t[ ";
    ss << spacing[0] << " " << spacing[1] << " ";
    ss << "]" << std::endl;

    ss << "Origin (mm): \t\t[ ";
    ss << origin[0] << " " << origin[1] << " ";
    ss << "]" << std::endl;

    ss << "Total elements: \t" << num_elements << std::endl; //Get the total number of pixels
    // ss << std::endl;

    return ss.str();
};

// ===========================================
// Overloading Functions
// ===========================================
template <typename pixel_type>
void image_base_2d<pixel_type>::operator = (image_base_2d<pixel_type> & input)
{
    // delete &data;
    data.reset();
    data = input.get_data();
    update(input);
};

// image_base_2d& operator + (image_base_2d & image2)
// {
//     int mat[3][3];
//     for(int i=0; i<3; i++)
//     {
//             for(int j=0; j<3; j++)
//             {
//                     mat[i][j]=a[i][j]+x.a[i][j];
//             }
//     }
// };