/*
* @Author: jose
* @Date:   2019-11-05 13:55:00
* @Last Modified by:   jose
* @Last Modified time: 2019-11-08 11:14:34
*/

#include "image_base.hpp"

// ===========================================
//      Functions of Class image_base_2d
// ===========================================


// ===========================================
// Create Functions
// ===========================================
// Constructor
image_base_2d::image_base_2d()
{
    init(0, 0);
    data.reset();
};

// Constructor
image_base_2d::image_base_2d(int w, int h)
{
    init(w, h);
    data.reset();
    data = std::shared_ptr<pixel_type[]>(new pixel_type[width*height]);
    // data = std::make_shared<float[]>(width*height); // dynamic allocation
};

// Constructor
image_base_2d::image_base_2d(std::shared_ptr<pixel_type[]> buffer, int w, int h)
{
    init(w, h);
    data.reset();
    data = buffer;
};

// Destructor
image_base_2d::~image_base_2d()
{
    data.reset();
};

// Empty init
void image_base_2d::init(int w, int h)
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
void image_base_2d::update(image_base_2d & input)
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
int image_base_2d::get_width()
{
    return width;
};

int image_base_2d::get_height()
{
    return height;
};

int image_base_2d::get_total_elements()
{
    return width*height;
};

std::vector<int> image_base_2d::get_size()
{
    return size;
};

std::vector<pixel_type> image_base_2d::get_spacing()
{
    return spacing;
};

std::vector<pixel_type> image_base_2d::get_origin()
{
    return origin;
};

std::shared_ptr<pixel_type[]> image_base_2d::get_data()
{
    return data;
};

int image_base_2d::get_ptr_count()
{
    return data.use_count();
};

// ===========================================
// Print Functions
// ===========================================
void image_base_2d::print(std::string msg)
{
    std::cout << image_base_2d::info(msg);
};

void image_base_2d::print_data()
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

void image_base_2d::print_ptr_count()
{
    std::cout << "internal image ptr count: ";
    std::cout << data.use_count() << std::endl; // print existing shared pointer
};

std::string image_base_2d::info(std::string msg)
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
void image_base_2d::operator = (image_base_2d & input)
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