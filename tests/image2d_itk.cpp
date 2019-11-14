/*
* @Author: jose
* @Date:   2019-11-07 10:12:34
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2019-11-12 15:50:39
*/

// std libs
#include <iostream>
#include <memory>

// itk
#include <itkImage.h>
#include <itkImageFileReader.h>

// local libs
#include "../inc/image_base.hpp"


int main(int argc, char *argv[])
{
    if( argc < 2 )
    {
        std::cerr << "Usage: " << std::endl;
        std::cerr << argv[0] << " file_name.ext" << std::endl;
        return EXIT_FAILURE;
    }

    // Definitions
    typedef itk::Image< unsigned short, 2 >     ImageType;
    typedef itk::ImageFileReader<ImageType>     ReaderType;

    // Objects
    ImageType::Pointer image_itk = ImageType::New();
    ReaderType::Pointer reader = ReaderType::New();

    // Set the image filename itk
    // reader->SetFileName("images/cameraman20x16.tif");
    reader->SetFileName(argv[1]);
    reader->Update();
    
    // Read the image from reader
    image_itk = reader->GetOutput();
    ImageType::RegionType region = image_itk->GetLargestPossibleRegion();
    ImageType::SizeType size = region.GetSize();
    std::unique_ptr<unsigned short> buffer(image_itk->GetBufferPointer()); // buffer = image->GetBufferPointer();

    // Print data
    std::cout << "Image description (itk):\n" << image_itk;
    std::cout << "\nImage size: " << size << std::endl;
    std::cout << "Image buffer address: " << image_itk->GetBufferPointer() << std::endl;
    std::cout << "Local pointer (unique) address: " << buffer.get() << std::endl;
    buffer.release(); // release the data to avoid double free error

    unsigned short * p = image_itk->GetBufferPointer();

    int width = size[0];
    int height = size[1];
    // Print the pixel values
    std::cout << "\nPixel values\n";
    for (int i=0; i<height; i++)
    {
        for (int j=0; j<width; j++)
        {
            std::cout << *(p++) << " ";
        }
        std::cout << "\n";
    };

    // Read the image with image_base_2d interface of itk
    image_base_2d<unsigned short> image1;
    image1.read(argv[1]);
    image1.print("Our Image");
    image1.print_data("Pixel values:");
    std::cout << "Image ptr count: " << image1.get_ptr_count() << std::endl;

    return 0;
};