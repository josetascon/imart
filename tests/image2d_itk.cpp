/*
* @Author: jose
* @Date:   2019-11-07 10:12:34
* @Last Modified by:   jose
* @Last Modified time: 2019-11-08 10:41:54
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

    // Print data
    // image
    
    // Read the image from reader
    image_itk = reader->GetOutput();

    std::cout << "Image description (itk):\n" << image_itk;

    ImageType::RegionType region = image_itk->GetLargestPossibleRegion();
    ImageType::SizeType size = region.GetSize();

    std::cout << "\nImage size: " << size << std::endl;
    std::cout << "Image buffer address: " << image_itk->GetBufferPointer() << std::endl;

    std::unique_ptr<unsigned short> buffer(image_itk->GetBufferPointer()); // buffer = image->GetBufferPointer();
    std::cout << "Local pointer (unique) address: " << buffer.get() << std::endl;
    buffer.release(); // release the data to avoid double free error

    unsigned short * p = image_itk->GetBufferPointer();

    int width = size[0];
    int height = size[1];
    // print the pixels
    std::cout << "\nPixel values\n";
    for (register int i=0; i<width; i++)
    {
        for (register int j=0; j<height; j++)
        {
            std::cout << *(p++) << " ";
        }
        std::cout << "\n";
    };

    return 0;
};