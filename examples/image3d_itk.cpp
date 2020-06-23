/*
* @Author: jose
* @Date:   2019-11-07 10:12:34
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-03-06 12:50:12
*/

// std libs
#include <iostream>
#include <memory>

// itk
#include <itkImage.h>
#include <itkImageFileReader.h>

// local libs
#include "../src/image.h"
#include "../src/utils/timer.h"


int main(int argc, char *argv[])
{
    if( argc < 2 )
    {
        std::cerr << "Usage: " << std::endl;
        std::cerr << argv[0] << " file_name.ext" << std::endl;
        return EXIT_FAILURE;
    }

    // ===========================================
    // ITK Image
    // ===========================================
    
    // Definitions
    typedef itk::Image< short, 3 >     ImageType;
    typedef itk::ImageFileReader<ImageType>     ReaderType;

    // Objects
    ImageType::Pointer image_itk = ImageType::New();
    ReaderType::Pointer reader = ReaderType::New();

    // Set the image filename itk
    // reader->SetFileName("images/cameraman20x16.tif");
    reader->SetFileName(argv[1]);

    timer tt1("ms");
    tt1.start();
    // Read the image from reader
    reader->Update();
    image_itk = reader->GetOutput();
    tt1.finish();
    
    // 
    ImageType::RegionType region = image_itk->GetLargestPossibleRegion();
    ImageType::SizeType size = region.GetSize();
    std::unique_ptr<short> buffer(image_itk->GetBufferPointer()); // buffer = image->GetBufferPointer();

    // Print data
    std::cout << "Image description (itk):\n" << image_itk;
    std::cout << "\nImage size: " << size << std::endl;
    std::cout << "Image buffer address: " << image_itk->GetBufferPointer() << std::endl;
    std::cout << "Local pointer (unique) address: " << buffer.get() << std::endl;
    buffer.release(); // release the data to avoid double free error

    short * p = image_itk->GetBufferPointer();

    int w = size[0];
    int h = size[1];
    int l = size[2];
    // Print the pixel values
    std::cout << "\nPixel values\n";
    for(int k = 0; k < l; k++)
    {
        std::cout << "[ ";
        for(int i = 0; i < h; i++)
        {
            for(int j=0; j < w; j++)
            {
                std::cout << p[j + i*w + k*w*h] << " "; // valgrind error solved
            };
            if(i < h-1){std::cout << std::endl << "  ";};
        };
        std::cout << "]" << std::endl;
    };
    std::cout << std::endl;

    // ===========================================
    // Project Image
    // ===========================================

    // 3D image
    image<short> image3(3);
    image3.read(argv[1]);

    // Read the image with image interface of itk
    timer tt2("ms");
    tt2.start();
    image3.read(argv[1]);
    tt2.finish();

    std::vector<double> direction;
    direction = image3.get_direction();

    for (int i; i < direction.size(); i++){ std::cout << direction[i]; };

    image3.print("Our Image");
    image3.print_data("Pixel values:");
    std::cout << "Image ptr count: " << image3.get_ptr_count() << std::endl;
    
    tt1.print("Reading image time itk: ");
    tt2.print("Reading image time ours: ");

    // image3.write("output.mha");

    return 0;
};