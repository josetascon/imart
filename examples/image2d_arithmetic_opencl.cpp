/*
* @Author: Jose Tascon
* @Date:   2019-11-07 10:13:08
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-06-18 15:24:24
*/


// std libs
#include <iostream>
#include <memory>
#include <vector>

// local libs
#include "../src/image.h"
#include "../src/vector_opencl.h"

using namespace imart;

int main()
{

    // ============================================
    //      Testing ImageBase2D basic operations
    // ============================================
    std::cout << std::endl;
    std::cout << "===================== ";
    std::cout << "Test class image, warmup addition";
    std::cout << " =====================";
    std::cout << std::endl;

    std::initializer_list<double> l1{1.1, 2.1, 3.1, 4.1, 5.1, 6.1};
    std::initializer_list<double> l2{2.1, 1.1, 0.1, -0.9, -1.9, -2.9};
    vector_opencl<double>::pointer buffer1 = std::make_shared<vector_opencl<double>>(l1);
    vector_opencl<double>::pointer buffer2 = std::make_shared<vector_opencl<double>>(l2);
    
    
    image_opencl<double> image1(buffer1, 3,2);
    image_opencl<double> image2(buffer2, 3,2);
    image_opencl<double> image3;

    std::cout << "Adding 2 images:" << std::endl;
    std::cout << "image3 = image1 + image2" << std::endl;
    image3 = image1 + image2;

    image3.print("image3 Info");
    
    image1.print_data("image1:");
    image2.print_data("image2:");
    image3.print_data("image3:");

    std::cout << "image1 ptr count: " << image1.get_ptr_count() << std::endl;
    std::cout << "image2 ptr count: " << image2.get_ptr_count() << std::endl;
    std::cout << "image3 ptr count: " << image3.get_ptr_count() << std::endl;

    std::cout << std::endl;
    std::cout << "===================== ";
    std::cout << "Test class image, scalar operations";
    std::cout << " =====================";
    std::cout << std::endl;

    image_opencl<int> image4(3,4);
    image_opencl<int> image5(3,4);
    image4.ones();

    image4.print_data("image4:");
    image5 = image4 + 3;
    image5.print_data("image5 = image4 + 3: ");
    image5 = 8 + image4;
    image5.print_data("image5 = 8 + image4: ");
    image5 = image5 - 6;
    image5.print_data("image5 = image5 - 6: ");
    image5 = 7 - image5;
    image5.print_data("image5 = 7 - image5: ");
    image5 = image5*6;
    image5.print_data("image5 = image5 * 6: ");
    image5 = 2*image5;
    image5.print_data("image5 = 2 * image5: ");
    image5 = image5/4;
    image5.print_data("image5 = image5 / 4: ");
    image5 = 36/image5;
    image5.print_data("image5 = 36/image5: ");

    image_opencl<float> image6(3,4);
    image6.ones();
    image6 = image6*3.0;
    image6.print_data("image6:");
    image6 = image6^2.0;
    image6.print_data("image6 = image6^2.0: ");

    image_opencl<float> image7(3,4);
    image_opencl<float> image8(3,4);
    image_opencl<float> image9(3,4);
    image7.assign(3.0);
    image8.assign(0.5);

    image7.print_data("image7:");
    image8.print_data("image8:");
    image9 = image7*3.0 + image8*2.0 + 1.2;
    image9.print_data("image9 = image7*3.0 + image8*2.0 + 1.2");
    

    /*
    std::cout << std::endl;
    std::cout << "===================== ";
    std::cout << "Test class image, arithmetic operations";
    std::cout << " =====================";
    std::cout << std::endl;
    
    // TODO **********
    std::cout << std::endl;
    std::cout << "===================== ";
    std::cout << "Test class image, matrix operations";
    std::cout << " =====================";
    std::cout << std::endl;

    ptr_vector_f buffer3 = std::make_shared<std::vector<float>>(6);
    *buffer3 = {1.1, 2.2, 3.1, 4.1, 2.2, 1.1};
    ptr_vector_f buffer4 = std::make_shared<std::vector<float>>(3);
    *buffer4 = {2.0, 0.5, 1.0};

    image<float> matrix1(buffer3, 3,2);
    image<float> matrix2(buffer4, 1,3);
    image<float> matrix3(1,2);

    matrix3 = matrix1._x_(matrix2);
    matrix1.print_data("matrix1:");
    matrix2.print_data("matrix2:");
    matrix3.print_data("matrix3 = matrix1 x matrix2:");
    // std::cout << "ptr: " << matrix3.get_data() << std::endl;

    ptr_vector_f buffer5 = std::make_shared<std::vector<float>>(6);
    *buffer5 = {1.1, 2.1, 3.1, 4.1, 2.1, 1.1};
    ptr_vector_f buffer6 = std::make_shared<std::vector<float>>(12);
    *buffer6 = {2.0, 1.5, 2.0, 1.0, 3.0, 1.0, 5.0, 2.5, 0.5, 1.0, 0.0, 0.5};

    image<float> matrix4(buffer5, 3,2);
    image<float> matrix5(buffer6, 4,3);
    image<float> matrix6(4,2);
    // matrix6.print();

    matrix6 = matrix4._x_(matrix5);
    matrix4.print_data("matrix4:");
    matrix5.print_data("matrix5:");
    matrix6.print_data("matrix6 = matrix4 x matrix5:");
    // std::cout << "ptr: " << matrix5.get_data() << std::endl;
    */

    return 0;
}