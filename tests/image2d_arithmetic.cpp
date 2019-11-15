/*
* @Author: Jose Tascon
* @Date:   2019-11-07 10:13:08
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2019-11-15 11:09:19
*/


// std libs
#include <iostream>
#include <memory>

// local libs
#include "../src/image_2d.h"

int main()
{
    // ============================================
    //      Testing ImageBase2D basic operations
    // ============================================
    std::cout << std::endl;
    std::cout << "===================== ";
    std::cout << "Test class image_2d, warmup addition";
    std::cout << " =====================";
    std::cout << std::endl;
    std::shared_ptr<float[]> buffer1(new float[6] {1.1, 2.1, 3.1, 4.1, 5.1, 6.1});
    std::shared_ptr<float[]> buffer2(new float[6] {2.1, 1.1, 0.1, -0.9, -1.9, -2.9});
    

    image_2d<float> image1(buffer1, 3,2);
    image_2d<float> image2(buffer2, 3,2);
    image_2d<float> image3;

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
    std::cout << "Test class image_2d, scalar operations";
    std::cout << " =====================";
    std::cout << std::endl;

    image_2d<int> image4(3,4);
    image_2d<int> image5(3,4);
    image4.ones();

    image4.print_data("image4 :");
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

    std::cout << std::endl;
    std::cout << "===================== ";
    std::cout << "Test class image_2d, arithmetic operations";
    std::cout << " =====================";
    std::cout << std::endl;
    

    std::cout << std::endl;
    std::cout << "===================== ";
    std::cout << "Test class image_2d, matrix operations";
    std::cout << " =====================";
    std::cout << std::endl;

    std::shared_ptr<float[]> buffer3(new float[6] {1.1, 2.2, 3.1, 4.1, 2.2, 1.1});
    std::shared_ptr<float[]> buffer4(new float[3] {2.0, 0.5, 1.0});
    image_2d<float> matrix1(buffer3, 3,2);
    image_2d<float> matrix2(buffer4, 1,3);
    image_2d<float> matrix3(1,2);

    matrix3 = matrix1._x_(matrix2);
    matrix1.print_data("matrix1:");
    matrix2.print_data("matrix2:");
    matrix3.print_data("matrix3 = matrix1 x matrix2:");
    // std::cout << "ptr: " << matrix3.get_data() << std::endl;

    std::shared_ptr<float[]> buffer5(new float[6] {1.1, 2.1, 3.1, 4.1, 2.1, 1.1});
    std::shared_ptr<float[]> buffer6(new float[12] {2.0, 1.5, 2.0, 1.0, 3.0, 1.0, 5.0, 2.5, 0.5, 1.0, 0.0, 0.5});
    image_2d<float> matrix4(buffer5, 3,2);
    image_2d<float> matrix5(buffer6, 4,3);
    image_2d<float> matrix6(4,2);
    // matrix6.print();

    matrix6 = matrix4._x_(matrix5);
    matrix4.print_data("matrix4:");
    matrix5.print_data("matrix5:");
    matrix6.print_data("matrix6 = matrix4 x matrix5:");
    // std::cout << "ptr: " << matrix5.get_data() << std::endl;



    return 0;
}