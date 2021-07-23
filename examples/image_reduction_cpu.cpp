/*
* @Author: Jose Tascon
* @Date:   2020-01-28 09:41:41
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2021-07-20 13:39:20
*/

// std libs
#include <iostream>
#include <memory>

// local libs
#include "../src/image.h"
#include "../src/utils/timer.h"

using namespace imart;
using type = double;

int main()
{
    std::cout << "===================== ";
    std::cout << "Test image, reductions functions";
    std::cout << " =====================";
    std::cout << std::endl;

    timer t1;

    // Image 1
    type min1, max1, sum1, prod1;
    image_cpu<type> image1(811,604);
    image1.ones();
    
    t1.start();
    sum1 = image1.sum();
    t1.finish();
    
    image1.print("Image1 in 2D");
    std::cout << "Image1 is filled with ones. Sum = " << sum1 << std::endl;
    std::cout << t1 ;
    if (sum1 == image1.get_total_elements()) { std::cout << "Reduction sum [OK]\n"; }
    else { std::cout<<  "Reduction sum [FAIL]\n"; };

    std::cout << "\nImage1 random initilization\n";
    image1.random();
    // image1.print_data();
    min1 = image1.min();
    max1 = image1.max();
    sum1 = image1.sum();
    int idx_min1 = image1.argmin();
    int idx_max1 = image1.argmax();
    // prod1 = image1.prod();
    std::cout << "Image1 min: " << min1 << std::endl;
    std::cout << "Image1 argmin: " << idx_min1 << std::endl;
    std::cout << "Image1 (argmin): " << image1[idx_min1] << std::endl;
    std::cout << "Image1 max: " << max1 << std::endl;
    std::cout << "Image1 argmax: " << idx_max1 << std::endl;
    std::cout << "Image1 (argmax): " << image1[idx_max1] << std::endl;
    std::cout << "Image1 sum: " << sum1 << std::endl;
    // std::cout << "Image1 prod: " << prod1 << std::endl;

    // Image 2
    type min2, max2, sum2, prod2;
    image_cpu<type> image2(201,125,153);
    image2.ones();
    // std::cout << image2(0,0,0);
    // std::cout << image2(200,200,100);
    
    t1.start();
    sum2 = image2.sum();
    t1.finish();

    image2.print("Image2 in 3D");
    std::cout << "Image2 is filled with ones. Sum = " << sum2 << std::endl;
    std::cout << t1;
    if (sum2 == image2.get_total_elements()) { std::cout << "Reduction sum [OK]\n"; }
    else { std::cout<<  "Reduction sum [FAIL]\n"; };

    std::cout << "\nImage2 random initilization\n";
    image2.random();
    min2 = image2.min();
    max2 = image2.max();
    sum2 = image2.sum();
    // prod2 = image2.prod();
    std::cout << "Image2 min: " << min2 << std::endl;
    std::cout << "Image2 max: " << max2 << std::endl;
    std::cout << "Image2 sum: " << sum2 << std::endl;
    // std::cout << "Image2 prod: " << prod2 << std::endl;
    std::cout << std::endl;

    std::cout << "===================== ";
    std::cout << "Test image, dot product function";
    std::cout << " =====================";
    std::cout << std::endl;

    image_cpu<type> image3(5,1);
    *image3.get_data() = {1.0, 2.0, 3.0, 1.0, 4.0};
    
    type x;
    image_cpu<type> result = image3*image3;
    x = result.sum();

    type dot;
    dot = image3.dot(image3);
    image3.print_data("Image3 values:");
    std::cout << "Self dot product of Image3: " << dot << std::endl;
    if (x == dot) { std::cout << "Dot product [OK]\n"; }
    else { std::cout<<  "Dot product [FAIL]\n"; };

    return 0;
};
