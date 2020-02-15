/*
* @Author: Jose Tascon
* @Date:   2019-11-07 10:13:08
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-02-15 19:06:34
*/


// std libs
#include <iostream>
#include <memory>
#include <vector>
#include <complex>

// local libs
#include "../src/image.h"

using type = double;

int main()
{
    image<type> img1(4,3);
    *img1.get_data() = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0};
    img1.print_data("image1");

    image<std::complex<type>> img2(4,3);
    image<std::complex<type>> img3(4,3);
    img2 = fft(img1);
    // img2.print();
    img2.print_data("image2 = fft(image1)");

    img3 = ifft(img2);
    // img3.print();
    img3.print_data("image3 = ifft(image2)");

    int w = 7;
    int h = 6;
    image<type> img10(w,h);
    type * p = img10.ptr();
    for(int i = 0; i < w*h; i++)
    {
        p[i] = i;
    };
    typename image_base<type>::vector grad(2);
    grad = gradient(img10);
    img10.print_data("i");
    grad[0]->print_data("didx");
    grad[1]->print_data("didy");

    // CONTINUE unpad*** Febrero 07 2020

    // // test padding
    // image<type> img20(3,2);
    // img20.ones();
    // image<type> a = pad(img20, std::vector<int>{2,1}, std::vector<int>{1,2});
    // img20.print_data("im to pad");
    // a.print_data("pad");

    return 0;
}