/*
* @Author: Jose Tascon
* @Date:   2020-09-09 15:41:57
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-09-22 17:51:19
*/

// std libs
#include <iostream>
#include <vector>

// local libs
#include "../src/image.h"
#include "../src/image_utils.h"
#include "../src/viewer.h"

using namespace imart;

int main(int argc, char *argv[])
{
    if( argc < 2 )
    {
        std::cerr << "Usage: " << std::endl;
        std::cerr << argv[0] << " input_file.mha" << std::endl;
        return EXIT_FAILURE;
    };

    using type = unsigned char;
    using image_type = image_cpu<float>;

    auto input = image_cpu<type>::new_pointer();
    // auto img = image_cpu<type>::new_pointer();
    input->read(argv[1]);
    // cast(*input, *img);

    auto input2 = image_cpu<type>::new_pointer();
    input2->read(argv[2]);

    auto inputfloat = image_cpu<float>::new_pointer();
    auto input2float = image_cpu<float>::new_pointer();
    cast(*input, *inputfloat);
    cast(*input2, *input2float);

    viewer<image_type> view;
    view.subplot(1,2);
    view.add_image(inputfloat);
    view.add_image(input2float);
    view.setup();

    for (size_t i = 0; i<10; i++)
    {
        view.visualize();
        sleep(1);
    };



    return EXIT_SUCCESS;
};