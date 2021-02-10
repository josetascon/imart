/*
* @Author: Jose Tascon
* @Date:   2020-09-09 15:41:57
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2021-02-01 23:48:32
*/

// std libs
#include <iostream>
#include <vector>

// local libs
#include "../src/utils/timer.h"
#include "../src/image.h"
#include "../src/image_utils.h"
#include "../src/viewer.h"

using namespace imart;

int main(int argc, char *argv[])
{
    if( argc != 3 )
    {
        std::cerr << "Usage: " << std::endl;
        std::cerr << argv[0] << " image1.ext image2.ext" << std::endl;
        return EXIT_FAILURE;
    };

    using type = unsigned char;
    using image_type = image_cpu<type>;

    auto input = image_cpu<type>::new_pointer();
    auto operation = image_cpu<type>::new_pointer();
    // auto img = image_cpu<type>::new_pointer();
    input->read(argv[1]);
    // cast(*input, *img);

    auto input2 = image_cpu<type>::new_pointer();
    input2->read(argv[2]);

    // auto inputfloat = image_cpu<float>::new_pointer();
    // auto input2float = image_cpu<float>::new_pointer();
    // cast(*input, *inputfloat);
    // cast(*input2, *input2float);

    timer t("ms");
    t.start();

    viewer<image_type> view;
    view.subplot(1,2);
    view.add_image(input);
    view.add_image(input2);
    view.setup();

    t.lap();
    printf("Setup time: %5.2f [ms]\n", t.get_elapsed());    

    for (size_t i = 0; i<10; i++)
    {
        view.visualize();
        t.lap();
        printf("Visualize time: %5.2f [ms]\n", t.get_elapsed());
        sleep(1);

        // input->equal(*input - 10);
        // view.update_image(input, 0);

        *input = (*input - 10);

        t.lap();
        printf("Sleep time: %5.2f [ms]\n", t.get_elapsed());
        view.update_image(input, 0);
        t.lap();
        printf("Update time: %5.2f [ms]\n", t.get_elapsed());
    };



    return EXIT_SUCCESS;
};