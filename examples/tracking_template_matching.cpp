/*
* @Author: Jose Tascon
* @Date:   2019-11-18 13:30:52
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2021-07-22 10:14:04
*/

// std libs
#include <iostream>
#include <memory>
#include <sstream>
#include <iomanip>
#include <filesystem>
#include <thread>

// local libs
#include "../src/image.h"
#include "../src/image_utils.h"
#include "../src/template_matching.h"

using namespace imart;

int main(int argc, char *argv[])
{
    // using type = float;
    using type = double;
    using image_type = image_cpu<type>;
    
    if( argc < 4 )
    {
        std::cerr << "Usage: " << std::endl;
        std::cerr << argv[0] << " <fixed> <template> <moving>" << std::endl;
        return EXIT_FAILURE;
    }
    // Files
    std::string file_fixed = argv[1];
    std::string file_template = argv[2];
    std::string file_moving = argv[3];

    // ============================================
    //      Testing tracking with two images
    // ============================================

    // Images
    auto img_fixed = image_cpu<type>::new_pointer();
    auto img_moving = image_cpu<type>::new_pointer();
    auto img_template = image_cpu<type>::new_pointer();
    
    // Read
    img_fixed->read(file_fixed);
    img_moving->read(file_moving);
    img_template->read(file_template);

    auto bbox_fixed = bounding_box(img_template);

    timer t("ms");
    t.start();
    auto tm = template_matching<type, vector_cpu<type>>::new_pointer(img_fixed, bbox_fixed);
    tm->set_slide(std::vector<int>{30,30});
    auto bbox_moving = tm->apply(img_moving);

    t.finish();
    printf("Total time: \t%5.2f [ms]\n", t.get_elapsed());
    
    printf("Template bounding box: (%d, %d, %d, %d)\n",
        bbox_fixed[0][0], bbox_fixed[0][1], bbox_fixed[1][0], bbox_fixed[1][1]);
    printf("Matched bounding box: (%d, %d, %d, %d)\n",
        bbox_moving[0][0], bbox_moving[0][1], bbox_moving[1][0], bbox_moving[1][1]);
    

    return 0;
};