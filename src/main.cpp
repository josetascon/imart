/*
* @Author: Jose Tascon
* @Date:   2019-11-15 11:26:18
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-09-03 10:41:03
*/

#include "../src/inherit.h"
#include "../src/object.h"
#include "../src/opencl_object.h"
#include "../src/cuda_object.h"
#include "../src/space_object.h"
#include "../src/data_object.h"
#include "../src/process_object.h"
#include "../src/pairwise_object.h"

#include "../src/cuda_buffer.h"

#include "../src/image.h"
#include "../src/image_utils.h"
#include "../src/grid.h"

#include "../src/transform.h"
#include "../src/affine.h"
#include "../src/dfield.h"

#include "../src/interpolator.h"
#include "../src/inearest.h"
#include "../src/ilinear.h"

#include "../src/metric.h"
#include "../src/ssd.h"
#include "../src/cc.h"
#include "../src/demons.h"

#include "../src/optimizer.h"
#include "../src/gradient_descent.h"

#include "../src/registration.h"

#include "../src/utils/timer.h"
#include "../src/utils/opencl_config.h"

;