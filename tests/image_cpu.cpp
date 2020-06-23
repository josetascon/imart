/*
* @Author: Jose Tascon
* @Date:   2020-06-06 00:00:00
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-06-17 10:00:37
*/

// std libs
#include <iostream>
#include <memory>
#include <typeinfo>

// gtest header
#include <gtest/gtest.h>

// local header
#include "../src/image.h"

using namespace imart;

template <typename T>
class test_image_cpu : public ::testing::Test{};

using test_types = ::testing::Types<float, double>;

// ============================================
//          Testing Constructors
// ============================================
TYPED_TEST_CASE(test_image_cpu, test_types);

TYPED_TEST(test_image_cpu, constructor)
{
    image_cpu<TypeParam> img1;
    ASSERT_TRUE("image" == img1.get_name());
    ASSERT_EQ( 0, img1.get_total_elements() );
    ASSERT_EQ( 0, img1.get_width() );
    ASSERT_EQ( 0, img1.get_height() );
    ASSERT_EQ( 1, img1.get_length() );
    ASSERT_EQ( 1, img1.get_channels() );

    image_cpu<TypeParam> img2(20,10);
    ASSERT_EQ( 20*10, img2.get_total_elements() );
    ASSERT_EQ( 20, img2.get_width() );
    ASSERT_EQ( 10, img2.get_height() );
    ASSERT_EQ( 1, img2.get_length() );
    ASSERT_EQ( 1, img2.get_channels() );

    image_cpu<TypeParam> img3(121,50,37);
    ASSERT_TRUE("image" == img3.get_name());
    ASSERT_EQ( 121*50*37, img3.get_total_elements() );
    ASSERT_EQ( 121, img3.get_width() );
    ASSERT_EQ( 50, img3.get_height() );
    ASSERT_EQ( 37, img3.get_length() );
    ASSERT_EQ( 1, img3.get_channels() );

    std::vector<int> size({20,30,15});
    image_cpu<TypeParam> img4(size, 1);
    ASSERT_EQ( 20*30*15, img4.get_total_elements() );
    ASSERT_EQ( 20, img4.get_width() );
    ASSERT_EQ( 30, img4.get_height() );
    ASSERT_EQ( 15, img4.get_length() );
    ASSERT_EQ( 1, img4.get_channels() );

    // list

    // image_cpu<TypeParam> img5(10,-5.2);
    // ASSERT_TRUE("image" == img5.get_name());
    // ASSERT_FLOAT_EQ( -5.2, img5[0] );
    // ASSERT_FLOAT_EQ( -5.2, img5[4] );
    // ASSERT_FLOAT_EQ( -5.2, img5[9] );

    // image_cpu<TypeParam> img6(img5);
    // ASSERT_TRUE("image" == img6.get_name());
    // ASSERT_FLOAT_EQ( -5.2, img5[0] );
    // ASSERT_FLOAT_EQ( -5.2, img5[4] );
    // ASSERT_FLOAT_EQ( -5.2, img5[9] );
}

// ============================================
//          Testing Pointers
// ============================================
TYPED_TEST(test_image_cpu, pointers)
{
    // using imgcpu_pointer = typename image_cpu<TypeParam>::pointer;
    // imgcpu_pointer img1 = image_cpu<TypeParam>::new_pointer();
    // ASSERT_TRUE("image" == img1->get_name());
    // ASSERT_EQ( 0, img1->get_total_elements() );
    // ASSERT_EQ( 0, img1->get_width() );
    // ASSERT_EQ( 0, img1->get_height() );
    // ASSERT_EQ( 1, img1->get_length() );
    // ASSERT_EQ( 1, img1->get_channels() );

    // imgcpu_pointer img2 = image_cpu<TypeParam>::new_pointer(7,8);
    // ASSERT_EQ( 7*8, img2->get_total_elements() );
    // ASSERT_EQ( 7, img1->get_width() );
    // ASSERT_EQ( 8, img1->get_height() );
    // ASSERT_EQ( 1, img1->get_length() );
    // ASSERT_EQ( 1, img1->get_channels() );

    // imgcpu_pointer img3 = image_cpu<TypeParam>::new_pointer();
    // ASSERT_TRUE("image" == img3->get_name());
    // ASSERT_EQ( 201, img3->get_total_elements() );

    // imgcpu_pointer img4 = image_cpu<TypeParam>::new_pointer(1,-7.0);
    // ASSERT_FLOAT_EQ( -7.0, img4->operator[](0) );

    // imgcpu_pointer img5 = image_cpu<TypeParam>::new_pointer(21,11.2231);
    // ASSERT_TRUE("image" == img5->get_name());
    // ASSERT_FLOAT_EQ( 11.2231, img5->operator[](0) );
    // ASSERT_FLOAT_EQ( 11.2231, img5->operator[](15) );
    // ASSERT_FLOAT_EQ( 11.2231, img5->operator[](20) );
}
/*
// ============================================
//          Testing Clone Methods
// ============================================
TYPED_TEST(test_image_cpu, clones)
{
    using imgcpu_pointer = typename image_cpu<TypeParam>::pointer;
    imgcpu_pointer img1 = image_cpu<TypeParam>::new_pointer(8,3.14159);
    imgcpu_pointer img2 = img1->clone();
    imgcpu_pointer img3 = img1->copy();
    imgcpu_pointer img4 = img1->mimic();

    ASSERT_TRUE("image" == img2->get_name());
    ASSERT_TRUE("image" == img3->get_name());
    ASSERT_TRUE("image" == img4->get_name());

    ASSERT_EQ( 8, img1->size() );
    ASSERT_EQ( 8, img2->size() );
    ASSERT_EQ( 8, img3->size() );
    ASSERT_EQ( 8, img4->size() );

    ASSERT_FALSE(img1.get() == img2.get());
    ASSERT_FALSE(img1.get() == img3.get());
    ASSERT_FALSE(img1.get() == img4.get());

    ASSERT_FLOAT_EQ( 3.14159, img1->operator[](0) );
    ASSERT_FLOAT_EQ( 3.14159, img1->operator[](7) );
    ASSERT_FLOAT_EQ( 3.14159, img2->operator[](0) );
    ASSERT_FLOAT_EQ( 3.14159, img2->operator[](7) );
    EXPECT_FLOAT_EQ( 0.0, img3->operator[](0) );
    EXPECT_FLOAT_EQ( 0.0, img3->operator[](7) );
    EXPECT_FLOAT_EQ( 0.0, img4->operator[](0) );
    EXPECT_FLOAT_EQ( 0.0, img4->operator[](7) );
}
*/
/*
// ============================================
//          Testing Initialization Methods
// ============================================
TYPED_TEST(test_image_cpu, initialization)
{
    int size = 50;
    using imgcpu_pointer = typename image_cpu<TypeParam>::pointer;
    imgcpu_pointer img1 = image_cpu<TypeParam>::new_pointer(size);
    imgcpu_pointer img2 = img1->mimic();
    imgcpu_pointer img3 = img1->mimic();

    img1->ones();
    img2->assign(TypeParam(3.0));
    img3->random();

    ASSERT_FLOAT_EQ( 1.0, img1->operator[](0) );
    ASSERT_FLOAT_EQ( 1.0, img1->operator[](int(size/2)) );
    ASSERT_FLOAT_EQ( 1.0, img1->operator[](size-4) );
    ASSERT_FLOAT_EQ( 3.0, img2->operator[](0) );
    ASSERT_FLOAT_EQ( 3.0, img2->operator[](int(size/2)) );
    ASSERT_FLOAT_EQ( 3.0, img2->operator[](size-4) );
    ASSERT_FALSE( img3->operator[](0) == img3->operator[](1) );
    ASSERT_FALSE( img3->operator[](0) == img3->operator[](int(size/2)) );
    ASSERT_FALSE( img3->operator[](int(size/2)) == img3->operator[](size-4) );

    img1->zeros();
    img2->assign(TypeParam(-2.431));

    ASSERT_FLOAT_EQ( 0.0, img1->operator[](0) );
    ASSERT_FLOAT_EQ( 0.0, img1->operator[](int(size/2)) );
    ASSERT_FLOAT_EQ( 0.0, img1->operator[](size-4) );
    ASSERT_FLOAT_EQ( -2.431, img2->operator[](0) );
    ASSERT_FLOAT_EQ( -2.431, img2->operator[](int(size/2)) );
    ASSERT_FLOAT_EQ( -2.431, img2->operator[](size-4) );
    ASSERT_FALSE( img3->operator[](0) == img3->operator[](1) );
    ASSERT_FALSE( img3->operator[](0) == img3->operator[](int(size/2)) );
    ASSERT_FALSE( img3->operator[](int(size/2)) == img3->operator[](size-4) );
}

// ============================================
//          Testing Vector Operations
// ============================================
TYPED_TEST(test_image_cpu, vector_operations)
{
    using imgcpu_pointer = typename image_cpu<TypeParam>::pointer;
    imgcpu_pointer img1 = image_cpu<TypeParam>::new_pointer(30,2.1);
    imgcpu_pointer img2 = image_cpu<TypeParam>::new_pointer(30,1.1);
    imgcpu_pointer img3;
    imgcpu_pointer img4 = image_cpu<TypeParam>::new_pointer(30,4.2);
    imgcpu_pointer img5;
    imgcpu_pointer img6 = image_cpu<TypeParam>::new_pointer(30,14.2123);
    imgcpu_pointer vec7;
    imgcpu_pointer vec8 = image_cpu<TypeParam>::new_pointer(30,-2*14.2123);
    imgcpu_pointer vec9;
    imgcpu_pointer img10 = image_cpu<TypeParam>::new_pointer(30,4.0);
    imgcpu_pointer img11;

    img3 = *img1 + *img2;
    ASSERT_TRUE("image" == img3->get_name());
    ASSERT_EQ( 30, img3->size() );
    ASSERT_FLOAT_EQ( 3.2, img3->operator[](0) );
    ASSERT_FLOAT_EQ( 3.2, img3->operator[](7) );
    ASSERT_FLOAT_EQ( 3.2, img3->operator[](18) );
    ASSERT_FLOAT_EQ( 3.2, img3->operator[](29) );

    img5 = *img3 - *img4;
    ASSERT_TRUE("image" == img5->get_name());
    ASSERT_EQ( 30, img5->size() );
    ASSERT_FLOAT_EQ( -1.0, img5->operator[](0) );
    ASSERT_FLOAT_EQ( -1.0, img5->operator[](9) );
    ASSERT_FLOAT_EQ( -1.0, img5->operator[](16) );
    ASSERT_FLOAT_EQ( -1.0, img5->operator[](29) );

    vec7 = (*img5) * (*img6);
    ASSERT_TRUE("image" == vec7->get_name());
    ASSERT_EQ( 30, vec7->size() );
    ASSERT_FLOAT_EQ( -14.2123, vec7->operator[](0) );
    ASSERT_FLOAT_EQ( -14.2123, vec7->operator[](3) );
    ASSERT_FLOAT_EQ( -14.2123, vec7->operator[](12) );
    ASSERT_FLOAT_EQ( -14.2123, vec7->operator[](29) );

    vec9 = (*vec7) / (*vec8);
    ASSERT_TRUE("image" == vec9->get_name());
    ASSERT_EQ( 30, vec9->size() );
    ASSERT_FLOAT_EQ( 0.5, vec9->operator[](0) );
    ASSERT_FLOAT_EQ( 0.5, vec9->operator[](4) );
    ASSERT_FLOAT_EQ( 0.5, vec9->operator[](21) );
    ASSERT_FLOAT_EQ( 0.5, vec9->operator[](29) );

    img11 = (*vec9) ^ (*img10);
    ASSERT_TRUE("image" == img11->get_name());
    ASSERT_EQ( 30, img11->size() );
    ASSERT_FLOAT_EQ( 0.0625, img11->operator[](0) );
    ASSERT_FLOAT_EQ( 0.0625, img11->operator[](5) );
    ASSERT_FLOAT_EQ( 0.0625, img11->operator[](23) );
    ASSERT_FLOAT_EQ( 0.0625, img11->operator[](29) );

    // Check vector values again after operations
    ASSERT_FLOAT_EQ( 2.1, img1->operator[](0) );
    ASSERT_FLOAT_EQ( 2.1, img1->operator[](29) );
    ASSERT_FLOAT_EQ( 1.1, img2->operator[](0) );
    ASSERT_FLOAT_EQ( 1.1, img2->operator[](29) );
    ASSERT_FLOAT_EQ( 3.2, img3->operator[](0) );
    ASSERT_FLOAT_EQ( 3.2, img3->operator[](29) );
    ASSERT_FLOAT_EQ( 4.2, img4->operator[](0) );
    ASSERT_FLOAT_EQ( 4.2, img4->operator[](29) );
    ASSERT_FLOAT_EQ( -1.0, img5->operator[](0) );
    ASSERT_FLOAT_EQ( -1.0, img5->operator[](29) );
    ASSERT_FLOAT_EQ( 14.2123, img6->operator[](0) );
    ASSERT_FLOAT_EQ( 14.2123, img6->operator[](29) );
    ASSERT_FLOAT_EQ( -14.2123, vec7->operator[](0) );
    ASSERT_FLOAT_EQ( -14.2123, vec7->operator[](29) );
    ASSERT_FLOAT_EQ( -2*14.2123, vec8->operator[](0) );
    ASSERT_FLOAT_EQ( -2*14.2123, vec8->operator[](29) );
    ASSERT_FLOAT_EQ( 0.5, vec9->operator[](0) );
    ASSERT_FLOAT_EQ( 0.5, vec9->operator[](29) );
    ASSERT_FLOAT_EQ( 4.0, img10->operator[](0) );
    ASSERT_FLOAT_EQ( 4.0, img10->operator[](29) );
}

// ============================================
//          Testing Scalar Operations
// ============================================
TYPED_TEST(test_image_cpu, scalar_operations)
{
    using imgcpu_pointer = typename image_cpu<TypeParam>::pointer;
    imgcpu_pointer img1 = image_cpu<TypeParam>::new_pointer(14,1.5);
    imgcpu_pointer img2 = image_cpu<TypeParam>::new_pointer(16,-4.632);
    imgcpu_pointer img3 = image_cpu<TypeParam>::new_pointer(18,-12.0);
    imgcpu_pointer img4 = image_cpu<TypeParam>::new_pointer(20,8.4);
    imgcpu_pointer img5 = image_cpu<TypeParam>::new_pointer(22,2.0);
    imgcpu_pointer img6;
    imgcpu_pointer vec7;
    imgcpu_pointer vec8;
    imgcpu_pointer vec9;
    imgcpu_pointer img10;

    // Right hand side
    img6 = *img1 + 4.5;
    ASSERT_TRUE("image" == img6->get_name());
    ASSERT_EQ( 14, img6->size() );
    ASSERT_FLOAT_EQ( 6.0, img6->operator[](0) );
    ASSERT_FLOAT_EQ( 6.0, img6->operator[](7) );
    ASSERT_FLOAT_EQ( 6.0, img6->operator[](13) );

    vec7 = *img2 - 1.111;
    ASSERT_TRUE("image" == vec7->get_name());
    ASSERT_EQ( 16, vec7->size() );
    ASSERT_FLOAT_EQ( -5.743, vec7->operator[](0) );
    ASSERT_FLOAT_EQ( -5.743, vec7->operator[](6) );
    ASSERT_FLOAT_EQ( -5.743, vec7->operator[](15) );

    vec8 = (*img3) * (-0.9);
    ASSERT_TRUE("image" == vec8->get_name());
    ASSERT_EQ( 18, vec8->size() );
    ASSERT_FLOAT_EQ( 10.8, vec8->operator[](0) );
    ASSERT_FLOAT_EQ( 10.8, vec8->operator[](13) );
    ASSERT_FLOAT_EQ( 10.8, vec8->operator[](17) );

    vec9 = (*img4) / (-4.0);
    ASSERT_TRUE("image" == vec9->get_name());
    ASSERT_EQ( 20, vec9->size() );
    ASSERT_FLOAT_EQ( -2.1, vec9->operator[](0) );
    ASSERT_FLOAT_EQ( -2.1, vec9->operator[](11) );
    ASSERT_FLOAT_EQ( -2.1, vec9->operator[](19) );

    img10 = (*img5) ^ (3.0);
    ASSERT_TRUE("image" == img10->get_name());
    ASSERT_EQ( 22, img10->size() );
    ASSERT_FLOAT_EQ( 8.0, img10->operator[](0) );
    ASSERT_FLOAT_EQ( 8.0, img10->operator[](8) );
    ASSERT_FLOAT_EQ( 8.0, img10->operator[](21) );

    // Left hand size
    img6 = TypeParam(4.5) + (*img1);
    ASSERT_TRUE("image" == img6->get_name());
    ASSERT_EQ( 14, img6->size() );
    ASSERT_FLOAT_EQ( 6.0, img6->operator[](0) );
    ASSERT_FLOAT_EQ( 6.0, img6->operator[](7) );
    ASSERT_FLOAT_EQ( 6.0, img6->operator[](13) );

    vec7 = TypeParam(1.111) - (*img2);
    ASSERT_TRUE("image" == vec7->get_name());
    ASSERT_EQ( 16, vec7->size() );
    ASSERT_FLOAT_EQ( 5.743, vec7->operator[](0) );
    ASSERT_FLOAT_EQ( 5.743, vec7->operator[](6) );
    ASSERT_FLOAT_EQ( 5.743, vec7->operator[](15) );

    vec8 = TypeParam(0.9)*(*img3);
    ASSERT_TRUE("image" == vec8->get_name());
    ASSERT_EQ( 18, vec8->size() );
    ASSERT_FLOAT_EQ( -10.8, vec8->operator[](0) );
    ASSERT_FLOAT_EQ( -10.8, vec8->operator[](13) );
    ASSERT_FLOAT_EQ( -10.8, vec8->operator[](17) );

    vec9 = TypeParam(-4.0)/(*img4);
    ASSERT_TRUE("image" == vec9->get_name());
    ASSERT_EQ( 20, vec9->size() );
    ASSERT_FLOAT_EQ( -0.47619051, vec9->operator[](0) );
    ASSERT_FLOAT_EQ( -0.47619051, vec9->operator[](11) );
    ASSERT_FLOAT_EQ( -0.47619051, vec9->operator[](19) );

    // Check original values after operations
    ASSERT_EQ( 14, img1->size() );
    ASSERT_FLOAT_EQ( 1.5, img1->operator[](0) );
    ASSERT_FLOAT_EQ( 1.5, img1->operator[](7) );
    ASSERT_FLOAT_EQ( 1.5, img1->operator[](13) );

    ASSERT_EQ( 16, img2->size() );
    ASSERT_FLOAT_EQ( -4.632, img2->operator[](0) );
    ASSERT_FLOAT_EQ( -4.632, img2->operator[](6) );
    ASSERT_FLOAT_EQ( -4.632, img2->operator[](15) );

    ASSERT_EQ( 18, img3->size() );
    ASSERT_FLOAT_EQ( -12.0, img3->operator[](0) );
    ASSERT_FLOAT_EQ( -12.0, img3->operator[](13) );
    ASSERT_FLOAT_EQ( -12.0, img3->operator[](17) );

    ASSERT_EQ( 20, vec9->size() );
    ASSERT_FLOAT_EQ( 8.4, img4->operator[](0) );
    ASSERT_FLOAT_EQ( 8.4, img4->operator[](11) );
    ASSERT_FLOAT_EQ( 8.4, img4->operator[](19) );

    ASSERT_EQ( 22, img5->size() );
    ASSERT_FLOAT_EQ( 2.0, img5->operator[](0) );
    ASSERT_FLOAT_EQ( 2.0, img5->operator[](8) );
    ASSERT_FLOAT_EQ( 2.0, img5->operator[](21) );
}

// ============================================
//          Testing Reduction Operations
// ============================================
TYPED_TEST(test_image_cpu, reduction_operations)
{
    int size = 32;
    using imgcpu_pointer = typename image_cpu<TypeParam>::pointer;
    imgcpu_pointer img1 = image_cpu<TypeParam>::new_pointer(size,0.0);
    imgcpu_pointer img2 = image_cpu<TypeParam>::new_pointer(size,0.0);

    TypeParam dot1 = 0; TypeParam dot2 = 0;
    TypeParam add1 = 0; TypeParam add2 = 0;
    
    for (int i=0; i<size; ++i)
    {
        img1->operator[](i) = TypeParam(i+100);
        img2->operator[](i) = TypeParam(-i*0.1);
        add1 += (i+100);
        dot1 += (i+100)*(i+100);
        add2 += -i*0.1;
        dot2 += (-i*0.1)*(-i*0.1);
    };

    ASSERT_FLOAT_EQ( 100, img1->min() );
    ASSERT_FLOAT_EQ( 131, img1->max() );
    ASSERT_FLOAT_EQ( add1, img1->sum() );
    ASSERT_FLOAT_EQ( dot1, img1->dot(*img1) );

    ASSERT_FLOAT_EQ( -3.1, img2->min() );
    ASSERT_FLOAT_EQ( 0.0, img2->max() );
    ASSERT_FLOAT_EQ( add2, img2->sum() );
    ASSERT_FLOAT_EQ( dot2, img2->dot(*img2) );
}

// ============================================
//          Testing Other Functions
// ============================================
TYPED_TEST(test_image_cpu, auxialiary_functions)
{
    int size = 32;
    using imgcpu_pointer = typename image_cpu<TypeParam>::pointer;

    // cast function
    imgcpu_pointer img1 = image_cpu<TypeParam>::new_pointer(size,-1.1);

    image_cpu<int>::pointer img2 = img1->template cast<int>();
    image_cpu<float>::pointer img3 = img1->template cast<float>();
    image_cpu<double>::pointer img4 = img1->template cast<double>();
    image_cpu<unsigned int>::pointer img5 = img1->template cast<unsigned int>();

    int vari = -1.0;
    float varf = -1.1;
    double vard = -1.1;
    unsigned int varui = static_cast<unsigned int>(varf);

    ASSERT_TRUE(typeid(int) == typeid(img2->operator[](10)));
    ASSERT_TRUE(typeid(float) == typeid(img3->operator[](10)));
    ASSERT_TRUE(typeid(double) == typeid(img4->operator[](10)));
    ASSERT_TRUE(typeid(unsigned int) == typeid(img5->operator[](10)));
    ASSERT_EQ( vari, img2->operator[](10) );
    ASSERT_EQ( varui, img5->operator[](10) );
    ASSERT_FLOAT_EQ( varf, img3->operator[](10) );
    ASSERT_FLOAT_EQ( vard, img4->operator[](10) );

    // normalize function
    imgcpu_pointer img11 = image_cpu<TypeParam>::new_pointer(size,0.0);
    imgcpu_pointer img12 = image_cpu<TypeParam>::new_pointer(size,0.0);
    imgcpu_pointer img13;
    imgcpu_pointer img14;

    for (int i=0; i<size; ++i)
    {
        img11->operator[](i) = TypeParam(-i-131.2);
        img12->operator[](i) = TypeParam(i+10.5);
    };

    img13 = img11->normalize();
    img14 = img12->normalize();

    ASSERT_FLOAT_EQ( -131.2, img11->operator[](0) );
    ASSERT_FLOAT_EQ( -131.2 - size + 1.0, img11->operator[](size-1) );
    ASSERT_FLOAT_EQ( 10.5, img12->operator[](0) );
    ASSERT_FLOAT_EQ( 10.5 + size - 1.0, img12->operator[](size-1) );
    ASSERT_FLOAT_EQ( 1.0, img13->operator[](0) );
    ASSERT_FLOAT_EQ( 1/double(size-1), img13->operator[](size-2) );
    ASSERT_FLOAT_EQ( 0.0, img13->operator[](size-1) );
    ASSERT_FLOAT_EQ( 0.0, img14->operator[](0) );
    ASSERT_FLOAT_EQ( 1/double(size-1), img14->operator[](1) );
    ASSERT_FLOAT_EQ( 1.0, img14->operator[](size-1) );
}
*/


int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}