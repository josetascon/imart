/*
* @Author: Jose Tascon
* @Date:   2020-06-06 00:00:00
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-08-31 19:29:15
*/

// std libs
#include <iostream>
#include <memory>
#include <typeinfo>

// gtest header
#include <gtest/gtest.h>

// local header
#include "../src/dfield.h"

using namespace imart;

template <typename T>
class test_dfield_cuda : public ::testing::Test{};

using test_types = ::testing::Types<float, double>;

// ============================================
//          Testing Constructors
// ============================================
TYPED_TEST_CASE(test_dfield_cuda, test_types);

TYPED_TEST(test_dfield_cuda, constructor)
{
    dfield_cuda<TypeParam> field1;
    ASSERT_TRUE("dfield" == field1.get_name());
    ASSERT_EQ( 2, field1.get_dimension() );
    ASSERT_EQ( 0, field1.get_size()[0] );
    ASSERT_EQ( 0, field1.get_size()[1] );
    ASSERT_FALSE( nullptr == field1.get_parameters(0).get() );
    ASSERT_TRUE( "image" == field1.get_parameters(0)->get_name() );
    ASSERT_FALSE( nullptr == field1.get_parameters(1).get() );
    ASSERT_TRUE( "image" == field1.get_parameters(1)->get_name() );
    
    dfield_cuda<TypeParam> field2(2);
    ASSERT_TRUE("dfield" == field2.get_name());
    ASSERT_EQ( 2, field2.get_dimension() );
    ASSERT_EQ( 0, field2.get_size()[0] );
    ASSERT_EQ( 0, field2.get_size()[1] );
    ASSERT_FALSE( nullptr == field2.get_parameters(0).get() );
    ASSERT_TRUE( "image" == field2.get_parameters(0)->get_name() );
    ASSERT_FALSE( nullptr == field2.get_parameters(1).get() );
    ASSERT_TRUE( "image" == field2.get_parameters(1)->get_name() );

    dfield_cuda<TypeParam> field3(3);
    ASSERT_TRUE("dfield" == field3.get_name());
    ASSERT_EQ( 3, field3.get_dimension() );
    ASSERT_EQ( 0, field3.get_size()[0] );
    ASSERT_EQ( 0, field3.get_size()[1] );
    ASSERT_EQ( 0, field3.get_size()[2] );
    ASSERT_FALSE( nullptr == field3.get_parameters(0).get() );
    ASSERT_TRUE( "image" == field3.get_parameters(0)->get_name() );
    ASSERT_FALSE( nullptr == field3.get_parameters(1).get() );
    ASSERT_TRUE( "image" == field3.get_parameters(1)->get_name() );
    ASSERT_FALSE( nullptr == field3.get_parameters(2).get() );
    ASSERT_TRUE( "image" == field3.get_parameters(2)->get_name() );
    
    dfield_cuda<TypeParam> field4(std::vector<int>{22,14});
    ASSERT_TRUE("dfield" == field4.get_name());
    ASSERT_EQ( 2, field4.get_dimension() );
    ASSERT_EQ( 22, field4.get_size()[0] );
    ASSERT_EQ( 14, field4.get_size()[1] );
    ASSERT_FALSE( nullptr == field4.get_parameters(0).get() );
    ASSERT_TRUE( "image" == field4.get_parameters(0)->get_name() );
    ASSERT_EQ( 22, field4.get_parameters(0)->get_width() );
    ASSERT_EQ( 14, field4.get_parameters(0)->get_height() );
    ASSERT_EQ( 1, field4.get_parameters(0)->get_length() );
    ASSERT_EQ( 1, field4.get_parameters(0)->get_channels() );
    ASSERT_FALSE( nullptr == field4.get_parameters(1).get() );
    ASSERT_TRUE( "image" == field4.get_parameters(1)->get_name() );
    ASSERT_EQ( 22, field4.get_parameters(1)->get_width() );
    ASSERT_EQ( 14, field4.get_parameters(1)->get_height() );
    ASSERT_EQ( 1, field4.get_parameters(1)->get_length() );
    ASSERT_EQ( 1, field4.get_parameters(1)->get_channels() );
    ASSERT_EQ( 0, field4.get_parameters(0)->operator[](0) ); // check zero initialization
    ASSERT_EQ( 0, field4.get_parameters(0)->operator[](5) );
    ASSERT_EQ( 0, field4.get_parameters(1)->operator[](0) );
    ASSERT_EQ( 0, field4.get_parameters(1)->operator[](5) );

    dfield_cuda<TypeParam> field5(std::vector<int>{8,54,31});
    ASSERT_TRUE("dfield" == field5.get_name());
    ASSERT_EQ( 3, field5.get_dimension() );
    ASSERT_EQ( 8, field5.get_size()[0] );
    ASSERT_EQ( 54, field5.get_size()[1] );
    ASSERT_EQ( 31, field5.get_size()[2] );
    ASSERT_FALSE( nullptr == field5.get_parameters(0).get() );
    ASSERT_TRUE( "image" == field5.get_parameters(0)->get_name() );
    ASSERT_EQ( 8, field5.get_parameters(0)->get_width() );
    ASSERT_EQ( 54, field5.get_parameters(0)->get_height() );
    ASSERT_EQ( 31, field5.get_parameters(0)->get_length() );
    ASSERT_EQ( 1, field5.get_parameters(0)->get_channels() );
    ASSERT_FALSE( nullptr == field5.get_parameters(1).get() );
    ASSERT_TRUE( "image" == field5.get_parameters(1)->get_name() );
    ASSERT_EQ( 8, field5.get_parameters(1)->get_width() );
    ASSERT_EQ( 54, field5.get_parameters(1)->get_height() );
    ASSERT_EQ( 31, field5.get_parameters(1)->get_length() );
    ASSERT_EQ( 1, field5.get_parameters(1)->get_channels() );
    ASSERT_FALSE( nullptr == field5.get_parameters(2).get() );
    ASSERT_TRUE( "image" == field5.get_parameters(2)->get_name() );
    ASSERT_EQ( 8, field5.get_parameters(2)->get_width() );
    ASSERT_EQ( 54, field5.get_parameters(2)->get_height() );
    ASSERT_EQ( 31, field5.get_parameters(2)->get_length() );
    ASSERT_EQ( 1, field5.get_parameters(2)->get_channels() );
    ASSERT_EQ( 0, field5.get_parameters(0)->operator[](0) ); // check zero initialization
    ASSERT_EQ( 0, field5.get_parameters(0)->operator[](5) );
    ASSERT_EQ( 0, field5.get_parameters(1)->operator[](0) );
    ASSERT_EQ( 0, field5.get_parameters(1)->operator[](5) );
    ASSERT_EQ( 0, field5.get_parameters(2)->operator[](0) );
    ASSERT_EQ( 0, field5.get_parameters(2)->operator[](5) );
    
    auto img1 = image_cuda<TypeParam>::new_pointer(18,13);
    dfield_cuda<TypeParam> field6(img1);
    ASSERT_TRUE("dfield" == field6.get_name());
    ASSERT_EQ( 2, field6.get_dimension() );
    ASSERT_EQ( 18, field6.get_size()[0] );
    ASSERT_EQ( 13, field6.get_size()[1] );
    ASSERT_FALSE( nullptr == field6.get_parameters(0).get() );
    ASSERT_TRUE( "image" == field6.get_parameters(0)->get_name() );
    ASSERT_EQ( 18, field6.get_parameters(0)->get_width() );
    ASSERT_EQ( 13, field6.get_parameters(0)->get_height() );
    ASSERT_EQ( 1, field6.get_parameters(0)->get_length() );
    ASSERT_EQ( 1, field6.get_parameters(0)->get_channels() );
    ASSERT_FALSE( nullptr == field6.get_parameters(1).get() );
    ASSERT_TRUE( "image" == field6.get_parameters(1)->get_name() );
    ASSERT_EQ( 18, field6.get_parameters(1)->get_width() );
    ASSERT_EQ( 13, field6.get_parameters(1)->get_height() );
    ASSERT_EQ( 1, field6.get_parameters(1)->get_length() );
    ASSERT_EQ( 1, field6.get_parameters(1)->get_channels() );
    ASSERT_EQ( 0, field6.get_parameters(0)->operator[](0) ); // check zero initialization
    ASSERT_EQ( 0, field6.get_parameters(0)->operator[](5) );
    ASSERT_EQ( 0, field6.get_parameters(1)->operator[](0) );
    ASSERT_EQ( 0, field6.get_parameters(1)->operator[](5) );

    auto img2 = image_cuda<TypeParam>::new_pointer(19,11);
    auto grid2 = image_cuda<TypeParam>::new_pointer(*img2);
    dfield_cuda<TypeParam> field7(grid2);
    ASSERT_TRUE("dfield" == field7.get_name());
    ASSERT_EQ( 2, field7.get_dimension() );
    ASSERT_EQ( 19, field7.get_size()[0] );
    ASSERT_EQ( 11, field7.get_size()[1] );
    ASSERT_FALSE( nullptr == field7.get_parameters(0).get() );
    ASSERT_TRUE( "image" == field7.get_parameters(0)->get_name() );
    ASSERT_EQ( 19, field7.get_parameters(0)->get_width() );
    ASSERT_EQ( 11, field7.get_parameters(0)->get_height() );
    ASSERT_EQ( 1, field7.get_parameters(0)->get_length() );
    ASSERT_EQ( 1, field7.get_parameters(0)->get_channels() );
    ASSERT_FALSE( nullptr == field7.get_parameters(1).get() );
    ASSERT_TRUE( "image" == field7.get_parameters(1)->get_name() );
    ASSERT_EQ( 19, field7.get_parameters(1)->get_width() );
    ASSERT_EQ( 11, field7.get_parameters(1)->get_height() );
    ASSERT_EQ( 1, field7.get_parameters(1)->get_length() );
    ASSERT_EQ( 1, field7.get_parameters(1)->get_channels() );
    ASSERT_EQ( 0, field7.get_parameters(0)->operator[](1) ); // check zero initialization
    ASSERT_EQ( 0, field7.get_parameters(0)->operator[](20) );
    ASSERT_EQ( 0, field7.get_parameters(1)->operator[](1) );
    ASSERT_EQ( 0, field7.get_parameters(1)->operator[](20) );

    auto img3 = image_cuda<TypeParam>::new_pointer(14,15,12);
    dfield_cuda<TypeParam> field8(img3);
    ASSERT_TRUE("dfield" == field8.get_name());
    ASSERT_EQ( 3, field8.get_dimension() );
    ASSERT_EQ( 14, field8.get_size()[0] );
    ASSERT_EQ( 15, field8.get_size()[1] );
    ASSERT_EQ( 12, field8.get_size()[2] );
    ASSERT_FALSE( nullptr == field8.get_parameters(0).get() );
    ASSERT_TRUE( "image" == field8.get_parameters(0)->get_name() );
    ASSERT_EQ( 14, field8.get_parameters(0)->get_width() );
    ASSERT_EQ( 15, field8.get_parameters(0)->get_height() );
    ASSERT_EQ( 12, field8.get_parameters(0)->get_length() );
    ASSERT_EQ( 1, field8.get_parameters(0)->get_channels() );
    ASSERT_FALSE( nullptr == field8.get_parameters(1).get() );
    ASSERT_TRUE( "image" == field8.get_parameters(1)->get_name() );
    ASSERT_EQ( 14, field8.get_parameters(1)->get_width() );
    ASSERT_EQ( 15, field8.get_parameters(1)->get_height() );
    ASSERT_EQ( 12, field8.get_parameters(1)->get_length() );
    ASSERT_EQ( 1, field8.get_parameters(1)->get_channels() );
    ASSERT_FALSE( nullptr == field8.get_parameters(2).get() );
    ASSERT_TRUE( "image" == field8.get_parameters(2)->get_name() );
    ASSERT_EQ( 14, field8.get_parameters(2)->get_width() );
    ASSERT_EQ( 15, field8.get_parameters(2)->get_height() );
    ASSERT_EQ( 12, field8.get_parameters(2)->get_length() );
    ASSERT_EQ( 1, field8.get_parameters(2)->get_channels() );
    ASSERT_EQ( 0, field8.get_parameters(0)->operator[](0) ); // check zero initialization
    ASSERT_EQ( 0, field8.get_parameters(0)->operator[](11) );
    ASSERT_EQ( 0, field8.get_parameters(1)->operator[](0) );
    ASSERT_EQ( 0, field8.get_parameters(1)->operator[](11) );
    ASSERT_EQ( 0, field8.get_parameters(2)->operator[](0) );
    ASSERT_EQ( 0, field8.get_parameters(2)->operator[](11) );

    auto img4 = image_cuda<TypeParam>::new_pointer(31,21,11);
    auto grid4 = grid_cuda<TypeParam>::new_pointer(*img4);
    dfield_cuda<TypeParam> field9(grid4);
    ASSERT_TRUE("dfield" == field9.get_name());
    ASSERT_EQ( 3, field9.get_dimension() );
    ASSERT_EQ( 31, field9.get_size()[0] );
    ASSERT_EQ( 21, field9.get_size()[1] );
    ASSERT_EQ( 11, field9.get_size()[2] );
    ASSERT_FALSE( nullptr == field9.get_parameters(0).get() );
    ASSERT_TRUE( "image" == field9.get_parameters(0)->get_name() );
    ASSERT_EQ( 31, field9.get_parameters(0)->get_width() );
    ASSERT_EQ( 21, field9.get_parameters(0)->get_height() );
    ASSERT_EQ( 11, field9.get_parameters(0)->get_length() );
    ASSERT_EQ( 1, field9.get_parameters(0)->get_channels() );
    ASSERT_FALSE( nullptr == field9.get_parameters(1).get() );
    ASSERT_TRUE( "image" == field9.get_parameters(1)->get_name() );
    ASSERT_EQ( 31, field9.get_parameters(1)->get_width() );
    ASSERT_EQ( 21, field9.get_parameters(1)->get_height() );
    ASSERT_EQ( 11, field9.get_parameters(1)->get_length() );
    ASSERT_EQ( 1, field9.get_parameters(1)->get_channels() );
    ASSERT_FALSE( nullptr == field9.get_parameters(2).get() );
    ASSERT_TRUE( "image" == field9.get_parameters(2)->get_name() );
    ASSERT_EQ( 31, field9.get_parameters(2)->get_width() );
    ASSERT_EQ( 21, field9.get_parameters(2)->get_height() );
    ASSERT_EQ( 11, field9.get_parameters(2)->get_length() );
    ASSERT_EQ( 1, field9.get_parameters(2)->get_channels() );
    ASSERT_EQ( 0, field9.get_parameters(0)->operator[](2) ); // check zero initialization
    ASSERT_EQ( 0, field9.get_parameters(0)->operator[](19) );
    ASSERT_EQ( 0, field9.get_parameters(1)->operator[](2) );
    ASSERT_EQ( 0, field9.get_parameters(1)->operator[](19) );
    ASSERT_EQ( 0, field9.get_parameters(2)->operator[](2) );
    ASSERT_EQ( 0, field9.get_parameters(2)->operator[](19) );

    dfield_cuda<TypeParam> field10(field9);
    ASSERT_TRUE("dfield" == field10.get_name());
    ASSERT_EQ( 3, field10.get_dimension() );
    ASSERT_EQ( 31, field10.get_size()[0] );
    ASSERT_EQ( 21, field10.get_size()[1] );
    ASSERT_EQ( 11, field10.get_size()[2] );
    ASSERT_FALSE( nullptr == field10.get_parameters(0).get() );
    ASSERT_FALSE( field9.get_parameters(0).get() == field10.get_parameters(0).get() );
    ASSERT_TRUE( "image" == field10.get_parameters(0)->get_name() );
    ASSERT_EQ( 31, field10.get_parameters(0)->get_width() );
    ASSERT_EQ( 21, field10.get_parameters(0)->get_height() );
    ASSERT_EQ( 11, field10.get_parameters(0)->get_length() );
    ASSERT_EQ( 1, field10.get_parameters(0)->get_channels() );
    ASSERT_FALSE( nullptr == field10.get_parameters(1).get() );
    ASSERT_FALSE( field9.get_parameters(1).get() == field10.get_parameters(1).get() );
    ASSERT_TRUE( "image" == field10.get_parameters(1)->get_name() );
    ASSERT_EQ( 31, field10.get_parameters(1)->get_width() );
    ASSERT_EQ( 21, field10.get_parameters(1)->get_height() );
    ASSERT_EQ( 11, field10.get_parameters(1)->get_length() );
    ASSERT_EQ( 1, field10.get_parameters(1)->get_channels() );
    ASSERT_FALSE( nullptr == field10.get_parameters(2).get() );
    ASSERT_FALSE( field9.get_parameters(2).get() == field10.get_parameters(2).get() );
    ASSERT_TRUE( "image" == field10.get_parameters(2)->get_name() );
    ASSERT_EQ( 31, field10.get_parameters(2)->get_width() );
    ASSERT_EQ( 21, field10.get_parameters(2)->get_height() );
    ASSERT_EQ( 11, field10.get_parameters(2)->get_length() );
    ASSERT_EQ( 1, field10.get_parameters(2)->get_channels() );
    ASSERT_EQ( 0, field10.get_parameters(0)->operator[](2) ); // check zero initialization
    ASSERT_EQ( 0, field10.get_parameters(0)->operator[](19) );
    ASSERT_EQ( 0, field10.get_parameters(1)->operator[](2) );
    ASSERT_EQ( 0, field10.get_parameters(1)->operator[](19) );
    ASSERT_EQ( 0, field10.get_parameters(2)->operator[](2) );
    ASSERT_EQ( 0, field10.get_parameters(2)->operator[](19) );
}

// ============================================
//          Testing Pointers
// ============================================
TYPED_TEST(test_dfield_cuda, pointers)
{
    using fieldcuda_pointer = typename dfield_cuda<TypeParam>::pointer;
    fieldcuda_pointer field1 = dfield_cuda<TypeParam>::new_pointer();
    ASSERT_TRUE("dfield" == field1->get_name());
    ASSERT_EQ( 2, field1->get_dimension() );
    ASSERT_EQ( 0, field1->get_size()[0] );
    ASSERT_EQ( 0, field1->get_size()[1] );
    ASSERT_FALSE( nullptr == field1->get_parameters(0).get() );
    ASSERT_TRUE( "image" == field1->get_parameters(0)->get_name() );
    ASSERT_FALSE( nullptr == field1->get_parameters(1).get() );
    ASSERT_TRUE( "image" == field1->get_parameters(1)->get_name() );

    fieldcuda_pointer field2 = dfield_cuda<TypeParam>::new_pointer(2);
    ASSERT_TRUE("dfield" == field2->get_name());
    ASSERT_EQ( 2, field2->get_dimension() );
    ASSERT_EQ( 0, field2->get_size()[0] );
    ASSERT_EQ( 0, field2->get_size()[1] );
    ASSERT_FALSE( nullptr == field2->get_parameters(0).get() );
    ASSERT_TRUE( "image" == field2->get_parameters(0)->get_name() );
    ASSERT_FALSE( nullptr == field2->get_parameters(1).get() );
    ASSERT_TRUE( "image" == field2->get_parameters(1)->get_name() );

    fieldcuda_pointer field3 = dfield_cuda<TypeParam>::new_pointer(3);
    ASSERT_TRUE("dfield" == field3->get_name());
    ASSERT_EQ( 3, field3->get_dimension() );
    ASSERT_EQ( 0, field3->get_size()[0] );
    ASSERT_EQ( 0, field3->get_size()[1] );
    ASSERT_EQ( 0, field3->get_size()[2] );
    ASSERT_FALSE( nullptr == field3->get_parameters(0).get() );
    ASSERT_TRUE( "image" == field3->get_parameters(0)->get_name() );
    ASSERT_FALSE( nullptr == field3->get_parameters(1).get() );
    ASSERT_TRUE( "image" == field3->get_parameters(1)->get_name() );
    ASSERT_FALSE( nullptr == field3->get_parameters(2).get() );
    ASSERT_TRUE( "image" == field3->get_parameters(2)->get_name() );

    fieldcuda_pointer field4 = dfield_cuda<TypeParam>::new_pointer(std::vector<int>{7,8});
    ASSERT_TRUE("dfield" == field4->get_name());
    ASSERT_EQ( 2, field4->get_dimension() );
    ASSERT_EQ( 7, field4->get_size()[0] );
    ASSERT_EQ( 8, field4->get_size()[1] );
    ASSERT_FALSE( nullptr == field4->get_parameters(0).get() );
    ASSERT_TRUE( "image" == field4->get_parameters(0)->get_name() );
    ASSERT_FALSE( nullptr == field4->get_parameters(1).get() );
    ASSERT_TRUE( "image" == field4->get_parameters(1)->get_name() );
    ASSERT_EQ( 7*8*1, field4->get_parameters(0)->get_total_elements() );
    ASSERT_EQ( 7, field4->get_parameters(0)->get_width() );
    ASSERT_EQ( 8, field4->get_parameters(0)->get_height() );
    ASSERT_EQ( 1, field4->get_parameters(0)->get_length() );
    ASSERT_EQ( 1, field4->get_parameters(0)->get_channels() );
    ASSERT_EQ( 7*8*1, field4->get_parameters(1)->get_total_elements() );
    ASSERT_EQ( 7, field4->get_parameters(1)->get_width() );
    ASSERT_EQ( 8, field4->get_parameters(1)->get_height() );
    ASSERT_EQ( 1, field4->get_parameters(1)->get_length() );
    ASSERT_EQ( 1, field4->get_parameters(1)->get_channels() );
    ASSERT_EQ( 0, field4->get_parameters(0)->operator[](3) ); // check zero initialization
    ASSERT_EQ( 0, field4->get_parameters(0)->operator[](25) );
    ASSERT_EQ( 0, field4->get_parameters(1)->operator[](3) );
    ASSERT_EQ( 0, field4->get_parameters(1)->operator[](25) );

    fieldcuda_pointer field5 = dfield_cuda<TypeParam>::new_pointer(std::vector<int>{9,5,3});
    ASSERT_TRUE("dfield" == field5->get_name());
    ASSERT_EQ( 3, field5->get_dimension() );
    ASSERT_EQ( 9, field5->get_size()[0] );
    ASSERT_EQ( 5, field5->get_size()[1] );
    ASSERT_EQ( 3, field5->get_size()[2] );
    ASSERT_FALSE( nullptr == field5->get_parameters(0).get() );
    ASSERT_TRUE( "image" == field5->get_parameters(0)->get_name() );
    ASSERT_FALSE( nullptr == field5->get_parameters(1).get() );
    ASSERT_TRUE( "image" == field5->get_parameters(1)->get_name() );
    ASSERT_FALSE( nullptr == field5->get_parameters(2).get() );
    ASSERT_TRUE( "image" == field5->get_parameters(2)->get_name() );
    ASSERT_EQ( 9*5*3, field5->get_parameters(0)->get_total_elements() );
    ASSERT_EQ( 9, field5->get_parameters(0)->get_width() );
    ASSERT_EQ( 5, field5->get_parameters(0)->get_height() );
    ASSERT_EQ( 3, field5->get_parameters(0)->get_length() );
    ASSERT_EQ( 1, field5->get_parameters(0)->get_channels() );
    ASSERT_EQ( 9*5*3, field5->get_parameters(1)->get_total_elements() );
    ASSERT_EQ( 9, field5->get_parameters(1)->get_width() );
    ASSERT_EQ( 5, field5->get_parameters(1)->get_height() );
    ASSERT_EQ( 3, field5->get_parameters(1)->get_length() );
    ASSERT_EQ( 1, field5->get_parameters(1)->get_channels() );
    ASSERT_EQ( 9*5*3, field5->get_parameters(2)->get_total_elements() );
    ASSERT_EQ( 9, field5->get_parameters(2)->get_width() );
    ASSERT_EQ( 5, field5->get_parameters(2)->get_height() );
    ASSERT_EQ( 3, field5->get_parameters(2)->get_length() );
    ASSERT_EQ( 1, field5->get_parameters(2)->get_channels() );
    ASSERT_EQ( 0, field5->get_parameters(0)->operator[](4) ); // check zero initialization
    ASSERT_EQ( 0, field5->get_parameters(0)->operator[](28) );
    ASSERT_EQ( 0, field5->get_parameters(1)->operator[](4) );
    ASSERT_EQ( 0, field5->get_parameters(1)->operator[](28) );
    ASSERT_EQ( 0, field5->get_parameters(2)->operator[](4) );
    ASSERT_EQ( 0, field5->get_parameters(2)->operator[](28) );
}

// ============================================
//          Testing Inherited Space
// ============================================
TYPED_TEST(test_dfield_cuda, inheritance)
{
    using fieldcuda_pointer = typename dfield_cuda<TypeParam>::pointer;
    fieldcuda_pointer field2 = dfield_cuda<TypeParam>::new_pointer(2);
    ASSERT_EQ( 2, field2->get_dimension() );
    ASSERT_EQ( 2, field2->get_size().size() );
    ASSERT_EQ( 0, field2->get_size()[0] );
    ASSERT_EQ( 0, field2->get_size()[1] );
    ASSERT_FLOAT_EQ( 1.0, field2->get_spacing()[0] );
    ASSERT_FLOAT_EQ( 1.0, field2->get_spacing()[1] );
    ASSERT_FLOAT_EQ( 0.0, field2->get_origin()[0] );
    ASSERT_FLOAT_EQ( 0.0, field2->get_origin()[1] );
    ASSERT_FLOAT_EQ( 1.0, field2->get_direction()[0] );
    ASSERT_FLOAT_EQ( 0.0, field2->get_direction()[1] );
    ASSERT_FLOAT_EQ( 0.0, field2->get_direction()[2] );
    ASSERT_FLOAT_EQ( 1.0, field2->get_direction()[3] );

    // change properties
    field2->set_spacing(std::vector<double>({2.2,1.1}));
    field2->set_origin(std::vector<double>({-10.2,8.1}));
    field2->set_direction(std::vector<double>({1.1,-0.1, 0.2, 0.9}));

    ASSERT_FLOAT_EQ( 2.2, field2->get_spacing()[0] );
    ASSERT_FLOAT_EQ( 1.1, field2->get_spacing()[1] );
    ASSERT_FLOAT_EQ( -10.2, field2->get_origin()[0] );
    ASSERT_FLOAT_EQ( 8.1, field2->get_origin()[1] );
    ASSERT_FLOAT_EQ( 1.1, field2->get_direction()[0] );
    ASSERT_FLOAT_EQ( -0.1, field2->get_direction()[1] );
    ASSERT_FLOAT_EQ( 0.2, field2->get_direction()[2] );
    ASSERT_FLOAT_EQ( 0.9, field2->get_direction()[3] );

    fieldcuda_pointer field3 = dfield_cuda<TypeParam>::new_pointer(3);
    ASSERT_EQ( 3, field3->get_dimension() );
    ASSERT_EQ( 3, field3->get_size().size() );
    ASSERT_EQ( 0, field3->get_size()[0] );
    ASSERT_EQ( 0, field3->get_size()[1] );
    ASSERT_EQ( 0, field3->get_size()[2] );
    ASSERT_FLOAT_EQ( 1.0, field3->get_spacing()[0] );
    ASSERT_FLOAT_EQ( 1.0, field3->get_spacing()[1] );
    ASSERT_FLOAT_EQ( 1.0, field3->get_spacing()[2] );
    ASSERT_FLOAT_EQ( 0.0, field3->get_origin()[0] );
    ASSERT_FLOAT_EQ( 0.0, field3->get_origin()[1] );
    ASSERT_FLOAT_EQ( 0.0, field3->get_origin()[2] );
    ASSERT_FLOAT_EQ( 1.0, field3->get_direction()[0] );
    ASSERT_FLOAT_EQ( 0.0, field3->get_direction()[1] );
    ASSERT_FLOAT_EQ( 0.0, field3->get_direction()[2] );
    ASSERT_FLOAT_EQ( 0.0, field3->get_direction()[3] );
    ASSERT_FLOAT_EQ( 1.0, field3->get_direction()[4] );
    ASSERT_FLOAT_EQ( 0.0, field3->get_direction()[5] );
    ASSERT_FLOAT_EQ( 0.0, field3->get_direction()[6] );
    ASSERT_FLOAT_EQ( 0.0, field3->get_direction()[7] );
    ASSERT_FLOAT_EQ( 1.0, field3->get_direction()[8] );

    // change properties
    field3->set_spacing(std::vector<double>({2.2,1.1,3.0}));
    field3->set_origin(std::vector<double>({-10.2,8.1,41.123}));
    field3->set_direction(std::vector<double>({1.1,-0.1, 0.2, 0.05, 0.9, -0.2, -0.12, 0.15, 1.0 }));

    ASSERT_FLOAT_EQ( 2.2, field3->get_spacing()[0] );
    ASSERT_FLOAT_EQ( 1.1, field3->get_spacing()[1] );
    ASSERT_FLOAT_EQ( 3.0, field3->get_spacing()[2] );
    ASSERT_FLOAT_EQ( -10.2, field3->get_origin()[0] );
    ASSERT_FLOAT_EQ( 8.1, field3->get_origin()[1] );
    ASSERT_FLOAT_EQ( 41.123, field3->get_origin()[2] );
    ASSERT_FLOAT_EQ( 1.1, field3->get_direction()[0] );
    ASSERT_FLOAT_EQ( -0.1, field3->get_direction()[1] );
    ASSERT_FLOAT_EQ( 0.2, field3->get_direction()[2] );
    ASSERT_FLOAT_EQ( 0.05, field3->get_direction()[3] );
    ASSERT_FLOAT_EQ( 0.9, field3->get_direction()[4] );
    ASSERT_FLOAT_EQ( -0.2, field3->get_direction()[5] );
    ASSERT_FLOAT_EQ( -0.12, field3->get_direction()[6] );
    ASSERT_FLOAT_EQ( 0.15, field3->get_direction()[7] );
    ASSERT_FLOAT_EQ( 1.0, field3->get_direction()[8] );
}

// ============================================
//          Testing Clone Methods
// ============================================
TYPED_TEST(test_dfield_cuda, clones)
{
    using fieldcuda_pointer = typename dfield_cuda<TypeParam>::pointer;
    fieldcuda_pointer field1 = dfield_cuda<TypeParam>::new_pointer(std::vector<int>{8,5});
    field1->get_parameters(0)->assign(3.14159);
    field1->get_parameters(1)->assign(2.71828);
    field1->set_spacing(std::vector<double>({3.2,1.1}));
    field1->set_origin(std::vector<double>({-0.5,0.5}));
    field1->set_direction(std::vector<double>({1.05,-0.2, -0.3, -0.9}));

    fieldcuda_pointer field2 = field1->clone();
    fieldcuda_pointer field3 = field1->copy();
    fieldcuda_pointer field4 = field1->mimic();

    ASSERT_TRUE("dfield" == field2->get_name());
    ASSERT_TRUE("dfield" == field3->get_name());
    ASSERT_TRUE("dfield" == field4->get_name());

    ASSERT_EQ( 40, field1->get_parameters(0)->get_total_elements() );
    ASSERT_EQ( 40, field2->get_parameters(0)->get_total_elements() );
    ASSERT_EQ( 40, field3->get_parameters(0)->get_total_elements() );
    ASSERT_EQ( 40, field4->get_parameters(0)->get_total_elements() );
    ASSERT_EQ( 8, field1->get_parameters(0)->get_width() );
    ASSERT_EQ( 8, field2->get_parameters(0)->get_width() );
    ASSERT_EQ( 8, field3->get_parameters(0)->get_width() );
    ASSERT_EQ( 8, field4->get_parameters(0)->get_width() );
    ASSERT_EQ( 5, field1->get_parameters(0)->get_height() );
    ASSERT_EQ( 5, field2->get_parameters(0)->get_height() );
    ASSERT_EQ( 5, field3->get_parameters(0)->get_height() );
    ASSERT_EQ( 5, field4->get_parameters(0)->get_height() );
    ASSERT_EQ( 1, field1->get_parameters(0)->get_length() );
    ASSERT_EQ( 1, field2->get_parameters(0)->get_length() );
    ASSERT_EQ( 1, field3->get_parameters(0)->get_length() );
    ASSERT_EQ( 1, field4->get_parameters(0)->get_length() );
    ASSERT_EQ( 1, field1->get_parameters(0)->get_channels() );
    ASSERT_EQ( 1, field2->get_parameters(0)->get_channels() );
    ASSERT_EQ( 1, field3->get_parameters(0)->get_channels() );
    ASSERT_EQ( 1, field4->get_parameters(0)->get_channels() );
    ASSERT_EQ( 40, field1->get_parameters(1)->get_total_elements() );
    ASSERT_EQ( 40, field2->get_parameters(1)->get_total_elements() );
    ASSERT_EQ( 40, field3->get_parameters(1)->get_total_elements() );
    ASSERT_EQ( 40, field4->get_parameters(1)->get_total_elements() );
    ASSERT_EQ( 8, field1->get_parameters(1)->get_width() );
    ASSERT_EQ( 8, field2->get_parameters(1)->get_width() );
    ASSERT_EQ( 8, field3->get_parameters(1)->get_width() );
    ASSERT_EQ( 8, field4->get_parameters(1)->get_width() );
    ASSERT_EQ( 5, field1->get_parameters(1)->get_height() );
    ASSERT_EQ( 5, field2->get_parameters(1)->get_height() );
    ASSERT_EQ( 5, field3->get_parameters(1)->get_height() );
    ASSERT_EQ( 5, field4->get_parameters(1)->get_height() );
    ASSERT_EQ( 1, field1->get_parameters(1)->get_length() );
    ASSERT_EQ( 1, field2->get_parameters(1)->get_length() );
    ASSERT_EQ( 1, field3->get_parameters(1)->get_length() );
    ASSERT_EQ( 1, field4->get_parameters(1)->get_length() );
    ASSERT_EQ( 1, field1->get_parameters(1)->get_channels() );
    ASSERT_EQ( 1, field2->get_parameters(1)->get_channels() );
    ASSERT_EQ( 1, field3->get_parameters(1)->get_channels() );
    ASSERT_EQ( 1, field4->get_parameters(1)->get_channels() );

    ASSERT_FLOAT_EQ( 3.2, field2->get_spacing()[0] );
    ASSERT_FLOAT_EQ( 1.1, field2->get_spacing()[1] );
    ASSERT_FLOAT_EQ( -0.5, field2->get_origin()[0] );
    ASSERT_FLOAT_EQ( 0.5, field2->get_origin()[1] );
    ASSERT_FLOAT_EQ( 1.05, field2->get_direction()[0] );
    ASSERT_FLOAT_EQ( -0.2, field2->get_direction()[1] );
    ASSERT_FLOAT_EQ( -0.3, field2->get_direction()[2] );
    ASSERT_FLOAT_EQ( -0.9, field2->get_direction()[3] );

    ASSERT_FLOAT_EQ( 3.2, field3->get_spacing()[0] );
    ASSERT_FLOAT_EQ( 1.1, field3->get_spacing()[1] );
    ASSERT_FLOAT_EQ( -0.5, field3->get_origin()[0] );
    ASSERT_FLOAT_EQ( 0.5, field3->get_origin()[1] );
    ASSERT_FLOAT_EQ( 1.05, field3->get_direction()[0] );
    ASSERT_FLOAT_EQ( -0.2, field3->get_direction()[1] );
    ASSERT_FLOAT_EQ( -0.3, field3->get_direction()[2] );
    ASSERT_FLOAT_EQ( -0.9, field3->get_direction()[3] );

    ASSERT_FLOAT_EQ( 3.2, field4->get_spacing()[0] );
    ASSERT_FLOAT_EQ( 1.1, field4->get_spacing()[1] );
    ASSERT_FLOAT_EQ( -0.5, field4->get_origin()[0] );
    ASSERT_FLOAT_EQ( 0.5, field4->get_origin()[1] );
    ASSERT_FLOAT_EQ( 1.05, field4->get_direction()[0] );
    ASSERT_FLOAT_EQ( -0.2, field4->get_direction()[1] );
    ASSERT_FLOAT_EQ( -0.3, field4->get_direction()[2] );
    ASSERT_FLOAT_EQ( -0.9, field4->get_direction()[3] );

    ASSERT_FALSE(field1.get() == field2.get());
    ASSERT_FALSE(field1.get() == field3.get());
    ASSERT_FALSE(field1.get() == field4.get());

    ASSERT_FALSE(field1->get_parameters_vector().get() == field2->get_parameters_vector().get());
    ASSERT_FALSE(field1->get_parameters_vector().get() == field3->get_parameters_vector().get());
    ASSERT_FALSE(field1->get_parameters_vector().get() == field4->get_parameters_vector().get());

    ASSERT_FALSE(field1->get_parameters(0)->get_data() == field2->get_parameters(0)->get_data());
    ASSERT_TRUE(field1->get_parameters(0)->get_data() == field3->get_parameters(0)->get_data());
    ASSERT_FALSE(field1->get_parameters(0)->get_data() == field4->get_parameters(0)->get_data());
    ASSERT_FALSE(field1->get_parameters(1)->get_data() == field2->get_parameters(1)->get_data());
    ASSERT_TRUE(field1->get_parameters(1)->get_data() == field3->get_parameters(1)->get_data());
    ASSERT_FALSE(field1->get_parameters(1)->get_data() == field4->get_parameters(1)->get_data());

    ASSERT_FLOAT_EQ( 3.14159, field1->get_parameters(0)->operator[](0) );
    ASSERT_FLOAT_EQ( 3.14159, field1->get_parameters(0)->operator[](7) );
    ASSERT_FLOAT_EQ( 3.14159, field1->get_parameters(0)->operator[](39) );
    ASSERT_FLOAT_EQ( 3.14159, field2->get_parameters(0)->operator[](0) );
    ASSERT_FLOAT_EQ( 3.14159, field2->get_parameters(0)->operator[](7) );
    ASSERT_FLOAT_EQ( 3.14159, field2->get_parameters(0)->operator[](39) );
    ASSERT_FLOAT_EQ( 3.14159, field3->get_parameters(0)->operator[](0) );
    ASSERT_FLOAT_EQ( 3.14159, field3->get_parameters(0)->operator[](7) );
    ASSERT_FLOAT_EQ( 3.14159, field3->get_parameters(0)->operator[](39) );
    EXPECT_FLOAT_EQ( 0.0, field4->get_parameters(0)->operator[](0) );
    EXPECT_FLOAT_EQ( 0.0, field4->get_parameters(0)->operator[](7) );
    EXPECT_FLOAT_EQ( 0.0, field4->get_parameters(0)->operator[](39) );
    ASSERT_FLOAT_EQ( 2.71828, field1->get_parameters(1)->operator[](0) );
    ASSERT_FLOAT_EQ( 2.71828, field1->get_parameters(1)->operator[](7) );
    ASSERT_FLOAT_EQ( 2.71828, field1->get_parameters(1)->operator[](39) );
    ASSERT_FLOAT_EQ( 2.71828, field2->get_parameters(1)->operator[](0) );
    ASSERT_FLOAT_EQ( 2.71828, field2->get_parameters(1)->operator[](7) );
    ASSERT_FLOAT_EQ( 2.71828, field2->get_parameters(1)->operator[](39) );
    ASSERT_FLOAT_EQ( 2.71828, field3->get_parameters(1)->operator[](0) );
    ASSERT_FLOAT_EQ( 2.71828, field3->get_parameters(1)->operator[](7) );
    ASSERT_FLOAT_EQ( 2.71828, field3->get_parameters(1)->operator[](39) );
    EXPECT_FLOAT_EQ( 0.0, field4->get_parameters(1)->operator[](0) );
    EXPECT_FLOAT_EQ( 0.0, field4->get_parameters(1)->operator[](7) );
    EXPECT_FLOAT_EQ( 0.0, field4->get_parameters(1)->operator[](39) );


    //transform clone dfield
    using transformcuda_pointer = typename transform_cuda<TypeParam>::pointer;
    transformcuda_pointer field5 = field1->clone();
    transformcuda_pointer field6 = field1->copy();
    transformcuda_pointer field7 = field1->mimic();

    ASSERT_TRUE("dfield" == field5->get_name());
    ASSERT_TRUE("dfield" == field6->get_name());
    ASSERT_TRUE("dfield" == field7->get_name());

    ASSERT_EQ( 40, field5->get_parameters(0)->get_total_elements() );
    ASSERT_EQ( 40, field6->get_parameters(0)->get_total_elements() );
    ASSERT_EQ( 40, field7->get_parameters(0)->get_total_elements() );
    ASSERT_EQ( 8, field5->get_parameters(0)->get_width() );
    ASSERT_EQ( 8, field6->get_parameters(0)->get_width() );
    ASSERT_EQ( 8, field7->get_parameters(0)->get_width() );
    ASSERT_EQ( 5, field5->get_parameters(0)->get_height() );
    ASSERT_EQ( 5, field6->get_parameters(0)->get_height() );
    ASSERT_EQ( 5, field7->get_parameters(0)->get_height() );
    ASSERT_EQ( 1, field5->get_parameters(0)->get_length() );
    ASSERT_EQ( 1, field6->get_parameters(0)->get_length() );
    ASSERT_EQ( 1, field7->get_parameters(0)->get_length() );
    ASSERT_EQ( 1, field5->get_parameters(0)->get_channels() );
    ASSERT_EQ( 1, field6->get_parameters(0)->get_channels() );
    ASSERT_EQ( 1, field7->get_parameters(0)->get_channels() );
    ASSERT_EQ( 40, field5->get_parameters(1)->get_total_elements() );
    ASSERT_EQ( 40, field6->get_parameters(1)->get_total_elements() );
    ASSERT_EQ( 40, field7->get_parameters(1)->get_total_elements() );
    ASSERT_EQ( 8, field5->get_parameters(1)->get_width() );
    ASSERT_EQ( 8, field6->get_parameters(1)->get_width() );
    ASSERT_EQ( 8, field7->get_parameters(1)->get_width() );
    ASSERT_EQ( 5, field5->get_parameters(1)->get_height() );
    ASSERT_EQ( 5, field6->get_parameters(1)->get_height() );
    ASSERT_EQ( 5, field7->get_parameters(1)->get_height() );
    ASSERT_EQ( 1, field5->get_parameters(1)->get_length() );
    ASSERT_EQ( 1, field6->get_parameters(1)->get_length() );
    ASSERT_EQ( 1, field7->get_parameters(1)->get_length() );
    ASSERT_EQ( 1, field5->get_parameters(1)->get_channels() );
    ASSERT_EQ( 1, field6->get_parameters(1)->get_channels() );
    ASSERT_EQ( 1, field7->get_parameters(1)->get_channels() );

    ASSERT_FLOAT_EQ( 3.2, field5->get_spacing()[0] );
    ASSERT_FLOAT_EQ( 1.1, field5->get_spacing()[1] );
    ASSERT_FLOAT_EQ( -0.5, field5->get_origin()[0] );
    ASSERT_FLOAT_EQ( 0.5, field5->get_origin()[1] );
    ASSERT_FLOAT_EQ( 1.05, field5->get_direction()[0] );
    ASSERT_FLOAT_EQ( -0.2, field5->get_direction()[1] );
    ASSERT_FLOAT_EQ( -0.3, field5->get_direction()[2] );
    ASSERT_FLOAT_EQ( -0.9, field5->get_direction()[3] );

    ASSERT_FLOAT_EQ( 3.2, field6->get_spacing()[0] );
    ASSERT_FLOAT_EQ( 1.1, field6->get_spacing()[1] );
    ASSERT_FLOAT_EQ( -0.5, field6->get_origin()[0] );
    ASSERT_FLOAT_EQ( 0.5, field6->get_origin()[1] );
    ASSERT_FLOAT_EQ( 1.05, field6->get_direction()[0] );
    ASSERT_FLOAT_EQ( -0.2, field6->get_direction()[1] );
    ASSERT_FLOAT_EQ( -0.3, field6->get_direction()[2] );
    ASSERT_FLOAT_EQ( -0.9, field6->get_direction()[3] );

    ASSERT_FLOAT_EQ( 3.2, field7->get_spacing()[0] );
    ASSERT_FLOAT_EQ( 1.1, field7->get_spacing()[1] );
    ASSERT_FLOAT_EQ( -0.5, field7->get_origin()[0] );
    ASSERT_FLOAT_EQ( 0.5, field7->get_origin()[1] );
    ASSERT_FLOAT_EQ( 1.05, field7->get_direction()[0] );
    ASSERT_FLOAT_EQ( -0.2, field7->get_direction()[1] );
    ASSERT_FLOAT_EQ( -0.3, field7->get_direction()[2] );
    ASSERT_FLOAT_EQ( -0.9, field7->get_direction()[3] );

    ASSERT_FALSE(field1.get() == field5.get());
    ASSERT_FALSE(field1.get() == field6.get());
    ASSERT_FALSE(field1.get() == field7.get());

    ASSERT_FALSE(field1->get_parameters_vector().get() == field5->get_parameters_vector().get());
    ASSERT_FALSE(field1->get_parameters_vector().get() == field6->get_parameters_vector().get());
    ASSERT_FALSE(field1->get_parameters_vector().get() == field7->get_parameters_vector().get());

    ASSERT_FALSE(field1->get_parameters(0)->get_data() == field5->get_parameters(0)->get_data());
    ASSERT_TRUE(field1->get_parameters(0)->get_data() == field6->get_parameters(0)->get_data());
    ASSERT_FALSE(field1->get_parameters(0)->get_data() == field7->get_parameters(0)->get_data());
    ASSERT_FALSE(field1->get_parameters(1)->get_data() == field5->get_parameters(1)->get_data());
    ASSERT_TRUE(field1->get_parameters(1)->get_data() == field6->get_parameters(1)->get_data());
    ASSERT_FALSE(field1->get_parameters(1)->get_data() == field7->get_parameters(1)->get_data());

    ASSERT_FLOAT_EQ( 3.14159, field5->get_parameters(0)->operator[](0) );
    ASSERT_FLOAT_EQ( 3.14159, field5->get_parameters(0)->operator[](7) );
    ASSERT_FLOAT_EQ( 3.14159, field5->get_parameters(0)->operator[](39) );
    ASSERT_FLOAT_EQ( 3.14159, field6->get_parameters(0)->operator[](0) );
    ASSERT_FLOAT_EQ( 3.14159, field6->get_parameters(0)->operator[](7) );
    ASSERT_FLOAT_EQ( 3.14159, field6->get_parameters(0)->operator[](39) );
    EXPECT_FLOAT_EQ( 0.0, field7->get_parameters(0)->operator[](0) );
    EXPECT_FLOAT_EQ( 0.0, field7->get_parameters(0)->operator[](7) );
    EXPECT_FLOAT_EQ( 0.0, field7->get_parameters(0)->operator[](39) );
    ASSERT_FLOAT_EQ( 2.71828, field5->get_parameters(1)->operator[](0) );
    ASSERT_FLOAT_EQ( 2.71828, field5->get_parameters(1)->operator[](7) );
    ASSERT_FLOAT_EQ( 2.71828, field5->get_parameters(1)->operator[](39) );
    ASSERT_FLOAT_EQ( 2.71828, field6->get_parameters(1)->operator[](0) );
    ASSERT_FLOAT_EQ( 2.71828, field6->get_parameters(1)->operator[](7) );
    ASSERT_FLOAT_EQ( 2.71828, field6->get_parameters(1)->operator[](39) );
    EXPECT_FLOAT_EQ( 0.0, field7->get_parameters(1)->operator[](0) );
    EXPECT_FLOAT_EQ( 0.0, field7->get_parameters(1)->operator[](7) );
    EXPECT_FLOAT_EQ( 0.0, field7->get_parameters(1)->operator[](39) );


    // second mimic
    transformcuda_pointer field8 = field6->clone();
    transformcuda_pointer field9 = field6->copy();
    transformcuda_pointer field10 = field6->mimic();

    ASSERT_TRUE("dfield" == field8->get_name());
    ASSERT_TRUE("dfield" == field9->get_name());
    ASSERT_TRUE("dfield" == field10->get_name());

    ASSERT_EQ( 40, field8->get_parameters(0)->get_total_elements() );
    ASSERT_EQ( 40, field9->get_parameters(0)->get_total_elements() );
    ASSERT_EQ( 40, field10->get_parameters(0)->get_total_elements() );
    ASSERT_EQ( 8, field8->get_parameters(0)->get_width() );
    ASSERT_EQ( 8, field9->get_parameters(0)->get_width() );
    ASSERT_EQ( 8, field10->get_parameters(0)->get_width() );
    ASSERT_EQ( 5, field8->get_parameters(0)->get_height() );
    ASSERT_EQ( 5, field9->get_parameters(0)->get_height() );
    ASSERT_EQ( 5, field10->get_parameters(0)->get_height() );
    ASSERT_EQ( 1, field8->get_parameters(0)->get_length() );
    ASSERT_EQ( 1, field9->get_parameters(0)->get_length() );
    ASSERT_EQ( 1, field10->get_parameters(0)->get_length() );
    ASSERT_EQ( 1, field8->get_parameters(0)->get_channels() );
    ASSERT_EQ( 1, field9->get_parameters(0)->get_channels() );
    ASSERT_EQ( 1, field10->get_parameters(0)->get_channels() );
    ASSERT_EQ( 40, field8->get_parameters(1)->get_total_elements() );
    ASSERT_EQ( 40, field9->get_parameters(1)->get_total_elements() );
    ASSERT_EQ( 40, field10->get_parameters(1)->get_total_elements() );
    ASSERT_EQ( 8, field8->get_parameters(1)->get_width() );
    ASSERT_EQ( 8, field9->get_parameters(1)->get_width() );
    ASSERT_EQ( 8, field10->get_parameters(1)->get_width() );
    ASSERT_EQ( 5, field8->get_parameters(1)->get_height() );
    ASSERT_EQ( 5, field9->get_parameters(1)->get_height() );
    ASSERT_EQ( 5, field10->get_parameters(1)->get_height() );
    ASSERT_EQ( 1, field8->get_parameters(1)->get_length() );
    ASSERT_EQ( 1, field9->get_parameters(1)->get_length() );
    ASSERT_EQ( 1, field10->get_parameters(1)->get_length() );
    ASSERT_EQ( 1, field8->get_parameters(1)->get_channels() );
    ASSERT_EQ( 1, field9->get_parameters(1)->get_channels() );
    ASSERT_EQ( 1, field10->get_parameters(1)->get_channels() );

    ASSERT_FLOAT_EQ( 3.2, field8->get_spacing()[0] );
    ASSERT_FLOAT_EQ( 1.1, field8->get_spacing()[1] );
    ASSERT_FLOAT_EQ( -0.5, field8->get_origin()[0] );
    ASSERT_FLOAT_EQ( 0.5, field8->get_origin()[1] );
    ASSERT_FLOAT_EQ( 1.05, field8->get_direction()[0] );
    ASSERT_FLOAT_EQ( -0.2, field8->get_direction()[1] );
    ASSERT_FLOAT_EQ( -0.3, field8->get_direction()[2] );
    ASSERT_FLOAT_EQ( -0.9, field8->get_direction()[3] );

    ASSERT_FLOAT_EQ( 3.2, field9->get_spacing()[0] );
    ASSERT_FLOAT_EQ( 1.1, field9->get_spacing()[1] );
    ASSERT_FLOAT_EQ( -0.5, field9->get_origin()[0] );
    ASSERT_FLOAT_EQ( 0.5, field9->get_origin()[1] );
    ASSERT_FLOAT_EQ( 1.05, field9->get_direction()[0] );
    ASSERT_FLOAT_EQ( -0.2, field9->get_direction()[1] );
    ASSERT_FLOAT_EQ( -0.3, field9->get_direction()[2] );
    ASSERT_FLOAT_EQ( -0.9, field9->get_direction()[3] );

    ASSERT_FLOAT_EQ( 3.2, field10->get_spacing()[0] );
    ASSERT_FLOAT_EQ( 1.1, field10->get_spacing()[1] );
    ASSERT_FLOAT_EQ( -0.5, field10->get_origin()[0] );
    ASSERT_FLOAT_EQ( 0.5, field10->get_origin()[1] );
    ASSERT_FLOAT_EQ( 1.05, field10->get_direction()[0] );
    ASSERT_FLOAT_EQ( -0.2, field10->get_direction()[1] );
    ASSERT_FLOAT_EQ( -0.3, field10->get_direction()[2] );
    ASSERT_FLOAT_EQ( -0.9, field10->get_direction()[3] );

    ASSERT_FALSE(field1.get() == field8.get());
    ASSERT_FALSE(field1.get() == field9.get());
    ASSERT_FALSE(field1.get() == field10.get());

    ASSERT_FALSE(field1->get_parameters_vector().get() == field8->get_parameters_vector().get());
    ASSERT_FALSE(field1->get_parameters_vector().get() == field9->get_parameters_vector().get());
    ASSERT_FALSE(field1->get_parameters_vector().get() == field10->get_parameters_vector().get());

    ASSERT_FALSE(field1->get_parameters(0)->get_data() == field8->get_parameters(0)->get_data());
    ASSERT_TRUE(field1->get_parameters(0)->get_data() == field9->get_parameters(0)->get_data());
    ASSERT_FALSE(field1->get_parameters(0)->get_data() == field10->get_parameters(0)->get_data());
    ASSERT_FALSE(field1->get_parameters(1)->get_data() == field8->get_parameters(1)->get_data());
    ASSERT_TRUE(field1->get_parameters(1)->get_data() == field9->get_parameters(1)->get_data());
    ASSERT_FALSE(field1->get_parameters(1)->get_data() == field10->get_parameters(1)->get_data());

    ASSERT_FLOAT_EQ( 3.14159, field8->get_parameters(0)->operator[](0) );
    ASSERT_FLOAT_EQ( 3.14159, field8->get_parameters(0)->operator[](7) );
    ASSERT_FLOAT_EQ( 3.14159, field8->get_parameters(0)->operator[](39) );
    ASSERT_FLOAT_EQ( 3.14159, field9->get_parameters(0)->operator[](0) );
    ASSERT_FLOAT_EQ( 3.14159, field9->get_parameters(0)->operator[](7) );
    ASSERT_FLOAT_EQ( 3.14159, field9->get_parameters(0)->operator[](39) );
    EXPECT_FLOAT_EQ( 0.0, field10->get_parameters(0)->operator[](0) );
    EXPECT_FLOAT_EQ( 0.0, field10->get_parameters(0)->operator[](7) );
    EXPECT_FLOAT_EQ( 0.0, field10->get_parameters(0)->operator[](39) );
    ASSERT_FLOAT_EQ( 2.71828, field8->get_parameters(1)->operator[](0) );
    ASSERT_FLOAT_EQ( 2.71828, field8->get_parameters(1)->operator[](7) );
    ASSERT_FLOAT_EQ( 2.71828, field8->get_parameters(1)->operator[](39) );
    ASSERT_FLOAT_EQ( 2.71828, field9->get_parameters(1)->operator[](0) );
    ASSERT_FLOAT_EQ( 2.71828, field9->get_parameters(1)->operator[](7) );
    ASSERT_FLOAT_EQ( 2.71828, field9->get_parameters(1)->operator[](39) );
    EXPECT_FLOAT_EQ( 0.0, field10->get_parameters(1)->operator[](0) );
    EXPECT_FLOAT_EQ( 0.0, field10->get_parameters(1)->operator[](7) );
    EXPECT_FLOAT_EQ( 0.0, field10->get_parameters(1)->operator[](39) );
}
/*

// ============================================
//          Testing Initialization Methods
// ============================================
TYPED_TEST(test_dfield_cuda, initialization)
{
    int size = 50;
    using fieldcuda_pointer = typename dfield_cuda<TypeParam>::pointer;
    fieldcuda_pointer field1 = dfield_cuda<TypeParam>::new_pointer(size);
    fieldcuda_pointer field2 = field1->mimic();
    fieldcuda_pointer field3 = field1->mimic();

    field1->ones();
    field2->assign(TypeParam(3.0));
    field3->random();

    ASSERT_FLOAT_EQ( 1.0, field1->operator[](0) );
    ASSERT_FLOAT_EQ( 1.0, field1->operator[](int(size/2)) );
    ASSERT_FLOAT_EQ( 1.0, field1->operator[](size-4) );
    ASSERT_FLOAT_EQ( 3.0, field2->operator[](0) );
    ASSERT_FLOAT_EQ( 3.0, field2->operator[](int(size/2)) );
    ASSERT_FLOAT_EQ( 3.0, field2->operator[](size-4) );
    ASSERT_FALSE( field3->operator[](0) == field3->operator[](1) );
    ASSERT_FALSE( field3->operator[](0) == field3->operator[](int(size/2)) );
    ASSERT_FALSE( field3->operator[](int(size/2)) == field3->operator[](size-4) );

    field1->zeros();
    field2->assign(TypeParam(-2.431));

    ASSERT_FLOAT_EQ( 0.0, field1->operator[](0) );
    ASSERT_FLOAT_EQ( 0.0, field1->operator[](int(size/2)) );
    ASSERT_FLOAT_EQ( 0.0, field1->operator[](size-4) );
    ASSERT_FLOAT_EQ( -2.431, field2->operator[](0) );
    ASSERT_FLOAT_EQ( -2.431, field2->operator[](int(size/2)) );
    ASSERT_FLOAT_EQ( -2.431, field2->operator[](size-4) );
    ASSERT_FALSE( field3->operator[](0) == field3->operator[](1) );
    ASSERT_FALSE( field3->operator[](0) == field3->operator[](int(size/2)) );
    ASSERT_FALSE( field3->operator[](int(size/2)) == field3->operator[](size-4) );
}
/*
// ============================================
//          Testing Vector Operations
// ============================================
TYPED_TEST(test_dfield_cuda, vector_operations)
{
    using fieldcuda_pointer = typename dfield_cuda<TypeParam>::pointer;
    fieldcuda_pointer field1 = dfield_cuda<TypeParam>::new_pointer(30,2.1);
    fieldcuda_pointer field2 = dfield_cuda<TypeParam>::new_pointer(30,1.1);
    fieldcuda_pointer field3;
    fieldcuda_pointer field4 = dfield_cuda<TypeParam>::new_pointer(30,4.2);
    fieldcuda_pointer field5;
    fieldcuda_pointer field6 = dfield_cuda<TypeParam>::new_pointer(30,14.2123);
    fieldcuda_pointer vec7;
    fieldcuda_pointer vec8 = dfield_cuda<TypeParam>::new_pointer(30,-2*14.2123);
    fieldcuda_pointer vec9;
    fieldcuda_pointer field10 = dfield_cuda<TypeParam>::new_pointer(30,4.0);
    fieldcuda_pointer field11;

    field3 = *field1 + *field2;
    ASSERT_TRUE("dfield" == field3->get_name());
    ASSERT_EQ( 30, field3->size() );
    ASSERT_FLOAT_EQ( 3.2, field3->operator[](0) );
    ASSERT_FLOAT_EQ( 3.2, field3->operator[](7) );
    ASSERT_FLOAT_EQ( 3.2, field3->operator[](18) );
    ASSERT_FLOAT_EQ( 3.2, field3->operator[](29) );

    field5 = *field3 - *field4;
    ASSERT_TRUE("dfield" == field5->get_name());
    ASSERT_EQ( 30, field5->size() );
    ASSERT_FLOAT_EQ( -1.0, field5->operator[](0) );
    ASSERT_FLOAT_EQ( -1.0, field5->operator[](9) );
    ASSERT_FLOAT_EQ( -1.0, field5->operator[](16) );
    ASSERT_FLOAT_EQ( -1.0, field5->operator[](29) );

    vec7 = (*field5) * (*field6);
    ASSERT_TRUE("dfield" == vec7->get_name());
    ASSERT_EQ( 30, vec7->size() );
    ASSERT_FLOAT_EQ( -14.2123, vec7->operator[](0) );
    ASSERT_FLOAT_EQ( -14.2123, vec7->operator[](3) );
    ASSERT_FLOAT_EQ( -14.2123, vec7->operator[](12) );
    ASSERT_FLOAT_EQ( -14.2123, vec7->operator[](29) );

    vec9 = (*vec7) / (*vec8);
    ASSERT_TRUE("dfield" == vec9->get_name());
    ASSERT_EQ( 30, vec9->size() );
    ASSERT_FLOAT_EQ( 0.5, vec9->operator[](0) );
    ASSERT_FLOAT_EQ( 0.5, vec9->operator[](4) );
    ASSERT_FLOAT_EQ( 0.5, vec9->operator[](21) );
    ASSERT_FLOAT_EQ( 0.5, vec9->operator[](29) );

    field11 = (*vec9) ^ (*field10);
    ASSERT_TRUE("dfield" == field11->get_name());
    ASSERT_EQ( 30, field11->size() );
    ASSERT_FLOAT_EQ( 0.0625, field11->operator[](0) );
    ASSERT_FLOAT_EQ( 0.0625, field11->operator[](5) );
    ASSERT_FLOAT_EQ( 0.0625, field11->operator[](23) );
    ASSERT_FLOAT_EQ( 0.0625, field11->operator[](29) );

    // Check vector values again after operations
    ASSERT_FLOAT_EQ( 2.1, field1->operator[](0) );
    ASSERT_FLOAT_EQ( 2.1, field1->operator[](29) );
    ASSERT_FLOAT_EQ( 1.1, field2->operator[](0) );
    ASSERT_FLOAT_EQ( 1.1, field2->operator[](29) );
    ASSERT_FLOAT_EQ( 3.2, field3->operator[](0) );
    ASSERT_FLOAT_EQ( 3.2, field3->operator[](29) );
    ASSERT_FLOAT_EQ( 4.2, field4->operator[](0) );
    ASSERT_FLOAT_EQ( 4.2, field4->operator[](29) );
    ASSERT_FLOAT_EQ( -1.0, field5->operator[](0) );
    ASSERT_FLOAT_EQ( -1.0, field5->operator[](29) );
    ASSERT_FLOAT_EQ( 14.2123, field6->operator[](0) );
    ASSERT_FLOAT_EQ( 14.2123, field6->operator[](29) );
    ASSERT_FLOAT_EQ( -14.2123, vec7->operator[](0) );
    ASSERT_FLOAT_EQ( -14.2123, vec7->operator[](29) );
    ASSERT_FLOAT_EQ( -2*14.2123, vec8->operator[](0) );
    ASSERT_FLOAT_EQ( -2*14.2123, vec8->operator[](29) );
    ASSERT_FLOAT_EQ( 0.5, vec9->operator[](0) );
    ASSERT_FLOAT_EQ( 0.5, vec9->operator[](29) );
    ASSERT_FLOAT_EQ( 4.0, field10->operator[](0) );
    ASSERT_FLOAT_EQ( 4.0, field10->operator[](29) );
}

// ============================================
//          Testing Scalar Operations
// ============================================
TYPED_TEST(test_dfield_cuda, scalar_operations)
{
    using fieldcuda_pointer = typename dfield_cuda<TypeParam>::pointer;
    fieldcuda_pointer field1 = dfield_cuda<TypeParam>::new_pointer(14,1.5);
    fieldcuda_pointer field2 = dfield_cuda<TypeParam>::new_pointer(16,-4.632);
    fieldcuda_pointer field3 = dfield_cuda<TypeParam>::new_pointer(18,-12.0);
    fieldcuda_pointer field4 = dfield_cuda<TypeParam>::new_pointer(20,8.4);
    fieldcuda_pointer field5 = dfield_cuda<TypeParam>::new_pointer(22,2.0);
    fieldcuda_pointer field6;
    fieldcuda_pointer vec7;
    fieldcuda_pointer vec8;
    fieldcuda_pointer vec9;
    fieldcuda_pointer field10;

    // Right hand side
    field6 = *field1 + 4.5;
    ASSERT_TRUE("dfield" == field6->get_name());
    ASSERT_EQ( 14, field6->size() );
    ASSERT_FLOAT_EQ( 6.0, field6->operator[](0) );
    ASSERT_FLOAT_EQ( 6.0, field6->operator[](7) );
    ASSERT_FLOAT_EQ( 6.0, field6->operator[](13) );

    vec7 = *field2 - 1.111;
    ASSERT_TRUE("dfield" == vec7->get_name());
    ASSERT_EQ( 16, vec7->size() );
    ASSERT_FLOAT_EQ( -5.743, vec7->operator[](0) );
    ASSERT_FLOAT_EQ( -5.743, vec7->operator[](6) );
    ASSERT_FLOAT_EQ( -5.743, vec7->operator[](15) );

    vec8 = (*field3) * (-0.9);
    ASSERT_TRUE("dfield" == vec8->get_name());
    ASSERT_EQ( 18, vec8->size() );
    ASSERT_FLOAT_EQ( 10.8, vec8->operator[](0) );
    ASSERT_FLOAT_EQ( 10.8, vec8->operator[](13) );
    ASSERT_FLOAT_EQ( 10.8, vec8->operator[](17) );

    vec9 = (*field4) / (-4.0);
    ASSERT_TRUE("dfield" == vec9->get_name());
    ASSERT_EQ( 20, vec9->size() );
    ASSERT_FLOAT_EQ( -2.1, vec9->operator[](0) );
    ASSERT_FLOAT_EQ( -2.1, vec9->operator[](11) );
    ASSERT_FLOAT_EQ( -2.1, vec9->operator[](19) );

    field10 = (*field5) ^ (3.0);
    ASSERT_TRUE("dfield" == field10->get_name());
    ASSERT_EQ( 22, field10->size() );
    ASSERT_FLOAT_EQ( 8.0, field10->operator[](0) );
    ASSERT_FLOAT_EQ( 8.0, field10->operator[](8) );
    ASSERT_FLOAT_EQ( 8.0, field10->operator[](21) );

    // Left hand size
    field6 = TypeParam(4.5) + (*field1);
    ASSERT_TRUE("dfield" == field6->get_name());
    ASSERT_EQ( 14, field6->size() );
    ASSERT_FLOAT_EQ( 6.0, field6->operator[](0) );
    ASSERT_FLOAT_EQ( 6.0, field6->operator[](7) );
    ASSERT_FLOAT_EQ( 6.0, field6->operator[](13) );

    vec7 = TypeParam(1.111) - (*field2);
    ASSERT_TRUE("dfield" == vec7->get_name());
    ASSERT_EQ( 16, vec7->size() );
    ASSERT_FLOAT_EQ( 5.743, vec7->operator[](0) );
    ASSERT_FLOAT_EQ( 5.743, vec7->operator[](6) );
    ASSERT_FLOAT_EQ( 5.743, vec7->operator[](15) );

    vec8 = TypeParam(0.9)*(*field3);
    ASSERT_TRUE("dfield" == vec8->get_name());
    ASSERT_EQ( 18, vec8->size() );
    ASSERT_FLOAT_EQ( -10.8, vec8->operator[](0) );
    ASSERT_FLOAT_EQ( -10.8, vec8->operator[](13) );
    ASSERT_FLOAT_EQ( -10.8, vec8->operator[](17) );

    vec9 = TypeParam(-4.0)/(*field4);
    ASSERT_TRUE("dfield" == vec9->get_name());
    ASSERT_EQ( 20, vec9->size() );
    ASSERT_FLOAT_EQ( -0.47619051, vec9->operator[](0) );
    ASSERT_FLOAT_EQ( -0.47619051, vec9->operator[](11) );
    ASSERT_FLOAT_EQ( -0.47619051, vec9->operator[](19) );

    // Check original values after operations
    ASSERT_EQ( 14, field1->size() );
    ASSERT_FLOAT_EQ( 1.5, field1->operator[](0) );
    ASSERT_FLOAT_EQ( 1.5, field1->operator[](7) );
    ASSERT_FLOAT_EQ( 1.5, field1->operator[](13) );

    ASSERT_EQ( 16, field2->size() );
    ASSERT_FLOAT_EQ( -4.632, field2->operator[](0) );
    ASSERT_FLOAT_EQ( -4.632, field2->operator[](6) );
    ASSERT_FLOAT_EQ( -4.632, field2->operator[](15) );

    ASSERT_EQ( 18, field3->size() );
    ASSERT_FLOAT_EQ( -12.0, field3->operator[](0) );
    ASSERT_FLOAT_EQ( -12.0, field3->operator[](13) );
    ASSERT_FLOAT_EQ( -12.0, field3->operator[](17) );

    ASSERT_EQ( 20, vec9->size() );
    ASSERT_FLOAT_EQ( 8.4, field4->operator[](0) );
    ASSERT_FLOAT_EQ( 8.4, field4->operator[](11) );
    ASSERT_FLOAT_EQ( 8.4, field4->operator[](19) );

    ASSERT_EQ( 22, field5->size() );
    ASSERT_FLOAT_EQ( 2.0, field5->operator[](0) );
    ASSERT_FLOAT_EQ( 2.0, field5->operator[](8) );
    ASSERT_FLOAT_EQ( 2.0, field5->operator[](21) );
}

// ============================================
//          Testing Reduction Operations
// ============================================
TYPED_TEST(test_dfield_cuda, reduction_operations)
{
    int size = 32;
    using fieldcuda_pointer = typename dfield_cuda<TypeParam>::pointer;
    fieldcuda_pointer field1 = dfield_cuda<TypeParam>::new_pointer(size,0.0);
    fieldcuda_pointer field2 = dfield_cuda<TypeParam>::new_pointer(size,0.0);

    TypeParam dot1 = 0; TypeParam dot2 = 0;
    TypeParam add1 = 0; TypeParam add2 = 0;
    
    for (int i=0; i<size; ++i)
    {
        field1->operator[](i) = TypeParam(i+100);
        field2->operator[](i) = TypeParam(-i*0.1);
        add1 += (i+100);
        dot1 += (i+100)*(i+100);
        add2 += -i*0.1;
        dot2 += (-i*0.1)*(-i*0.1);
    };

    ASSERT_FLOAT_EQ( 100, field1->min() );
    ASSERT_FLOAT_EQ( 131, field1->max() );
    ASSERT_FLOAT_EQ( add1, field1->sum() );
    ASSERT_FLOAT_EQ( dot1, field1->dot(*field1) );

    ASSERT_FLOAT_EQ( -3.1, field2->min() );
    ASSERT_FLOAT_EQ( 0.0, field2->max() );
    ASSERT_FLOAT_EQ( add2, field2->sum() );
    ASSERT_FLOAT_EQ( dot2, field2->dot(*field2) );
}

// ============================================
//          Testing Other Functions
// ============================================
TYPED_TEST(test_dfield_cuda, auxialiary_functions)
{
    int size = 32;
    using fieldcuda_pointer = typename dfield_cuda<TypeParam>::pointer;

    // cast function
    fieldcuda_pointer field1 = dfield_cuda<TypeParam>::new_pointer(size,-1.1);

    dfield_cuda<int>::pointer field2 = field1->template cast<int>();
    dfield_cuda<float>::pointer field3 = field1->template cast<float>();
    dfield_cuda<double>::pointer field4 = field1->template cast<double>();
    dfield_cuda<unsigned int>::pointer field5 = field1->template cast<unsigned int>();

    int vari = -1.0;
    float varf = -1.1;
    double vard = -1.1;
    unsigned int varui = static_cast<unsigned int>(varf);

    ASSERT_TRUE(typeid(int) == typeid(field2->operator[](10)));
    ASSERT_TRUE(typeid(float) == typeid(field3->operator[](10)));
    ASSERT_TRUE(typeid(double) == typeid(field4->operator[](10)));
    ASSERT_TRUE(typeid(unsigned int) == typeid(field5->operator[](10)));
    ASSERT_EQ( vari, field2->operator[](10) );
    ASSERT_EQ( varui, field5->operator[](10) );
    ASSERT_FLOAT_EQ( varf, field3->operator[](10) );
    ASSERT_FLOAT_EQ( vard, field4->operator[](10) );

    // normalize function
    fieldcuda_pointer field11 = dfield_cuda<TypeParam>::new_pointer(size,0.0);
    fieldcuda_pointer field12 = dfield_cuda<TypeParam>::new_pointer(size,0.0);
    fieldcuda_pointer field13;
    fieldcuda_pointer field14;

    for (int i=0; i<size; ++i)
    {
        field11->operator[](i) = TypeParam(-i-131.2);
        field12->operator[](i) = TypeParam(i+10.5);
    };

    field13 = field11->normalize();
    field14 = field12->normalize();

    ASSERT_FLOAT_EQ( -131.2, field11->operator[](0) );
    ASSERT_FLOAT_EQ( -131.2 - size + 1.0, field11->operator[](size-1) );
    ASSERT_FLOAT_EQ( 10.5, field12->operator[](0) );
    ASSERT_FLOAT_EQ( 10.5 + size - 1.0, field12->operator[](size-1) );
    ASSERT_FLOAT_EQ( 1.0, field13->operator[](0) );
    ASSERT_FLOAT_EQ( 1/double(size-1), field13->operator[](size-2) );
    ASSERT_FLOAT_EQ( 0.0, field13->operator[](size-1) );
    ASSERT_FLOAT_EQ( 0.0, field14->operator[](0) );
    ASSERT_FLOAT_EQ( 1/double(size-1), field14->operator[](1) );
    ASSERT_FLOAT_EQ( 1.0, field14->operator[](size-1) );
}
*/


int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}