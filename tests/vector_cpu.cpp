/*
* @Author: Jose Tascon
* @Date:   2020-06-06 00:00:00
* @Last Modified by:   Jose Tascon
* @Last Modified time: 2020-09-07 16:55:57
*/

// std libs
#include <iostream>
#include <memory>
#include <typeinfo>

// gtest header
#include <gtest/gtest.h>

// local header
#include "../src/vector_cpu.h"

using namespace imart;

template <typename T>
class test_vector_cpu : public ::testing::Test{};

using test_types = ::testing::Types<float, double>;

// ============================================
//          Testing Constructors
// ============================================
TYPED_TEST_CASE(test_vector_cpu, test_types);

TYPED_TEST(test_vector_cpu, constructor)
{
    vector_cpu<TypeParam> vec1;
    ASSERT_TRUE("vector_cpu" == vec1.get_name());
    ASSERT_EQ( 0, vec1.size() );

    vector_cpu<TypeParam> vec2(0);
    ASSERT_EQ( 0, vec2.size() );

    vector_cpu<TypeParam> vec3(121);
    ASSERT_TRUE("vector_cpu" == vec3.get_name());
    ASSERT_EQ( 121, vec3.size() );

    vector_cpu<TypeParam> vec4(1,3.0);
    ASSERT_FLOAT_EQ( 3.0, vec4[0] );

    vector_cpu<TypeParam> vec5(10,-5.2);
    ASSERT_TRUE("vector_cpu" == vec5.get_name());
    ASSERT_FLOAT_EQ( -5.2, vec5[0] );
    ASSERT_FLOAT_EQ( -5.2, vec5[4] );
    ASSERT_FLOAT_EQ( -5.2, vec5[9] );

    vector_cpu<TypeParam> vec6(vec5);
    ASSERT_TRUE("vector_cpu" == vec6.get_name());
    ASSERT_FLOAT_EQ( -5.2, vec5[0] );
    ASSERT_FLOAT_EQ( -5.2, vec5[4] );
    ASSERT_FLOAT_EQ( -5.2, vec5[9] );

    vector_cpu<TypeParam> vec7{0.2, -0.2, 5.1, 8.234, -1.0};
    ASSERT_TRUE("vector_cpu" == vec7.get_name());
    ASSERT_EQ( 5, vec7.size() );
    ASSERT_FLOAT_EQ( 0.2, vec7[0] );
    ASSERT_FLOAT_EQ( -0.2, vec7[1] );
    ASSERT_FLOAT_EQ( 5.1, vec7[2] );
    ASSERT_FLOAT_EQ( 8.234, vec7[3] );
    ASSERT_FLOAT_EQ( -1.0, vec7[4] );
}

// ============================================
//          Testing Pointers
// ============================================
TYPED_TEST(test_vector_cpu, pointers)
{
    using vcpu_pointer = typename vector_cpu<TypeParam>::pointer;
    vcpu_pointer vec1 = vector_cpu<TypeParam>::new_pointer();
    ASSERT_TRUE("vector_cpu" == vec1->get_name());
    ASSERT_EQ( 0, vec1->size() );

    vcpu_pointer vec2 = vector_cpu<TypeParam>::new_pointer(0);
    ASSERT_EQ( 0, vec2->size() );

    vcpu_pointer vec3 = vector_cpu<TypeParam>::new_pointer(201);
    ASSERT_TRUE("vector_cpu" == vec3->get_name());
    ASSERT_EQ( 201, vec3->size() );

    vcpu_pointer vec4 = vector_cpu<TypeParam>::new_pointer(1,-7.0);
    ASSERT_FLOAT_EQ( -7.0, vec4->operator[](0) );

    vcpu_pointer vec5 = vector_cpu<TypeParam>::new_pointer(21,11.2231);
    ASSERT_TRUE("vector_cpu" == vec5->get_name());
    ASSERT_FLOAT_EQ( 11.2231, vec5->operator[](0) );
    ASSERT_FLOAT_EQ( 11.2231, vec5->operator[](15) );
    ASSERT_FLOAT_EQ( 11.2231, vec5->operator[](20) );
}

// ============================================
//          Testing Clone Methods
// ============================================
TYPED_TEST(test_vector_cpu, clones)
{
    using vcpu_pointer = typename vector_cpu<TypeParam>::pointer;
    vcpu_pointer vec1 = vector_cpu<TypeParam>::new_pointer(8,3.14159);
    vcpu_pointer vec2 = vec1->clone();
    vcpu_pointer vec3 = vec1->copy();
    vcpu_pointer vec4 = vec1->mimic();

    ASSERT_TRUE("vector_cpu" == vec2->get_name());
    ASSERT_TRUE("vector_cpu" == vec3->get_name());
    ASSERT_TRUE("vector_cpu" == vec4->get_name());

    ASSERT_EQ( 8, vec1->size() );
    ASSERT_EQ( 8, vec2->size() );
    ASSERT_EQ( 8, vec3->size() );
    ASSERT_EQ( 8, vec4->size() );

    ASSERT_FALSE(vec1.get() == vec2.get());
    ASSERT_FALSE(vec1.get() == vec3.get());
    ASSERT_FALSE(vec1.get() == vec4.get());

    ASSERT_FLOAT_EQ( 3.14159, vec1->operator[](0) );
    ASSERT_FLOAT_EQ( 3.14159, vec1->operator[](7) );
    ASSERT_FLOAT_EQ( 3.14159, vec2->operator[](0) );
    ASSERT_FLOAT_EQ( 3.14159, vec2->operator[](7) );
    EXPECT_FLOAT_EQ( 0.0, vec3->operator[](0) );
    EXPECT_FLOAT_EQ( 0.0, vec3->operator[](7) );
    EXPECT_FLOAT_EQ( 0.0, vec4->operator[](0) );
    EXPECT_FLOAT_EQ( 0.0, vec4->operator[](7) );
}

// ============================================
//          Testing Initialization Methods
// ============================================
TYPED_TEST(test_vector_cpu, initialization)
{
    int size = 50;
    using vcpu_pointer = typename vector_cpu<TypeParam>::pointer;
    vcpu_pointer vec1 = vector_cpu<TypeParam>::new_pointer(size);
    vcpu_pointer vec2 = vec1->mimic();
    vcpu_pointer vec3 = vec1->mimic();

    vec1->ones();
    vec2->assign(TypeParam(3.0));
    vec3->random();

    ASSERT_FLOAT_EQ( 1.0, vec1->operator[](0) );
    ASSERT_FLOAT_EQ( 1.0, vec1->operator[](int(size/2)) );
    ASSERT_FLOAT_EQ( 1.0, vec1->operator[](size-4) );
    ASSERT_FLOAT_EQ( 3.0, vec2->operator[](0) );
    ASSERT_FLOAT_EQ( 3.0, vec2->operator[](int(size/2)) );
    ASSERT_FLOAT_EQ( 3.0, vec2->operator[](size-4) );
    ASSERT_FALSE( vec3->operator[](0) == vec3->operator[](1) );
    ASSERT_FALSE( vec3->operator[](0) == vec3->operator[](int(size/2)) );
    ASSERT_FALSE( vec3->operator[](int(size/2)) == vec3->operator[](size-4) );

    vec1->zeros();
    vec2->assign(TypeParam(-2.431));

    ASSERT_FLOAT_EQ( 0.0, vec1->operator[](0) );
    ASSERT_FLOAT_EQ( 0.0, vec1->operator[](int(size/2)) );
    ASSERT_FLOAT_EQ( 0.0, vec1->operator[](size-4) );
    ASSERT_FLOAT_EQ( -2.431, vec2->operator[](0) );
    ASSERT_FLOAT_EQ( -2.431, vec2->operator[](int(size/2)) );
    ASSERT_FLOAT_EQ( -2.431, vec2->operator[](size-4) );
    ASSERT_FALSE( vec3->operator[](0) == vec3->operator[](1) );
    ASSERT_FALSE( vec3->operator[](0) == vec3->operator[](int(size/2)) );
    ASSERT_FALSE( vec3->operator[](int(size/2)) == vec3->operator[](size-4) );
}

// ============================================
//          Testing Vector Operations
// ============================================
TYPED_TEST(test_vector_cpu, vector_operations)
{
    using vcpu_pointer = typename vector_cpu<TypeParam>::pointer;
    vcpu_pointer vec1 = vector_cpu<TypeParam>::new_pointer(30,2.1);
    vcpu_pointer vec2 = vector_cpu<TypeParam>::new_pointer(30,1.1);
    vcpu_pointer vec3;
    vcpu_pointer vec4 = vector_cpu<TypeParam>::new_pointer(30,4.2);
    vcpu_pointer vec5;
    vcpu_pointer vec6 = vector_cpu<TypeParam>::new_pointer(30,14.2123);
    vcpu_pointer vec7;
    vcpu_pointer vec8 = vector_cpu<TypeParam>::new_pointer(30,-2*14.2123);
    vcpu_pointer vec9;
    vcpu_pointer vec10 = vector_cpu<TypeParam>::new_pointer(30,4.0);
    vcpu_pointer vec11;

    vec3 = *vec1 + *vec2;
    ASSERT_TRUE("vector_cpu" == vec3->get_name());
    ASSERT_EQ( 30, vec3->size() );
    ASSERT_FLOAT_EQ( 3.2, vec3->operator[](0) );
    ASSERT_FLOAT_EQ( 3.2, vec3->operator[](7) );
    ASSERT_FLOAT_EQ( 3.2, vec3->operator[](18) );
    ASSERT_FLOAT_EQ( 3.2, vec3->operator[](29) );

    vec5 = *vec3 - *vec4;
    ASSERT_TRUE("vector_cpu" == vec5->get_name());
    ASSERT_EQ( 30, vec5->size() );
    ASSERT_FLOAT_EQ( -1.0, vec5->operator[](0) );
    ASSERT_FLOAT_EQ( -1.0, vec5->operator[](9) );
    ASSERT_FLOAT_EQ( -1.0, vec5->operator[](16) );
    ASSERT_FLOAT_EQ( -1.0, vec5->operator[](29) );

    vec7 = (*vec5) * (*vec6);
    ASSERT_TRUE("vector_cpu" == vec7->get_name());
    ASSERT_EQ( 30, vec7->size() );
    ASSERT_FLOAT_EQ( -14.2123, vec7->operator[](0) );
    ASSERT_FLOAT_EQ( -14.2123, vec7->operator[](3) );
    ASSERT_FLOAT_EQ( -14.2123, vec7->operator[](12) );
    ASSERT_FLOAT_EQ( -14.2123, vec7->operator[](29) );

    vec9 = (*vec7) / (*vec8);
    ASSERT_TRUE("vector_cpu" == vec9->get_name());
    ASSERT_EQ( 30, vec9->size() );
    ASSERT_FLOAT_EQ( 0.5, vec9->operator[](0) );
    ASSERT_FLOAT_EQ( 0.5, vec9->operator[](4) );
    ASSERT_FLOAT_EQ( 0.5, vec9->operator[](21) );
    ASSERT_FLOAT_EQ( 0.5, vec9->operator[](29) );

    vec11 = (*vec9) ^ (*vec10);
    ASSERT_TRUE("vector_cpu" == vec11->get_name());
    ASSERT_EQ( 30, vec11->size() );
    ASSERT_FLOAT_EQ( 0.0625, vec11->operator[](0) );
    ASSERT_FLOAT_EQ( 0.0625, vec11->operator[](5) );
    ASSERT_FLOAT_EQ( 0.0625, vec11->operator[](23) );
    ASSERT_FLOAT_EQ( 0.0625, vec11->operator[](29) );

    // Check vector values again after operations
    ASSERT_FLOAT_EQ( 2.1, vec1->operator[](0) );
    ASSERT_FLOAT_EQ( 2.1, vec1->operator[](29) );
    ASSERT_FLOAT_EQ( 1.1, vec2->operator[](0) );
    ASSERT_FLOAT_EQ( 1.1, vec2->operator[](29) );
    ASSERT_FLOAT_EQ( 3.2, vec3->operator[](0) );
    ASSERT_FLOAT_EQ( 3.2, vec3->operator[](29) );
    ASSERT_FLOAT_EQ( 4.2, vec4->operator[](0) );
    ASSERT_FLOAT_EQ( 4.2, vec4->operator[](29) );
    ASSERT_FLOAT_EQ( -1.0, vec5->operator[](0) );
    ASSERT_FLOAT_EQ( -1.0, vec5->operator[](29) );
    ASSERT_FLOAT_EQ( 14.2123, vec6->operator[](0) );
    ASSERT_FLOAT_EQ( 14.2123, vec6->operator[](29) );
    ASSERT_FLOAT_EQ( -14.2123, vec7->operator[](0) );
    ASSERT_FLOAT_EQ( -14.2123, vec7->operator[](29) );
    ASSERT_FLOAT_EQ( -2*14.2123, vec8->operator[](0) );
    ASSERT_FLOAT_EQ( -2*14.2123, vec8->operator[](29) );
    ASSERT_FLOAT_EQ( 0.5, vec9->operator[](0) );
    ASSERT_FLOAT_EQ( 0.5, vec9->operator[](29) );
    ASSERT_FLOAT_EQ( 4.0, vec10->operator[](0) );
    ASSERT_FLOAT_EQ( 4.0, vec10->operator[](29) );
}

// ============================================
//          Testing Scalar Operations
// ============================================
TYPED_TEST(test_vector_cpu, scalar_operations)
{
    using vcpu_pointer = typename vector_cpu<TypeParam>::pointer;
    vcpu_pointer vec1 = vector_cpu<TypeParam>::new_pointer(14,1.5);
    vcpu_pointer vec2 = vector_cpu<TypeParam>::new_pointer(16,-4.632);
    vcpu_pointer vec3 = vector_cpu<TypeParam>::new_pointer(18,-12.0);
    vcpu_pointer vec4 = vector_cpu<TypeParam>::new_pointer(20,8.4);
    vcpu_pointer vec5 = vector_cpu<TypeParam>::new_pointer(22,2.0);
    vcpu_pointer vec6;
    vcpu_pointer vec7;
    vcpu_pointer vec8;
    vcpu_pointer vec9;
    vcpu_pointer vec10;

    // Right hand side
    vec6 = *vec1 + 4.5;
    ASSERT_TRUE("vector_cpu" == vec6->get_name());
    ASSERT_EQ( 14, vec6->size() );
    ASSERT_FLOAT_EQ( 6.0, vec6->operator[](0) );
    ASSERT_FLOAT_EQ( 6.0, vec6->operator[](7) );
    ASSERT_FLOAT_EQ( 6.0, vec6->operator[](13) );

    vec7 = *vec2 - 1.111;
    ASSERT_TRUE("vector_cpu" == vec7->get_name());
    ASSERT_EQ( 16, vec7->size() );
    ASSERT_FLOAT_EQ( -5.743, vec7->operator[](0) );
    ASSERT_FLOAT_EQ( -5.743, vec7->operator[](6) );
    ASSERT_FLOAT_EQ( -5.743, vec7->operator[](15) );

    vec8 = (*vec3) * (-0.9);
    ASSERT_TRUE("vector_cpu" == vec8->get_name());
    ASSERT_EQ( 18, vec8->size() );
    ASSERT_FLOAT_EQ( 10.8, vec8->operator[](0) );
    ASSERT_FLOAT_EQ( 10.8, vec8->operator[](13) );
    ASSERT_FLOAT_EQ( 10.8, vec8->operator[](17) );

    vec9 = (*vec4) / (-4.0);
    ASSERT_TRUE("vector_cpu" == vec9->get_name());
    ASSERT_EQ( 20, vec9->size() );
    ASSERT_FLOAT_EQ( -2.1, vec9->operator[](0) );
    ASSERT_FLOAT_EQ( -2.1, vec9->operator[](11) );
    ASSERT_FLOAT_EQ( -2.1, vec9->operator[](19) );

    vec10 = (*vec5) ^ (3.0);
    ASSERT_TRUE("vector_cpu" == vec10->get_name());
    ASSERT_EQ( 22, vec10->size() );
    ASSERT_FLOAT_EQ( 8.0, vec10->operator[](0) );
    ASSERT_FLOAT_EQ( 8.0, vec10->operator[](8) );
    ASSERT_FLOAT_EQ( 8.0, vec10->operator[](21) );

    // Left hand size
    vec6 = TypeParam(4.5) + (*vec1);
    ASSERT_TRUE("vector_cpu" == vec6->get_name());
    ASSERT_EQ( 14, vec6->size() );
    ASSERT_FLOAT_EQ( 6.0, vec6->operator[](0) );
    ASSERT_FLOAT_EQ( 6.0, vec6->operator[](7) );
    ASSERT_FLOAT_EQ( 6.0, vec6->operator[](13) );

    vec7 = TypeParam(1.111) - (*vec2);
    ASSERT_TRUE("vector_cpu" == vec7->get_name());
    ASSERT_EQ( 16, vec7->size() );
    ASSERT_FLOAT_EQ( 5.743, vec7->operator[](0) );
    ASSERT_FLOAT_EQ( 5.743, vec7->operator[](6) );
    ASSERT_FLOAT_EQ( 5.743, vec7->operator[](15) );

    vec8 = TypeParam(0.9)*(*vec3);
    ASSERT_TRUE("vector_cpu" == vec8->get_name());
    ASSERT_EQ( 18, vec8->size() );
    ASSERT_FLOAT_EQ( -10.8, vec8->operator[](0) );
    ASSERT_FLOAT_EQ( -10.8, vec8->operator[](13) );
    ASSERT_FLOAT_EQ( -10.8, vec8->operator[](17) );

    vec9 = TypeParam(-4.0)/(*vec4);
    ASSERT_TRUE("vector_cpu" == vec9->get_name());
    ASSERT_EQ( 20, vec9->size() );
    ASSERT_FLOAT_EQ( -0.47619051, vec9->operator[](0) );
    ASSERT_FLOAT_EQ( -0.47619051, vec9->operator[](11) );
    ASSERT_FLOAT_EQ( -0.47619051, vec9->operator[](19) );

    // Check original values after operations
    ASSERT_EQ( 14, vec1->size() );
    ASSERT_FLOAT_EQ( 1.5, vec1->operator[](0) );
    ASSERT_FLOAT_EQ( 1.5, vec1->operator[](7) );
    ASSERT_FLOAT_EQ( 1.5, vec1->operator[](13) );

    ASSERT_EQ( 16, vec2->size() );
    ASSERT_FLOAT_EQ( -4.632, vec2->operator[](0) );
    ASSERT_FLOAT_EQ( -4.632, vec2->operator[](6) );
    ASSERT_FLOAT_EQ( -4.632, vec2->operator[](15) );

    ASSERT_EQ( 18, vec3->size() );
    ASSERT_FLOAT_EQ( -12.0, vec3->operator[](0) );
    ASSERT_FLOAT_EQ( -12.0, vec3->operator[](13) );
    ASSERT_FLOAT_EQ( -12.0, vec3->operator[](17) );

    ASSERT_EQ( 20, vec9->size() );
    ASSERT_FLOAT_EQ( 8.4, vec4->operator[](0) );
    ASSERT_FLOAT_EQ( 8.4, vec4->operator[](11) );
    ASSERT_FLOAT_EQ( 8.4, vec4->operator[](19) );

    ASSERT_EQ( 22, vec5->size() );
    ASSERT_FLOAT_EQ( 2.0, vec5->operator[](0) );
    ASSERT_FLOAT_EQ( 2.0, vec5->operator[](8) );
    ASSERT_FLOAT_EQ( 2.0, vec5->operator[](21) );
}

// ============================================
//          Testing Reduction Operations
// ============================================
TYPED_TEST(test_vector_cpu, reduction_operations)
{
    int size = 32;
    using vcpu_pointer = typename vector_cpu<TypeParam>::pointer;
    vcpu_pointer vec1 = vector_cpu<TypeParam>::new_pointer(size,0.0);
    vcpu_pointer vec2 = vector_cpu<TypeParam>::new_pointer(size,0.0);

    TypeParam dot1 = 0; TypeParam dot2 = 0;
    TypeParam add1 = 0; TypeParam add2 = 0;
    
    for (int i=0; i<size; ++i)
    {
        vec1->operator[](i) = TypeParam(i+100);
        vec2->operator[](i) = TypeParam(-i*0.1);
        add1 += (i+100);
        dot1 += (i+100)*(i+100);
        add2 += -i*0.1;
        dot2 += (-i*0.1)*(-i*0.1);
    };

    ASSERT_FLOAT_EQ( 100, vec1->min() );
    ASSERT_FLOAT_EQ( 131, vec1->max() );
    ASSERT_FLOAT_EQ( add1, vec1->sum() );
    ASSERT_FLOAT_EQ( dot1, vec1->dot(*vec1) );

    ASSERT_FLOAT_EQ( -3.1, vec2->min() );
    ASSERT_FLOAT_EQ( 0.0, vec2->max() );
    ASSERT_FLOAT_EQ( add2, vec2->sum() );
    ASSERT_FLOAT_EQ( dot2, vec2->dot(*vec2) );
}

// ============================================
//          Testing Other Functions
// ============================================
TYPED_TEST(test_vector_cpu, auxialiary_functions)
{
    int size = 32;
    using vcpu_pointer = typename vector_cpu<TypeParam>::pointer;

    // cast function
    vcpu_pointer vec1 = vector_cpu<TypeParam>::new_pointer(size,-1.1);

    vector_cpu<int>::pointer vec2 = vec1->template cast<int>();
    vector_cpu<float>::pointer vec3 = vec1->template cast<float>();
    vector_cpu<double>::pointer vec4 = vec1->template cast<double>();
    vector_cpu<unsigned int>::pointer vec5 = vec1->template cast<unsigned int>();

    int vari = -1.0;
    float varf = -1.1;
    double vard = -1.1;
    // unsigned int varui = static_cast<unsigned int>(varf);

    ASSERT_TRUE(typeid(int) == typeid(vec2->operator[](10)));
    ASSERT_TRUE(typeid(float) == typeid(vec3->operator[](10)));
    ASSERT_TRUE(typeid(double) == typeid(vec4->operator[](10)));
    ASSERT_TRUE(typeid(unsigned int) == typeid(vec5->operator[](10)));
    ASSERT_EQ( vari, vec2->operator[](10) );
    // ASSERT_EQ( varui, vec5->operator[](10) );
    ASSERT_FLOAT_EQ( varf, vec3->operator[](10) );
    ASSERT_FLOAT_EQ( vard, vec4->operator[](10) );

    // normalize function
    vcpu_pointer vec11 = vector_cpu<TypeParam>::new_pointer(size,0.0);
    vcpu_pointer vec12 = vector_cpu<TypeParam>::new_pointer(size,0.0);
    vcpu_pointer vec13;
    vcpu_pointer vec14;

    for (int i=0; i<size; ++i)
    {
        vec11->operator[](i) = TypeParam(-i-131.2);
        vec12->operator[](i) = TypeParam(i+10.5);
    };

    vec13 = vec11->normalize();
    vec14 = vec12->normalize();

    ASSERT_FLOAT_EQ( -131.2, vec11->operator[](0) );
    ASSERT_FLOAT_EQ( -131.2 - size + 1.0, vec11->operator[](size-1) );
    ASSERT_FLOAT_EQ( 10.5, vec12->operator[](0) );
    ASSERT_FLOAT_EQ( 10.5 + size - 1.0, vec12->operator[](size-1) );
    ASSERT_FLOAT_EQ( 1.0, vec13->operator[](0) );
    ASSERT_FLOAT_EQ( 1/double(size-1), vec13->operator[](size-2) );
    ASSERT_FLOAT_EQ( 0.0, vec13->operator[](size-1) );
    ASSERT_FLOAT_EQ( 0.0, vec14->operator[](0) );
    ASSERT_FLOAT_EQ( 1/double(size-1), vec14->operator[](1) );
    ASSERT_FLOAT_EQ( 1.0, vec14->operator[](size-1) );
}

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}