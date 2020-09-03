/*
* @Author: jose
* @Date:   2020-08-24 00:00:00
* @Last Modified by:   jose
* @Last Modified time: 2020-08-24 00:00:00
*/

#ifndef __CUDA_OBJECT_H__
#define __CUDA_OBJECT_H__

// std libs
#include <iostream>     // std::cout
#include <vector>       // std::vector
#include <cassert>      // assert
#include <cmath>        // math

// local libs
#include "object.h"

namespace imart
{

class cuda_object: public inherit<cuda_object, object>
{
public:
    // Type definitions
    using self    = cuda_object;
    using pointer = std::shared_ptr<self>;
    using vector  = std::vector<self::pointer>;

protected:
    // ===========================================
    // Internal Variables
    // ===========================================
    int _error_;
    int threads;
    std::vector<int> grid;  // grid with blocks
    std::vector<int> block; // threads per block

    int blocksy;
    float inv_blocksy;

    bool status_init;
    bool status_three_dim;

    // ===========================================
    // Functions
    // ===========================================
    void init();
    void check_error(int err);

public:
    // ===========================================
    // Create Functions
    // ===========================================
    // Constructors
    cuda_object();
    // Destructor
    // ~cuda_object();

    // ===========================================
    // Get Functions
    // ===========================================
    int get_threads() const;
    std::vector<int> get_grid() const;
    std::vector<int> get_block() const;
    int get_blocksy() const;
    float get_blocksy_inverse() const;

    // ===========================================
    // Set Functions
    // ===========================================
    void set_threads(int num_threads);
    void set_grid(std::vector<int> grid_size);
    void set_block(std::vector<int> block_size);

    // ===========================================
    // Print Functions
    // ===========================================
    virtual std::string info(std::string msg);

    // ===========================================
    // Functions
    // ===========================================
    // compute
    void setup(int max_size);
    void setup(std::vector<int> & dims);

    template<typename Func, typename ...Args>
    void execute( Func const & function , Args const &... args );
};

// ===========================================
//      Functions of Class metric
// ===========================================

// ===========================================
// Constructor Functions
// ===========================================
// Constructor
cuda_object::cuda_object()
{   
    this->class_name = "cuda_object";
    init();
};

void cuda_object::init()
{
    // Set internal variables
    status_init = false;
    status_three_dim = false;
    _error_ = 0;

    // Init block and grid
    threads = 256;

    blocksy = 0;
    inv_blocksy = 0;

    // Initialization succeded
    status_init = true;

    // cuda_check_gpu();
    // cuda_print();
};

// Constructor
// cuda_object::~cuda_object()
// {   
//     // Clean and free memory
//     cl::flush();
//     _error_ = cl::finish();
//     imart_assert_cl(_error_, "Error during finish");
// };

// ===========================================
// Get Functions
// ===========================================
int cuda_object::get_threads() const
{
    return threads;
};

std::vector<int> cuda_object::get_grid() const
{
    return grid;
};

std::vector<int> cuda_object::get_block() const
{
    return block;
};

int cuda_object::get_blocksy() const
{
    return blocksy;
};

float cuda_object::get_blocksy_inverse() const
{
    return inv_blocksy;
};

// ===========================================
// Set Functions
// ===========================================
void cuda_object::set_threads(int num_threads)
{
    threads = num_threads;
};

void cuda_object::set_grid(std::vector<int> grid_size)
{
    grid = grid_size;
};

void cuda_object::set_block(std::vector<int> block_size)
{
    block = block_size;
};

void cuda_object::setup(int max_size)
{
    if (max_size < 1) return; // finish if size is zero
    
    // define grid and block size
    block = std::vector<int>{threads, 0, 0};
    grid = std::vector<int>{1 + (max_size - 1)/threads, 0, 0};
    // grid = std::vector<int>{(max_size + threads - 1)/threads, 0, 0};

};

// ===========================================
// Print Functions
// ===========================================
std::string cuda_object::info(std::string msg)
{
    std::stringstream ss;
    std::string title = "CUDA Object Information";
    if (msg != "") { title = msg; };
    // Summary of the optimizer information
    ss << object::info(title);
    // if(status_init)
    // if(status_program)
    // if(status_kernel)
    return ss.str();
};

// ===========================================
// Functions
// ===========================================
void cuda_object::setup(std::vector<int> & dims)
{
    assert(dims.size() > 0 && dims.size() < 4);
    // std::cout << "Run Kernel" << std::endl;
    if (dims.size() == 1)
    {
        if (dims[0] < 1) return; // finish if size is zero
        // define grid and block size
        threads = 256;
        block = std::vector<int>{threads, 0, 0};
        grid = std::vector<int>{1 + (dims[0] - 1)/threads, 0, 0};
    }
    else if (dims.size() == 2)
    {
        if (dims[0] < 1 or dims[1] < 1) return; // finish if size is zero
        // define grid and block size
        threads = 32;
        block = std::vector<int>{threads, threads, 0};
        grid = std::vector<int>{1 + (dims[0] - 1)/threads, 1 + (dims[1] - 1)/threads, 0};

    }
    else if (dims.size() == 3)
    {
        if (dims[0] < 1 or dims[1] < 1 or dims[2] < 1) return; // finish if size is zero
        // define grid and block size
        threads = 8;
        block = std::vector<int>{threads, threads, threads};
        grid = std::vector<int>{1 + (dims[0] - 1)/threads, (1 + (dims[1] - 1)/threads), (1 + (dims[2] - 1)/threads)};

        blocksy = 1 + (dims[1] - 1)/threads;
        inv_blocksy = 1/float(blocksy);
        status_three_dim = true;
    }
    else ;    
};

template<typename Func, typename ...Args>
void cuda_object::execute( Func const & function , Args const &... args )
{
    // std::cout << "Run kernel with arguments: " << std::endl;
    if constexpr (sizeof...(args) > 0) 
    {
        // option 1
        // cuda_launcher launch(grid, block);
        // launch.execute(function, args ...);
        
        // option 2
        // function(args...);

        // option 3
        function(grid, block, args...);

        // if (status_three_dim)
        //     function(grid, block, args...);
        // else ; 
        //     function(grid, block, args..., blocksy, inv_blocksy);
        // status_three_dim = false;
    };
};

}; //end namespace

#endif