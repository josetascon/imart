/*
* @Author: jose
* @Date:   2019-11-20 00:00:00
* @Last Modified by:   jose
* @Last Modified time: 2019-11-20 00:00:00
*/

#ifndef __PROCESS_OBJECT_H__
#define __PROCESS_OBJECT_H__

// std libs
#include <iostream>     // std::cout
#include <cassert>      // assert

// local libs
#include "inherit.h"
#include "object.h"

namespace imart
{

// Class object
class process_object: public inherit<process_object, object>
{
public:
    //Type definitions
    using self    = process_object;
    using pointer = std::shared_ptr<self>;
    using vector  = std::vector<self::pointer>;

protected:
    // ===========================================
    // Internal Variables
    // ===========================================
    int num_inputs;                         // number of inputs
    int num_outputs;                        // number of outputs
    bool ready;                             // ready to run
    std::vector<object::pointer> vinput;    // vector with pointers of inputs
    std::vector<object::pointer> voutput;   // vector with pointers of outputs

    // ===========================================
    // Functions
    // ===========================================
    virtual void init(int inputs, int outputs);     // init default properties
    virtual std::string info(std::string msg);
    virtual void execute_();

    // ===========================================
    // Set Functions
    // ===========================================
    void set_total_inputs(int num);
    void set_total_outputs(int num);
    template<typename ... Args> void setup_input(Args... args);
    template<typename ... Args> void setup_output(Args... args);

public:
    // ===========================================
    // Constructor Functions
    // ===========================================
    process_object();                                 // constructor empty
    process_object(int inputs, int outputs);          // constructor with number of inputs and outputs
    process_object(const process_object & input);     // constructor with same type    
    ~process_object();                                // destructor empty

    // ===========================================
    // Create Functions
    // ===========================================
    virtual void clone_(const process_object & input);// copy everything
    virtual void copy_ (const process_object & input);// share data
    virtual void mimic_(const process_object & input);// copy meta data

    // ===========================================
    // Get Functions
    // ===========================================
    int get_total_inputs() const;
    int get_total_outputs() const;
    std::vector<object::pointer> get_input() const;
    std::vector<object::pointer> get_output() const;

    // ===========================================
    // Set Functions
    // ===========================================
    virtual void set_input();                       // to be defined 
    virtual void set_output();                      // to be defined

    // ===========================================
    // Functions
    // ===========================================
    void execute();
};


// ===========================================
//          Functions of Class object
// ===========================================

// ===========================================
// Create Functions
// ===========================================
//! Constructor empty
process_object::process_object()
{
    this->class_name = "process_object";
    init(1,1);
};

//! Constructor with number of dimensions
process_object::process_object(int inputs, int outputs)
{
    this->class_name = "process_object";
    init(inputs, outputs);
};

//! Constructor to clone
process_object::process_object(const process_object & input)
{
    clone_(input);                // call the virtual
};

// Destructor
process_object::~process_object()
{
    ;
};

// ===========================================
// Create Functions
// ===========================================
// Initialization function
void process_object::init(int inputs, int outputs)
{
    // Attributes initialization
    set_total_inputs(inputs);
    set_total_inputs(outputs);
};

// template <typename type>
// void process_object::copy_properties(const process_object & input)
// {

// };

void process_object::clone_(const process_object & input)
{
    copy_(input);
};

void process_object::copy_(const process_object & input)
{
    mimic_(input);
    vinput = input.get_input();
    voutput = input.get_output();
};

void process_object::mimic_(const process_object & input)
{
    num_inputs = input.get_total_inputs();
    num_outputs = input.get_total_outputs();
    ready = false;
};

// ===========================================
// Get Functions
// ===========================================
int process_object::get_total_inputs() const
{
    return num_inputs;
};

int process_object::get_total_outputs() const
{
    return num_outputs;
};

std::vector<object::pointer> process_object::get_input() const
{
    return vinput;
};

std::vector<object::pointer> process_object::get_output() const
{
    return voutput;
};

// ===========================================
// Set Functions
// ===========================================
void process_object::set_total_inputs(int num)
{
    num_inputs = num;
};

void process_object::set_total_outputs(int num)
{
    num_outputs = num;
};

template<typename ... Args>
void process_object::setup_input(Args... args)
{
    std::vector<object::pointer> in({args ...});
    assert(num_inputs == in.size());
    vinput = in;
};

template<typename ... Args>
void process_object::setup_output(Args... args)
{
    std::vector<object::pointer> out({args ...});
    assert(num_outputs == out.size());
    voutput = out;
};

// ===========================================
// Print Functions
// ===========================================
std::string process_object::info(std::string msg)
{
    std::stringstream ss;
    std::string title = "Process Object Information";
    if (msg != "") { title = msg; };

    // Summary of the object information
    // ss << object::info(title);
    ss << "Inputs #: \t\t\t" << num_inputs << std::endl;
    for(int i = 0; i < vinput.size(); i++) 
    {
        ss << "\tInput " << i+1 << ":\t\t";
        ss << vinput[i]->get_name() << "\n";
    };
    ss << "Outputs #: \t\t\t" << num_outputs << std::endl;
    for(int i = 0; i < voutput.size(); i++) 
    {
        ss << "\tOutput " << i+1 << ":\t\t";
        ss << voutput[i]->get_name() << "\n";
    };
    return ss.str();
};

// ===========================================
// Functions
// ===========================================
void process_object::execute()
{
    execute_();
}

}; //end namespace

#endif