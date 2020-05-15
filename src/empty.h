/*
* @Author: jose
* @Date:   2020-01-27 00:00:00
* @Last Modified by:   jose
* @Last Modified time: 2020-01-27 00:00:00
*/

#ifndef __CLASS_H__
#define __CLASS_H__

#include "image.h"

namespace imart
{

template <typename pixel_type>
class myclass: public object<pixel_type>
{
public:
    //Type definitions
    ;
protected:
    // ===========================================
    // Internal Variables
    // ===========================================

    // ===========================================
    // Functions
    // ===========================================
    void init();

public:
    // ===========================================
    // Create Functions
    // ===========================================
    // Constructors
    

    // ===========================================
    // Get Functions
    // ===========================================
    int get_abc() const;

    // ===========================================
    // Print Functions
    // ===========================================
    virtual std::string info(std::string msg);
    virtual std::string info_data(std::string msg);

    // ===========================================
    // Functions
    // ===========================================
    // compute
    void execute();
};






// ===========================================
//      Functions of Class metric
// ===========================================


// ===========================================
// Create Functions
// ===========================================
// Constructor
template <typename pixel_type>
myclass<pixel_type>::myclass()
{
    init();
};

template <typename pixel_type>
void myclass<pixel_type>::init()
{
    ;
};


template <typename pixel_type>
std::string myclass::<pixel_type>::info(std::string msg)
{
    std::stringstream ss;
    std::string title = "Optimizer Information";
    if (msg != "") { title = msg; };
    // Summary of the optimizer information
    ss << object<pixel_type>::info(title);
};

}; //end namespace

#endif