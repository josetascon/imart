/*
* @Author: jose
* @Date:   2020-08-20 00:00:00
* @Last Modified by:   jose
* @Last Modified time: 2020-08-20 00:00:00
*/

#ifndef __TRACKER_H__
#define __TRACKER_H__

#include "process_object.h"
#include "image.h"
#include "transform.h"
#include "interpolator.h"
// #include "viewer.h"

namespace imart
{

// Class tracker
template <typename type, typename container=vector_cpu<type>>
class tracker: public inherit<tracker<type,container>, process_object>
{

public:
    //Type definitions
    using self    = tracker;
    using pointer = std::shared_ptr<self>;
    using vector  = std::vector<self::pointer>;
    
protected:
    // ===========================================
    // Internal Variables
    // ===========================================
    typename image<type,container>::pointer fixed;
    typename image<type,container>::pointer moving;
    typename image<type,container>::pointer current;
    std::vector<std::vector<int>> box_fixed;
    std::vector<std::vector<int>> box_moving;
    std::vector<std::vector<int>> box_current;
    bool current_mode;
    int count_current;

    // ===========================================
    // Functions
    // ===========================================
    void init(int d);
    void init( typename image<type,container>::pointer fixed_image, 
               typename image<type,container>::pointer moving_image,
               std::vector<std::vector<int>> bbox_fixed );

public:
    // ===========================================
    // Create Functions
    // ===========================================
    // Constructors
    tracker();
    tracker(int d);
    tracker(typename image<type,container>::pointer fixed_image, 
            std::vector<std::vector<int>> bbox_fixed);
    
    // ===========================================
    // Get Functions
    // ===========================================
    typename image<type,container>::pointer get_fixed() const;
    typename image<type,container>::pointer get_moving() const;
    typename image<type,container>::pointer get_current() const;
    std::vector<std::vector<int>> get_box_fixed() const;
    std::vector<std::vector<int>> get_box_moving() const;
    std::vector<std::vector<int>> get_box_current() const;
    bool get_current_mode() const;
    
    // ===========================================
    // Set Functions
    // ===========================================
    virtual void set_fixed(typename image<type,container>::pointer fixed_image);
    virtual void set_moving(typename image<type,container>::pointer moving_image);
    virtual void set_current(typename image<type,container>::pointer current_image);
    virtual void set_box_fixed(std::vector<std::vector<int>> bbox_fixed);
    virtual void set_box_moving(std::vector<std::vector<int>> bbox_moving);
    virtual void set_box_current(std::vector<std::vector<int>> bbox_current);
    virtual void set_current_mode(bool mode);
    
    // ===========================================
    // Print Functions
    // ===========================================
    std::string info(std::string msg);

    // ===========================================
    // Functions
    // ===========================================
    virtual std::vector<std::vector<int>> apply(typename image<type,container>::pointer moving_image);
    virtual std::vector<std::vector<int>> operator () (typename image<type,container>::pointer moving_image);
};


// ===========================================
//      Functions of Class tracker
// ===========================================

// ===========================================
// Create Functions
// ===========================================
// Constructor empty
template <typename type, typename container>
tracker<type,container>::tracker()
{
    this->class_name = "tracker";
    init(2);
};

template <typename type, typename container>
tracker<type,container>::tracker(int d)
{
    this->class_name = "tracker";
    init(d);
};

template <typename type, typename container>
tracker<type,container>::tracker( typename image<type,container>::pointer fixed_image, 
                                std::vector<std::vector<int>> bbox_fixed)
{
    this->class_name = "tracker";
    auto moving_image = image<type,container>::new_pointer(fixed_image->get_dimension());
    init(fixed_image, moving_image, bbox_fixed);
};

template <typename type, typename container>
void tracker<type,container>::init(int d)
{
    // std::cout << "init tracker" << std::endl;
    auto fixed_image = image<type,container>::new_pointer(d);
    auto moving_image = image<type,container>::new_pointer(d);
    std::vector<std::vector<int>> box = std::vector<std::vector<int>>{std::vector<int>(d,0),std::vector<int>(d,0)};
    
    init(fixed_image, moving_image, box);
    // std::cout << "end init tracker" << std::endl;
}

template <typename type, typename container>
void tracker<type,container>::init( typename image<type,container>::pointer fixed_image, 
                                    typename image<type,container>::pointer moving_image,
                                    std::vector<std::vector<int>> bbox_fixed )
{
    // std::cout << "tracker init" << std::endl;
    assert(fixed_image->get_dimension() == moving_image->get_dimension());
    int d = fixed_image->get_dimension();

    current_mode = false;

    fixed = fixed_image;
    moving = moving_image;
    current = fixed_image->clone();
    count_current = 0;
    box_fixed = bbox_fixed;
    box_moving = std::vector<std::vector<int>>{std::vector<int>(d,0),std::vector<int>(d,0)};
    box_current = bbox_fixed;
    
    this->set_total_inputs(2);      //process_object::init
    this->set_total_outputs(0);     //process_object::init
    this->setup_input(fixed, moving);
    this->setup_output();           // no defined output, this object is a base class
    // std::cout << "tracker init end" << std::endl;
    return;
};

// ===========================================
// Get Functions
// ===========================================
template <typename type, typename container>
typename image<type,container>::pointer tracker<type,container>::get_fixed() const
{
    return fixed;
};

template <typename type, typename container>
typename image<type,container>::pointer tracker<type,container>::get_moving() const
{
    return moving;
};

template <typename type, typename container>
typename image<type,container>::pointer tracker<type,container>::get_current() const
{
    return current;
};

template <typename type, typename container>
std::vector<std::vector<int>> tracker<type,container>::get_box_fixed() const
{
    return box_fixed;
};

template <typename type, typename container>
std::vector<std::vector<int>> tracker<type,container>::get_box_moving() const
{
    return box_moving;
};

template <typename type, typename container>
std::vector<std::vector<int>> tracker<type,container>::get_box_current() const
{
    return box_current;
};

template <typename type, typename container>
bool tracker<type,container>::get_current_mode() const
{
    return current_mode;
};

// ===========================================
// Set Functions
// ===========================================
template <typename type, typename container>
void tracker<type,container>::set_fixed(typename image<type,container>::pointer fixed_image)
{
    fixed = fixed_image;
    this->set_input(0, fixed); // update process input
};

template <typename type, typename container>
void tracker<type,container>::set_moving(typename image<type,container>::pointer moving_image)
{
    moving = moving_image;
    this->set_input(1, moving); // update process input
};

template <typename type, typename container>
void tracker<type,container>::set_current(typename image<type,container>::pointer current_image)
{
    current = current_image;
};

template <typename type, typename container>
void tracker<type,container>::set_box_fixed(std::vector<std::vector<int>> bbox_fixed)
{
    box_fixed = bbox_fixed;
};

template <typename type, typename container>
void tracker<type,container>::set_box_moving(std::vector<std::vector<int>> bbox_moving)
{
    box_moving = bbox_moving;
};

template <typename type, typename container>
void tracker<type,container>::set_box_current(std::vector<std::vector<int>> bbox_current)
{
    box_current = bbox_current;
};

template <typename type, typename container>
void tracker<type,container>::set_current_mode(bool mode)
{
    current_mode = mode;
};

// ===========================================
// Print Functions
// ===========================================
template <typename type, typename container>
std::string tracker<type,container>::info(std::string msg)
{
    std::stringstream ss;
    std::string title = "Tracker Object Information";
    if (msg != "") { title = msg; };

    ss << object::info(title);
    ss << process_object::info("");
    return ss.str();
};

// ===========================================
// Functions
// ===========================================
template <typename type, typename container>
std::vector<std::vector<int>> tracker<type,container>::apply(typename image<type,container>::pointer moving_image)
{
    return box_moving;
};

template <typename type, typename container>
std::vector<std::vector<int>> tracker<type,container>::operator() (typename image<type,container>::pointer moving_image)
{
    return apply(moving_image);
};

}; //end namespace

#endif