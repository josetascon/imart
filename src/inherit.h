/*
* @Author: jose
* @Date:   2020-06-03 00:00:00
* @Last Modified by:   jose
* @Last Modified time: 2020-06-03 00:00:00
*/

#ifndef __INHERIT_H__
#define __INHERIT_H__

// std libs
#include <iostream>     // std::cout
#include <memory>       // std::shared_ptr

namespace imart
{

//! Class as base of everything
class base
{
public:
    virtual ~base() = default;

    std::shared_ptr<base> clone() const
    {
        return std::shared_ptr<base>(this->clone_method());
    };

    std::shared_ptr<base> copy() const
    {
        return std::shared_ptr<base>(this->copy_method());
    };

    std::shared_ptr<base> mimic() const
    {
        return std::shared_ptr<base>(this->mimic_method());
    };

    private:
    virtual base * clone_method() const = 0;
    virtual base * copy_method() const = 0;
    virtual base * mimic_method() const = 0;
};

//! Class template to allow virtual inheritance
template <typename t>
class virtualize : virtual public t
{
   using t::t;
};

//! Class template to allow multiple inheritance and/or virtual inheritance
template <typename self, typename ... parent>
class inherit_multiple : public parent ...
{
public:

    std::shared_ptr<self> clone() const
    {
        return std::shared_ptr<self>(static_cast<self *>(this->clone_method()));
    };

    std::shared_ptr<self> copy() const
    {
        return std::shared_ptr<self>(static_cast<self *>(this->copy_method()));
    };

    std::shared_ptr<self> mimic() const
    {
        return std::shared_ptr<self>(static_cast<self *>(this->mimic_method()));
    };

    template<typename... ARGS>
    static std::shared_ptr<self> new_pointer(const ARGS&... args)
    {
        return std::make_shared<self>(args...);
    };

protected:
    virtual std::string info(std::string msg){ return msg;}; // defined to avoid multiple versions

private:
    virtual inherit_multiple * clone_method() const override
    {
        // std::cout << "clone method inherit" << std::endl;
        // self * p = new self(static_cast<const self & >(*this));
        self * p = new self();
        p->clone_(static_cast<const self & >(*this));
        return p;
    };
    virtual inherit_multiple * copy_method() const override
    {
        // std::cout << "copy method inherit" << std::endl;
        self * p = new self();
        p->copy_(static_cast<const self & >(*this));
        return p;
    };
    virtual inherit_multiple * mimic_method() const override
    {
        // std::cout << "mimic method inherit" << std::endl;
        self * p = new self();
        p->mimic_(static_cast<const self & >(*this));
        return p;
    };
};

//! Class template to allow single inheritance and transfer parent class constructor
template <typename self, typename parent>
class inherit : public parent
{
public:
    using parent::parent;

    std::shared_ptr<self> clone() const
    {
        return std::shared_ptr<self>(static_cast<self *>(this->clone_method()));
    };

    std::shared_ptr<self> copy() const
    {
        return std::shared_ptr<self>(static_cast<self *>(this->copy_method()));
    };

    std::shared_ptr<self> mimic() const
    {
        return std::shared_ptr<self>(static_cast<self *>(this->mimic_method()));
    };

    template<typename... ARGS>
    static std::shared_ptr<self> new_pointer(const ARGS&... args)
    {
        return std::make_shared<self>(args...);
    };

private:
    virtual inherit * clone_method() const override
    {
        // std::cout << "clone method inherit" << std::endl;
        // self * p = new self(static_cast<const self & >(*this));
        self * p = new self();
        p->clone_(static_cast<const self & >(*this));
        return p;
    };
    virtual inherit * copy_method() const override
    {
        // std::cout << "copy method inherit" << std::endl;
        self * p = new self();
        p->copy_(static_cast<const self & >(*this));
        return p;
    };
    virtual inherit * mimic_method() const override
    {
        // std::cout << "mimic method inherit" << std::endl;
        self * p = new self();
        p->mimic_(static_cast<const self & >(*this));
        return p;
    };
};

}; //end namespace

#endif