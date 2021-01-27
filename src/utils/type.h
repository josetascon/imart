/*
* @Author: jose
* @Date:   2020-10-08 00:00:00
* @Last Modified by:   jose
* @Last Modified time: 2020-10-08 00:00:00
*/

#ifndef __UTILS_TYPE_H__
#define __UTILS_TYPE_H__

// std libs
#include <iostream>     // std::cout
#include <string>       // std::string
#include <typeinfo>     // operator typeids

template<typename type>
std::string string_type(){ return typeid(type).name();};
template<>
std::string string_type<unsigned char>(){ return "unsigned char";};
template<>
std::string string_type<char>(){ return "char";};
template<>
std::string string_type<unsigned short>(){ return "unsigned short";};
template<>
std::string string_type<short>(){ return "short";};
template<>
std::string string_type<unsigned int>(){ return "unsigned int";};
template<>
std::string string_type<int>(){ return "int";};
template<>
std::string string_type<float>(){ return "float";};
template<>
std::string string_type<double>(){ return "double";};

#endif