/*
* @Author: jose
* @Date:   2020-10-08 00:00:00
* @Last Modified by:   jose
* @Last Modified time: 2020-10-08 00:00:00
*/

#ifndef __UTILS_ASSERTION_H__
#define __UTILS_ASSERTION_H__

// std libs
#include <iostream>     // std::cout
#include <string>       // std::string
#include <cassert>      // assert

// boost libs
#include <boost/stacktrace.hpp>

# define imart_assert(status, msg) \
    imart_assert_error(status, __FILE__, __LINE__, msg)

void imart_assert_error(bool status, const char* file, int line, const char* msg)
{
    if (not status)
    {
        std::cerr << "\n******* IMART Error *******"
                  << "\n[Error] Information:\t" << msg
                  << "\n[Error] Error code:\t" << status
                  // << "\n[Error] Description:\t" << imart_error_code(status)
                  << "\n[Error] File:\t\t" << file
                  << "\n[Error] Line:\t\t" << line// << std::endl;
                  << "\n[Error] Backtrace:\n" << boost::stacktrace::stacktrace() << std::endl;
        assert(status);
        exit(0);
    };
};

std::string imart_error_code(int error)
{
    // Print error message
    std::string error_code;
    switch(error)
    {
    case -1:
        error_code = "imart fault";
        break;
    default:
        error_code = "imart assert code " + std::to_string(error);
    }
    return error_code;
};

#endif