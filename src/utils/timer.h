/*
* @Author: jose
* @Date:   2019-11-15 13:55:00
* @Last Modified by:   jose
* @Last Modified time: 2019-11-15 13:55:00
*/

#ifndef __TIMER_H__
#define __TIMER_H__

// std libs
#include <iostream>     // std::cout
#include <sstream>      // stringstream
#include <ctime>
#include <cmath>

// class that works as a stopwatch
class timer
{

private:
    clock_t t1; // start clock
    clock_t t2; // since start
    clock_t t3; // since last lap
    bool laps;
    float lap_time;
    float total_time;
    float scale;
    std::string state;
    std::string units;

    void warning(); // warning message

public:
    // ===========================================
    // Create Functions
    // ===========================================
    timer();
    timer(std::string units);
    ~timer();
    // ===========================================
    // Get Functions
    // ===========================================
    clock_t get_clock();
    float get_elapsed();
    float get_time(); // get final time (total)
    std::string get_state();
    std::string get_units();

    void set_units(std::string units = "ms");

    // ===========================================
    // Print Functions
    // ===========================================
    void print(std::string msg = "");
    std::string info(std::string msg);
    friend std::ostream & operator << (std::ostream & os, timer & tt);

    // ===========================================
    // Functions
    // ===========================================
    void start(); // run
    void lap(bool show=true);
    void finish(bool show=false); // reset
};








// ===========================================
//          Functions of Class timer
// ===========================================
timer::timer(): timer("s") { };

timer::timer(std::string units)
{
    scale = 1.0;
    state = "zero";
    bool laps = false;
    timer::set_units(units);
}

timer::~timer()
{
    ;
};

// ===========================================
// Get Functions
// ===========================================
clock_t timer::get_clock()
{
    return t1;
};

float timer::get_elapsed()
{
    return lap_time;
};

float timer::get_time()
{
    return total_time;
};

std::string timer::get_state()
{
    return state;
};

std::string timer::get_units()
{
    return units;
};

void timer::set_units(std::string units)
{
    if (units == "ms")
    {
        this->units = units;
        scale = pow(10.0,3);
    }
    else if (units == "us")
    {
        this->units = units;
        scale = pow(10.0,6);
    }
    else if (units == "ns")
    {
        this->units = units;
        scale = pow(10.0,9);
    }
    else
    {
        this->units = "s";
        scale = 1.0;
    };
};

// ===========================================
// Print Functions
// ===========================================
void timer::print(std::string msg)
{
    std::cout << timer::info(msg);
};

std::string timer::info(std::string msg)
{
    std::stringstream ss;
    std::string title = "[timer]";
    if (msg != "") { title = msg; };
    // more info here
    ss << title << "\t\t";
    
    if (laps)
    {
        ss << "lap time: " << lap_time;
        ss <<" [" << units << "]\t||\t";
    }
    ss << "total time: " << total_time;
    ss <<" [" << units << "]";
    ss << std::endl;
    
    return ss.str();
};

std::ostream & operator << (std::ostream & os, timer & tt)
{
    os << tt.info("");
    return os;
};

void timer::warning()
{
    std::cout << "[warning]\t\t timer not initialized";
    std::cout << std::endl;
};

// ===========================================
// Functions
// ===========================================
void timer::start()
{
    if (state == "zero" || state == "reset" )
    {
        clock_t taux = clock();
        t1 = taux;
        t2 = taux;
        t3 = taux;
        state = "running";
        laps = false;
    }
    else { timer::warning(); };
};

void timer::lap(bool show)
{
    if (state == "running")
    {
        clock_t taux = clock();
        t2 = taux - t1;
        t3 = taux - t3;

        total_time = ((float)t2)*scale/CLOCKS_PER_SEC;
        lap_time = ((float)t3)*scale/CLOCKS_PER_SEC;

        t3 = taux;
        laps = true;
        
        if (show){ timer::print(); };
    }
    else { timer::warning(); };
};

void timer::finish(bool show)
{
    if (state == "running")
    {
        clock_t taux = clock();
        t2 = taux - t1;
        t3 = taux - t3;

        total_time = ((float)t2)*scale/CLOCKS_PER_SEC;
        lap_time = ((float)t3)*scale/CLOCKS_PER_SEC;

        if (show){ timer::print(); };

        // Reset all clocks
        t1 = t2 = t3 = 0;
        state = "reset";
        // timer::start(); // ??
    }
    else { this->warning(); };
};

#endif