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
#include <chrono>
#include <ctime>
#include <cmath>

// class that works as a stopwatch
class timer
{

private:
    std::clock_t cpu_t1; // start clock
    std::clock_t cpu_t2; // since start
    std::clock_t cpu_t3; // since last lap
    std::chrono::system_clock::time_point t1; // start clock
    std::chrono::system_clock::time_point t2; // since start
    std::chrono::system_clock::time_point t3; // since last lap

    bool laps;
    double cpu_lap_time;
    double cpu_total_time;
    double lap_time;
    double total_time;
    double scale;
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
    std::chrono::system_clock::time_point get_clock();
    double get_elapsed();
    double get_time(); // get final time (total)
    std::clock_t get_clock_cpu();
    double get_elapsed_cpu();
    double get_time_cpu(); // get final time (total)
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
    void lap(std::string msg="");
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
std::clock_t timer::get_clock_cpu()
{
    return cpu_t1;
};

double timer::get_elapsed_cpu()
{
    return cpu_lap_time;
};

double timer::get_time_cpu()
{
    return cpu_total_time;
};

std::chrono::system_clock::time_point timer::get_clock()
{
    return t1;
};

double timer::get_elapsed()
{
    return lap_time;
};

double timer::get_time()
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
    std::string title = "[Timer]";
    if (msg != "") { title = msg; };
    // more info here
    ss << title;
    
    if (laps)
    {
        ss << "[Lap time][" << lap_time;
        ss <<" " << units << "]\t";
    }
    ss << "[Elapsed time][" << total_time;
    ss <<" " << units << "]";
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
        std::clock_t taux = clock();
        auto wall_taux = std::chrono::system_clock::now();
        cpu_t1 = cpu_t2 = cpu_t3 = taux;
        t1 = t2 = t3 = wall_taux;
        state = "running";
        laps = false;
    }
    else { timer::warning(); };
};

void timer::lap(std::string msg)
{
    if (state == "running")
    {
        // measure clocks
        std::clock_t taux = clock();
        auto wall_taux = std::chrono::system_clock::now();

        // cpu clock
        cpu_t2 = taux - cpu_t1;
        cpu_t3 = taux - cpu_t3;
        cpu_total_time = ((double)cpu_t2)*scale/CLOCKS_PER_SEC;
        cpu_lap_time = ((double)cpu_t3)*scale/CLOCKS_PER_SEC;

        //wall clock
        // t2 = wall_taux - t1;
        // t3 = wall_taux - t3;
        total_time = scale*(std::chrono::duration<double>(wall_taux - t1).count());
        lap_time = scale*(std::chrono::duration<double>(wall_taux - t3).count());

        cpu_t3 = taux;
        t3 = wall_taux;
        laps = true;
        
        if (msg != ""){ timer::print(msg); };
    }
    else { timer::warning(); };
};

void timer::finish(bool show)
{
    if (state == "running")
    {
        // measure clocks
        std::clock_t taux = clock();
        auto wall_taux = std::chrono::system_clock::now();

        // cpu clock
        cpu_t2 = taux - cpu_t1;
        cpu_t3 = taux - cpu_t3;
        cpu_total_time = ((double)cpu_t2)*scale/CLOCKS_PER_SEC;
        cpu_lap_time = ((double)cpu_t3)*scale/CLOCKS_PER_SEC;

        //wall clock
        // t2 = wall_taux - t1;
        // t3 = wall_taux - t3;
        total_time = scale*(std::chrono::duration<double>(wall_taux - t1).count());
        lap_time = scale*(std::chrono::duration<double>(wall_taux - t3).count());

        laps = false;

        if (show){ timer::print(); };

        // Reset all clocks
        cpu_t1 = cpu_t2 = cpu_t3 = 0;
        t1 = t2 = t3 = std::chrono::system_clock::time_point();
        state = "reset";
        // timer::start(); // ??
    }
    else { this->warning(); };
};

#endif