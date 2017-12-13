#pragma once
#include "initial_condition.hpp"
#include <iostream>
#include <math.h>
using namespace std;

class phys_quantity{
private:
  double density;
  double velocity;
  double pressure;
  double sound_speed;
public:
  phys_quantity(){
    density = 0.0;
    velocity = 0.0;
    pressure = 0.0;
    sound_speed = 0.0;
  }
  double get_density()  {return density;}
  double get_velocity() {return velocity;}
  double get_pressure() {return pressure;}
  double get_sound()    {return sound_speed;}
  double get_mach()     {return velocity/sound_speed;}

  double set_density  (double new_dens)   {density = new_dens;}
  double set_velocity (double new_vel)    {velocity = new_vel;}
  double set_pressure (double new_press)  {pressure = new_press;}
  double set_sound    (double new_sound)  {sound_speed = new_sound;}

  void calc_sound () {
    sound_speed = sqrt(GAMMA*pressure/density);
  }
  void set_condition(){
    double dens, vel, press;
    cout << "  'density' 'velocity' 'pressure' = ";
    for( ;!(cin >> density >> velocity >> pressure); ){
      cout << "Loading Error!!" << endl;
    }
    calc_sound();
  }

};
