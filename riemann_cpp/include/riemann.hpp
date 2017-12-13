#pragma once

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
  
}
