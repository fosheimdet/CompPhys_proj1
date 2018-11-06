#pragma once
#include "CelestialBody.h"
void euler(CelestialBody& theBody, double t_max, int N);
void velocityVerlet(CelestialBody& theBody, double t_max, int N, bool relativistic); 