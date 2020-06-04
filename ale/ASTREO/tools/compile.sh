#!/bin/bash

# compile viewsat
g++ -o viewsat.x viewsat.cpp -lm -L/usr/local/lib -lEopSpw -last2Body -lastMath -lcoordFK5 -lastUtils -lSGP4 -lastIOD -lastTime  -I/usr/local/include/vallado

# compile photo
g++ -o photosat.x photosat.cpp -lm -L/usr/local/lib -lEopSpw -last2Body -lastMath -lcoordFK5 -lastUtils -lSGP4 -lastIOD -lastTime  -I/usr/local/include/vallado -I/usr/local/include/astro
g++ -o photostar.x photostar.cpp -lm -L/usr/local/lib -lEopSpw -last2Body -lastMath -lcoordFK5 -lastUtils -lSGP4 -lastIOD -lastTime  -I/usr/local/include/vallado -I/usr/local/include/astro

# compile sorting
g++ -o sort.x sort.cpp
