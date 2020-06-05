/*
 * MIT License
 *
 * Copyright © 2017 S5Lab
 * Created by Leonardo Parisi (leonardo.parisi[at]gmail.com)
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef _H_ASTRO_H
#define _H_ASTRO_H

// Vallado
#include <vallado/astUtils.h>
#include <vallado/astTime.h>
#include <vallado/EopSpw.h>
#include <vallado/coordFK5.h>
#include <vallado/SGP4.h>
#include <vallado/astTime.h>
#include <vallado/astMath.h>
#include <vallado/ast2Body.h>
#include <vallado/astIOD.h>

// Utils
#include "./utils/define.h"
#include "./utils/angle.hpp"
#include "./utils/utils.hpp"
#include "./utils/date.hpp"
#include "./utils/tle.hpp"
#include "./utils/rk4.hpp"
#include "./utils/quaternion.hpp"
#include "./utils/attitude.hpp"
#include "./utils/curl.hpp"

// Converter
#include "./converter/iau80.hpp"
#include "./converter/eopc.hpp"
#include "./converter/ecef.hpp"
#include "./converter/teme.hpp"
#include "./converter/eci.hpp"
#include "./converter/lla2ecef.hpp"
#include "./converter/eci2jnow.hpp"
#include "./converter/converter.h"
#include "./converter/effect.hpp"
#include "./utils/photo.hpp" // Questo non adrebbe qua

// Propagator
#include "./propagator/sgp4.hpp"

// Corps
#include "./corps/observatory.hpp"
#include "./corps/satellite.hpp"
#include "./corps/sun.hpp"



#endif /* _H_ASTRO_H */
