/*
 * GNU GENERAL PUBLIC LICENSE
 *
 * Copyright (C) 2017
 * Created by Leonardo Parisi (leonardo.parisi[at]gmail.com)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
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
#include "./utils/poe.hpp"

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
