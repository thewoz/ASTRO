# ASTRO


## File list

###### Nel folder Converter:

* ecef.hpp      -- Funzioni di converzione da ECEF a TEME e J2000
* eci.hpp        -- Funzioni di conversione da J2000 a TEME e ECEF
* teme.hpp     -- Funzioni di converzione da TEME a ECEF e j2000
* ll2ecef.hpp  -- Funzioni di converzione da Latitudine e Longitudine in ECEF
* iau80.hpp    -- Funzioni per la gestione dei parametri iau80 che servono nelle funzioni di converzione
* eopc.hpp     -- Funzioni per la gestione dei parametri eopc che servono nelle funzioni di converzione


###### Nel folder Propagator:

* SGP4.hpp  -- Funzioni per la propagazione con SGP4


###### Nel folder Utils:

* Attitude.hpp -- Funzioni per propagare l'assetto
* RK4.hpp  -- Funzione per itegrare tramite Runge-Kutta
* Cubicspline.hpp -- Funzione per fittare con una spline
* Angle.hpp  -- Funzioni per convertire tra angoli e radianti
* Curl.hpp -- Funzioni per scaricare file dal web
* Quaternion.hpp -- Funzioni per gestire i quaterioni
* TLE.hpp -- Funzioni per gestire i TLE
* Date.hpp -- Funzioni per gestire le date
* Lagrange.hpp -- Funzione per fittare tramite Lagrange
* Utils,hpp -- Funzioni per calcolarsi l'elevazione del satellite , l'angolo theta con il sole e il fattore di attenuazione atmosferico


###### Nel folder Corps:

* Satellite.hpp -- Classe che gestisce il satellite
* Observatory.hpp -- Classe che gestisce l'osservatorio
* Sun.hpp -- Classe che gestisce il sole


## FIXME

* Controllare le funzioni che rigurdano la conversione da Longitudine e Latidutine in ECEF:
  Ce ne sta una in lla2Ecef.h (quella che uso) una in vallado lla2ecef() e una in Observatory.hpp (_convert())
* Controllare che gli step di integrazione del orbita del satellite del sole e della stazione a terra siano corretti
* Controllare la posizione del sole in che sistema di riferimento si trova
* Controllare perche nell'ossservatorio non chiamo site di vallado

## TODO

* Finire di implementare le funzioni di calcolo di RA DEC del satellite e del sole
* Finire di implementare le funzioni di calcolo di illuminazione
* Implementare le funzioni di scarico del TLE


## NOTE

* La funzione invjday() di SPG4 è stata modifica in Vallado
* La funzione convtime in Vallado è stata modificata


## INSTALL

- Prima va installata la libreria di Vallado 
- Poi fare il make install della libreria Astro
