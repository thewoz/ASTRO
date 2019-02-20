# ASTRO


Nel folder Converter:

ecef.hpp      -- Funzioni di converzione da ECEF a TEME e J2000
eci.hpp        -- Funzioni di conversione da J2000 a TEME e ECEF
teme.hpp     -- Funzioni di converzione da TEME a ECEF e j2000
ll2ecef.hpp  -- Funzioni di converzione da Latitudine e Longitudine in ECEF
iau80.hpp    -- Funzioni per la gestione dei parametri iau80 che servono nelle funzioni di converzione
eopc.hpp     -- Funzioni per la gestione dei parametri eopc che servono nelle funzioni di converzione


Nel folder Propagator:

SGP4.hpp  -- Funzioni per la propagazione con SGP4


Nel folder Utils:

Attitude.hpp -- Funzioni per propagare l'assetto
RK4.hpp  -- Funzione per itegrare tramite Runge-Kutta
Cubicspline.hpp -- Funzione per fittare con una spline
Angle.hpp  -- Funzioni per convertire tra angoli e radianti
Curl.hpp -- Funzioni per scaricare file dal web
Quaternion.hpp -- Funzioni per gestire i quaterioni
TLE.hpp -- Funzioni per gestire i TLE
Date.hpp -- Funzioni per gestire le date
Lagrange.hpp -- Funzione per fittare tramite Lagrange
Utils,hpp -- Funzioni per calcolarsi l'elevazione del satellite , l'angolo theta con il sole e il fattore di attenuazione atmosferico


Nel folder Corps:

Satellite.hpp -- Classe che gestisce il satellite
Observatory.hpp -- Classe che gestisce l'osservatorio
Sun.hpp -- Classe che gestisce il sole


TO FIX:
1) Controllare le funzioni che riguardano la posizione del sole
2) Controllare le funzioni che rigurdano la conversione da Longitudine e Latidutine in ECEF:
    Ce ne sta una in lla2Ecef.h (quella che uso) una in vallado lla2ecef() e una in Observatory.hpp (_convert())
3) Controllare che gli step di integrazione del orbita del satellite del sole e della stazione a terra siano corretti


TODO:
1) Finire di implementare le funzioni di calcolo di RA DEC del satellite e del sole
2) Finire di implementare le funzioni di calcolo di illuminazione
3) Implementare le funzioni di scarico del TLE


REMEMBER:
1) La funzione invjday() di SPG4 è stata modifica in Vallado
2) La funzione convtime in Vallado è stata modificata
