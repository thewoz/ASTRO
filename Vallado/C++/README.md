I sorgenti del codice di vallado sono stati scaricati dal sito CelesTrak sono aggiornati alla versione rilasciata il 18 Aprile 2018

Modifiche dei sorgenti che non alterano il funzionamento del codice originale:

1) Sono stati cancellati tutti i caratteri non ASCII
2) Sono stati tolti gli include locali
3) Sono state modificate tutte le strcpy_s in strcpy
4) Sono state modificate tutte le scanf_s in scanf
5) Sono stati correti i commenti lasciati aperti
6) Sono stati corretti gli errori di tipo nelle scanf e printf
7) Sono stati modificati i pi con M_PI
8) Sono stati modificati gli include passando ad esempio da math.h a cmath
9) Sono stati tolti la maggior parte dei warning in compilazione
10) È stata aggiunta una versione con meno argomenti di convtime( ) in astTime
11) Sono stati eliminati i file non utili

Modifiche che posso alterare il funzionamento del codice rispetto all'originale:

1) Modificata convtime( ) in astTime.h allineandola con quella della versione in Matlab
2) Modificata invjday( ) in SPG4.h rendendola uguale a quella in astTime.h

È stato aggiunto il file astUtils.h con al suo interno le funzioni (prese dalla versione di Matlab):

1) sight( )
2) light( )


All'interno si trova sia il codice in C++ come sopra modificato sia la verione in Matlab.
È stato aggiunto il progetto in xcode per compilare il codice C++ 

Nel Makefile sono presenti i seguenti comandi:

1) make (che compila solamente il codice)
2) make install che compila e installa il codice
3) make uninstall che lo disistalla


Di fatto per installarlo basta fare "make install" da linea di comando

Per utilizzarlo in altri progetti ad esempio:

Per gli header in compilazione -I/usr/local/include/  nel codice #include <vallado/SGP4.h>)
Per le librerie in compilazione -L/usr/local/lib -l/usr/local/include/vallado/ -lSGP4

