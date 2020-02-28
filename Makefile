INCLUDE_INSTALL_PATH=/usr/local/include
LIBRARY_INSTALL_PATH=/usr/local/lib

LIBRARY_NAME=astro


all:


setup:


development:
	@rm -rf $(INCLUDE_INSTALL_PATH)/$(LIBRARY_NAME)
	@ln -s $(shell pwd)/include $(INCLUDE_INSTALL_PATH)/$(LIBRARY_NAME)
	@cp astro.pc /usr/local/lib/pkgconfig/

install:
	rm -rf $(INCLUDE_INSTALL_PATH)/$(LIBRARY_NAME)
	mkdir -p $(INCLUDE_INSTALL_PATH)/$(LIBRARY_NAME)
	cp -R ./include/* $(INCLUDE_INSTALL_PATH)/$(LIBRARY_NAME)
	@cp astro.pc /usr/local/lib/pkgconfig/


uninstall:
	@rm -rf $(INCLUDE_INSTALL_PATH)/$(LIBRARY_NAME)


clean:
	@rm -rf ./bin


test:
	@mkdir -p ./bin
	g++ -march=native -Os -std=c++17 -I./include -DTESTEOPC     ./main.cpp -L/usr/local/lib/ -lcoordFK5 -lastTime -lastMath -lastIOD -lSGP4 -lastUtils -last2Body -lEopSpw -lcurl -o ./bin/testEopc
	g++ -march=native -Os -std=c++17 -I./include -DTESTDATE     ./main.cpp -L/usr/local/lib/ -lcoordFK5 -lastTime -lastMath -lastIOD -lSGP4 -lastUtils -last2Body -lEopSpw -lcurl -o ./bin/testDate
	g++ -march=native -Os -std=c++17 -I./include -DTESTSATELLITE ./main.cpp -L/usr/local/lib/ -lcoordFK5 -lastTime -lastMath -lastIOD -lSGP4 -lastUtils -last2Body -lEopSpw -lcurl -o ./bin/testSatellite
	g++ -march=native -Os -std=c++17 -I./include -DTESTSUN      ./main.cpp -L/usr/local/lib/ -lcoordFK5 -lastTime -lastMath -lastIOD -lSGP4 -lastUtils -last2Body -lEopSpw -lcurl -o ./bin/testSun
	g++ -march=native -Os -std=c++17 -I./include -DTESTSTATION  ./main.cpp -L/usr/local/lib/ -lcoordFK5 -lastTime -lastMath -lastIOD -lSGP4 -lastUtils -last2Body -lEopSpw -lcurl -o ./bin/testStation
	g++ -march=native -Os -std=c++17 -I./include -DTESTATTITUDE ./main.cpp -L/usr/local/lib/ -lcoordFK5 -lastTime -lastMath -lastIOD -lSGP4 -lastUtils -last2Body -lEopSpw -lcurl -o ./bin/testAttitude
	g++ -march=native -Os -std=c++17 -I./include -DTESTCONVERT  ./main.cpp -L/usr/local/lib/ -lcoordFK5 -lastTime -lastMath -lastIOD -lSGP4 -lastUtils -last2Body -lEopSpw -lcurl -o ./bin/testConvert








