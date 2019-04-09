INCLUDE_INSTALL_PATH=/usr/local/include
LIBRARY_INSTALL_PATH=/usr/local/lib

LIBRARY_NAME=astro


all:


setup:


development:
	@rm -rf $(INCLUDE_INSTALL_PATH)/$(LIBRARY_NAME)
	@ln -s $(shell pwd)/include $(INCLUDE_INSTALL_PATH)/$(LIBRARY_NAME)

install:
	rm -rf $(INCLUDE_INSTALL_PATH)/$(LIBRARY_NAME)
	mkdir -p $(INCLUDE_INSTALL_PATH)/$(LIBRARY_NAME)
	cp -R ./include/* $(INCLUDE_INSTALL_PATH)/$(LIBRARY_NAME)


uninstall:
	@rm -rf $(INCLUDE_INSTALL_PATH)/$(LIBRARY_NAME)


clean:
	@rm -rf ./bin


test:
	@mkdir -p ./bin
	g++ -march=native -Os -std=c++17 -I./include -I/opt//local/include/ -DTESTEOPC     ./main.cpp -L/usr/local/lib/ -lcoordFK5 -lastTime -lastMath -lastIOD -lSGP4 -lastUtils -last2Body -lEopSpw -lcurl -o ./bin/testEopc
	g++ -march=native -Os -std=c++17 -I./include -I/opt//local/include/  -DTESTDATE     ./main.cpp -L/usr/local/lib/ -lcoordFK5 -lastTime -lastMath -lastIOD -lSGP4 -lastUtils -last2Body -lEopSpw -lcurl -o ./bin/testDate
	g++ -march=native -Os -std=c++17 -I./include -I/opt//local/include/  -DTESTSATELLITE ./main.cpp -L/usr/local/lib/ -lcoordFK5 -lastTime -lastMath -lastIOD -lSGP4 -lastUtils -last2Body -lEopSpw -lcurl -o ./bin/testSatellite
	g++ -march=native -Os -std=c++17 -I./include -I/opt//local/include/  -DTESTSUN      ./main.cpp -L/usr/local/lib/ -lcoordFK5 -lastTime -lastMath -lastIOD -lSGP4 -lastUtils -last2Body -lEopSpw -lcurl -o ./bin/testSun
	g++ -march=native -Os -std=c++17 -I./include -I/opt//local/include/  -DTESTSTATION  ./main.cpp -L/usr/local/lib/ -lcoordFK5 -lastTime -lastMath -lastIOD -lSGP4 -lastUtils -last2Body -lEopSpw -lcurl -o ./bin/testStation
	g++ -march=native -Os -std=c++17 -I./include -I/opt//local/include/  -DTESTATTITUDE ./main.cpp -L/usr/local/lib/ -lcoordFK5 -lastTime -lastMath -lastIOD -lSGP4 -lastUtils -last2Body -lEopSpw -lcurl -o ./bin/testAttitude




