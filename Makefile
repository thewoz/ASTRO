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


#test:
#	@mkdir -p ./bin
#	g++ -march=native -Os -std=c++11 -I./include -DTESTEOPC     ./main.cpp -L/usr/local/lib/vallado/                      -lastTime                   -lcurl -o ./bin/testEopc
#	g++ -march=native -Os -std=c++11 -I./include -TESTSATELLITE ./main.cpp -L/usr/local/lib/vallado/ -last2Body -lastMath -lastTime -lcoordFK5 -lSGP4 -lcurl -o ./bin/testSatellite
#	g++ -march=native -Os -std=c++11 -I./include -DTESTSUN      ./main.cpp -L/usr/local/lib/vallado/ -last2Body -lastMath -lastTime -lcoordFK5 -lSGP4 -lcurl -o ./bin/testSun
#	g++ -march=native -Os -std=c++11 -I./include -DTESTSTATION  ./main.cpp -L/usr/local/lib/vallado/ -last2Body -lastMath -lastTime -lcoordFK5 -lSGP4 -lcurl -o ./bin/testStation




