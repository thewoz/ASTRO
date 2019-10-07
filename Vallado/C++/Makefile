INCLUDE_INSTALL_PATH=/usr/local/include
LIBRARY_INSTALL_PATH=/usr/local/lib

LIBRARY_NAME=vallado


all: setup SGP4 astTime astMath ast2Body astUtils coordFK5 astIOD EopSpw


setup: clear uninstall
	@mkdir -p ./build/include
	@mkdir -p ./build/lib


install: setup SGP4 astTime astMath ast2Body astUtils coordFK5 astIOD EopSpw
	@mkdir -p $(INCLUDE_INSTALL_PATH)/$(LIBRARY_NAME)/
	@mkdir -p $(LIBRARY_INSTALL_PATH)/
	@cp ./build/include/* $(INCLUDE_INSTALL_PATH)/$(LIBRARY_NAME)/
	@cp ./build/lib/*     $(LIBRARY_INSTALL_PATH)/
	@cp -r ./data         $(INCLUDE_INSTALL_PATH)/$(LIBRARY_NAME)/
	@rm -r ./build


uninstall:
	@rm -rf $(INCLUDE_INSTALL_PATH)/$(LIBRARY_NAME)/
	@rm -rf $(LIBRARY_INSTALL_PATH)/$(LIBRARY_NAME)/


clear:
	@rm -rf $(INCLUDE_INSTALL_PATH)/$(LIBRARY_NAME)
	@rm -rf $(LIBRARY_INSTALL_PATH)/$(LIBRARY_NAME)


# process SGP4
SGP4: setup
	g++ -march=native -Os -c ./src/SGP4/SGP4.cpp -o ./build/SGP4.o
	@ar rcs ./build/lib/libSGP4.a ./build/SGP4.o
	@rm ./build/SGP4.o
	@cp ./src/SGP4/SGP4.h ./build/include


# process astTime
astTime: setup
	g++ -march=native -Os -c ./src/astTime/astTime.cpp -o ./build/astTime.o
	@ar rcs ./build/lib/libastTime.a ./build/astTime.o
	@rm ./build/astTime.o
	@cp ./src/astTime/astTime.h ./build/include


# process astMath
astMath: setup
	g++ -march=native -Os -c ./src/astMath/astMath.cpp -o ./build/astMath.o
	@ar rcs ./build/lib/libastMath.a ./build/astMath.o
	@rm ./build/astMath.o
	@cp ./src/astMath/astMath.h ./build/include


# process ast2Body
ast2Body: setup astTime astMath
	g++ -march=native -Os -c -I./src/astMath -I./src/astTime ./src/ast2Body/ast2Body.cpp -o ./build/ast2Body.o
	@ar rcs ./build/lib/libast2Body.a ./build/ast2Body.o
	@rm ./build/ast2Body.o
	@cp ./src/ast2Body/ast2Body.h ./build/include


# process astUtils
astUtils: setup astMath ast2Body
	g++ -march=native -Os -c -I./src/astMath -I./src/astTime -I./src/ast2Body ./src/astUtils/astUtils.cpp -o ./build/astUtils.o
	@ar rcs ./build/lib/libastUtils.a ./build/astUtils.o
	@rm ./build/astUtils.o
	@cp ./src/astUtils/astUtils.h ./build/include


# process coordFK5
coordFK5: setup astTime astMath ast2Body
	@g++ -march=native -Os -c -I./src/astMath -I./src/astTime -I./src/ast2Body ./src/coordFK5/coordFK5.cpp -o ./build/coordFK5.o
	@ar rcs ./build/lib/libcoordFK5.a ./build/coordFK5.o
	@rm ./build/coordFK5.o
	@cp ./src/coordFK5/coordFK5.h ./build/include


# process astIOD
astIOD: setup astTime astMath ast2Body astUtils
	g++ -march=native -Os -c -I./src/astMath -I./src/astTime -I./src/ast2Body -I./src/astUtils ./src/astIOD/astIOD.cpp -o ./build/astIOD.o
	@ar rcs ./build/lib/libastIOD.a ./build/astIOD.o
	@rm ./build/astIOD.o
	@cp ./src/astIOD/astIOD.h ./build/include


# process EopSpw
EopSpw: setup astTime astMath ast2Body
	g++ -march=native -Os -c -I./src/astMath -I./src/astTime -I./src/ast2Body ./src/EopSpw/EopSpw.cpp -o ./build/EopSpw.o
	@ar rcs ./build/lib/libEopSpw.a ./build/EopSpw.o
	@rm ./build/EopSpw.o
	@cp ./src/EopSpw/EopSpw.h ./build/include
