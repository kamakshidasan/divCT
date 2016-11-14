all : serialct stitch

serialct :
	g++ -O3 -std=c++0x ct_serial.cpp -o serialct
	
stitch :
	gcc -O3 -fopenmp stitch.c -o stitch
