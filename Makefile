all: PBWT rPBWT

PBWT: PBWT.cpp
	g++ -std=c++17 -Wshadow -Wall -o PBWT PBWT.cpp -O2 -Wno-unused-result

rPBWT: rPBWT.cpp
	g++ -std=c++17 -Wshadow -Wall -o rPBWT rPBWT.cpp -O2 -Wno-unused-result

clean:
	rm -f PBWT rPBWT

.PHONY: all clean 
