all: test_main.o VelocityGrid.o CollisionNodes.o KorobovGenerator.o UtilsRandom.o
	g++ test_main.o VelocityGrid.o CollisionNodes.o KorobovGenerator.o UtilsRandom.o -o test_main
test_main.o: test_main.cpp
	g++ -c test_main.cpp
CollisionNodes.o: CollisionNodes.cpp
	g++ -c CollisionNodes.cpp
KorobovGenerator.o: KorobovGenerator.cpp
	g++ -c KorobovGenerator.cpp
VelocityGrid.o: VelocityGrid.cpp
	g++ -c VelocityGrid.cpp
UtilsRandom.o: UtilsRandom.cpp
	g++ -c UtilsRandom.cpp
clean:
	rm *.o test_main
