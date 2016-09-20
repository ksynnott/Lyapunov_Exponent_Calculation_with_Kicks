CC = g++ 

# compiler flags:
#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings
CFLAGS  = -g -Wall

# the build target executable:
TARGET = AMproj
XHEADERS = Attach.h

OBJECTS = $(TARGET).o GramSchmidt.o RungeKutta.o Lyapunov.o

all: $(TARGET)

$(TARGET): $(OBJECTS) 
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJECTS)
	
$(TARGET).o: $(TARGET).cpp $(XHEADERS)
	$(CC) $(CFLAGS) -c $(TARGET).cpp
	
GramSchmidt.o: GramSchmidt.cpp GramSchmidt.h $(XHEADERS)	
	$(CC) $(CFLAGS) -c GramSchmidt.cpp
	
RungeKutta.o: RungeKutta.cpp RungeKutta.h $(XHEADERS)	
	$(CC) $(CFLAGS) -c RungeKutta.cpp
	
	
Lyapunov.o: Lyapunov.cpp Lyapunov.h $(XHEADERS)	
	$(CC) $(CFLAGS) -c Lyapunov.cpp
	
run:
	./$(TARGET)

clean:
	$(RM) $(TARGET) $(OBJECTS)
	
plot:
	python PlotLypunov.py
	
RMPlots:
	rm *.txt
	
MaxPlot:
	python MaxPlotter.py

LyapunovPlot:
	python PlotLypunov.py
	
Cplot:
	python ConvergencePlot.py 
	
runplot:
	./$(TARGET)
	python plotter.py