#include <stdio.h>
#include <inttypes.h>
#include <Windows.h>

/* Include files needed to use VnSensor. */
#include "vn/sensors.h"

void asciiOrBinaryAsyncMessageReceived(void *userData, VnUartPacket *packet, size_t runningIndex);
int processErrorReceived(char* errorMessage, VnError errorCode);

FILE *gnuplotPipe;
char * commandsForGnuplot = "plot 'data.temp' using 1:2 title 'x' with lines,\
								'data.temp' using 1:3 title 'y' with lines,\
								'data.temp' using 1:4 title 'z' with lines";
FILE * temp;

int main(void)
{
	VnSensor vs;
	char modelNumber[30];
	char strConversions[50];
	vec3f ypr;
	YawPitchRollMagneticAccelerationAndAngularRatesRegister reg;
	VpeBasicControlRegister vpeReg;
	uint32_t oldHz, newHz;
	VnAsciiAsync asyncType;
	BinaryOutputRegister bor;
	VnError error;

	// Open a pipe to gnuplot
	temp = fopen("data.temp", "w");

    gnuplotPipe = _popen("gnuplot -persist", "w");
    if (gnuplotPipe == NULL) {
        fprintf(stderr, "Error opening pipe to gnuplot\n");
        return 1;
    }

	const char SENSOR_PORT[] = "COM3";
	const uint32_t SENSOR_BAUDRATE = 115200;

	/* We first need to initialize our VnSensor structure. */
	VnSensor_initialize(&vs);

	/* Now connect to our sensor. */
	if ((error = VnSensor_connect(&vs, SENSOR_PORT, SENSOR_BAUDRATE)) != E_NONE)
		return processErrorReceived("Error connecting to sensor.", error);

	/* Let's query the sensor's model number. */
	if ((error = VnSensor_readModelNumber(&vs, modelNumber, sizeof(modelNumber))) != E_NONE)
	{
		//return processErrorReceived("Error reading model number.", error);
		printf("Error reading model number.\n");
		goto finish;
	}
		
	printf("Model Number: %s\n", modelNumber);

	/* Get some orientation data from the sensor. */
	if ((error = VnSensor_readYawPitchRoll(&vs, &ypr)) != E_NONE)
	{
		//return processErrorReceived("Error reading yaw pitch roll.", error);
		printf("Error reading yaw pitch roll.\n");
		goto finish;
	}
		
	str_vec3f(strConversions, ypr);
	printf("Current YPR: %s\n", strConversions);

	/* Get some orientation and IMU data. */
	if ((error = VnSensor_readYawPitchRollMagneticAccelerationAndAngularRates(&vs, &reg)) != E_NONE)
		return processErrorReceived("Error reading orientation and IMU data.", error);
	str_vec3f(strConversions, reg.yawPitchRoll);
	printf("Current YPR: %s\n", strConversions);
	str_vec3f(strConversions, reg.mag);
	printf("Current Magnetic: %s\n", strConversions);
	str_vec3f(strConversions, reg.accel);
	printf("Current Acceleration: %s\n", strConversions);
	str_vec3f(strConversions, reg.gyro);
	printf("Current Angular Rates: %s\n", strConversions);

	
	if ((error = VnSensor_writeAsyncDataOutputFrequency(&vs, 40, true)) != E_NONE)
	{
		//return processErrorReceived("Error writing async data output frequency.", error);
		printf("Error writing async data output frequency.\n");
		goto finish;
	}
		
	if ((error = VnSensor_readAsyncDataOutputFrequency(&vs, &newHz)) != E_NONE)
	{
		//return processErrorReceived("Error reading async data output frequency.", error);
		printf("Error reading async data output frequency.\n");
		goto finish;
	}
		
	//printf("Old Async Frequency: %d Hz\n", oldHz);
	printf("New Async Frequency: %d Hz\n", newHz);

	BinaryOutputRegister_initialize(
		&bor,ASYNCMODE_PORT1, 20,
		COMMONGROUP_TIMESTARTUP | COMMONGROUP_YAWPITCHROLL,	/* Note use of binary OR to configure flags. */
		TIMEGROUP_NONE, IMUGROUP_NONE, GPSGROUP_NONE, ATTITUDEGROUP_NONE, INSGROUP_NONE, GPSGROUP_NONE);

	if ((error = VnSensor_writeBinaryOutput1(&vs, &bor, true)) != E_NONE)
	{
		//return processErrorReceived("Error writing binary output 1.", error);
		printf("Error writing binary output 1.\n");
		goto finish;
	}

	VnSensor_registerAsyncPacketReceivedHandler(&vs, asciiOrBinaryAsyncMessageReceived, NULL);

	printf("Starting sleep...\n");
	VnThread_sleepSec(50);

	VnSensor_unregisterAsyncPacketReceivedHandler(&vs);
	
	finish:

	 _pclose(gnuplotPipe);

	/* Now disconnect from the sensor since we are finished. */
	if ((error = VnSensor_disconnect(&vs)) != E_NONE)
		return processErrorReceived("Error disconnecting from sensor.", error);

	return 0;
}

void asciiOrBinaryAsyncMessageReceived(void *userData, VnUartPacket *packet, size_t runningIndex)
{
	vec3f ypr;
	char strConversions[50];
	static int cnt = 0;

	/* Silence 'unreferenced formal parameters' warning in Visual Studio. */
	(userData);
	(runningIndex);

	if (VnUartPacket_type(packet) == PACKETTYPE_ASCII && VnUartPacket_determineAsciiAsyncType(packet) == VNYPR)
	{
		VnUartPacket_parseVNYPR(packet, &ypr);
		str_vec3f(strConversions, ypr);
		printf("ASCII Async YPR: %s\n", strConversions);

		return;
	}

	if (VnUartPacket_type(packet) == PACKETTYPE_BINARY)
	{
		uint64_t timeStartup;

		/* First make sure we have a binary packet type we expect since there
		 * are many types of binary output types that can be configured. */
		if (!VnUartPacket_isCompatible(packet,
			COMMONGROUP_TIMESTARTUP | COMMONGROUP_YAWPITCHROLL,
			TIMEGROUP_NONE,
			IMUGROUP_NONE,
			GPSGROUP_NONE,
			ATTITUDEGROUP_NONE,
			INSGROUP_NONE,
      GPSGROUP_NONE))
			/* Not the type of binary packet we are expecting. */
			return;

		/* Ok, we have our expected binary output packet. Since there are many
		 * ways to configure the binary data output, the burden is on the user
		 * to correctly parse the binary packet. However, we can make use of
		 * the parsing convenience methods provided by the Packet structure.
		 * When using these convenience methods, you have to extract them in
		 * the order they are organized in the binary packet per the User Manual. */
		timeStartup = VnUartPacket_extractUint64(packet);
		ypr = VnUartPacket_extractVec3f(packet);

		str_vec3f(strConversions, ypr);
		cnt++;
		//printf("Binary Async TimeStartup: %" PRIu64 "\n", timeStartup);
		printf("%d - Binary Async YPR: %s\n", cnt, strConversions);
		
		// -------------------------------------------------------------

		//fprintf(gnuplotPipe, "plot '-' using 1:1 with lines\n"); // Specifying x and y columns
		static float xmin = 0; // Adjust as needed
   	    static float xmax = 100; // Adjust as needed
   		//if(cnt==1)
			//fprintf(gnuplotPipe, "set xrange [%f:%f]\n", xmin, xmax);	
		
		//fprintf(gnuplotPipe,"plot '-'\n");
		//fprintf(gnuplotPipe, "%d %f \n", cnt, ypr.c[0]);
        //fprintf(gnuplotPipe, "e\n");
        //fflush(gnuplotPipe);

		fprintf(temp, "%d %f %f %f \n", cnt, ypr.c[0], ypr.c[1], ypr.c[2]); 				fflush(temp);
		fprintf(gnuplotPipe, "%s \n", commandsForGnuplot); 								fflush(gnuplotPipe);
		//fprintf(gnuplotPipe, "%s \n", commandsForGnuplot[1]); 								fflush(gnuplotPipe);
		//fprintf(gnuplotPipe, "%s \n", commandsForGnuplot[2]); 								fflush(gnuplotPipe);
		fprintf(gnuplotPipe, "e\n"); 														fflush(gnuplotPipe);

		//Sleep(100);

		// ------------------------------------------------------------

		return;
	}
}

int processErrorReceived(char* errorMessage, VnError errorCode)
{
	char errorCodeStr[100];
	strFromVnError(errorCodeStr, errorCode);
	printf("%s\nERROR: %s\n", errorMessage, errorCodeStr);
	return -1;
}
