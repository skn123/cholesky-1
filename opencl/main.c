#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <time.h>
#include <malloc.h>
#include "spd_matrix.h"
#include "cholesky.h"
#include <CL/cl.h>

void CL_CALLBACK onOpenCLError(const char *errinfo,  const void *private_info,
                               size_t cb, void *user_data)
{
    printf("Error while creating context or working in this context : %s", errinfo);
}

typedef struct DeviceDesc{
	cl_device_id    deviceId;
	cl_device_type  deviceType;
	char*           deviceTypeString;
	char*           deviceName;
} DeviceDesc;

int main(int argc, char *argv[])
{
    srand( time( NULL ) );
	if (!argv[1])
	{
	    printf("Specify matrix dimension.\n");
	    exit(-1);
	}
	int dimension = atoi(argv[1]);

    //miejsce na to co zwracaja funkcje opencl
    cl_int result;

    //platformy
    cl_uint             numEntries = 1;
    cl_platform_id*     platforms;
    cl_uint             numPlatforms;
    //urzadzenia
    cl_uint             maxDevices = 1;
    cl_device_id*       deviceIDs;
    cl_uint             numDevices;
    //kontekst
    cl_context_properties*  properties = 0;
    cl_uint                 usedDevices = 1;
    //kolejka polecen z wlaczona opcja profilera
    cl_command_queue_properties commandQueueProperties = CL_QUEUE_PROFILING_ENABLE;

    //dostepne platformy
    platforms = (cl_platform_id*)malloc(sizeof(cl_platform_id)*numEntries);
    result = clGetPlatformIDs(numEntries, platforms, &numPlatforms);
    if(result != CL_SUCCESS) exit(1);

    //dostepne urzadzenia
    deviceIDs = (cl_device_id*)malloc(maxDevices*sizeof(cl_device_id));
    result = clGetDeviceIDs(platforms[0], CL_DEVICE_TYPE_GPU, maxDevices, deviceIDs, &numDevices);
    if(result != CL_SUCCESS) exit(2);

    //tworzymy kontekst
    //We use the clCreateContext function
    cl_context context = clCreateContext(properties, usedDevices, deviceIDs, &onOpenCLError, NULL, &result);
    if(result != CL_SUCCESS) exit(3);

    //tworzymy kolejke polecen
    cl_command_queue commands = clCreateCommandQueue(context, deviceIDs[0], commandQueueProperties, &result);
    if(result != CL_SUCCESS) exit(4);

    //wyciaga informacje o urzadzeniu
    DeviceDesc device_desc;
    device_desc.deviceId = deviceIDs[0];
    size_t actualSize;
    result = clGetDeviceInfo(deviceIDs[0], CL_DEVICE_TYPE, sizeof(cl_device_type), &device_desc.deviceType, &actualSize);

    //sprawdza typ urzadzenia
    switch(device_desc.deviceType)
    {
         case CL_DEVICE_TYPE_CPU:
            device_desc.deviceTypeString = "Processor";
            break;
         case CL_DEVICE_TYPE_GPU:
            device_desc.deviceTypeString = "Graphics card";
            break;
         case CL_DEVICE_TYPE_ACCELERATOR:
            device_desc.deviceTypeString = "Accelerator";
            break;
         default:
            device_desc.deviceTypeString = "NONE";
            break;
    }

    //Wyciaga nazwe urzadzenia
    size_t deviceNameLength = 4096;
    char* tempDeviceName = (char*)malloc(4096);
    result |= clGetDeviceInfo(deviceIDs[0], CL_DEVICE_NAME, deviceNameLength, tempDeviceName, &actualSize);
    if(result == CL_SUCCESS){
        device_desc.deviceName = (char*)malloc(actualSize);
        memcpy(device_desc.deviceName, tempDeviceName, actualSize);
        free(tempDeviceName);
    }
    if(result != CL_SUCCESS)
    {
        printf("Error while getting device info\n");
        exit(0);
    }

    printf("%s: %s\n", device_desc.deviceTypeString,device_desc.deviceName);

    //ladujemy kod kernela z pliku
    char *kernels;
    long input_file_size;
    FILE *input_file = fopen("kernels.c", "rb");
    fseek(input_file, 0, SEEK_END);
    input_file_size = ftell(input_file);
    rewind(input_file);
    kernels = malloc(input_file_size * (sizeof(char)));
    fread(kernels, sizeof(char), input_file_size, input_file);
    fclose(input_file);

    cl_program program = clCreateProgramWithSource(context, 1, (const char **)&kernels, NULL, &result);
    if(result != CL_SUCCESS) exit(5);

    //kompilujemy kernel
    result = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
    if(result != CL_SUCCESS)
    {
        size_t length;
        char buffer[2048];
        clGetProgramBuildInfo(program, deviceIDs[0], CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &length);
        printf("--- Build log ---\n %s\n", buffer);
        exit(6);
    }


    //wyciaga funkcje kernela ze skompilowanego kodu
    cl_kernel choldc_gpu = clCreateKernel(program, "choldc_gpu", &result);
    if(result != CL_SUCCESS) exit(7);

    printf("Dimension: %d\n", dimension);
    float norm1, norm2, norm3;
    float** A = dmatrix(1, dimension, 1, dimension);
    float** A_clone = dmatrix(1, dimension, 1, dimension);
    float** L = dmatrix(1, dimension, 1, dimension);
    float** L_t = dmatrix(1, dimension, 1, dimension);

    //Generowanie macierzy SPD
	A = generate_random_matrix(A, dimension);
    A = construct_symetric_matrix(A, dimension);
	A = matrix_positive_definite(A, dimension);
    norm1 = frobenius_norm(A, dimension);

    size_t numberOfValues = dimension * dimension;
    size_t sizeOfBuffers = numberOfValues * sizeof(float);
    float* inputDoubles = (float*)malloc(sizeOfBuffers);
    inputDoubles = convert_to_array(A, dimension);

    //tworzy bufor do wczytywania danych na karte
    cl_mem inputBuffer = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeOfBuffers, NULL, &result);
    if(result != CL_SUCCESS) exit(8);

    //kopiujemy dane na karte
    cl_bool     blockingWrite = CL_TRUE;
    size_t      offset = 0;
    cl_event    dataInputCopyEvent;
    cl_event*   eventsToWait = NULL;
    cl_uint     numEvents = 0;

    result = clEnqueueWriteBuffer(commands, inputBuffer, blockingWrite, offset, sizeOfBuffers, inputDoubles, numEvents, eventsToWait, &dataInputCopyEvent);
    if(result != CL_SUCCESS) exit(9);

    //tworzymy bufor do zczytywania danych z karty (takiej samej wielkosci jak input, bo macierz L jest taka sama)
    cl_mem outputBuffer = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeOfBuffers, NULL, &result);
    if(result != CL_SUCCESS) exit(10);

    //ustawiamy argumenty dla funkcji choldc_gpu
    result = 0;
    result |= clSetKernelArg(choldc_gpu, 0, sizeof(cl_mem), &inputBuffer);
    result |= clSetKernelArg(choldc_gpu, 1, sizeof(cl_mem), &outputBuffer);
    result |= clSetKernelArg(choldc_gpu, 2, sizeof(size_t), &dimension);
    if(result != CL_SUCCESS) exit(11);

    //wykonujemy funkcje kernal
    cl_uint   workDim = 1;
    size_t*   globalWorkOffset = NULL;
    size_t    globalWorkSize =  dimension;
    size_t    localWorkSize = dimension;
    cl_event  kernelExecEvent;
              eventsToWait = NULL;
              numEvents = 0;
//    printf("%d\n", CL_KERNEL_WORK_GROUP_SIZE);
    //pobiera info o ilosci watkow
    result = clGetKernelWorkGroupInfo(choldc_gpu, deviceIDs[0], CL_KERNEL_WORK_GROUP_SIZE, sizeof(localWorkSize), &localWorkSize, NULL);
//    if(localWorkSize > dimension) localWorkSize = dimension;
    if(result != CL_SUCCESS) exit(12);

    result = clEnqueueNDRangeKernel(commands, choldc_gpu, workDim, globalWorkOffset, &globalWorkSize, NULL, numEvents, eventsToWait, &kernelExecEvent);
    clWaitForEvents(1, &kernelExecEvent);
    cl_ulong start = 0, end = 0;
    clGetEventProfilingInfo(kernelExecEvent, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, NULL);
    clGetEventProfilingInfo(kernelExecEvent, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &end, NULL);
    cl_double g_NDRangePureExecTimeMs = (cl_double)(end - start)*(cl_double)(1e-09);
    printf("Time method 1: % 20.16lf\n", g_NDRangePureExecTimeMs);
    if(result != CL_SUCCESS)
    {
        printf("%d\n", (unsigned int) result);
        exit(13);
    }
    //czeka az sie wszystko wykona
    clFinish(commands);

    //wyciagamy wyniki
    cl_bool     blockingRead = CL_TRUE;
                offset = 0;
    float*     resultArray;
    cl_event    readResultsEvent;
                eventsToWait = NULL;
                numEvents = 0;

    resultArray = (float*)malloc(numberOfValues * sizeof(float));

    //kopiuje dane z bufora wyjsciowego
    clEnqueueReadBuffer(commands, outputBuffer, blockingRead, offset, sizeOfBuffers, resultArray, numEvents, eventsToWait, &readResultsEvent);

    //wypisuje wynik
    L = convert_to_matrix(resultArray, dimension);
    L_t = clone_matrix(L, dimension);
    L_t = transpose_matrix(L_t, dimension);
    A_clone = multiply(L, L_t, A_clone, dimension);
    norm2 = frobenius_norm(A_clone, dimension);
    printf("Error method 1: % 20.16lf\n", fabs(norm1 - norm2));
//    print_matrix(L, dimension);
    print_matrix_to_file(L, dimension);
//    print_matrix_to_file(L_t, dimension);





//    //Faktoryzacja Choleskyego metoda 1.
//    L = choldc(A, L, dimension);
//    L_t = clone_matrix(L, dimension);
//    L_t = transpose_matrix(L_t, dimension);
//
//    //Oblicza norme nowej macierzy
//    A_clone = multiply(L, L_t, A_clone, dimension);
//    norm2 = frobenius_norm(A_clone, dimension);
//
//    //Faktoryzacja Choleskyego metoda 2.
//    L = choldc2(A, L, dimension);
//    L_t = clone_matrix(L, dimension);
//    L_t = transpose_matrix(L_t, dimension);
//
//    //Oblicza norme nowej macierzy
//    A_clone = multiply(L, L_t, A_clone, dimension);
//    norm3 = frobenius_norm(A_clone, dimension);
//
//    printf("Error method 1: % 20.16lf\n", fabs(norm1 - norm2));
//    printf("Error method 2: % 20.16lf\n", fabs(norm1 - norm3));




    free(platforms);
    free(deviceIDs);
    free(inputDoubles);
    free(resultArray);

return 0;
}
