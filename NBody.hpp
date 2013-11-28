#ifndef NBODY_H_
#define NBODY_H_

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "CLUtil.hpp"

#define GROUP_SIZE 128

#define SAMPLE_VERSION "AMD-APP-SDK-v2.9.214.1"

using namespace appsdk;

class NBody
{
        cl_double setupTime;
        cl_double kernelTime;               /**< time taken to run kernel and read result back */

        cl_float delT;                      /**< dT (timestep) */
        cl_float espSqr;                    /**< Softening Factor*/
        cl_float* initPos;                  /**< initial position */
        cl_float* initVel;                  /**< initial velocity */
        cl_float* vel;                      /**< Output velocity */
        cl_context context;                 /**< CL context */
        cl_device_id *devices;              /**< CL device list */
        cl_mem particlePos[2];              // positions of particles
        cl_mem particleVel[2];              // velocity of particles
        int currentPosBufferIndex;
        unsigned int currentIterationCL;    // current iteration in the simulation
        float* mappedPosBuffer;             // mapped pointer of the position buffer
        int mappedPosBufferIndex;
        cl_command_queue commandQueue;      /**< CL command queue */
        cl_program program;                 /**< CL program */
        cl_kernel kernel;                   /**< CL kernel */
        size_t groupSize;                   /**< Work-Group size */

        int iterations;
        SDKDeviceInfo
        deviceInfo;                /**< Structure to store device information*/
        KernelWorkGroupInfo
        kernelInfo;          /**< Structure to store kernel related info */

        int fpsTimer;
        int timerNumFrames;

        SDKTimer *sampleTimer;      /**< SDKTimer object */

    private:

        float random(float randMax, float randMin);
        float randomGauss(float randMax, float randMin);

    public:

        CLCommandArgs   *sampleArgs;   /**< CLCommand argument class */

        cl_int numParticles;
        bool    isFirstLuanch;
        cl_event glEvent;
        /**
        * Constructor
        * Initialize member variables
        */
        explicit NBody()
            : setupTime(0),
              kernelTime(0),
              delT(60.0f * 60.0f * 24.0f * 365.25f * 1000.0f * 1000.0f),
              espSqr(1e25f),
              initPos(NULL),
              initVel(NULL),
              vel(NULL),
              devices(NULL),
              groupSize(GROUP_SIZE),
              iterations(1),
              currentPosBufferIndex(0),
              currentIterationCL(0),
              mappedPosBuffer(NULL),
              fpsTimer(0),
              timerNumFrames(0),
              isFirstLuanch(true),
              glEvent(NULL)
        {
            sampleArgs = new CLCommandArgs();
            sampleTimer = new SDKTimer();
            sampleArgs->sampleVerStr = SAMPLE_VERSION;
            numParticles = 10000;
        }

        ~NBody();

        /**
        * Allocate and initialize host memory array with random values
        * @return SDK_SUCCESS on success and SDK_FAILURE on failure
        */
        int setupNBody();

        void generateGalaxy(int count, int offset, float branches, float pitch, float cX, float cY, float cZ, float vX, float vY, float vZ);

        /**
         * Override from SDKSample, Generate binary image of given kernel
         * and exit application
        * @return SDK_SUCCESS on success and SDK_FAILURE on failure
         */
        int genBinaryImage();

        /**
        * OpenCL related initialisations.
        * Set up Context, Device list, Command Queue, Memory buffers
        * Build CL kernel program executable
        * @return SDK_SUCCESS on success and SDK_FAILURE on failure
        */
        int setupCL();

        /**
        * Set values for kernels' arguments
        * @return SDK_SUCCESS on success and SDK_FAILURE on failure
        */
        int setupCLKernels();

        /**
        * Enqueue calls to the kernels
        * on to the command queue, wait till end of kernel execution.
        * Get kernel start and end time if timing is enabled
        * @return SDK_SUCCESS on success and SDK_FAILURE on failure
        */
        int runCLKernels();

        /**
        * Reference CPU implementation of Binomial Option
        * for performance comparison
        */
        void nBodyCPUReference(float* currentPos, float* currentVel
                               , float* newPos, float* newVel);


        float* getMappedParticlePositions();
        void releaseMappedParticlePositions();

        /**
        * Override from SDKSample. Print sample stats.
        */
        void printStats();

        /**
        * Override from SDKSample. Initialize
        * command line parser, add custom options
        * @return SDK_SUCCESS on success and SDK_FAILURE on failure
        */
        int initialize();

        /**
        * Override from SDKSample, adjust width and height
        * of execution domain, perform all sample setup
        * @return SDK_SUCCESS on success and SDK_FAILURE on failure
        */
        int setup();

        /**
        * Override from SDKSample
        * Run OpenCL NBody
        * @return SDK_SUCCESS on success and SDK_FAILURE on failure
        */
        int run();

        /**
        * Override from SDKSample
        * Cleanup memory allocations
        * @return SDK_SUCCESS on success and SDK_FAILURE on failure
        */
        int cleanup();

        /**
        * Override from SDKSample
        * Verify against reference implementation
        * @return SDK_SUCCESS on success and SDK_FAILURE on failure
        */
        int verifyResults();


        // init the timer for FPS calculation
        void initFPSTimer()
        {
            timerNumFrames = 0;
            fpsTimer = sampleTimer->createTimer();
            sampleTimer->resetTimer(fpsTimer);
            sampleTimer->startTimer(fpsTimer);
        };

        // calculate FPS
        double getFPS()
        {
            sampleTimer->stopTimer(fpsTimer);
            double elapsedTime = sampleTimer->readTimer(fpsTimer);
            double fps = timerNumFrames/elapsedTime;
            timerNumFrames = 0;
            sampleTimer->resetTimer(fpsTimer);
            sampleTimer->startTimer(fpsTimer);
            return fps;
        };
};

#endif // NBODY_H_
