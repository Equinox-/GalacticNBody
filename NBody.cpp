#include "NBody.hpp"
#include <GL/glut.h>
#include <cmath>
#include <malloc.h>
#include <math.h>

int numBodies;      /**< No. of particles*/
cl_float* pos;      /**< Output position */
void* me;           /**< Pointing to NBody class */
cl_bool display;

float
NBody::random(float randMin, float randMax)
{
    float result;
    result =(float)rand() / (float)RAND_MAX;

    return ((1.0f - result) * randMin + result *randMax);
}

float
NBody::randomGauss(float randMin, float randMax)
{
    float center = (randMax+randMin)/2.0f;
    float result = (2.0f * sqrt((float)rand() / (float)RAND_MAX)) - 1.0f;
    return center + (randMax - center) * result;
}

inline float min(float a, float b) {
    return a<b?a:b;
}

// G = 6.67384e-11f                         m^3 * kg^-1 * s^-2
// G = 6.67384e-11 / (9.4605284E15 ^ 3)     ly^3 * kg^-1 * s^-2
// G = 6.67384e-11 / (9.4605284E15 ^ 3) * 1.9891E30
// G = 1.567783995250E-28
float blackHoleMass = 2.6e6f;
float galacticRadius = 100e3f;
float diskRadius = 1e3;
float galacticBuldge = 10e3f;
void NBody::generateGalaxy(int count, int offset, float branches, float pitch, float cX, float cY, float cZ, float vX, float vY, float vZ) {
    float totalMass = blackHoleMass;
    for(int i = offset + 1; i < offset + count; ++i)
    {
        int index = i << 2;
        int vIndex = index;

        // First 3 values are position in x,y and z direction
        /*float angle = random(0, 6.28);
        float branch = floor(random(0, 1)*branches);
        angle -= (6.28 / branches * branch);
        float magnitude = branch / branches * galacticRadius; // Choose an arm
        magnitude += (galacticRadius*(angle / 6.28f)); // Add curve stuff
        magnitude += randomGauss(-magnitude/10, magnitude/10); // Add some randomness

       

        angle += randomGauss(pow(normalMag,2)-1.75f, 1.75f-pow(normalMag,2)); // More random*/

        float branch = floor(random(0, 1) * branches);
        float angle = 6.28f - (((float)(i-offset)) * 6.28 / count);//pow(random(0, 1), 0.75) * 6.28;
        float magnitude = (6.5f - angle) / 6.28 * galacticRadius;
        angle += (6.28f / branches) * branch;

        float normalMag = magnitude/galacticRadius;
        angle += random(pow(normalMag,5)-1.25f, 1.25f-pow(normalMag,5));

        float mass = randomGauss(0.5, (2.0f - normalMag) * 500.0f);

        initPos[index] = magnitude * cos(angle);
        initPos[index+1] = magnitude * sin(angle);
	    float lpitch = random(-pow(1.0f-normalMag, 5.0f), pow(1.0f-normalMag, 5.0f));
        initPos[index+2] = galacticBuldge * lpitch + (diskRadius * randomGauss(-1.0f,1.0f));

        float vMag = (float) sqrt(1.567783995250E-28 * (totalMass + (2.0E6 * magnitude)) / magnitude);//*/pow((blackHoleMass+(1e5f * normalMag)) * 6.67384e-11f * 1.989e30f / 1e18f / magnitude, 0.55f) * (1.0f - fabs(lpitch));
        initVel[vIndex] = sin(angle) * -vMag;
        initVel[vIndex+1] = cos(angle) * vMag;
        initVel[vIndex+2] = 0;

        // Mass value
        initPos[index + 3] = mass;

        totalMass += mass;
    }

    initPos[offset*4] = initPos[offset*4+1] = initPos[offset*4+2] = 0;
    initPos[offset*4+3] = blackHoleMass;
    initVel[offset*4] = initVel[offset*4+1] = initVel[offset*4+2] = initVel[offset*4+3] = 0;
    for(int i = offset; i < offset + count; ++i)
    {
        int index = 4 * i;
        int vIndex = index;
        float oldX = initPos[index];
        float oldZ = initPos[index+2];

        float oldVX = initVel[vIndex];
        float oldVZ = initVel[vIndex+2];

        initPos[index] = oldX * cos(pitch) + oldZ * sin(pitch) + cX;
        initPos[index+1] += cY;
        initPos[index+2] = oldZ * cos(pitch) + oldX * sin (pitch) + cZ;

        initVel[vIndex] = oldVX * cos(pitch) + oldVZ * sin(pitch) + vX;
        initVel[vIndex+1] += vY;
        initVel[vIndex+2] = oldVZ * cos(pitch) + oldVX * sin (pitch) + vZ;
    }
}

int
NBody::setupNBody()
{
    // make sure numParticles is multiple of group size
    numParticles = (cl_int)(((size_t)numParticles < groupSize) ? groupSize :
                            numParticles);
    numParticles = (cl_int)((numParticles / groupSize) * groupSize);

    numBodies = numParticles;

    initPos = (cl_float*)malloc(numBodies * sizeof(cl_float4));
    initVel = (cl_float*)malloc(numBodies * sizeof(cl_float4));
    CHECK_ALLOCATION(initPos, "Failed to allocate host memory. (initPos)");

    for(int i = 0; i < numBodies; ++i)
    {
        int index = i * 4;
        initPos[index+3] = 1e30;

        float mag = (i == 0 ? 1.0f : random(1.0f, galacticRadius));
        float angle = random(0, 6.28f);
        float angle2 = random(0, 6.28f);
        initPos[index] = mag * cos(angle2) * cos(angle);
        initPos[index+1] = mag * cos(angle2) * sin(angle);
        initPos[index+2] = mag * sin(angle2);

        float velMag = 0.00001f * mag / galacticRadius;
        initVel[index] = velMag * -sin(angle);
        initVel[index+1] = velMag * cos(angle);
        initVel[index+2] = 0;
        initVel[index+3] = 0;
    }

    // initialization of inputs
    int countPer = numBodies;
    generateGalaxy(countPer, 0, 3, 0, 0, 0, 0, 0, 0, 0);
//    generateGalaxy(countPer, countPer, 3, 4.25f, 75000, 0, 75000, -.01f/315576000000.0f, 0, -.01f/315576000000.0f);
    centerObject = 0;
    /*
    for(int i = 0; i < numBodies; ++i)
    {
        int index = 4 * i;
        int vIndex = index;

        // First 3 values are position in x,y and z direction
        float angle = random(0, 6.28);
        float branch = floor(random(0, branches));
        angle -= (6.28 / branches * branch);
        float magnitude = floor(branch) / branches * 40e3f; // Choose an arm
        magnitude += (40e3f*(angle / 6.28f)); // Add curve stuff
        magnitude += randomGauss(-magnitude/10, magnitude/10); // Add some randomness

        angle += randomGauss((magnitude/61e3f)-1.25f, 1.25f-(magnitude/61e3f)); // More random

        float mass = random(0.5, 100);

        initPos[index] = magnitude * cos(angle);
        initPos[index+1] = magnitude * sin(angle);
        initPos[index+2] = random(0, 10 * sqrt(abs(61e3 - magnitude)));
        if (random(-1,1) > 0) {
            initPos[index+2] *= -1.0;
        }

        float vMag = 1.0867674e-12 * pow(magnitude, 0.33);
        initVel[vIndex] = sin(angle) * vMag;
        initVel[vIndex+1] = cos(angle) * -vMag;
        initVel[vIndex+2] = 0;

        // Mass value
        initPos[index + 3] = mass;
    }

    initPos[0] = initPos[1] = initPos[2] = 0;
    initPos[3] = 40000;
    initVel[0] = initVel[1] = initVel[2] = 0;*/
    return SDK_SUCCESS;
}

int
NBody::genBinaryImage()
{
    bifData binaryData;
    binaryData.kernelName = std::string("NBody_Kernels.cl");
    binaryData.flagsStr = std::string("");
    if(sampleArgs->isComplierFlagsSpecified())
    {
        binaryData.flagsFileName = std::string(sampleArgs->flags.c_str());
    }

    binaryData.binaryName = std::string(sampleArgs->dumpBinary.c_str());
    int status = generateBinaryImage(binaryData);
    return status;
}


int
NBody::setupCL()
{
    cl_int status = CL_SUCCESS;

    cl_device_type dType;

    if(sampleArgs->deviceType.compare("cpu") == 0)
    {
        dType = CL_DEVICE_TYPE_CPU;
    }
    else //deviceType = "gpu"
    {
        dType = CL_DEVICE_TYPE_GPU;
        if(sampleArgs->isThereGPU() == false)
        {
            std::cout << "GPU not found. Falling back to CPU device" << std::endl;
            dType = CL_DEVICE_TYPE_CPU;
        }
    }

    /*
     * Have a look at the available platforms and pick either
     * the AMD one if available or a reasonable default.
     */
    cl_platform_id platform = NULL;
    int retValue = getPlatform(platform, sampleArgs->platformId,
                               sampleArgs->isPlatformEnabled());
    CHECK_ERROR(retValue, SDK_SUCCESS, "getPlatform() failed");

    // Display available devices.
    retValue = displayDevices(platform, dType);
    CHECK_ERROR(retValue, SDK_SUCCESS, "displayDevices() failed");


    /*
     * If we could find our platform, use it. Otherwise use just available platform.
     */
    cl_context_properties cps[3] =
    {
        CL_CONTEXT_PLATFORM,
        (cl_context_properties)platform,
        0
    };

    context = clCreateContextFromType(
                  cps,
                  dType,
                  NULL,
                  NULL,
                  &status);
    CHECK_OPENCL_ERROR( status, "clCreateContextFromType failed.");

    // getting device on which to run the sample
    status = getDevices(context, &devices, sampleArgs->deviceId,
                        sampleArgs->isDeviceIdEnabled());
    CHECK_ERROR(status, SDK_SUCCESS, "getDevices() failed");

    {
        // The block is to move the declaration of prop closer to its use
        cl_command_queue_properties prop = 0;
        commandQueue = clCreateCommandQueue(
                           context,
                           devices[sampleArgs->deviceId],
                           prop,
                           &status);
        CHECK_OPENCL_ERROR( status, "clCreateCommandQueue failed.");
    }

    //Set device info of given cl_device_id
    retValue = deviceInfo.setDeviceInfo(devices[sampleArgs->deviceId]);
    CHECK_ERROR(retValue, SDK_SUCCESS, "SDKDeviceInfo::setDeviceInfo() failed");

    /*
    * Create and initialize memory objects
    */
    size_t bufferSize = numBodies * sizeof(cl_float4);
    for (int i = 0; i < 2; i++)
    {
        particlePos[i] = clCreateBuffer(context, CL_MEM_READ_WRITE, bufferSize, 0,
                                        &status);
        CHECK_OPENCL_ERROR(status, "clCreateBuffer failed. (particlePos)");
        particleVel[i] = clCreateBuffer(context, CL_MEM_READ_WRITE, bufferSize, 0,
                                        &status);
        CHECK_OPENCL_ERROR(status, "clCreateBuffer failed. (particleVel)");
    }

    // Initialize position buffer
    status = clEnqueueWriteBuffer(commandQueue,particlePos[0],CL_TRUE,0,bufferSize,
                                  initPos,0,0,NULL);
    CHECK_OPENCL_ERROR(status, "clEnqueueWriteBuffer failed. ");

    // Initialize the velocity buffer
    status = clEnqueueWriteBuffer(commandQueue, particleVel[0], CL_TRUE, 0, bufferSize, initVel, 0, 0, NULL);
    CHECK_OPENCL_ERROR(status, "clEnqueueWriteBuffer failed. ");

    status = clFlush(commandQueue);
    CHECK_OPENCL_ERROR(status, "clFlush failed. ");

    // create a CL program using the kernel source
    buildProgramData buildData;
    buildData.kernelName = std::string("NBody_Kernels.cl");
    buildData.devices = devices;
    buildData.deviceId = sampleArgs->deviceId;
    buildData.flagsStr = std::string("");
    if(sampleArgs->isLoadBinaryEnabled())
    {
        buildData.binaryName = std::string(sampleArgs->loadBinary.c_str());
    }

    if(sampleArgs->isComplierFlagsSpecified())
    {
        buildData.flagsFileName = std::string(sampleArgs->flags.c_str());
    }

    retValue = buildOpenCLProgram(program, context, buildData);
    CHECK_ERROR(retValue, SDK_SUCCESS, "buildOpenCLProgram() failed");

    // get a kernel object handle for a kernel with the given name
    kernel = clCreateKernel(program,"nbody_sim",&status);
    CHECK_OPENCL_ERROR(status, "clCreateKernel failed.");

    return SDK_SUCCESS;
}


int
NBody::setupCLKernels()
{
    cl_int status;

    // Set appropriate arguments to the kernel

    // numBodies
    status = clSetKernelArg(
                 kernel,
                 2,
                 sizeof(cl_int),
                 (void *)&numBodies);
    CHECK_OPENCL_ERROR(status, "clSetKernelArg failed. (numBodies)");

    // time step
    status = clSetKernelArg(
                 kernel,
                 3,
                 sizeof(cl_float),
                 (void *)&delT);
    CHECK_OPENCL_ERROR(status, "clSetKernelArg failed. (delT)");

    // upward Pseudoprobability
    status = clSetKernelArg(
                 kernel,
                 4,
                 sizeof(cl_float),
                 (void *)&espSqr);
    CHECK_OPENCL_ERROR(status, "clSetKernelArg failed. (espSqr)");


    return SDK_SUCCESS;
}


int NBody::runCLKernels()
{
    cl_int status;

    int currentBuffer = currentPosBufferIndex;
    int nextBuffer = (currentPosBufferIndex+1)%2;

    /*
    * Enqueue a kernel run call.
    */
    size_t globalThreads[] = {numBodies};
    size_t localThreads[] = {groupSize};

    // Particle positions
    status = clSetKernelArg(kernel,0,sizeof(cl_mem),
                            (void*) (particlePos+currentBuffer));
    CHECK_OPENCL_ERROR(status, "clSetKernelArg failed. (updatedPos)");

    // Particle velocity
    status = clSetKernelArg(kernel,1,sizeof(cl_mem),
                            (void *) (particleVel+currentBuffer));
    CHECK_OPENCL_ERROR(status, "clSetKernelArg failed. (updatedVel)");

    // Particle positions
    status = clSetKernelArg(kernel,5,sizeof(cl_mem),
                            (void*) (particlePos+nextBuffer));
    CHECK_OPENCL_ERROR(status, "clSetKernelArg failed. (unewPos)");

    // Particle velocity
    status = clSetKernelArg(kernel,6,sizeof(cl_mem),
                            (void*) (particleVel+nextBuffer));
    CHECK_OPENCL_ERROR(status, "clSetKernelArg failed. (newVel)");

    status = clEnqueueNDRangeKernel(commandQueue,kernel,1,NULL,globalThreads,
                                    localThreads,0,NULL,NULL);
    CHECK_OPENCL_ERROR(status, "clEnqueueNDRangeKernel failed.");

    status = clFlush(commandQueue);
    CHECK_OPENCL_ERROR(status, "clFlush failed.");

    currentPosBufferIndex = nextBuffer;
    timerNumFrames++;
    return SDK_SUCCESS;
}

float* NBody::getMappedParticlePositions()
{
    cl_int status;
    mappedPosBufferIndex = currentPosBufferIndex;
    mappedPosBuffer = (float*) clEnqueueMapBuffer(commandQueue,
                      particlePos[mappedPosBufferIndex], CL_TRUE, CL_MAP_READ
                      , 0, numBodies*4*sizeof(float), 0, NULL, NULL, &status);
    return mappedPosBuffer;
}

void NBody::releaseMappedParticlePositions()
{
    if (mappedPosBuffer)
    {
        cl_int status = clEnqueueUnmapMemObject(commandQueue,
                                                particlePos[mappedPosBufferIndex], mappedPosBuffer, 0, NULL, NULL);
        mappedPosBuffer = NULL;
        clFlush(commandQueue);
    }
}

int
NBody::initialize()
{
    // Call base class Initialize to get default configuration
    int status = 0;
    if (sampleArgs->initialize() != SDK_SUCCESS)
    {
        return SDK_FAILURE;
    }

    Option *num_particles = new Option;
    CHECK_ALLOCATION(num_particles,
                     "error. Failed to allocate memory (num_particles)\n");

    num_particles->_sVersion = "x";
    num_particles->_lVersion = "particles";
    num_particles->_description = "Number of particles";
    num_particles->_type = CA_ARG_INT;
    num_particles->_value = &numParticles;

    sampleArgs->AddOption(num_particles);
    delete num_particles;

    Option *num_iterations = new Option;
    CHECK_ALLOCATION(num_iterations,
                     "error. Failed to allocate memory (num_iterations)\n");

    num_iterations->_sVersion = "i";
    num_iterations->_lVersion = "iterations";
    num_iterations->_description = "Number of iterations";
    num_iterations->_type = CA_ARG_INT;
    num_iterations->_value = &iterations;

    sampleArgs->AddOption(num_iterations);
    delete num_iterations;

    return SDK_SUCCESS;
}

int
NBody::setup()
{
    int status = 0;
    if(setupNBody() != SDK_SUCCESS)
    {
        return SDK_FAILURE;
    }

    int timer = sampleTimer->createTimer();
    sampleTimer->resetTimer(timer);
    sampleTimer->startTimer(timer);

    status = setupCL();
    if(status != SDK_SUCCESS)
    {
        return SDK_FAILURE;
    }

    sampleTimer->stopTimer(timer);
    // Compute setup time
    setupTime = (double)(sampleTimer->readTimer(timer));

    display = !sampleArgs->quiet && !sampleArgs->verify;

    return SDK_SUCCESS;
}

/**
* @brief Initialize GL
*/
void
GLInit()
{
    glClearColor(0.0 ,0.0, 0.0, 0.0);
    glClear(GL_COLOR_BUFFER_BIT);
    glClear(GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
}

/**
* @brief Glut Idle function
*/
void
idle()
{
    glutPostRedisplay();
}

/**
* @brief Glut reshape func
*
* @param w numParticles of OpenGL window
* @param h height of OpenGL window
*/
void
reShape(int w,int h)
{
    glViewport(0, 0, w, h);

    glViewport(0, 0, w, h);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluPerspective(45.0f, w/h, 1.0f, 1e27f);
    gluLookAt (0.0, 0.0, -2.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0);
}

int yaw = 0, pitch = 25;
float scale = 50000;// * 9E15;//9.4605284e9f;
void displayfunc()
{
    static int numFrames = 0;

    glClearColor(0.0 ,0.0, 0.0, 0.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glPointSize(1.0);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE);
    glEnable(GL_BLEND);
    glDepthMask(GL_FALSE);

    glColor3f(1.0f, 0.5f, 0.5f);

    NBody *nb = (NBody *)me;
    if (nb->isFirstLuanch)
    {
        nb->runCLKernels();
        nb->isFirstLuanch = false;
        return;
    }


    int numBodies = nb->numParticles;
    float* pos = nb->getMappedParticlePositions();
    nb->runCLKernels();

    glLoadIdentity();
    glRotatef(pitch, 1, 0, 0);
    glRotatef(yaw, 0, 1, 0);
    glRotatef(90, 1, 0, 0);

    float bDist = sqrt(pow(pos[nb->centerObject<<2]-pos[0],2) + pow(pos[(nb->centerObject<<2)+1]-pos[1],2) + pow(pos[(nb->centerObject<<2)+2] - pos[2], 2));
    // printf("Scale %f->%f\n", bDist, scale);

    glTranslated((pos[nb->centerObject<<2]+pos[0])/(2*-scale), (pos[(nb->centerObject<<2)+1]+pos[1])/(2*-scale), (pos[(nb->centerObject<<2)+2]+pos[2])/(2*-scale));

    glBegin(GL_POINTS);
    for(int i = 0; i < numBodies; ++i,pos+=4)
    {
        glPointSize(1.0);
        if (*(pos+3) > 0.0f) {
            if (*(pos+3) > 1E6) {
                glEnd();
                glPointSize(10.0);
                glBegin(GL_POINTS);
                glColor3f(1.0f, 0.5f, 0.5f);
                glVertex4f(*pos,*(pos+1),*(pos+2),scale);
                glEnd();
                glPointSize(1.0);
                glBegin(GL_POINTS);
            } else {

                glColor3f(1.0f, 0.5f, 0.5f);
                glVertex4f(*pos,*(pos+1),*(pos+2),scale);
            }
        } else {
            glColor3f(0.0f,1.0f,0.0f);
            glVertex4f(*pos,*(pos+1),*(pos+2),scale);
        }
    }
    glEnd();
    nb->releaseMappedParticlePositions();

    glFlush();
    glutSwapBuffers();

    numFrames++;
    if (numFrames >= 100)
    {
        char buf[256];
        sprintf(buf, "OpenCL N-body: %.02f FPS, %d bodies"
                , (float)nb->getFPS(), nb->numParticles);
        glutSetWindowTitle(buf);
        numFrames = 0;
    }
}

// keyboard function
void
keyboardFunc(unsigned char key, int mouseX, int mouseY)
{
    switch(key)
    {
        // If the user hits escape or Q, then exit

        // ESCAPE_KEY = 27
    case 27:
    {
        if(((NBody*)me)->cleanup() != SDK_SUCCESS)
        {
            exit(1);
        }
        else
        {
            exit(0);
        }
    }
    case 'q':
        scale /= 1.50f;
        break;
    case 'e':
        scale *= 1.50f;
        break;
    case 'w':
        pitch++;
        break;
    case 's':
        pitch--;
        break;
    case 'a':
        yaw++;
        break;
    case 'd':
        yaw--;
        break;
    default:
        break;
    }
}

int
NBody::cleanup()
{
    // Releases OpenCL resources (Context, Memory etc.)
    cl_int status;

    status = clReleaseKernel(kernel);
    CHECK_OPENCL_ERROR(status, "clReleaseKernel failed.(kernel)");

    status = clReleaseProgram(program);
    CHECK_OPENCL_ERROR(status, "clReleaseProgram failed.(program)");

    for (int i = 0; i < 2; i++)
    {
        status = clReleaseMemObject(particlePos[i]);
        CHECK_OPENCL_ERROR(status, "clReleaseMemObject failed.(particlePos)");
        status = clReleaseMemObject(particleVel[i]);
        CHECK_OPENCL_ERROR(status, "clReleaseMemObject failed.(particleVel)");
    }

    status = clReleaseCommandQueue(commandQueue);
    CHECK_OPENCL_ERROR(status, "clReleaseCommandQueue failed.(commandQueue)");

    status = clReleaseContext(context);
    CHECK_OPENCL_ERROR(status, "clReleaseContext failed.(context)");

    return SDK_SUCCESS;
}

NBody::~NBody()
{
    if (this->glEvent)
    {
        clReleaseEvent(this->glEvent);
    }
    // release program resources
    FREE(initPos);

    FREE(initVel);

#if defined (_WIN32)
    ALIGNED_FREE(pos);
#else
    FREE(pos);
#endif

#if defined (_WIN32)
    ALIGNED_FREE(vel);
#else
    FREE(vel);
#endif

    FREE(devices);
}


int
main(int argc, char * argv[])
{
    int status = 0;
    NBody clNBody;
    me = &clNBody;

    if(clNBody.initialize() != SDK_SUCCESS)
    {
        return SDK_FAILURE;
    }

    if (clNBody.sampleArgs->parseCommandLine(argc, argv) != SDK_SUCCESS)
    {
        return SDK_FAILURE;
    }

    if(clNBody.sampleArgs->isDumpBinaryEnabled())
    {
        return clNBody.genBinaryImage();
    }

    status = clNBody.setup();
    if(status != SDK_SUCCESS)
    {
        return SDK_FAILURE;
    }

    if(clNBody.setupCLKernels() != SDK_SUCCESS)
    {
        return SDK_FAILURE;
    }

    if(display)
    {
        // Run in  graphical window if requested
        glutInit(&argc, argv);
        glutInitWindowPosition(100,1280);
        glutInitWindowSize(1280,720);
        glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE );
        glutCreateWindow("OpenCL N-Body");
        GLInit();
        glutDisplayFunc(displayfunc);
        glutReshapeFunc(reShape);
        glutIdleFunc(idle);
        glutKeyboardFunc(keyboardFunc);
        sleep(10);
        glutMainLoop();
    }

    status = clNBody.cleanup();
    CHECK_ERROR(status, SDK_SUCCESS, "Sample CleanUP Failed");

    return SDK_SUCCESS;
}
