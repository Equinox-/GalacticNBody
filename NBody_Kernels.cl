#define UNROLL_FACTOR  8
__kernel 
void nbody_sim(__global float4* pos, __global float4* vel
		,int numBodies ,float deltaTime, float epsSqr
		,__global float4* newPosition, __global float4* newVelocity) {

    unsigned int gid = get_global_id(0);
    float4 myPosRoot = pos[gid];
    float4 oldVel = vel[gid];

    float4 acc[5];
    acc[0] = (float4)0.0f;
    acc[1] = (float4)0.0f;
    acc[2] = (float4)0.0f;
    acc[3] = (float4)0.0f;
    acc[4] = (float4)0.0f;


    int i = 0;
    int k = 1;
    for (; k < 5; ++k) {
        float mult = deltaTime * abs(3 - k) / 2.0f;
        float4 myPos;
        myPos.xyz = myPosRoot.xyz + oldVel.xyz * mult + acc[k-1].xyz * 0.5f * mult * mult;
        for (; (i+UNROLL_FACTOR) < numBodies; ) {
    #pragma unroll UNROLL_FACTOR
            for(int j = 0; j < UNROLL_FACTOR; j++,i++) {
                float4 p = pos[i];
                float4 r;
                r.xyz = p.xyz - myPos.xyz;
                float distSqr = r.x * r.x  +  r.y * r.y  +  r.z * r.z;

                float invDist = 1.0f / sqrt(distSqr + epsSqr);
                float invDistCube = p.w * invDist * invDist * invDist;
                float s = invDistCube;

                // accumulate effect of all particles
                acc[k].xyz += s * r.xyz;
            }
        }
        for (; i < numBodies; i++) {
            float4 p = pos[i];

            float4 r;
            r.xyz = p.xyz - myPos.xyz;
            float distSqr = r.x * r.x  +  r.y * r.y  +  r.z * r.z;

            float invDist = 1.0f / sqrt(distSqr + epsSqr);
            float invDistCube = p.w * invDist * invDist * invDist;
            float s = invDistCube;

            // accumulate effect of all particles
            acc[k].xyz += s * r.xyz;
        }
        acc[k].xyz *= 6.67384e-11f / 1e18f;
    }

    acc[0] = (float4)0.0f;
    for (k = 1; k<5; ++k) {
        acc[0].xyz += acc[k].xyz * (abs(3 - k) / 6.0f);
    }

    // updated position and velocity
    float4 newPos;
    newPos.xyz = myPosRoot.xyz + oldVel.xyz * deltaTime + acc[0].xyz * 0.5f * deltaTime * deltaTime;
    newPos.w = myPosRoot.w;

    float4 newVel;
    newVel.xyz = oldVel.xyz + acc[0].xyz * deltaTime;
    newVel.w = oldVel.w;

    // write to global memory
    newPosition[gid] = newPos;
    newVelocity[gid] = newVel;
}