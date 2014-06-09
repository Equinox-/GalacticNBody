#define UNROLL_FACTOR  8
__kernel 
void nbody_sim(__global float4* pos, __global float4* vel
		,int numBodies ,float deltaTime, float epsSqr
		,__global float4* newPosition, __global float4* newVelocity) {

    unsigned int gid = get_global_id(0);
    float4 myPosRoot = pos[gid];
    float4 oldVel = vel[gid];

    if (myPosRoot.w <= 0.0f) {
        newPosition[gid] = myPosRoot;
        newVelocity[gid] = oldVel;
        return;
    }

    float4 acc[5];
    acc[0] = (float4)0.0f;
    acc[1] = (float4)0.0f;
    acc[2] = (float4)0.0f;
    acc[3] = (float4)0.0f;
    acc[4] = (float4)0.0f;


    int i = 0;
    int k = 1;

    int lostMass = 0;
    
    for (; k < 5; ++k) {
        float mult = deltaTime * (float)abs(3 - k) / 2.0f;
        float4 myPos;
        myPos.xyz = myPosRoot.xyz + oldVel.xyz * mult + acc[k-1].xyz * 0.5f * mult * mult;
        for (; (i+UNROLL_FACTOR) < numBodies; ) {
    #pragma unroll UNROLL_FACTOR
            for(int j = 0; j < UNROLL_FACTOR; j++,i++) {
                float4 p = pos[i];
                float4 r;
                r.xyz = p.xyz - myPos.xyz;
                float distSqr = r.x * r.x  +  r.y * r.y  +  r.z * r.z;
                float dist = sqrt(distSqr + epsSqr);
                float invDist = 1.0f / dist;
                float invDistCube = (p.w + ((p.w > 1.0e6f)*(5.0E6f * dist))) * invDist * invDist * invDist;
                float s = invDistCube;

                // accumulate effect of all particles
                acc[k].xyz += (gid != i) * s * r.xyz;
                if (gid != i && distSqr < 1E2) {
                    acc[k].w += p.w;
                    if (p.w > 1.0e6f) {
                        lostMass = 1;
                    }
                }
            }
        }
        // for (; i < numBodies; i++) {
        //     float4 p = pos[i];

        //     float4 r;
        //     r.xyz = p.xyz - myPos.xyz;
        //     float distSqr = r.x * r.x  +  r.y * r.y  +  r.z * r.z;

        //     float invDist = 1.0f / sqrt(distSqr + epsSqr);
        //     float invDistCube = p.w * invDist * invDist * invDist;
        //     float s = invDistCube;

        //     // accumulate effect of all particles
        //     acc[k].xyz += (gid != i) * s * r.xyz;
        //     if (gid != i && distSqr < 2) {
        //         acc[k].w += p.w;
        //         if (p.w > 1.0e6f) {
        //             lostMass = 1;
        //         }
        //     }
        // }
        acc[k].xyz *= 1.567783995250E-28f;
    }

    acc[0] = (float4)0.0f;
    for (k = 1; k<5; ++k) {
        acc[0].xyz += acc[k].xyz * ((float)abs(3 - k) / 6.0f);
        acc[0].w += acc[k].w * ((float)abs(3-k) / 6.0f);
    }

    float4 newVel;
    newVel.xyz = oldVel.xyz + acc[0].xyz * deltaTime;
    newVel.w = oldVel.w;

    if (myPosRoot.w > 1.0e6f) {
        acc[0].x = acc[0].y = acc[0].z = 0;
    }

    // updated position and velocity
    float4 newPos;
    newPos.xyz = myPosRoot.xyz + oldVel.xyz * deltaTime + acc[0].xyz * 0.5f * deltaTime * deltaTime;
    newPos.w = myPosRoot.w;

    if (newPos.w > 1.0e6f && acc[0].w > 0.0f) {
        newPos.w += acc[0].w;
    }


    if (lostMass) {
        newPos.w = 0;
    }

    // write to global memory
    newPosition[gid] = newPos;
    newVelocity[gid] = newVel;
}